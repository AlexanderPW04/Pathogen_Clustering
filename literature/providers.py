"""Replaceable extraction providers with one interface.

- ClaudeProvider: the Anthropic Messages API with structured JSON outputs.
- MockProvider: deterministic fixture responses for tests and the offline demo.
- LocalModelProvider: a documented HTTP contract for a future locally hosted
  model (for example one served from ASU SOL); not trained or tested here.

Every provider returns a ProviderResult whose data is checked against the same
JSON Schema in Python, so the rest of the pipeline does not change when the
provider does.
"""

import copy
import hashlib
import json
from pathlib import Path
from tempfile import gettempdir

import requests

from literature.config import FALLBACK_WORST_CASE, MODEL_PRICES
from literature.schema import PROMPT_VERSION, schema_problems

FALLBACK_BETA = "server-side-fallback-2026-07-01"


class ProviderError(RuntimeError):
    """A request failed; `fatal` means the run should stop instead of continuing."""

    def __init__(self, message, fatal=True):
        super().__init__(message)
        self.fatal = fatal


def approximate_tokens(*texts):
    """Rough count (about 3 characters per token) when no counting API is available."""

    return sum(len(text) for text in texts) // 3 + 50


def result(data=None, input_tokens=0, output_tokens=0, cost=0.0, model="", error=None, cached=False):
    return {"data": data, "input_tokens": input_tokens, "output_tokens": output_tokens,
            "cost_usd": cost, "model": model, "error": error, "cached": cached}


def checked(data, schema):
    """Return (data, error) after validating structured output again in Python."""

    problems = schema_problems(data, schema)
    if problems:
        return None, "The response did not match the expected structure: " + "; ".join(problems[:3])
    return data, None


class ClaudeProvider:
    name = "Claude API (Anthropic Messages API, structured outputs)"
    is_live = True

    def __init__(self, settings, client=None):
        self.settings = settings
        self.model = settings["model"]
        self._client = client

    def client(self):
        if self._client is None:
            import anthropic  # imported only when a live request is approved
            # The SDK retries connection errors, 408, 409, 429, and 5xx twice
            # with exponential backoff; invalid requests are not retried.
            self._client = anthropic.Anthropic(timeout=120.0, max_retries=2)
        return self._client

    def params(self, system, user, schema, max_tokens):
        return {"model": self.model, "max_tokens": max_tokens, "system": system,
                "messages": [{"role": "user", "content": user}],
                "output_config": {"format": {"type": "json_schema", "schema": schema},
                                  "effort": self.settings["effort"]}}

    def count_tokens(self, system, user, schema):
        """Input tokens from the free token-counting endpoint, or an approximation."""

        import anthropic
        try:
            counted = self.client().messages.count_tokens(
                model=self.model, system=system, messages=[{"role": "user", "content": user}],
                output_config={"format": {"type": "json_schema", "schema": schema}})
            return counted.input_tokens, "Counted with the Anthropic token-counting API"
        except anthropic.APIError:
            return approximate_tokens(system, user), "Approximate (token counting unavailable)"

    def run(self, task, system, user, schema, max_tokens, article_id=""):
        import anthropic
        params = self.params(system, user, schema, max_tokens)
        try:
            if self.settings["refusal_fallback"]:
                # If a safety classifier declines, the API re-runs the request on
                # Anthropic's recommended fallback model inside the same call.
                response = self.client().beta.messages.create(**params, betas=[FALLBACK_BETA],
                                                              fallbacks="default")
            else:
                response = self.client().messages.create(**params)
        except (anthropic.AuthenticationError, anthropic.PermissionDeniedError) as error:
            raise ProviderError(f"Claude API rejected the credentials (HTTP {error.status_code}).") from None
        except anthropic.NotFoundError:
            raise ProviderError(f"Model {self.model} was not found for this API key.") from None
        except anthropic.BadRequestError as error:
            raise ProviderError(f"Claude API rejected the request as invalid (HTTP {error.status_code}).") from None
        except anthropic.RateLimitError:
            raise ProviderError("Claude API rate limit persisted after retries.") from None
        except anthropic.APIStatusError as error:
            raise ProviderError(f"Claude API error (HTTP {error.status_code}) after retries.") from None
        except anthropic.APIConnectionError:
            raise ProviderError("Could not connect to the Claude API.") from None

        input_tokens, output_tokens, cost = usage_cost(response, self.model)
        served = getattr(response, "model", self.model)
        if response.stop_reason == "refusal":
            return result(None, input_tokens, output_tokens, cost, served,
                          "The model declined this request (refusal); nothing was extracted.")
        if response.stop_reason == "max_tokens":
            return result(None, input_tokens, output_tokens, cost, served,
                          "The response reached the output-token limit and was discarded.")
        text = next((block.text for block in response.content if getattr(block, "type", "") == "text"), "")
        try:
            data = json.loads(text)
        except ValueError:
            return result(None, input_tokens, output_tokens, cost, served, "The response was not valid JSON.")
        data, error = checked(data, schema)
        return result(data, input_tokens, output_tokens, cost, served, error)


def price_for(model):
    return MODEL_PRICES.get(model, FALLBACK_WORST_CASE)


def usage_cost(response, model):
    """(input tokens, output tokens, estimated USD) from a response's usage."""

    usage = response.usage
    hops = getattr(usage, "iterations", None) or [usage]
    input_total = output_total = 0
    cost = 0.0
    for hop in hops:
        hop_model = getattr(hop, "model", None) or getattr(response, "model", None) or model
        input_price, output_price = price_for(hop_model)
        plain = getattr(hop, "input_tokens", 0) or 0
        written = getattr(hop, "cache_creation_input_tokens", 0) or 0
        read = getattr(hop, "cache_read_input_tokens", 0) or 0
        output = getattr(hop, "output_tokens", 0) or 0
        input_total += plain + written + read
        output_total += output
        # Cache writes cost 1.25x input; cache reads are charged at 0.1x (an upper bound).
        cost += (plain + 1.25 * written + 0.1 * read) * input_price / 1e6 + output * output_price / 1e6
    return input_total, output_total, cost


class MockProvider:
    """Returns fixture responses keyed by stage and article. Not an AI model."""

    name = "Offline mock provider (fixture responses, not AI output)"
    is_live = False
    model = "mock-fixture-v1"

    def __init__(self, path):
        self.responses = json.loads(Path(path).read_text(encoding="utf-8"))

    def count_tokens(self, system, user, schema):
        return approximate_tokens(system, user), "Approximate (characters / 3)"

    def run(self, task, system, user, schema, max_tokens, article_id=""):
        if task == "association":
            data = self.association(user, self.responses.get("association", {}).get(article_id, []))
        else:
            empty = {"facts": [], "context_statements": [], "insufficient_evidence_samples": [], "notes": ""}
            data = copy.deepcopy(self.responses.get("extraction", {}).get(article_id, empty))
        data, error = checked(data, schema)
        return result(data, approximate_tokens(system, user), 0, 0.0, self.model, error)


    @staticmethod
    def association(user, rules):
        """Decisions for the facts in the request; facts without a fixture rule are accepted."""

        facts = json.loads(user.split("<facts>", 1)[1].split("</facts>", 1)[0])
        decisions = []
        for fact in facts:
            rule = next((rule for rule in rules if rule["Sample_ID"] == fact["Sample_ID"]
                         and rule["Field"] == fact["Field"]), None)
            decisions.append({"Fact_ID": fact["Fact_ID"],
                              "Linkage": rule["Linkage"] if rule else "Specific to this sample",
                              "Explanation": rule["Explanation"] if rule else
                              "Mock default: the fixture records no concern for this fact."})
        return {"decisions": decisions}


class LocalModelProvider:
    """Future locally hosted model behind an HTTP endpoint (LOCAL_MODEL_URL).

    Request JSON: {"contract": PROMPT_VERSION, "task", "system", "user",
    "schema", "max_output_tokens"}. Response JSON: {"data": <object matching
    schema>, "model": "<name and version>", "usage": {"input_tokens",
    "output_tokens"}}. Cost is reported as zero (local compute).
    """

    name = "Locally hosted model (HTTP contract)"
    is_live = True

    def __init__(self, url, session=None):
        if not url.startswith(("http://", "https://")):
            raise ValueError("LOCAL_MODEL_URL must be an http(s) URL.")
        self.url = url
        self.model = "local-model"
        self.session = session or requests.Session()

    def count_tokens(self, system, user, schema):
        return approximate_tokens(system, user), "Approximate (characters / 3)"

    def run(self, task, system, user, schema, max_tokens, article_id=""):
        payload = {"contract": PROMPT_VERSION, "task": task, "system": system, "user": user,
                   "schema": schema, "max_output_tokens": max_tokens}
        try:
            response = self.session.post(self.url, json=payload, timeout=180)
            response.raise_for_status()
            body = response.json()
        except (requests.RequestException, ValueError):
            raise ProviderError("The local model service did not return a usable response.") from None
        usage = body.get("usage", {}) if isinstance(body, dict) else {}
        self.model = str(body.get("model", self.model)) if isinstance(body, dict) else self.model
        data, error = checked(body.get("data") if isinstance(body, dict) else None, schema)
        return result(data, int(usage.get("input_tokens", 0) or 0), int(usage.get("output_tokens", 0) or 0),
                      0.0, self.model, error)


class ExtractionCache:
    """JSON files keyed by provider model, prompt version, stage, and exact request text."""

    def __init__(self, folder=None):
        self.folder = Path(folder) if folder else Path(gettempdir()) / "genomic_surveillance_ai_cache"

    def key(self, model, task, system, user, provider=""):
        text = json.dumps([provider, model, PROMPT_VERSION, task, system, user], ensure_ascii=False)
        return hashlib.sha256(text.encode("utf-8")).hexdigest()

    def get(self, key):
        path = self.folder / f"{key}.json"
        try:
            stored = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            return None
        return result(stored["data"], stored["input_tokens"], stored["output_tokens"], 0.0,
                      stored["model"], cached=True)

    def put(self, key, value):
        if value["data"] is None:
            return  # failures are never cached
        try:
            self.folder.mkdir(parents=True, exist_ok=True)
            stored = {name: value[name] for name in ["data", "input_tokens", "output_tokens", "model"]}
            (self.folder / f"{key}.json").write_text(json.dumps(stored), encoding="utf-8")
        except OSError:
            pass  # caching is an optimization only
