"""Settings for literature enrichment, read only from the server environment.

Live (billable) Claude requests stay off unless AI_ENRICHMENT_ENABLED is true
and ANTHROPIC_API_KEY is set. The key itself is never returned, logged, or
sent to the browser; only "configured" / "not configured" is reported.
"""

import os
from pathlib import Path

ENV_FILE = Path(__file__).resolve().parents[1] / ".env"  # gitignored, optional

DEFAULT_MODEL = "claude-opus-5-5"

# Official first-party prices in USD per million tokens, checked on
# https://platform.claude.com/docs/en/about-claude/pricing on 2026-11-02.
# Prices change: set AI_PRICE_INPUT_PER_MTOK / AI_PRICE_OUTPUT_PER_MTOK to
# override, and re-check the page before relying on an estimate.
PRICES_CHECKED = "2026-11-02"
MODEL_PRICES = {
    "claude-opus-5-5": (4.00, 20.00),
    "claude-sonnet-5-5": (2.00, 10.00),
    "claude-haiku-5-5": (0.10, 0.50),  # prompts up to 100,000 tokens
    "claude-opus-4-8": (5.00, 25.00),
}
# A declined request can be re-run on another model (server-side fallback).
# Budget checks assume the most expensive likely fallback, Claude Opus 4.8.
FALLBACK_WORST_CASE = (5.00, 25.00)
FALLBACK_MODELS = {"claude-opus-5-5", "claude-sonnet-5-5", "claude-opus-5", "claude-fable-5-1"}

TRUE_WORDS = {"1", "true", "yes", "on"}


def load_env_file(path=ENV_FILE):
    """Read KEY=VALUE lines from a local .env file without overriding real settings."""

    if not path.is_file():
        return
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        key = key.strip()
        if key and key not in os.environ:
            os.environ[key] = value.strip().strip('"').strip("'")


def number_setting(name, default, minimum=0):
    """A non-negative number from the environment, or the default when invalid."""

    try:
        value = float(os.getenv(name, default))
    except ValueError:
        return float(default)
    return value if value >= minimum else float(default)


def ai_settings():
    """Current settings as a dictionary (never includes the API key)."""

    load_env_file()
    model = os.getenv("ANTHROPIC_MODEL", "").strip() or DEFAULT_MODEL
    prices = MODEL_PRICES.get(model)
    price_source = f"official price list checked {PRICES_CHECKED}" if prices else None
    if os.getenv("AI_PRICE_INPUT_PER_MTOK") and os.getenv("AI_PRICE_OUTPUT_PER_MTOK"):
        prices = (number_setting("AI_PRICE_INPUT_PER_MTOK", 0),
                  number_setting("AI_PRICE_OUTPUT_PER_MTOK", 0))
        price_source = "AI_PRICE_*_PER_MTOK environment override"
    fallback = (os.getenv("AI_REFUSAL_FALLBACK", "true").strip().casefold() in TRUE_WORDS
                and model in FALLBACK_MODELS)
    return {
        "enabled": os.getenv("AI_ENRICHMENT_ENABLED", "").strip().casefold() in TRUE_WORDS,
        "key_configured": bool(os.getenv("ANTHROPIC_API_KEY", "").strip()),
        "model": model,
        "effort": os.getenv("AI_EFFORT", "low").strip() or "low",
        "prices": prices,
        "price_source": price_source,
        "refusal_fallback": fallback,
        "max_cost_usd": number_setting("AI_MAX_COST_USD_PER_RUN", 1.00),
        "max_cost_usd_per_day": number_setting("AI_MAX_COST_USD_PER_DAY", 3.00),
        "max_articles": int(number_setting("AI_MAX_ARTICLES_PER_RUN", 5, 1)),
        "max_requests": int(number_setting("AI_MAX_REQUESTS_PER_RUN", 10, 1)),
        "max_input_tokens": int(number_setting("AI_MAX_INPUT_TOKENS_PER_RUN", 60000, 1)),
        "max_output_tokens": int(number_setting("AI_MAX_OUTPUT_TOKENS_PER_REQUEST", 4000, 256)),
        "access_code_set": bool(os.getenv("AI_ACCESS_CODE", "").strip()),
        "local_model_url": os.getenv("LOCAL_MODEL_URL", "").strip(),
        "cache_dir": os.getenv("AI_CACHE_DIR", "").strip(),
    }


def ai_status():
    """A short, credential-free description for the diagnostics panel."""

    settings = ai_settings()
    if not settings["enabled"]:
        state = "Disabled"
        message = "AI enrichment is off (AI_ENRICHMENT_ENABLED is not true). No paid requests are possible."
    elif not settings["key_configured"]:
        state = "Not configured"
        message = "AI enrichment is enabled but ANTHROPIC_API_KEY is not set on this server."
    elif settings["prices"] is None:
        state = "Not configured"
        message = (f"No price is known for {settings['model']}; set AI_PRICE_INPUT_PER_MTOK and "
                   "AI_PRICE_OUTPUT_PER_MTOK so spending can be bounded.")
    else:
        state = "Configured"
        message = (f"Live requests allowed after explicit approval; model {settings['model']}, "
                   f"limit ${settings['max_cost_usd']:.2f} per run.")
    return {"status": state, "message": message, "model": settings["model"]}


def ncbi_settings():
    """Optional NCBI contact details; requests work without them at 3 per second."""

    load_env_file()
    return {"api_key": os.getenv("NCBI_API_KEY", "").strip(),
            "email": os.getenv("NCBI_EMAIL", "").strip(),
            "max_requests": int(number_setting("LITERATURE_MAX_REQUESTS", 80, 1)),
            "max_samples": int(number_setting("LITERATURE_MAX_SAMPLES", 25, 1))}
