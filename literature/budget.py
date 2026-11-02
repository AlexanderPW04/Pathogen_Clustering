"""Per-run and per-day limits on requests, tokens, and estimated spending.

Before each request the worst case is checked: the counted (or estimated)
input tokens plus the full output-token allowance, priced at the configured
model's rate (or the fallback model's rate when that could be higher). A
request that could cross a limit is not sent.
"""

from datetime import date
from threading import Lock

from literature.config import FALLBACK_WORST_CASE

DAILY_LEDGER = {"date": None, "spent": 0.0}
DAILY_LOCK = Lock()


def spent_today():
    with DAILY_LOCK:
        if DAILY_LEDGER["date"] != date.today():
            DAILY_LEDGER["date"] = date.today()
            DAILY_LEDGER["spent"] = 0.0
        return DAILY_LEDGER["spent"]


def add_daily_spending(amount):
    spent_today()
    with DAILY_LOCK:
        DAILY_LEDGER["spent"] += amount


class RunBudget:
    def __init__(self, settings, live):
        self.live = live
        self.max_requests = settings["max_requests"]
        self.max_input_tokens = settings["max_input_tokens"]
        self.max_output_tokens = settings["max_output_tokens"]
        self.max_cost = settings["max_cost_usd"]
        self.max_daily_cost = settings["max_cost_usd_per_day"]
        prices = settings["prices"] or FALLBACK_WORST_CASE
        if settings.get("refusal_fallback"):
            prices = (max(prices[0], FALLBACK_WORST_CASE[0]), max(prices[1], FALLBACK_WORST_CASE[1]))
        self.prices = prices if live else (0.0, 0.0)
        self.requests = 0
        self.cached = 0
        self.input_tokens = 0
        self.output_tokens = 0
        self.cost = 0.0
        self.stopped = None

    def worst_case_cost(self, input_tokens):
        return (input_tokens * self.prices[0] + self.max_output_tokens * self.prices[1]) / 1e6

    def refusal(self, input_tokens):
        """Why the next request must not be sent, or None if it fits every limit."""

        if self.requests >= self.max_requests:
            return f"request limit reached ({self.max_requests} per run)"
        if self.input_tokens + input_tokens > self.max_input_tokens:
            return f"input-token limit would be exceeded ({self.max_input_tokens:,} per run)"
        worst = self.worst_case_cost(input_tokens)
        if self.cost + worst > self.max_cost:
            return f"spending limit could be exceeded (${self.max_cost:.2f} per run)"
        if self.live and spent_today() + worst > self.max_daily_cost:
            return f"daily spending limit could be exceeded (${self.max_daily_cost:.2f} per day)"
        return None

    def record(self, outcome):
        if outcome["cached"]:
            self.cached += 1
            return
        self.requests += 1
        self.input_tokens += outcome["input_tokens"]
        self.output_tokens += outcome["output_tokens"]
        self.cost += outcome["cost_usd"]
        if self.live:
            add_daily_spending(outcome["cost_usd"])

    def summary(self):
        return {"Requests sent": self.requests, "Cached results reused": self.cached,
                "Input tokens": self.input_tokens, "Output tokens": self.output_tokens,
                "Estimated cost (USD)": round(self.cost, 4), "Stopped because": self.stopped or ""}
