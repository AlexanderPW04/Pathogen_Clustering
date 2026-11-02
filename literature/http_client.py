"""A small polite HTTP client for NCBI and PubTator (free public services).

- Rate limit: one limiter per service, shared by every client in this server
  process, so two searches at once still stay under the service's limit
  (E-utilities: 3 requests/second, 10 with NCBI_API_KEY; PubTator3: its API
  page asks for no more than three requests per second).
- Request budget: a cap on requests per run; retries count toward it.
- Retries: only for temporary failures (connection errors, HTTP 429 and 5xx),
  at most `retries` times, waiting the server's Retry-After time when given
  (up to MAX_RETRY_AFTER seconds) or 1 s, 2 s, 4 s otherwise.
- Cache: successful responses are kept for an hour and shared between runs,
  so repeating a search does not ask the service again.
- Timeouts: separate connection and read timeouts.
It never follows URLs supplied by article text and never logs credentials.
"""

import time
from email.utils import parsedate_to_datetime
from threading import Lock

import requests

RETRY_STATUS = {429, 500, 502, 503, 504}
MAX_RETRY_AFTER = 30  # seconds; a longer requested wait ends the attempt instead
CACHE_SECONDS = 3600
CACHE_LIMIT = 500


class RequestLimitReached(RuntimeError):
    """The per-run request cap was reached; no further requests are sent."""


class ServiceUnavailable(RuntimeError):
    """A service could not be reached or returned an error after retries."""


class RateLimiter:
    """Keeps at least `interval` seconds between the starts of requests to one service."""

    def __init__(self, interval):
        self.interval = interval
        self.last = 0.0
        self.lock = Lock()

    def wait(self, interval=None):
        with self.lock:
            delay = (interval if interval is not None else self.interval) - (time.monotonic() - self.last)
            if delay > 0:
                time.sleep(delay)
            self.last = time.monotonic()


class ResponseCache:
    """Successful responses for a limited time, safe to share between threads."""

    def __init__(self, seconds=CACHE_SECONDS, limit=CACHE_LIMIT):
        self.seconds, self.limit = seconds, limit
        self.items = {}
        self.lock = Lock()

    def get(self, key):
        with self.lock:
            item = self.items.get(key)
            if item is None or time.monotonic() - item[0] > self.seconds:
                self.items.pop(key, None)
                return None
            return item[1]

    def put(self, key, value):
        with self.lock:
            if len(self.items) >= self.limit:
                del self.items[min(self.items, key=lambda name: self.items[name][0])]
            self.items[key] = (time.monotonic(), value)


# Shared by all clients that use the real network (one entry per service).
SHARED_LIMITERS = {}
SHARED_CACHE = ResponseCache()
SHARED_LOCK = Lock()


def shared_limiter(name, interval):
    with SHARED_LOCK:
        limiter = SHARED_LIMITERS.setdefault(name, RateLimiter(interval))
        limiter.interval = min(limiter.interval, interval) if limiter.interval else interval
        return limiter


def retry_after_seconds(response):
    """Seconds requested by a Retry-After header (number or HTTP date), or None."""

    value = (getattr(response, "headers", None) or {}).get("Retry-After")
    if not value:
        return None
    try:
        return max(0.0, float(value))
    except ValueError:
        pass
    try:
        when = parsedate_to_datetime(value)
        return max(0.0, when.timestamp() - time.time())
    except (TypeError, ValueError):
        return None


def cache_key(url, params, expect):
    items = []
    for name, value in sorted((params or {}).items()):
        items.append((name, tuple(value) if isinstance(value, (list, tuple)) else value))
    return (url, tuple(items), expect)


class ApiClient:
    def __init__(self, name, allowed_prefixes, min_interval=0.34, max_requests=80,
                 timeout=(5, 30), retries=2, session=None):
        self.name = name
        self.allowed_prefixes = tuple(allowed_prefixes)
        self.min_interval = min_interval
        self.max_requests = max_requests
        self.timeout = timeout
        self.retries = retries
        # A supplied session (tests, fixtures) gets its own limiter and cache,
        # so recorded responses never mix with real ones.
        shared = session is None
        self.session = session or requests.Session()
        self.limiter = shared_limiter(name, min_interval) if shared else RateLimiter(min_interval)
        self.cache = SHARED_CACHE if shared else ResponseCache()
        self.requests_sent = 0
        self.cache_hits = 0
        self.failures = []
        self._lock = Lock()

    def get(self, url, params=None, expect="json"):
        """GET a URL from an allowed service; return parsed JSON or text."""

        if not url.startswith(self.allowed_prefixes):
            raise ValueError(f"{self.name}: refusing an unexpected URL.")
        key = cache_key(url, params, expect)
        cached = self.cache.get(key)
        if cached is not None:
            self.cache_hits += 1
            return cached
        problem, wait = "", None
        for attempt in range(self.retries + 1):
            with self._lock:
                if self.requests_sent >= self.max_requests:
                    raise RequestLimitReached(
                        f"{self.name}: the limit of {self.max_requests} requests for this run was reached.")
                self.requests_sent += 1
            self.limiter.wait(None if wait is None else max(wait, self.min_interval))
            wait = None
            try:
                response = self.session.get(url, params=params, timeout=self.timeout)
            except requests.RequestException as error:
                problem = f"{self.name} could not be reached ({type(error).__name__})."
            else:
                if response.status_code == 200:
                    try:
                        value = response.json() if expect == "json" else response.text
                    except ValueError:
                        problem = f"{self.name} returned a response that could not be read."
                        break
                    self.cache.put(key, value)
                    return value
                problem = (f"{self.name} rate limit reached (HTTP 429)" if response.status_code == 429
                           else f"{self.name} returned HTTP {response.status_code}")
                problem += f" after {attempt + 1} attempt(s)."
                if response.status_code not in RETRY_STATUS:
                    break  # invalid requests are not retried
                wait = retry_after_seconds(response)
                if wait is not None and wait > MAX_RETRY_AFTER:
                    problem = (f"{self.name} asked to wait {wait:.0f} s before retrying (HTTP "
                               f"{response.status_code}); the search stopped instead.")
                    break
            if attempt < self.retries and wait is None:
                wait = min(2 ** attempt, 4)  # exponential backoff: 1 s, 2 s, 4 s
        self.failures.append(problem)
        raise ServiceUnavailable(problem)
