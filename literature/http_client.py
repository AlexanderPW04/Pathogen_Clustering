"""A small polite HTTP client for NCBI and PubTator (free public services).

It adds timeouts, a minimum delay between requests (rate limit), a cap on the
number of requests per run, bounded retries for temporary failures only, and
an in-memory response cache. It never follows URLs supplied by article text.
"""

import time
from threading import Lock

import requests

RETRY_STATUS = {429, 500, 502, 503, 504}


class RequestLimitReached(RuntimeError):
    """The per-run request cap was reached; no further requests are sent."""


class ServiceUnavailable(RuntimeError):
    """A service could not be reached or returned an error after retries."""


class ApiClient:
    def __init__(self, name, allowed_prefixes, min_interval=0.34, max_requests=80,
                 timeout=15, retries=2, session=None):
        self.name = name
        self.allowed_prefixes = tuple(allowed_prefixes)
        self.min_interval = min_interval
        self.max_requests = max_requests
        self.timeout = timeout
        self.retries = retries
        self.session = session or requests.Session()
        self.requests_sent = 0
        self.cache = {}
        self.failures = []
        self._last = 0.0
        self._lock = Lock()

    def get(self, url, params=None, expect="json"):
        """GET a URL from an allowed service; return parsed JSON or text."""

        if not url.startswith(self.allowed_prefixes):
            raise ValueError(f"{self.name}: refusing an unexpected URL.")
        key = (url, tuple(sorted((params or {}).items())), expect)
        if key in self.cache:
            return self.cache[key]
        for attempt in range(self.retries + 1):
            with self._lock:
                if self.requests_sent >= self.max_requests:
                    raise RequestLimitReached(
                        f"{self.name}: the limit of {self.max_requests} requests for this run was reached.")
                wait = self.min_interval - (time.monotonic() - self._last)
                if wait > 0:
                    time.sleep(wait)
                self.requests_sent += 1
                self._last = time.monotonic()
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
                    self.cache[key] = value
                    return value
                problem = f"{self.name} returned HTTP {response.status_code}."
                if response.status_code not in RETRY_STATUS:
                    break  # invalid requests are not retried
            if attempt < self.retries:
                time.sleep(min(2 ** attempt, 4))  # bounded backoff: 1 s, 2 s
        self.failures.append(problem)
        raise ServiceUnavailable(problem)
