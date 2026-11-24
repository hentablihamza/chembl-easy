"""A small, robust HTTP client for the ChEMBL REST API.

Solves the problems students hit first when they script against ChEMBL:

* **Timeouts and transient failures** (HTTP 429/502/503/504, dropped
  connections) are retried with exponential backoff, honouring
  ``Retry-After`` when the server sends it.
* **Rate limiting** - a minimum interval between requests keeps scripts
  polite and avoids being throttled.
* **Re-downloading the same data** - every successful JSON response can be
  cached on disk, so re-running a notebook is instant and works offline.
* **Pagination** - :meth:`ChEMBLClient.paginate` walks ``limit``/``offset``
  pages (1000 records per page, the API maximum) and yields records.
* **Cryptic server errors** - a 500 caused by an unparsable query molecule is
  turned into a clear :class:`InvalidQueryError` instead of being retried.
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
import time
from pathlib import Path
from typing import Any, Callable, Iterator, Mapping

import requests

from chembl_easy.version import __version__

logger = logging.getLogger(__name__)

BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"
MAX_PAGE_SIZE = 1000
RETRY_STATUS = frozenset({429, 500, 502, 503, 504})
DEFAULT_CACHE_DIR = Path(
    os.environ.get("CHEMBL_EASY_CACHE", Path.home() / ".cache" / "chembl-easy")
)


class ChEMBLError(RuntimeError):
    """Raised when the ChEMBL API cannot return a usable response."""


class NotFoundError(ChEMBLError):
    """Raised for HTTP 404 (unknown ChEMBL ID or malformed URL)."""


class InvalidQueryError(ChEMBLError):
    """Raised when ChEMBL rejects the query itself (e.g. an unparsable SMILES)."""


def _cache_key(url: str, params: Mapping[str, Any] | None) -> str:
    items = sorted((str(k), str(v)) for k, v in (params or {}).items())
    raw = url + "?" + json.dumps(items)
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


class ChEMBLClient:
    """Retrying, caching, rate-limited client for ``/chembl/api/data``.

    Parameters
    ----------
    base_url:
        API root. Point it at a local mirror if you run one.
    cache_dir:
        Directory for the on-disk JSON cache. ``None`` disables caching.
    cache_ttl:
        Maximum age of a cache entry in seconds. ``None`` means entries never
        expire (fine for a fixed ChEMBL release).
    timeout:
        Per-request timeout in seconds. Activity pages of 1000 records can
        take 5-20 s, so the default is generous.
    max_retries:
        Number of retries after the first attempt for retryable failures.
    backoff:
        Base of the exponential backoff (``backoff * 2**attempt`` seconds).
    min_interval:
        Minimum seconds between two network requests (simple rate limiter).
    session:
        Any object with a ``requests``-like ``get`` method. Injected in tests.
    sleep:
        Sleep function, injectable so tests do not actually wait.
    """

    def __init__(
        self,
        base_url: str = BASE_URL,
        cache_dir: str | os.PathLike[str] | None = DEFAULT_CACHE_DIR,
        cache_ttl: float | None = None,
        timeout: float = 90.0,
        max_retries: int = 5,
        backoff: float = 1.5,
        min_interval: float = 0.1,
        session: Any | None = None,
        sleep: Callable[[float], None] = time.sleep,
    ) -> None:
        self.base_url = base_url.rstrip("/")
        self.cache_dir = Path(cache_dir) if cache_dir is not None else None
        self.cache_ttl = cache_ttl
        self.timeout = timeout
        self.max_retries = max_retries
        self.backoff = backoff
        self.min_interval = min_interval
        self._sleep = sleep
        self._last_request = 0.0
        self.n_requests = 0
        self.n_cache_hits = 0
        if session is None:
            session = requests.Session()
            session.headers.update(
                {
                    "Accept": "application/json",
                    "User-Agent": f"chembl-easy/{__version__} (+https://github.com/hentablihamza/chembl-easy)",
                }
            )
        self.session = session

    # ------------------------------------------------------------------ cache
    def _cache_path(self, key: str) -> Path | None:
        if self.cache_dir is None:
            return None
        return self.cache_dir / key[:2] / f"{key}.json"

    def _read_cache(self, key: str) -> dict | None:
        path = self._cache_path(key)
        if path is None or not path.exists():
            return None
        if self.cache_ttl is not None and time.time() - path.stat().st_mtime > self.cache_ttl:
            return None
        try:
            return json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            return None

    def _write_cache(self, key: str, data: dict) -> None:
        path = self._cache_path(key)
        if path is None:
            return
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_suffix(".tmp")
        tmp.write_text(json.dumps(data), encoding="utf-8")
        tmp.replace(path)

    def clear_cache(self) -> int:
        """Delete every cached response. Returns the number of files removed."""
        if self.cache_dir is None or not self.cache_dir.exists():
            return 0
        n = 0
        for path in self.cache_dir.rglob("*.json"):
            path.unlink()
            n += 1
        return n

    # ---------------------------------------------------------------- network
    def _throttle(self) -> None:
        wait = self.min_interval - (time.monotonic() - self._last_request)
        if wait > 0:
            self._sleep(wait)
        self._last_request = time.monotonic()

    def _retry_delay(self, attempt: int, response: Any | None) -> float:
        if response is not None:
            retry_after = getattr(response, "headers", {}).get("Retry-After")
            if retry_after:
                try:
                    return max(float(retry_after), 0.0)
                except ValueError:
                    pass
        return self.backoff * (2**attempt)

    def url(self, endpoint: str) -> str:
        endpoint = endpoint.strip("/")
        if not endpoint.endswith(".json"):
            endpoint += ".json"
        return f"{self.base_url}/{endpoint}"

    def get(self, endpoint: str, params: Mapping[str, Any] | None = None, use_cache: bool = True) -> dict:
        """GET ``endpoint`` (e.g. ``"activity"`` or ``"molecule/CHEMBL25"``) as JSON.

        Query values are always sent as URL parameters, never pasted into the
        path, so SMILES containing ``#``, ``/`` or ``\\`` are encoded safely.
        """
        url = self.url(endpoint)
        params = dict(params or {})
        key = _cache_key(url, params)
        if use_cache:
            cached = self._read_cache(key)
            if cached is not None:
                self.n_cache_hits += 1
                return cached

        last_error: str = ""
        for attempt in range(self.max_retries + 1):
            self._throttle()
            response = None
            try:
                self.n_requests += 1
                response = self.session.get(url, params=params, timeout=self.timeout)
            except (requests.ConnectionError, requests.Timeout) as exc:
                last_error = f"{type(exc).__name__}: {exc}"
            else:
                status = response.status_code
                if status == 200:
                    try:
                        data = response.json()
                    except ValueError as exc:
                        last_error = f"invalid JSON: {exc}"
                    else:
                        if use_cache:
                            self._write_cache(key, data)
                        return data
                elif status == 404:
                    raise NotFoundError(f"404 Not Found: {response.url if hasattr(response, 'url') else url}")
                elif status == 500 and "standardiz" in (response.text or "").lower():
                    raise InvalidQueryError(
                        "ChEMBL could not parse/standardize the query molecule. "
                        "Check the SMILES (chembl_easy.validate_smiles helps)."
                    )
                elif status in RETRY_STATUS:
                    last_error = f"HTTP {status}"
                else:
                    raise ChEMBLError(f"HTTP {status} for {url}: {(response.text or '')[:300]}")

            if attempt < self.max_retries:
                delay = self._retry_delay(attempt, response)
                logger.warning(
                    "ChEMBL request failed (%s); retry %d/%d in %.1f s",
                    last_error, attempt + 1, self.max_retries, delay,
                )
                self._sleep(delay)
        raise ChEMBLError(f"Giving up on {url} after {self.max_retries + 1} attempts ({last_error})")

    def count(self, endpoint: str, params: Mapping[str, Any] | None = None) -> int:
        """Total number of records matching ``params`` (one tiny request)."""
        params = dict(params or {})
        params.update({"limit": 1, "offset": 0})
        data = self.get(endpoint, params)
        return int(data["page_meta"]["total_count"])

    def paginate(
        self,
        endpoint: str,
        params: Mapping[str, Any] | None = None,
        key: str | None = None,
        page_size: int = MAX_PAGE_SIZE,
        max_records: int | None = None,
        progress: Callable[[int, int], None] | None = None,
    ) -> Iterator[dict]:
        """Yield every record of a list endpoint, page by page.

        ``key`` is the name of the list in the JSON payload (``"activities"``,
        ``"molecules"``, ``"targets"`` ...); it defaults to the first list
        found next to ``page_meta``.
        """
        page_size = max(1, min(int(page_size), MAX_PAGE_SIZE))
        if max_records is not None:
            page_size = min(page_size, max(1, max_records))
        base = dict(params or {})
        offset = 0
        yielded = 0
        while True:
            page = self.get(endpoint, {**base, "limit": page_size, "offset": offset})
            meta = page.get("page_meta", {})
            if key is None:
                key = next(k for k, v in page.items() if k != "page_meta" and isinstance(v, list))
            records = page.get(key, [])
            total = int(meta.get("total_count") or 0)
            for rec in records:
                yield rec
                yielded += 1
                if max_records is not None and yielded >= max_records:
                    return
            offset += len(records)
            if progress is not None:
                progress(min(offset, total), total)
            if not records or not meta.get("next") or offset >= total:
                return

    def status(self) -> dict:
        """ChEMBL release and database counts (never cached)."""
        return self.get("status", use_cache=False)
