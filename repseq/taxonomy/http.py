"""Shared HTTP plumbing for the NCBI / UniProt clients: retry + failure log.

Two jobs, both born from one bug (2026-11-12, Orthoherpesviridae pass 1):
a transient Entrez error (rate limit, timeout, 5xx) on a metadata lookup was
swallowed by ``except Exception: return None`` — one attempt, no retry, no
count — so 9 of 5,425 records silently lost their taxonomy, formed an
``Unknown`` group and contributed 6 duplicate representatives.

* :func:`get_with_retry` — every client request goes through it, so a
  transient error gets :data:`MAX_TRIES` attempts with a growing wait.
* :func:`record_failure` / :func:`failures` — a process-wide log of the
  lookups that STILL failed after the retries. The clients keep returning
  ``None`` (their callers' contract is unchanged) but the failure is now on
  record, so the CLI can warn, write ``{prefix}_lookup_failures.tsv`` and
  say so in the Methods prose. Process-wide rather than per-client because
  the failures arise at several pipeline steps (metadata resolution, GenBank
  record fetches for protein QC / HMM QC / segment metadata) and the report
  is written by code that never sees the client object.

A failed lookup is never cached, so a re-run retries exactly those records.
"""
from __future__ import annotations

import logging
import threading
import time
from typing import Any, Callable, Optional

import requests

logger = logging.getLogger(__name__)

#: Attempts per request, the first included.
MAX_TRIES = 3
#: HTTP statuses worth another attempt: rate limit + server-side trouble.
RETRY_STATUS = frozenset({429, 500, 502, 503, 504})
#: Longest single wait, also the cap on a server's ``Retry-After``.
_MAX_WAIT = 60.0

# Indirection so tests can run the retry path without really waiting.
_sleep = time.sleep

#: Stage names used in the failure log (and the TSV's ``stage`` column).
STAGE_METADATA = "metadata"        # organism / taxid / host ... by accession
STAGE_LINEAGE = "taxonomy"         # a record that lost its ranked lineage
# Internal: a lineage fetch that failed, keyed by TAXID. Many records share
# a taxid and another thread may fetch it successfully a second later, which
# clears this entry — so it cannot be what the report lists. The resolver
# turns it into a per-accession STAGE_LINEAGE entry at the moment a record
# is actually resolved without its lineage. Stages starting with "_" are
# left out of :func:`failures`.
STAGE_LINEAGE_TAXID = "_taxid"
STAGE_RECORD = "genbank_record"    # GenBank record: CDS proteins, source
STAGE_SEQUENCE = "sequence"        # nucleotide body (replay)

_lock = threading.Lock()
_failures: dict[tuple[str, str], str] = {}

# Give-up state for "there is no network at all". Retrying is right for a
# hiccup and ruinous for an outage: 2 s + 4 s of waiting on EVERY one of
# 120,000 per-record lookups turns a run that used to fall through in
# seconds into two days of a frozen progress bar. So once this many
# requests in a row have failed to even reach the service, stop retrying
# (one attempt each) until one gets through again.
_OFFLINE_AFTER = 5
_consecutive_unreachable = 0
_offline_announced = False


def _is_unreachable(exc: BaseException) -> bool:
    """Could not reach the service at all (as opposed to an HTTP answer)."""
    return isinstance(exc, (requests.ConnectionError, requests.Timeout))


def _note_outcome(exc: Optional[BaseException]) -> None:
    global _consecutive_unreachable, _offline_announced
    with _lock:
        if exc is not None and _is_unreachable(exc):
            _consecutive_unreachable += 1
            announce = (
                _consecutive_unreachable >= _OFFLINE_AFTER
                and not _offline_announced
            )
            if announce:
                _offline_announced = True
        else:
            # Any answer from the service, even an error status, means the
            # network is there.
            _consecutive_unreachable = 0
            _offline_announced = False
            announce = False
    if announce:
        logger.warning(
            "%d database requests in a row could not reach the service — "
            "the network or the service looks unavailable. Retries are "
            "suspended until a request gets through; lookups will fail "
            "fast and be listed as failures. If you are working offline on "
            "purpose, stop and re-run with --no-resolve.",
            _OFFLINE_AFTER,
        )


def _retries_suspended() -> bool:
    with _lock:
        return _consecutive_unreachable >= _OFFLINE_AFTER


def _is_transient(exc: BaseException) -> bool:
    if isinstance(exc, requests.HTTPError):
        resp = exc.response
        return resp is not None and resp.status_code in RETRY_STATUS
    return isinstance(
        exc,
        (
            requests.ConnectionError,
            requests.Timeout,
            requests.exceptions.ChunkedEncodingError,
        ),
    )


def _wait_seconds(exc: BaseException, attempt: int) -> float:
    """Server's ``Retry-After`` when it gives one, else 2 s, 4 s, ..."""
    resp = getattr(exc, "response", None)
    if resp is not None:
        try:
            retry_after = float(resp.headers.get("Retry-After", ""))
            if retry_after > 0:
                return min(retry_after, _MAX_WAIT)
        except (TypeError, ValueError, AttributeError):
            pass
    return min(2.0 ** attempt, _MAX_WAIT)


def get_with_retry(
    url: str,
    *,
    params: Optional[dict] = None,
    timeout: float = 30,
    throttle: Optional[Callable[[], None]] = None,
) -> Any:
    """``requests.get`` + ``raise_for_status`` with retries on transient errors.

    ``throttle`` (the client's rate limiter) runs before EVERY attempt, so a
    retry still respects the service's request spacing. A non-transient
    error (404, 400, a parse problem upstream) is raised at once; a
    transient one is raised only after :data:`MAX_TRIES` attempts.
    """
    for attempt in range(1, MAX_TRIES + 1):
        if throttle is not None:
            throttle()
        try:
            resp = requests.get(url, params=params, timeout=timeout)
            resp.raise_for_status()
            _note_outcome(None)
            return resp
        except Exception as exc:
            _note_outcome(exc)
            if (
                attempt >= MAX_TRIES
                or not _is_transient(exc)
                or (_is_unreachable(exc) and _retries_suspended())
            ):
                raise
            wait = _wait_seconds(exc, attempt)
            logger.info(
                "Request to %s failed (%s); attempt %d of %d, retrying in "
                "%.0f s.", url, describe_error(exc), attempt, MAX_TRIES, wait,
            )
            _sleep(wait)
    raise AssertionError("unreachable")  # pragma: no cover


def describe_error(exc: BaseException) -> str:
    """A short plain-English reason for the failure log."""
    if isinstance(exc, requests.HTTPError) and exc.response is not None:
        code = exc.response.status_code
        what = {
            429: "rate limit reached",
            500: "server error",
            502: "server error",
            503: "service unavailable",
            504: "server timed out",
        }.get(code, "request refused")
        return f"HTTP {code} ({what})"
    if isinstance(exc, requests.Timeout):
        return "no answer within the time limit"
    if isinstance(exc, requests.ConnectionError):
        return "could not connect"
    return f"{type(exc).__name__}: {exc}"


def record_failure(stage: str, key: Any, exc_or_msg: Any) -> str:
    """Log that the lookup of ``key`` at ``stage`` failed; return the reason."""
    msg = (
        exc_or_msg if isinstance(exc_or_msg, str) else describe_error(exc_or_msg)
    )
    with _lock:
        _failures[(stage, str(key))] = msg
    return msg


def clear_failure(stage: str, key: Any) -> None:
    """A later attempt for ``key`` succeeded — it is no longer a failure."""
    with _lock:
        _failures.pop((stage, str(key)), None)


def failure_for(stage: str, key: Any) -> Optional[str]:
    with _lock:
        return _failures.get((stage, str(key)))


def failures() -> list[dict[str, str]]:
    """Every reportable failure as ``{stage, id, error}``, in a stable order
    (internal ``_``-prefixed stages are left out)."""
    with _lock:
        items = sorted(_failures.items())
    return [
        {"stage": stage, "id": key, "error": msg}
        for (stage, key), msg in items
        if not stage.startswith("_")
    ]


def is_network_failure(exc: BaseException) -> bool:
    """Is this a failure of the LOOKUP rather than an answer about the record?

    A 4xx other than 429 is the service answering "no" (400 = an id it does
    not recognise, 404 = no such record) — that is a fact about the record,
    not something a re-run will change, so it is not reported as a failure.
    Everything else (rate limit, 5xx, timeout, no connection, or an
    unexpected error while reading the answer) is.
    """
    if isinstance(exc, requests.HTTPError) and exc.response is not None:
        code = exc.response.status_code
        return code == 429 or code >= 500
    return True


def reset_failures() -> None:
    """Start a fresh log (called once at the start of metadata resolution)."""
    global _consecutive_unreachable, _offline_announced
    with _lock:
        _failures.clear()
        _consecutive_unreachable = 0
        _offline_announced = False
