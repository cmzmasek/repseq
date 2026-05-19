"""Top-level cache-aware batch entry point for hmmscan.

``scan_proteins`` is the single function the pipeline calls. It:
    1. Resolves and indexes the database (with auto-hmmpress).
    2. Computes a db-signature for cache invalidation.
    3. Dedupes queries by AA sequence (many CDSes share identical AA
       across closely related isolates).
    4. Splits cached vs uncached, batches the uncached set into ONE
       hmmscan call, writes results back to the cache.
    5. Returns ``{query_id: [hit_dict, ...]}``.

``passes_cutoffs`` applies the E-value/GA + relative-length gate so
callers don't need to know the cutoff layout themselves.
"""
from __future__ import annotations

import hashlib
from typing import Any, Optional

from .database import (
    db_signature,
    ensure_pressed,
    parse_ga_cutoffs,
    resolve_database_path,
)
from .errors import HMMScanError
from .hmmscan import is_available, scan

CACHE_SOURCE = "hmmscan"


def _cache_key(protein_seq: str, db_sig: str) -> str:
    h = hashlib.sha256(protein_seq.encode("utf-8")).hexdigest()
    return f"{h}:{db_sig}"


def scan_proteins(
    proteins: dict[str, str],
    cfg: dict[str, Any],
    cache: Optional[Any] = None,
) -> dict[str, list[dict]]:
    """Run hmmscan over a batch of protein queries with persistent caching.

    Args:
        proteins: ``{query_id: protein_aa_sequence}``. ``query_id`` must
            be unique and free of whitespace; caller's responsibility.
        cfg: full repseq config dict.
        cache: optional ``TaxonomyCache``. When provided, cached hits
            are reused and new hits written back. When ``None``, no
            caching (every call hits hmmscan).

    Returns:
        ``{query_id: [hit_dict, ...]}``. Queries with no hits are absent
        from the dict.

    Raises:
        HMMDatabaseError: database resolution or indexing failed.
        HMMScanError: hmmscan invocation failed.
    """
    if not proteins:
        return {}
    if not is_available():
        raise HMMScanError("hmmscan is not on PATH")

    hcfg = cfg.get("hmm", {}) or {}
    db_path = resolve_database_path(hcfg.get("database"))
    ensure_pressed(db_path)
    sig = db_signature(db_path)

    results: dict[str, list[dict]] = {}
    # Dedup by sequence: identical AA from many CDSes only scanned once.
    seq_to_qids: dict[str, list[str]] = {}
    for qid, seq in proteins.items():
        if not seq:
            continue
        if cache is not None:
            entry = cache.get(CACHE_SOURCE, _cache_key(seq, sig))
            if entry is not None:
                results[qid] = entry.get("hits", [])
                continue
        seq_to_qids.setdefault(seq, []).append(qid)

    if seq_to_qids:
        # Deterministic per-unique-sequence scan id so we can round-trip
        # to the caller's qids after hmmscan returns.
        ordered_seqs = sorted(seq_to_qids.keys())
        scan_queries = {f"S{idx:07d}": seq for idx, seq in enumerate(ordered_seqs)}
        sid_to_seq = {sid: seq for sid, seq in scan_queries.items()}

        threads = hcfg.get("threads")
        if threads is None:
            threads = cfg.get("threads", 1)
        raw_hits = scan(db_path, scan_queries, threads=int(threads or 1))

        for sid, seq in sid_to_seq.items():
            hits = raw_hits.get(sid, [])
            if cache is not None:
                cache.set(CACHE_SOURCE, _cache_key(seq, sig), {"hits": hits})
            for qid in seq_to_qids[seq]:
                results[qid] = hits

    return results


def get_ga_cutoffs(cfg: dict[str, Any]) -> dict[str, Optional[float]]:
    """Parse GA cutoffs from the configured database (call once per run)."""
    db_path = resolve_database_path(cfg.get("hmm", {}).get("database"))
    return parse_ga_cutoffs(db_path)


def passes_cutoffs(
    hit: dict,
    ga_cutoffs: dict[str, Optional[float]],
    default_evalue: float,
    relative_length_cutoff: float,
    use_ga_when_available: bool,
) -> bool:
    """Apply E-value/GA AND coverage gates to one hit.

    Similarity gate:
        - If ``use_ga_when_available`` and the target HMM has a curated
          GA bit-score, require ``dom_score >= GA``.
        - Otherwise require ``dom_evalue <= default_evalue``.
    Coverage gate:
        - ``ali_span / hmm_len >= relative_length_cutoff``. The HMM
          model length is the denominator (per design): short HMM on
          long CDS is a valid domain hit; long HMM with short alignment
          is the reject case.
    """
    ga = ga_cutoffs.get(hit["target"]) if use_ga_when_available else None
    if ga is not None:
        similarity_pass = hit["dom_score"] >= ga
    else:
        similarity_pass = hit["dom_evalue"] <= default_evalue
    hmm_len = max(int(hit.get("hmm_len", 0)), 1)
    coverage = hit["ali_span"] / hmm_len
    length_pass = coverage >= relative_length_cutoff
    return similarity_pass and length_pass


def coverage_of(hit: dict) -> float:
    """Convenience: ali_span / hmm_len, clipped to [0, 1]."""
    hmm_len = max(int(hit.get("hmm_len", 0)), 1)
    return min(1.0, hit["ali_span"] / hmm_len)
