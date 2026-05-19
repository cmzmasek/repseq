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


# ---------------------------------------------------------------------------
# Token notation: single HMM ("Name") or multidomain ("A--B--C").
# ---------------------------------------------------------------------------
#
# In a multidomain token the HMMs are written in C-to-N order:
#     "HMM1--HMM2"     means HMM1 lies C-terminal to HMM2 on the protein.
#     "HMM1--HMM2--HMM3" means HMM1 is most C-terminal, HMM3 most N-terminal.
# This is the opposite of the natural N-to-C reading direction most molecular
# biology uses; the project designer fixed C-to-N as the convention in v0.14.0
# and it is documented in default_config.yaml, README, and CLAUDE.md. Do NOT
# silently flip the direction — that would silently invert real configs.
#
# A CDS satisfies a token when:
#   - single token: that HMM has at least one passing hit on the CDS;
#   - multidomain: every HMM in the token has a passing hit AND those hits
#     appear in C-to-N order along the CDS (each named domain starts strictly
#     after the previous one ends — non-overlapping). Extra domains anywhere
#     on the CDS are fine ("HMMX--A--B" CDS satisfies "A--B" because A still
#     lies C-terminal to B).
#
# The token operations stay free of config-shape knowledge so they can be
# unit-tested in isolation.

TOKEN_SEPARATOR = "--"


def parse_hmm_token(token: str) -> list[str]:
    """Split a token into its ordered list of HMM names.

    Single HMM ``"Name"`` → ``["Name"]``. Multidomain ``"A--B--C"`` →
    ``["A", "B", "C"]`` in declared (C-to-N) order. Empty / whitespace-only
    components raise ``ValueError`` so misformatted tokens like ``"A----B"``
    or ``" --A"`` don't silently become single-HMM tokens.
    """
    if not isinstance(token, str):
        raise ValueError(f"HMM token must be a string, got {type(token).__name__}")
    raw = token.strip()
    if not raw:
        raise ValueError("HMM token cannot be empty")
    parts = [p.strip() for p in raw.split(TOKEN_SEPARATOR)]
    if any(not p for p in parts):
        raise ValueError(
            f"HMM token {token!r} has an empty component "
            f"(use 'A{TOKEN_SEPARATOR}B' with no surrounding whitespace)"
        )
    return parts


def cds_satisfies_token(
    hits: list[dict],
    token_hmms: list[str],
) -> Optional[float]:
    """Does this CDS's hit list satisfy the token? Returns the worst (largest)
    domain E-value across the satisfying hits on success, or ``None`` on
    failure.

    Hits are expected to carry a pre-computed ``passing`` flag (set by the
    pipeline's ``_run_hmm_qc`` step) — only passing hits are considered.
    The returned E-value is the worst-of-set so the caller can rank
    candidate CDSes "best E across the weakest required domain wins."

    Ordering rule for multidomain tokens: the named domains must appear in
    C-to-N order along the CDS (token's first HMM is C-terminal). For each
    consecutive pair (Hi, Hi+1) we require ``Hi.ali_from > Hi+1.ali_to`` —
    strict non-overlap, with Hi strictly C-terminal to Hi+1. If a domain has
    multiple passing hits, we pick the best (lowest E-value) hit that still
    satisfies the order constraint via a greedy left-to-right walk; the
    walk only fails if no consistent assignment exists. Extra hits to HMMs
    not named in the token are ignored.
    """
    if not token_hmms:
        return None
    # Index passing hits by target name.
    by_target: dict[str, list[dict]] = {}
    for h in hits or []:
        if not h.get("passing"):
            continue
        by_target.setdefault(h["target"], []).append(h)
    # Every named HMM must have at least one passing hit.
    for name in token_hmms:
        if name not in by_target:
            return None
    if len(token_hmms) == 1:
        # Single-HMM token: pick the best hit and return its E-value.
        best = min(by_target[token_hmms[0]], key=lambda h: h["dom_evalue"])
        return float(best["dom_evalue"])
    # Multidomain: assign one hit per named HMM s.t. the C-to-N order
    # holds (hit_i.ali_from > hit_{i+1}.ali_to for all consecutive i, i+1).
    # Greedy: walk the token N→C (reverse: most-N-terminal first), at each
    # step pick the *most-N-terminal* passing hit whose ali_to is strictly
    # less than the next required ali_from.
    reversed_tokens = list(reversed(token_hmms))
    chosen: list[dict] = []
    prev_to = -1  # ali_to of the previously-placed (more N-terminal) hit
    for name in reversed_tokens:
        # Candidate hits whose ali_from is strictly C-terminal to prev_to.
        candidates = [h for h in by_target[name] if h["ali_from"] > prev_to]
        if not candidates:
            return None
        # Pick the candidate with the smallest ali_to (leaves the most room
        # for subsequent more-C-terminal hits).
        pick = min(candidates, key=lambda h: h["ali_to"])
        chosen.append(pick)
        prev_to = pick["ali_to"]
    return float(max(h["dom_evalue"] for h in chosen))
