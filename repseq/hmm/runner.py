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
    db_path = resolve_database_path(
        hcfg.get("database"), cache_dir=_cache_dir(cfg)
    )
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


def _cache_dir(cfg: dict[str, Any]) -> Optional[Path]:
    """Cache directory for combined directory-databases (None → temp)."""
    cd = cfg.get("cache_dir")
    return Path(cd).expanduser() if cd else None


def get_ga_cutoffs(cfg: dict[str, Any]) -> dict[str, Optional[float]]:
    """Parse GA cutoffs from the configured database (call once per run)."""
    db_path = resolve_database_path(
        cfg.get("hmm", {}).get("database"), cache_dir=_cache_dir(cfg)
    )
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
        - ``hmm_span / hmm_len >= relative_length_cutoff`` where
          ``hmm_span = hmm_to - hmm_from + 1`` is the span on the HMM
          MODEL (domtblout hmm-from/hmm-to), not the alignment span on
          the protein (ali-from/ali-to). Coverage asks "how much of the
          profile did this hit match", so both numerator and denominator
          live in HMM-model coordinates. Using the protein span as the
          numerator diverges whenever the CDS carries insertions or
          deletions relative to the model (insertions inflate it, indels
          deflate it). Per design: a short HMM matching part of a long
          CDS is a valid hit; a long HMM that aligns only a short stretch
          of its model is the reject case.
    """
    ga = ga_cutoffs.get(hit["target"]) if use_ga_when_available else None
    if ga is not None:
        similarity_pass = hit["dom_score"] >= ga
    else:
        similarity_pass = hit["dom_evalue"] <= default_evalue
    length_pass = coverage_of(hit) >= relative_length_cutoff
    return similarity_pass and length_pass


def coverage_of(hit: dict) -> float:
    """HMM-model coverage: ``(hmm_to - hmm_from + 1) / hmm_len``, clipped to [0, 1].

    The numerator is the span on the HMM MODEL, not the alignment span on
    the protein (``ali_span``, which stays in protein coordinates and is
    used for domain-ordering checks). The model span is the correct
    measure of "how much of the profile this hit covered".
    """
    hmm_len = max(int(hit.get("hmm_len", 0)), 1)
    hmm_span = max(0, int(hit.get("hmm_to", 0)) - int(hit.get("hmm_from", 0)) + 1)
    return min(1.0, hmm_span / hmm_len)


# ---------------------------------------------------------------------------
# Token notation: single HMM ("Name") or multidomain ("A--B--C").
# ---------------------------------------------------------------------------
#
# In a multidomain token the HMMs are written in N-to-C order — the natural
# molecular-biology reading direction (left-to-right = N-terminus to
# C-terminus):
#     "HMM1--HMM2"     means HMM1 lies N-terminal to HMM2 on the protein.
#     "HMM1--HMM2--HMM3" means HMM1 is most N-terminal, HMM3 most C-terminal.
#
# A CDS satisfies a token when:
#   - single token: that HMM has at least one passing hit on the CDS;
#   - multidomain: every HMM in the token has a passing hit AND those hits
#     appear in N-to-C order along the CDS (each named domain starts strictly
#     after the previous one ends — non-overlapping). Extra domains anywhere
#     on the CDS are fine ("A--B--HMMX" CDS satisfies "A--B" because A still
#     lies N-terminal to B).
#
# The token operations stay free of config-shape knowledge so they can be
# unit-tested in isolation.

TOKEN_SEPARATOR = "--"


def parse_hmm_token(token: str) -> list[str]:
    """Split a token into its ordered list of HMM names.

    Single HMM ``"Name"`` → ``["Name"]``. Multidomain ``"A--B--C"`` →
    ``["A", "B", "C"]`` in declared (N-to-C) order. Empty / whitespace-only
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
    *,
    overlap_tolerance: int = 0,
) -> Optional[float]:
    """Does this CDS's hit list satisfy the token? Returns the worst (largest)
    domain E-value across the satisfying hits on success, or ``None`` on
    failure.

    Hits are expected to carry a pre-computed ``passing`` flag (set by the
    pipeline's ``_run_hmm_qc`` step) — only passing hits are considered.
    The returned E-value is the worst-of-set so the caller can rank
    candidate CDSes "best E across the weakest required domain wins."

    Ordering rule for multidomain tokens — named domains must appear in
    N-to-C order along the CDS (token's first HMM is N-terminal). For each
    consecutive pair (Hi, Hi+1) we require:

    * ``Hi.ali_to <= Hi+1.ali_from + overlap_tolerance`` — boundaries may
      overlap by up to ``overlap_tolerance`` residues (default 0 → strict
      non-overlap, the pre-v0.22 behaviour). The tolerance accommodates
      Pfam profile boundaries that don't quite align to the molecular
      cleavage site (e.g. coronavirus S1/S2 across the furin site).
    * ``Hi+1.ali_from > Hi.ali_from`` AND ``Hi+1.ali_to > Hi.ali_to`` —
      strict forward progression of both endpoints. This rejects full
      containment and prevents a single fused region from masquerading
      as two distinct domains by matching the two profiles at the same
      span, regardless of how generous ``overlap_tolerance`` is.

    If a domain has multiple passing hits, we pick the *most-N-terminal*
    hit that still satisfies the order constraint via a greedy left-to-
    right walk; the walk only fails if no consistent assignment exists.
    Extra hits to HMMs not named in the token are ignored.
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
        # overlap_tolerance is irrelevant — there are no neighbours to overlap.
        best = min(by_target[token_hmms[0]], key=lambda h: h["dom_evalue"])
        return float(best["dom_evalue"])
    # Multidomain: walk the token N→C, at each step picking the most-N-
    # terminal hit consistent with the prior pick under the tolerance rule.
    chosen: list[dict] = []
    prev_from = -1  # ali_from of the previously-placed (more N-terminal) hit
    prev_to = -1    # ali_to of the same
    for name in token_hmms:
        # Candidates that progress forward on both endpoints AND don't
        # overlap the prior hit by more than the tolerance.
        candidates = [
            h for h in by_target[name]
            if h["ali_from"] > prev_from
            and h["ali_to"] > prev_to
            and h["ali_to"] - h["ali_from"] >= 0  # well-formed; defensive
            and (prev_to < 0 or h["ali_from"] >= prev_to - overlap_tolerance + 1)
        ]
        if not candidates:
            return None
        # Pick the candidate with the smallest ali_to (leaves the most room
        # for subsequent more-C-terminal hits).
        pick = min(candidates, key=lambda h: h["ali_to"])
        chosen.append(pick)
        prev_from = pick["ali_from"]
        prev_to = pick["ali_to"]
    return float(max(h["dom_evalue"] for h in chosen))
