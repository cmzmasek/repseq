"""Per-isolate taxonomy-consistency QC for segmented viruses.

Drops any isolate whose segments disagree on the taxonomy rank named
in ``segmented.taxonomy_consistency.rank`` (default ``species``).

The default is *on*: reassortment between distinct parent species is
real biology for many segmented viruses (peribunyaviruses,
orthomyxoviruses, …), so isolates with mixed-species segments would
otherwise pass through unflagged and contaminate the representative
set. A user who wants to keep them can set
``segmented.taxonomy_consistency.enabled: false``.
"""

from __future__ import annotations

import re
from collections import defaultdict
from typing import Optional

from ..models import Sequence


# Collapse internal whitespace runs to a single space before comparing
# species labels — protects against trailing/leading space and oddities
# like "Bunyamwera  virus" (double space) from inconsistent NCBI
# annotation.
_WS_RUN = re.compile(r"\s+")


def _norm_label(value: Optional[str]) -> Optional[str]:
    """Lower-case + whitespace-collapse a taxonomy label for comparison.

    Returns ``None`` for empty / unset values so the caller can treat
    "missing" distinctly from "different". Anything else round-trips
    through ``str.lower()`` and ``re.sub`` so "Bunyamwera virus" and
    "bunyamwera  virus " compare equal.
    """
    if value is None:
        return None
    s = _WS_RUN.sub(" ", str(value).strip()).lower()
    return s or None


def _rank_value(seq: Sequence, rank: str) -> Optional[str]:
    """Look up ``rank`` on a Sequence's resolved taxonomy, normalised."""
    if seq.taxonomy is None:
        return None
    raw = seq.taxonomy.get_rank(rank)
    return _norm_label(raw)


def filter_taxonomy_consistent_isolates(
    sequences: list[Sequence],
    *,
    rank: str = "species",
    policy=None,
    protected_out: Optional[list[tuple[str, str, str]]] = None,
) -> tuple[list[Sequence], list[tuple[str, str]]]:
    """Drop isolates whose segments disagree on the given taxonomy rank.

    Groups input sequences by ``seq.isolate_id``; for each group,
    compares the *populated* values at the requested rank. Missing
    values (resolver couldn't fill the rank, or the rank simply isn't
    in the lineage) are ignored — an isolate is only dropped when two
    or more segments carry distinct populated labels. Single-segment
    groups and groups with no populated labels are always kept (no
    mismatch detectable).

    Sequences without an isolate_id (UniProt input, missing accession,
    or runs where ``_populate_genbank_isolate_segment`` was skipped)
    bypass this filter entirely — the regex fallback in
    ``filter_complete_isolates`` will catch them later, and judging
    "consistency" without a known isolate grouping isn't meaningful.

    Returns:
        ``(kept_sequences, removed)`` where ``removed`` is a list of
        ``(accession_or_id, reason)`` pairs the caller can feed into
        ``QCReport.add_removed``. Order of ``kept_sequences`` matches
        the input order (this matters: downstream steps assume input
        order is stable for deterministic behaviour).
    """
    groups: dict[str, list[Sequence]] = defaultdict(list)
    floaters: list[Sequence] = []
    for seq in sequences:
        if seq.isolate_id:
            groups[seq.isolate_id].append(seq)
        else:
            floaters.append(seq)

    drop_ids: set[int] = set()
    removed: list[tuple[str, str]] = []
    reason = f"taxonomy_mismatch:{rank}"

    for iso_id, segs in groups.items():
        if len(segs) < 2:
            continue
        labels = {lbl for lbl in (_rank_value(s, rank) for s in segs) if lbl}
        if len(labels) <= 1:
            continue
        # Force-keep whitelist: naming any one segment protects the whole
        # isolate from being dropped. Record each segment as protected (one
        # row per segment, mirroring the per-segment removal granularity).
        if policy is not None and policy.protects_any(segs, "taxonomy_consistency"):
            if protected_out is not None:
                for seq in segs:
                    protected_out.append(
                        (seq.id, "taxonomy_consistency", reason)
                    )
            continue
        # Mismatch — drop every segment of this isolate. We record the
        # accession (matches the column header on _qc_removed.tsv) and
        # fall back to seq.id if no accession is set.
        for seq in segs:
            drop_ids.add(id(seq))
            removed.append((seq.accession or seq.id, reason))

    kept = [seq for seq in sequences if id(seq) not in drop_ids]
    return kept, removed
