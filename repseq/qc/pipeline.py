"""QC pipeline: duplicate removal, length filter, ambiguous chars, annotation filter."""

from __future__ import annotations

import hashlib
import re
import statistics
from typing import Any

from ..models import QCReport, Sequence, SequenceType


# ---------------------------------------------------------------------------
# Step 1 – Exact duplicate removal
# ---------------------------------------------------------------------------

def remove_duplicates(sequences: list[Sequence], report: QCReport) -> list[Sequence]:
    """Remove exact-duplicate sequences, keeping the highest-quality copy.

    Among byte-identical sequences the survivor is chosen by
    RefSeq > reviewed-UniProt > first-seen, so a curated record is never
    discarded in favour of an arbitrary earlier duplicate. (Length is not a
    tiebreaker here — exact duplicates are the same length by definition.)
    """
    groups: dict[str, list[Sequence]] = {}
    order: list[str] = []
    for seq in sequences:
        h = hashlib.md5(seq.sequence.encode()).hexdigest()
        if h not in groups:
            groups[h] = []
            order.append(h)
        groups[h].append(seq)

    kept: list[Sequence] = []
    for h in order:
        members = groups[h]
        if len(members) == 1:
            kept.append(members[0])
            continue
        # max() is stable, so ties resolve to the first-seen member.
        best = max(members, key=lambda s: (s.is_refseq, s.is_reviewed))
        kept.append(best)
        for seq in members:
            if seq is best:
                continue
            seq.qc_passed = False
            seq.qc_fail_reason = f"exact_duplicate_of:{best.id}"
            report.removed_duplicates += 1
            report.add_removed(seq.id, seq.qc_fail_reason)
    return kept


# ---------------------------------------------------------------------------
# Step 2 – Length filter
# ---------------------------------------------------------------------------

def length_filter(sequences: list[Sequence], cfg: dict[str, Any], report: QCReport) -> list[Sequence]:
    """Filter sequences by length (median-percent or min/max)."""
    if not sequences:
        return sequences

    mode = cfg.get("mode", "median_percent")
    kept: list[Sequence] = []

    if mode == "median_percent":
        lengths = [s.length for s in sequences]
        median = statistics.median(lengths)
        # min_percent == 0 (or absent) disables the lower bound; max_percent
        # is an optional upper cap so grossly oversized records (mis-joined
        # genomes, contaminants) can also be dropped.
        min_pct = cfg.get("min_percent", 50)
        max_pct = cfg.get("max_percent")
        min_len = int(median * min_pct / 100.0) if min_pct else None
        max_len = int(median * max_pct / 100.0) if max_pct else None
    else:  # min_max
        min_len = cfg.get("min_length")
        max_len = cfg.get("max_length")

    for seq in sequences:
        fail = False
        reason = None
        if min_len is not None and seq.length < min_len:
            fail = True
            reason = f"length_too_short:{seq.length}<{min_len}"
        elif max_len is not None and seq.length > max_len:
            fail = True
            reason = f"length_too_long:{seq.length}>{max_len}"

        if fail:
            seq.qc_passed = False
            seq.qc_fail_reason = reason
            report.removed_length += 1
            report.add_removed(seq.id, reason)
        else:
            kept.append(seq)

    return kept


# ---------------------------------------------------------------------------
# Step 3 – Ambiguous character filter
# ---------------------------------------------------------------------------

def ambiguous_filter(
    sequences: list[Sequence], threshold: float, report: QCReport
) -> list[Sequence]:
    """Remove sequences with ambiguous character fraction above threshold."""
    kept: list[Sequence] = []
    for seq in sequences:
        frac = seq.ambiguous_fraction
        if frac > threshold:
            reason = f"ambiguous_fraction:{frac:.3f}>{threshold}"
            seq.qc_passed = False
            seq.qc_fail_reason = reason
            report.removed_ambiguous += 1
            report.add_removed(seq.id, reason)
        else:
            kept.append(seq)
    return kept


# ---------------------------------------------------------------------------
# Step 4 – Annotation keyword filter
# ---------------------------------------------------------------------------

def _build_keyword_pattern(keywords: list[str]) -> re.Pattern:
    """Build a case-insensitive pattern matching any keyword.

    A word boundary is anchored only on a side where the keyword's edge
    character is itself a word character. Blindly wrapping every keyword in
    ``\\b...\\b`` breaks keywords with a non-word edge: ``\\b`` after the
    ``:`` in ``MAG:`` demands a following word character, so ``MAG: Genus``
    (the actual NCBI title format) would never match.
    """
    def _is_word(ch: str) -> bool:
        return ch.isalnum() or ch == "_"

    parts: list[str] = []
    for kw in keywords:
        if not kw:
            continue
        esc = re.escape(kw)
        left = r"\b" if _is_word(kw[0]) else ""
        right = r"\b" if _is_word(kw[-1]) else ""
        parts.append(f"{left}{esc}{right}")
    if not parts:
        return re.compile(r"(?!)")  # matches nothing
    return re.compile("(" + "|".join(parts) + ")", re.IGNORECASE)


def annotation_filter(
    sequences: list[Sequence], cfg: dict[str, Any], report: QCReport
) -> list[Sequence]:
    """Remove sequences whose annotation/description matches any keyword.

    Matching targets the parsed ``description`` rather than the raw header,
    so structured header fields (e.g. a UniProt ``OS=`` organism name) can't
    trigger a false positive on words like "partial". Falls back to the full
    header only when no description was parsed.
    """
    if not cfg.get("enabled", True):
        return sequences

    keywords = cfg.get("keywords", [])
    if not keywords:
        return sequences

    pattern = _build_keyword_pattern(keywords)
    kept: list[Sequence] = []

    for seq in sequences:
        m = pattern.search(seq.description or seq.header)
        if m:
            reason = f"annotation_keyword:{m.group(1)}"
            seq.qc_passed = False
            seq.qc_fail_reason = reason
            report.removed_annotation += 1
            report.add_removed(seq.id, reason)
        else:
            kept.append(seq)

    return kept


# ---------------------------------------------------------------------------
# Public pipeline entry point
# ---------------------------------------------------------------------------

def run_qc(sequences: list[Sequence], cfg: dict[str, Any]) -> tuple[list[Sequence], QCReport]:
    """Run the full QC pipeline and return (passed_sequences, report)."""
    report = QCReport()
    report.total_input = len(sequences)

    qc_cfg = cfg.get("qc", {})

    if qc_cfg.get("remove_duplicates", True):
        sequences = remove_duplicates(sequences, report)

    # In segmented mode the input is a mixed pool of individual segments with
    # wildly different lengths (e.g. influenza PB2 ~2300 nt vs NS ~890 nt). A
    # whole-pool median — or any single absolute bound — is not meaningful and
    # would wrongly drop the short segments, leaving every isolate "incomplete"
    # at the completeness step. Per-segment bounds are applied later via
    # segmented.viruses.<v>.segment_lengths instead, so skip the global filter.
    if cfg.get("segmented", {}).get("enabled"):
        report.length_filter_skipped = True
    else:
        sequences = length_filter(sequences, qc_cfg.get("length_filter", {}), report)

    thresh = qc_cfg.get("ambiguous_threshold", 0.05)
    sequences = ambiguous_filter(sequences, thresh, report)

    sequences = annotation_filter(
        sequences, qc_cfg.get("annotation_filter", {}), report
    )

    report.passed = len(sequences)
    return sequences, report
