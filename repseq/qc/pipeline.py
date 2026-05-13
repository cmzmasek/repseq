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
    """Remove exact duplicate sequences, keeping the first occurrence."""
    seen: dict[str, str] = {}  # hash -> first seq id
    kept: list[Sequence] = []
    for seq in sequences:
        h = hashlib.md5(seq.sequence.encode()).hexdigest()
        if h in seen:
            seq.qc_passed = False
            seq.qc_fail_reason = f"exact_duplicate_of:{seen[h]}"
            report.removed_duplicates += 1
            report.add_removed(seq.id, seq.qc_fail_reason)
        else:
            seen[h] = seq.id
            kept.append(seq)
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
        pct = cfg.get("min_percent", 50) / 100.0
        lengths = [s.length for s in sequences]
        median = statistics.median(lengths)
        min_len = int(median * pct)
        max_len = None
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
    """Build a case-insensitive pattern matching any keyword as a whole word."""
    escaped = [re.escape(kw) for kw in keywords]
    return re.compile(r"\b(" + "|".join(escaped) + r")\b", re.IGNORECASE)


def annotation_filter(
    sequences: list[Sequence], cfg: dict[str, Any], report: QCReport
) -> list[Sequence]:
    """Remove sequences whose header matches any annotation keyword."""
    if not cfg.get("enabled", True):
        return sequences

    keywords = cfg.get("keywords", [])
    if not keywords:
        return sequences

    pattern = _build_keyword_pattern(keywords)
    kept: list[Sequence] = []

    for seq in sequences:
        m = pattern.search(seq.header)
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

    sequences = length_filter(sequences, qc_cfg.get("length_filter", {}), report)

    thresh = qc_cfg.get("ambiguous_threshold", 0.05)
    sequences = ambiguous_filter(sequences, thresh, report)

    sequences = annotation_filter(
        sequences, qc_cfg.get("annotation_filter", {}), report
    )

    report.passed = len(sequences)
    return sequences, report
