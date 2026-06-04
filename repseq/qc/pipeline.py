"""QC pipeline: duplicate removal, length filter, ambiguous chars, annotation filter."""

from __future__ import annotations

import hashlib
import re
from typing import Any

from ..models import QCReport, Sequence, SequenceType
from ..overrides import ProtectionPolicy, protected_keep


# ---------------------------------------------------------------------------
# Step 1 – Exact duplicate removal
# ---------------------------------------------------------------------------

def remove_duplicates(
    sequences: list[Sequence],
    report: QCReport,
    *,
    policy: ProtectionPolicy | None = None,
) -> list[Sequence]:
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
            reason = f"exact_duplicate_of:{best.id}"
            if protected_keep(seq, "duplicates", reason, policy, report):
                kept.append(seq)
                continue
            seq.qc_passed = False
            seq.qc_fail_reason = reason
            report.removed_duplicates += 1
            report.add_removed(seq.id, reason)
    return kept


# ---------------------------------------------------------------------------
# Step 2 – Length filter
# ---------------------------------------------------------------------------

def genome_length_filter(
    sequences: list[Sequence],
    cfg: dict[str, Any],
    report: QCReport,
    *,
    policy: ProtectionPolicy | None = None,
) -> list[Sequence]:
    """Drop whole sequences outside absolute nucleotide-length bounds.

    ``cfg`` is the ``qc.genome_length_filter`` block: ``{min, max}`` (either
    optional). Non-segmented mode only — the caller (:func:`run_qc`) is
    responsible for not invoking this in segmented mode, where per-segment
    bounds apply instead. There is no median/relative mode: the bounds are
    absolute, so the user must know roughly how long the genome is. With
    neither bound set this is a no-op.
    """
    if not sequences:
        return sequences

    min_len = cfg.get("min")
    max_len = cfg.get("max")
    if min_len is None and max_len is None:
        return sequences

    kept: list[Sequence] = []
    for seq in sequences:
        reason = None
        if min_len is not None and seq.length < min_len:
            reason = f"length_too_short:{seq.length}<{min_len}"
        elif max_len is not None and seq.length > max_len:
            reason = f"length_too_long:{seq.length}>{max_len}"

        if reason:
            if protected_keep(seq, "length", reason, policy, report):
                kept.append(seq)
                continue
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
    sequences: list[Sequence],
    threshold: float,
    report: QCReport,
    *,
    policy: ProtectionPolicy | None = None,
) -> list[Sequence]:
    """Remove sequences with ambiguous character fraction above threshold."""
    kept: list[Sequence] = []
    for seq in sequences:
        frac = seq.ambiguous_fraction
        if frac > threshold:
            reason = f"ambiguous_fraction:{frac:.3f}>{threshold}"
            if protected_keep(seq, "ambiguous", reason, policy, report):
                kept.append(seq)
                continue
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
    sequences: list[Sequence],
    cfg: dict[str, Any],
    report: QCReport,
    *,
    policy: ProtectionPolicy | None = None,
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
            if protected_keep(seq, "annotation", reason, policy, report):
                kept.append(seq)
                continue
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
    segmented = bool(cfg.get("segmented", {}).get("enabled"))
    # Force-keep whitelist: named sequences bypass the QC removal stages
    # they would fail (overrides.protect_qc). Built once from cfg and
    # threaded into each stage; inactive (a no-op) when not configured.
    policy = ProtectionPolicy.from_cfg(cfg)

    if qc_cfg.get("remove_duplicates", True):
        # In segmented mode a single segment can be byte-identical between two
        # otherwise distinct isolates; dropping it here would leave one isolate
        # missing that segment and the whole isolate would be discarded as
        # incomplete. Dedup is instead applied to the concatenated per-isolate
        # sequences after build_concatenated_sequences (see cli._handle_segmented).
        if segmented:
            report.dedup_skipped = True
        else:
            sequences = remove_duplicates(sequences, report, policy=policy)

    # The whole-genome length filter is non-segmented-only and opt-in. In
    # segmented mode the input is a mixed pool of segments with very
    # different lengths (e.g. influenza PB2 ~2300 nt vs NS ~890 nt), so a
    # single whole-sequence bound is meaningless — per-segment bounds apply
    # later via segmented.viruses.<v>.segment_lengths instead. When the
    # filter is disabled (the default) it likewise doesn't run. Either way
    # length_filter_skipped records that the genome filter did not fire.
    glf = qc_cfg.get("genome_length_filter", {}) or {}
    if segmented or not glf.get("enabled", False):
        report.length_filter_skipped = True
    else:
        sequences = genome_length_filter(sequences, glf, report, policy=policy)

    thresh = qc_cfg.get("ambiguous_threshold", 0.05)
    sequences = ambiguous_filter(sequences, thresh, report, policy=policy)

    sequences = annotation_filter(
        sequences, qc_cfg.get("annotation_filter", {}), report, policy=policy
    )

    report.passed = len(sequences)
    return sequences, report
