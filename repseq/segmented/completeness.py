"""Segmented virus completeness filter and sequence concatenation."""

from __future__ import annotations

import re
from collections import defaultdict
from typing import Any, Optional

from ..models import QCReport, Sequence


# ---------------------------------------------------------------------------
# Isolate grouping
# ---------------------------------------------------------------------------

def extract_isolate_id(seq: Sequence, isolate_regex: str) -> Optional[str]:
    """Extract isolate identifier from header using a regex pattern.

    The regex should contain a named group 'isolate' or use group 1.
    Falls back to seq.isolate_id if already set.
    """
    if seq.isolate_id:
        return seq.isolate_id

    pattern = re.compile(isolate_regex, re.IGNORECASE)
    target = seq.header

    m = pattern.search(target)
    if not m:
        return None

    try:
        return m.group("isolate")
    except IndexError:
        pass

    if m.lastindex and m.lastindex >= 1:
        return m.group(1)

    return m.group(0)


def _normalise_isolate_id(raw: str) -> str:
    """Normalise an isolate ID for comparison (lowercase, collapse whitespace)."""
    return re.sub(r"\s+", " ", raw.strip().lower())


# ---------------------------------------------------------------------------
# Segment identification
# ---------------------------------------------------------------------------

def identify_segment(seq: Sequence, segment_names: list[str], segment_regex: Optional[str] = None) -> Optional[str]:
    """Return the segment name/number for a sequence.

    Checks seq.segment first, then scans the header using a regex or
    by searching for the segment names directly.
    """
    if seq.segment:
        # Normalise numeric segment references (e.g. "4" -> segment_names[3])
        if seq.segment.isdigit():
            idx = int(seq.segment) - 1
            if 0 <= idx < len(segment_names):
                return segment_names[idx]
        # Check if it's already a known segment name
        for name in segment_names:
            if seq.segment.upper() == name.upper():
                return name
        return seq.segment

    # Search header for segment names
    header_upper = seq.header.upper()

    if segment_regex:
        m = re.search(segment_regex, seq.header, re.IGNORECASE)
        if m:
            try:
                return m.group("segment")
            except IndexError:
                if m.lastindex:
                    raw = m.group(1)
                    if raw.isdigit():
                        idx = int(raw) - 1
                        if 0 <= idx < len(segment_names):
                            return segment_names[idx]
                    return raw

    # Direct name search
    for name in segment_names:
        pattern = re.compile(r"\b" + re.escape(name) + r"\b", re.IGNORECASE)
        if pattern.search(seq.header):
            return name

    # Numeric: "segment N"
    m = re.search(r"\bsegment\s+(\d+)\b", seq.header, re.IGNORECASE)
    if m:
        idx = int(m.group(1)) - 1
        if 0 <= idx < len(segment_names):
            return segment_names[idx]

    return None


# ---------------------------------------------------------------------------
# Completeness filter
# ---------------------------------------------------------------------------

def filter_complete_isolates(
    sequences: list[Sequence],
    virus_cfg: dict[str, Any],
    report: QCReport,
) -> tuple[list[Sequence], dict[str, list[Sequence]]]:
    """Keep only sequences from isolates that have all expected segments.

    Returns:
        (kept_sequences, isolate_map)
        where isolate_map maps isolate_id -> [Sequence, ...] in segment order.
    """
    expected_segments: list[str] = virus_cfg["segments"]
    isolate_regex: str = virus_cfg["isolate_regex"]
    segment_regex: Optional[str] = virus_cfg.get("segment_regex")

    # Group sequences by isolate
    isolate_map: dict[str, dict[str, Sequence]] = defaultdict(dict)
    unresolved: list[Sequence] = []

    for seq in sequences:
        isolate_raw = extract_isolate_id(seq, isolate_regex)
        if isolate_raw is None:
            unresolved.append(seq)
            continue
        seq.isolate_id = isolate_raw

        seg = identify_segment(seq, expected_segments, segment_regex)
        if seg is None:
            unresolved.append(seq)
            continue

        isolate_key = _normalise_isolate_id(isolate_raw)
        if seg in isolate_map[isolate_key]:
            # Keep the longer sequence if duplicated segment
            existing = isolate_map[isolate_key][seg]
            if seq.length > existing.length:
                isolate_map[isolate_key][seg] = seq
        else:
            isolate_map[isolate_key][seg] = seq

    # Keep only complete isolates
    kept_sequences: list[Sequence] = []
    complete_isolates: dict[str, list[Sequence]] = {}

    for isolate_key, seg_map in isolate_map.items():
        present = set(seg_map.keys())
        expected = set(expected_segments)
        if expected.issubset(present):
            ordered = [seg_map[s] for s in expected_segments]
            complete_isolates[isolate_key] = ordered
            kept_sequences.extend(ordered)
        else:
            missing = expected - present
            for seq in seg_map.values():
                reason = f"incomplete_isolate:missing_segments:{','.join(sorted(missing))}"
                seq.qc_passed = False
                seq.qc_fail_reason = reason
                report.removed_incomplete_isolates += 1
                report.add_removed(seq.id, reason)

    # Unresolved sequences are excluded with a warning
    for seq in unresolved:
        reason = "segmented_filter:could_not_identify_isolate_or_segment"
        seq.qc_passed = False
        seq.qc_fail_reason = reason
        report.removed_incomplete_isolates += 1
        report.add_removed(seq.id, reason)

    return kept_sequences, complete_isolates


# ---------------------------------------------------------------------------
# Sequence concatenation
# ---------------------------------------------------------------------------

def concatenate_isolate(
    segments: list[Sequence],
    isolate_id: str,
) -> Sequence:
    """Concatenate segment sequences for one isolate into a single Sequence."""
    combined_seq = "".join(s.sequence for s in segments)
    representative = segments[0]

    concat_header = (
        f"CONCAT|{isolate_id}|"
        + "|".join(s.accession or s.id for s in segments)
    )

    from ..models import SequenceType
    return Sequence(
        id=f"CONCAT|{isolate_id}",
        header=concat_header,
        sequence=combined_seq,
        seq_type=representative.seq_type,
        source=representative.source,
        accession=representative.accession,
        organism=representative.organism,
        description=representative.description,
        strain=representative.strain,
        host=representative.host,
        collection_date=representative.collection_date,
        country=representative.country,
        isolate_id=isolate_id,
        is_refseq=all(s.is_refseq for s in segments),
        is_reviewed=all(s.is_reviewed for s in segments),
        taxonomy=representative.taxonomy,
    )


def build_concatenated_sequences(
    complete_isolates: dict[str, list[Sequence]],
) -> list[Sequence]:
    """Return one concatenated Sequence per complete isolate."""
    return [
        concatenate_isolate(segs, isolate_id)
        for isolate_id, segs in complete_isolates.items()
    ]
