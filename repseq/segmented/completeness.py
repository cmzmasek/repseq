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
    """Normalise an isolate ID for comparison and downstream identification.

    Lowercases and replaces any whitespace run with ``_``; also replaces
    the pipe character (our own CONCAT separator) with ``_``. The result
    is used both as the dict key when grouping segments and — propagated
    through ``concatenate_isolate`` — as the ``seq.id`` of the resulting
    CONCAT sequence. Two downstream parsers depend on these constraints:

      * MMseqs2 takes the first whitespace-delimited token of a FASTA
        header as the sequence identifier, so whitespace inside seq.id
        silently breaks the cluster-TSV round-trip.
      * The segmented FASTA writer and the protein-fasta report read
        the isolate back out of seq.id via ``split("|")[1]``; a stray
        pipe inside the captured name would truncate the recovered id.
    """
    out = re.sub(r"\s+", "_", raw.strip().lower())
    out = out.replace("|", "_")
    return out


# ---------------------------------------------------------------------------
# Segment identification
# ---------------------------------------------------------------------------

def _build_alias_map(
    segment_names: list[str],
    segment_aliases: Optional[dict[str, list[str]]],
) -> dict[str, str]:
    """Return a lower-cased alias → canonical-name map.

    Each canonical name maps to itself; each alias maps to its canonical.
    Used by identify_segment to recognise synonyms like "large segment" → "L".
    """
    out: dict[str, str] = {}
    for canonical in segment_names:
        out[canonical.lower()] = canonical
    if segment_aliases:
        for canonical, syns in segment_aliases.items():
            for syn in syns:
                if syn:
                    out[syn.lower()] = canonical
    return out


def identify_segment(
    seq: Sequence,
    segment_names: list[str],
    segment_regex: Optional[str] = None,
    segment_aliases: Optional[dict[str, list[str]]] = None,
) -> Optional[str]:
    """Return the canonical segment name for a sequence.

    Resolution order:
      1. seq.segment as a numeric index (e.g. "4" → segment_names[3])
      2. seq.segment as a canonical name or alias (case-insensitive)
      3. custom segment_regex applied to the header
      4. Word-boundary search across canonical names AND aliases in the
         header, longest term first (so "large segment" wins over "large")
      5. "segment N" pattern in the header, mapped via the canonical list

    Returns None if nothing matches.
    """
    alias_to_canonical = _build_alias_map(segment_names, segment_aliases)

    if seq.segment:
        if seq.segment.isdigit():
            idx = int(seq.segment) - 1
            if 0 <= idx < len(segment_names):
                return segment_names[idx]
        canon = alias_to_canonical.get(seq.segment.lower())
        if canon:
            return canon
        return seq.segment

    if segment_regex:
        m = re.search(segment_regex, seq.header, re.IGNORECASE)
        if m:
            raw: Optional[str] = None
            try:
                raw = m.group("segment")
            except IndexError:
                if m.lastindex:
                    raw = m.group(1)
            if raw:
                if raw.isdigit():
                    idx = int(raw) - 1
                    if 0 <= idx < len(segment_names):
                        return segment_names[idx]
                return alias_to_canonical.get(raw.lower(), raw)

    # Word-boundary search over canonical names + aliases. Longest term
    # first so multi-word aliases (e.g. "large segment") win over short
    # substrings ("large"). Single-character terms (e.g. the "M" segment)
    # are excluded here: \bM\b matches any stray standalone "M" in a
    # header and would mis-assign the segment. Such segments must instead
    # be resolved via the numeric index, an explicit segment_regex, the
    # seq.segment field, or a multi-character alias.
    candidates = sorted(
        (t for t in alias_to_canonical if len(t) >= 2), key=len, reverse=True
    )
    for term in candidates:
        if re.search(r"\b" + re.escape(term) + r"\b", seq.header, re.IGNORECASE):
            return alias_to_canonical[term]

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
    segment_aliases: Optional[dict[str, list[str]]] = virus_cfg.get("segment_aliases")

    # Group sequences by isolate
    isolate_map: dict[str, dict[str, Sequence]] = defaultdict(dict)
    unresolved: list[Sequence] = []

    for seq in sequences:
        isolate_raw = extract_isolate_id(seq, isolate_regex)
        if isolate_raw is None:
            unresolved.append(seq)
            continue
        seq.isolate_id = isolate_raw

        seg = identify_segment(seq, expected_segments, segment_regex, segment_aliases)
        if seg is None:
            unresolved.append(seq)
            continue

        # Persist the resolved canonical segment so downstream output
        # (e.g. {prefix}_isolate_proteins.tsv) sees it.
        seq.segment = seg

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
        # A concatenated isolate has no single accession — its identity is
        # the isolate_id. Reporting segment 0's accession here would be
        # misleading; per-segment accessions remain in the header and in
        # the per-segment output files.
        accession=None,
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


def segment_length_filter(
    complete_isolates: dict[str, list[Sequence]],
    segment_names: list[str],
    segment_lengths: dict[str, dict[str, int]],
    report: QCReport,
) -> dict[str, list[Sequence]]:
    """Drop complete isolates where any segment falls outside configured length bounds.

    segment_lengths maps segment name → {min: N, max: M}; either bound is optional.
    Dropped sequences are recorded in the QC report under removed_length.
    """
    kept: dict[str, list[Sequence]] = {}
    for isolate_key, segs in complete_isolates.items():
        fail_reason: Optional[str] = None
        for seg_name, seq in zip(segment_names, segs):
            bounds = segment_lengths.get(seg_name)
            if not bounds:
                continue
            mn = bounds.get("min")
            mx = bounds.get("max")
            if mn is not None and seq.length < mn:
                fail_reason = f"segment_length:{seg_name}:{seq.length}<{mn}"
                break
            if mx is not None and seq.length > mx:
                fail_reason = f"segment_length:{seg_name}:{seq.length}>{mx}"
                break
        if fail_reason:
            for seq in segs:
                seq.qc_passed = False
                seq.qc_fail_reason = fail_reason
                report.removed_length += 1
                report.add_removed(seq.id, fail_reason)
        else:
            kept[isolate_key] = segs
    return kept


def build_concatenated_sequences(
    complete_isolates: dict[str, list[Sequence]],
) -> list[Sequence]:
    """Return one concatenated Sequence per complete isolate."""
    return [
        concatenate_isolate(segs, isolate_id)
        for isolate_id, segs in complete_isolates.items()
    ]
