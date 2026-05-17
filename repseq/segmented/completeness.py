"""Segmented virus completeness filter and sequence concatenation."""

from __future__ import annotations

import re
from collections import defaultdict
from typing import Any, Optional

from ..clustering.marker import select_marker_protein
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

def _first_non_empty(segments: list[Sequence], attr: str) -> Optional[Any]:
    """Return the first non-empty value of ``attr`` across ``segments``.

    Segments may disagree on per-isolate metadata (host, strain, date,
    country) when one segment was sequenced by a different submitter
    than the others; taking the first non-empty value is the least
    surprising default for the common case of identical metadata, and
    avoids dropping a usable value because segment 0 happens to be
    blank.
    """
    for seg in segments:
        v = getattr(seg, attr, None)
        if v not in (None, ""):
            return v
    return None


def concatenate_isolate(
    segments: list[Sequence],
    isolate_id: str,
    protein_sequence: Optional[str] = None,
) -> Sequence:
    """Concatenate segment sequences for one isolate into a single Sequence.

    Per-isolate metadata (organism, description, strain, host,
    collection_date, country, taxonomy) is taken as the first non-empty
    value across segments — segment 0 is the usual source when metadata
    is consistent, but a blank field on segment 0 falls through to the
    next segment that has one. ``protein_sequence`` is the
    in-segment-order concatenation of the isolate's marker proteins
    (see ``build_concatenated_sequences``); populated when
    protein-alphabet clustering is active, ``None`` otherwise.
    """
    combined_seq = "".join(s.sequence for s in segments)
    representative = segments[0]

    concat_header = (
        f"CONCAT|{isolate_id}|"
        + "|".join(s.accession or s.id for s in segments)
    )

    taxonomy = next(
        (s.taxonomy for s in segments if s.taxonomy is not None),
        None,
    )

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
        organism=_first_non_empty(segments, "organism"),
        description=_first_non_empty(segments, "description"),
        strain=_first_non_empty(segments, "strain"),
        host=_first_non_empty(segments, "host"),
        collection_date=_first_non_empty(segments, "collection_date"),
        country=_first_non_empty(segments, "country"),
        isolate_id=isolate_id,
        is_refseq=all(s.is_refseq for s in segments),
        is_reviewed=all(s.is_reviewed for s in segments),
        taxonomy=taxonomy,
        protein_sequence=protein_sequence,
        # Hand the per-segment Sequence objects to downstream output
        # (phyloXML multi-<sequence> emission). Each carries its own
        # accession + .proteins list, so the writer can list every
        # underlying nuc accession and protein without re-fetching.
        concat_segments=list(segments),
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
        fail_segment: Optional[str] = None
        fail_direction: Optional[str] = None
        for seg_name, seq in zip(segment_names, segs):
            bounds = segment_lengths.get(seg_name)
            if not bounds:
                continue
            mn = bounds.get("min")
            mx = bounds.get("max")
            if mn is not None and seq.length < mn:
                fail_reason = f"segment_length:{seg_name}:{seq.length}<{mn}"
                fail_segment = seg_name
                fail_direction = "too_short"
                break
            if mx is not None and seq.length > mx:
                fail_reason = f"segment_length:{seg_name}:{seq.length}>{mx}"
                fail_segment = seg_name
                fail_direction = "too_long"
                break
        if fail_reason:
            seg_counts = report.removed_length_by_segment.setdefault(
                fail_segment, {"too_short": 0, "too_long": 0}
            )
            seg_counts[fail_direction] += 1
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
    segment_names: Optional[list[str]] = None,
    cluster_protein: Optional[dict[str, list[str]]] = None,
    require_protein: bool = False,
    report: Optional[QCReport] = None,
) -> list[Sequence]:
    """Return one concatenated Sequence per complete isolate.

    When ``require_protein`` is True (protein-alphabet clustering), the
    marker protein for each segment is selected via
    ``select_marker_protein`` and concatenated in ``segment_names`` order
    onto the result's ``protein_sequence``. Isolates whose any segment
    has no qualifying marker are dropped and counted under
    ``report.removed_incomplete_isolates`` with reason
    ``incomplete_isolate:missing_marker_protein:<segment>``.

    ``cluster_protein`` maps segment_name → alias list; absent / empty
    entries fall through to "longest CDS" selection.
    """
    out: list[Sequence] = []
    cp = cluster_protein or {}
    for isolate_id, segs in complete_isolates.items():
        protein_concat: Optional[str] = None
        marker_ids: Optional[list[str]] = None
        if require_protein:
            if segment_names is None:
                segment_names = [seg.segment or "" for seg in segs]
            seg_by_name = {s.segment: s for s in segs if s.segment}
            parts: list[str] = []
            ids: list[str] = []
            missing_marker: Optional[str] = None
            for seg_name in segment_names:
                seg = seg_by_name.get(seg_name)
                if seg is None:
                    missing_marker = seg_name
                    break
                aliases = cp.get(seg_name) or []
                marker = select_marker_protein(seg.proteins, aliases)
                if marker is None:
                    missing_marker = seg_name
                    break
                parts.append(marker["sequence"])
                if marker.get("protein_id"):
                    ids.append(marker["protein_id"])
            if missing_marker is not None:
                if report is not None:
                    reason = (
                        f"incomplete_isolate:missing_marker_protein:{missing_marker}"
                    )
                    for seq in segs:
                        seq.qc_passed = False
                        seq.qc_fail_reason = reason
                        report.removed_incomplete_isolates += 1
                        report.add_removed(seq.id, reason)
                continue
            protein_concat = "".join(parts)
            marker_ids = ids or None
        concat = concatenate_isolate(segs, isolate_id, protein_sequence=protein_concat)
        concat.marker_protein_ids = marker_ids
        out.append(concat)
    return out
