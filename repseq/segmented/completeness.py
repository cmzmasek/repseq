"""Segmented virus completeness filter and sequence concatenation."""

from __future__ import annotations

import re
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Any, Iterator, Optional

from ..clustering.marker import (
    MarkerFailure,
    _format_failure_reason,
    describe_length_miss,
    select_marker_protein,
)
from ..models import QCReport, Sequence
from ..overrides import _norm_id


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

    Delegates to :func:`repseq.overrides._norm_id`, the single source of
    truth for id normalisation, so the grouping key / CONCAT ``seq.id`` and
    the force_select / protect_qc / exclude matcher cannot drift apart (a
    space-form override id must bind both the raw per-segment ``isolate_id``
    and this underscored CONCAT id). ``_norm_id`` returns ``None`` for empty
    input; an empty isolate name normalises to ``""`` as before.
    """
    return _norm_id(raw) or ""


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
      0. seq.segment equal to a canonical name (case-insensitive) — so a
         config naming its segments "2", "3" maps GenBank "2" to "2", not
         to the second name
      1. seq.segment as a numeric index (e.g. "4" → segment_names[3])
      2. seq.segment as a canonical name or alias (case-insensitive)
      3. custom segment_regex applied to the header
      4. Word-boundary search across canonical names AND aliases in the
         header, longest term first (so "large segment" wins over "large")
      5. "segment N" pattern in the header, mapped via the canonical list

    Returns None if nothing matches.
    """
    alias_to_canonical = _build_alias_map(segment_names, segment_aliases)
    canonical_lc = {name.lower(): name for name in segment_names}

    if seq.segment:
        exact = canonical_lc.get(seq.segment.lower())
        if exact is not None:
            return exact
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
                if raw.lower() in canonical_lc:
                    return canonical_lc[raw.lower()]
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
# Isolate assembly keys — shared by assembly and the isolate-level QC gates
# ---------------------------------------------------------------------------

def resolve_record_key(
    seq: Sequence, virus_cfg: dict[str, Any],
) -> tuple[Optional[str], Optional[str]]:
    """``(isolate_key, canonical_segment)`` for one record, exactly as
    :func:`filter_complete_isolates` will resolve them. Does not mutate.

    ``isolate_key`` is the normalised isolate id (GenBank ``/isolate`` or
    ``/strain``, else ``isolate_regex`` on the header); ``None`` when neither
    yields one. ``canonical_segment`` maps GenBank's ``/segment`` (e.g.
    ``"5"``) or an alias onto the configured segment name (e.g. ``"NP"``).
    """
    regex = virus_cfg.get("isolate_regex")
    if not seq.isolate_id and not regex:
        # No regex configured (validation requires one; hand-built cfgs may
        # lack it): an empty pattern would match every header.
        return None, None
    isolate_raw = extract_isolate_id(seq, regex or "")
    if isolate_raw is None:
        return None, None
    seg = identify_segment(
        seq,
        list(virus_cfg.get("segments") or []),
        virus_cfg.get("segment_regex"),
        virus_cfg.get("segment_aliases"),
    )
    return _normalise_isolate_id(isolate_raw), seg


def prefer_record(candidate: Sequence, existing: Sequence) -> bool:
    """The duplicate-segment rule: when an isolate has two records for one
    segment, assembly keeps the one with the longer nucleotide sequence (the
    first one seen wins a tie)."""
    return candidate.length > existing.length


@dataclass
class IsolateRecords:
    """One isolate as assembly will see it.

    ``records`` is every record of the isolate in input order, paired with
    its canonical segment (``None`` when unidentifiable). ``kept`` maps each
    segment to the ONE record assembly will use — a duplicate that loses to
    :func:`prefer_record` is discarded there, so a QC gate must judge only
    ``kept`` or a partial duplicate would remove an isolate that also has a
    good record for that segment.
    """

    records: list[tuple[Sequence, Optional[str]]] = field(default_factory=list)
    kept: dict[str, Sequence] = field(default_factory=dict)

    def judged(self) -> Iterator[tuple[str, Sequence]]:
        """``(segment, record)`` for each kept record, in input order."""
        for seq, seg in self.records:
            if seg is not None and self.kept.get(seg) is seq:
                yield seg, seq

    def sequences(self) -> list[Sequence]:
        return [seq for seq, _ in self.records]

    def has_all(self, expected: list[str]) -> bool:
        return set(expected) <= set(self.kept)


def group_isolate_records(
    sequences: list[Sequence], virus_cfg: dict[str, Any],
) -> dict[str, IsolateRecords]:
    """Group records into isolates the way :func:`filter_complete_isolates`
    does — same isolate key, same canonical segment, same duplicate rule —
    so the isolate-level QC gates (HMM identity, marker length) judge
    exactly the records assembly will keep. Records whose isolate cannot be
    resolved are left out; assembly drops them as unidentifiable.
    """
    groups: dict[str, IsolateRecords] = {}
    for seq in sequences:
        key, seg = resolve_record_key(seq, virus_cfg)
        if key is None:
            continue
        grp = groups.get(key)
        if grp is None:
            grp = groups[key] = IsolateRecords()
        grp.records.append((seq, seg))
        if seg is None:
            continue
        existing = grp.kept.get(seg)
        if existing is None or prefer_record(seq, existing):
            grp.kept[seg] = seq
    return groups


# ---------------------------------------------------------------------------
# Completeness filter
# ---------------------------------------------------------------------------

def detect_strain_collisions(
    sequences: list[Sequence],
) -> dict[tuple[str, str], list[str]]:
    """Find isolates whose grouping key came from ``/strain`` *and* have
    more than one distinct accession for the same segment — the
    over-merge signature created by the strain-as-isolate-id fallback.

    Returns ``{(isolate_id, segment): [accession, accession, ...]}`` for
    every colliding (isolate_id, segment) pair, with the accession list
    de-duplicated and sorted for stable output. An empty dict means
    nothing to warn about.

    Only sequences with ``isolate_id_source == "strain"`` are considered
    — ``/isolate``-derived ids are submitter-asserted unique per
    biological sample and so are not flagged here, and ``regex``-derived
    ids reflect a header convention the user explicitly chose. Missing
    ``isolate_id`` or ``segment`` (i.e. unresolvable records) are
    skipped: they'll fail the completeness filter for an unrelated
    reason and don't belong in this detector's output.
    """
    by_iso_seg: dict[tuple[str, str], set[str]] = defaultdict(set)
    for seq in sequences:
        if seq.isolate_id_source != "strain":
            continue
        if not seq.isolate_id or not seq.segment:
            continue
        by_iso_seg[(seq.isolate_id, seq.segment)].add(
            seq.accession or seq.id
        )
    return {
        key: sorted(accs)
        for key, accs in by_iso_seg.items()
        if len(accs) > 1
    }


def filter_complete_isolates(
    sequences: list[Sequence],
    virus_cfg: dict[str, Any],
    report: QCReport,
    extra_segments_action: str = "warn",
) -> tuple[list[Sequence], dict[str, list[Sequence]], dict[str, list[str]]]:
    """Keep only sequences from isolates that have all expected segments.

    ``extra_segments_action`` controls what happens when an isolate's
    seg_map contains segment names *outside* the expected list (e.g. a
    fourth segment named "X" for an L/M/S virus, or a non-canonical
    identifier that ``identify_segment`` returned unchanged). On
    ``"warn"`` (default) the isolate keeps its expected segments and
    extras are silently pruned from the concat, but the detection is
    surfaced via the returned ``extras_by_isolate`` dict so the caller
    can emit a warning. On ``"drop"`` the whole isolate is removed,
    every segment lands in ``_qc_removed.tsv`` with reason
    ``extra_segments:<comma-joined extras>``, and
    ``report.removed_extra_segments`` is bumped (one per dropped
    isolate; units are isolates, not segments).

    Returns:
        (kept_sequences, isolate_map, extras_by_isolate)
        where isolate_map maps isolate_id -> [Sequence, ...] in segment
        order and ``extras_by_isolate`` maps isolate_id -> sorted list of
        the extra (non-expected) segment names that were detected. The
        latter is populated regardless of action so the caller can warn
        even when no records were dropped.
    """
    expected_segments: list[str] = virus_cfg["segments"]
    isolate_regex: str = virus_cfg["isolate_regex"]
    segment_regex: Optional[str] = virus_cfg.get("segment_regex")
    segment_aliases: Optional[dict[str, list[str]]] = virus_cfg.get("segment_aliases")

    # Group sequences by isolate
    isolate_map: dict[str, dict[str, Sequence]] = defaultdict(dict)
    unresolved: list[Sequence] = []

    for seq in sequences:
        # Remember whether seq.isolate_id was already populated upstream
        # (GenBank source feature in cli._populate_genbank_isolate_segment)
        # so we only tag the regex fallback when it actually fires.
        had_isolate_id = seq.isolate_id is not None
        isolate_raw = extract_isolate_id(seq, isolate_regex)
        if isolate_raw is None:
            unresolved.append(seq)
            continue
        seq.isolate_id = isolate_raw
        if not had_isolate_id and seq.isolate_id_source is None:
            seq.isolate_id_source = "regex"

        seg = identify_segment(seq, expected_segments, segment_regex, segment_aliases)
        if seg is None:
            unresolved.append(seq)
            continue

        # Persist the resolved canonical segment so downstream output
        # (e.g. {prefix}_isolate_proteins.tsv) sees it.
        seq.segment = seg

        isolate_key = _normalise_isolate_id(isolate_raw)
        if seg in isolate_map[isolate_key]:
            # Duplicated segment: keep the record prefer_record picks (the
            # QC gates judge the same record via group_isolate_records).
            existing = isolate_map[isolate_key][seg]
            if prefer_record(seq, existing):
                isolate_map[isolate_key][seg] = seq
        else:
            isolate_map[isolate_key][seg] = seq

    # Detect isolates whose seg_map carries names outside the expected
    # list (e.g. a fourth segment, or a non-canonical identifier that
    # identify_segment returned unchanged). This runs BEFORE the
    # completeness check so the "drop" action gets to remove the
    # isolate before incomplete_isolate accounting kicks in.
    expected_set = set(expected_segments)
    extras_by_isolate: dict[str, list[str]] = {}
    for isolate_key, seg_map in isolate_map.items():
        extras = sorted(set(seg_map.keys()) - expected_set)
        if extras:
            extras_by_isolate[isolate_key] = extras

    dropped_for_extras: set[str] = set()
    if extras_by_isolate and extra_segments_action == "drop":
        for isolate_key, extras in extras_by_isolate.items():
            seg_map = isolate_map[isolate_key]
            reason = f"extra_segments:{','.join(extras)}"
            for seq in seg_map.values():
                seq.qc_passed = False
                seq.qc_fail_reason = reason
                report.add_removed(seq.accession or seq.id, reason)
            report.removed_extra_segments += 1
            dropped_for_extras.add(isolate_key)

    # Keep only complete isolates
    kept_sequences: list[Sequence] = []
    complete_isolates: dict[str, list[Sequence]] = {}

    for isolate_key, seg_map in isolate_map.items():
        if isolate_key in dropped_for_extras:
            continue
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

    return kept_sequences, complete_isolates, extras_by_isolate


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
        subtype=_first_non_empty(segments, "subtype"),
        collection_date=_first_non_empty(segments, "collection_date"),
        country=_first_non_empty(segments, "country"),
        isolate_id=isolate_id,
        isolate_id_source=_first_non_empty(segments, "isolate_id_source"),
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


def resolve_segment_marker(
    seg_name: str,
    segment_markers: Optional[dict],
    cluster_protein_per_seg: Optional[dict],
) -> tuple[Optional[str], list[str], list[str]]:
    """``(name, hmm_tokens, aliases)`` for one segment's marker.

    ``segment_markers`` wins over the legacy ``cluster_protein`` form. When
    ``segment_markers[seg]`` exists but declares no ``hmms:`` (alias-only),
    the ``hmms:`` of any ``cluster_protein[seg]`` dict entries are pulled in
    (and their aliases merged), so a config that split aliases and HMMs
    across the two maps still gets its HMM gate. For the legacy form,
    plain-string entries are aliases. The single resolver for clustering,
    the QC gates and the output pickers.
    """
    segment_markers = segment_markers or {}
    cluster_protein_per_seg = cluster_protein_per_seg or {}
    sm_name: Optional[str] = None
    sm_tokens: list[str] = []
    sm_aliases: list[str] = []
    if seg_name in segment_markers:
        spec = segment_markers[seg_name] or {}
        sm_name = spec.get("name")
        sm_tokens = list(spec.get("hmms") or [])
        sm_aliases = list(spec.get("aliases") or [])
        if sm_tokens:
            return sm_name, sm_tokens, sm_aliases
        # Fall through to cluster_protein for HMMs (and any extra aliases).
    if seg_name in cluster_protein_per_seg:
        tokens: list[str] = list(sm_tokens)
        aliases: list[str] = list(sm_aliases)
        name: Optional[str] = sm_name
        for entry in cluster_protein_per_seg[seg_name] or []:
            if isinstance(entry, str):
                aliases.append(entry)
            elif isinstance(entry, dict):
                name = name or entry.get("name")
                tokens.extend(entry.get("hmms") or [])
                aliases.extend(entry.get("aliases") or [])
        return name, tokens, aliases
    return sm_name, sm_tokens, sm_aliases


def segment_marker_specs(
    seg_name: str,
    segment_markers: Optional[dict[str, dict]],
    cluster_protein: Optional[dict[str, Any]],
) -> list:
    """The marker spec list ``select_marker_protein`` gets for one segment.

    ``segment_markers[seg]`` wins: one spec named after the segment,
    carrying its ``length_aa``, with tokens and aliases from
    :func:`resolve_segment_marker` (so an alias-only entry picks up
    ``cluster_protein[seg]``'s ``hmms:``, exactly as the HMM gate and the
    output pickers resolve it). Otherwise the legacy ``cluster_protein[seg]``
    list. Shared by isolate assembly and the marker-length QC stage so both
    judge a segment by the same spec.
    """
    sm = segment_markers or {}
    if seg_name in sm:
        spec = sm[seg_name] or {}
        _name, tokens, aliases = resolve_segment_marker(seg_name, sm, cluster_protein)
        return [{
            "name": seg_name,
            "aliases": aliases,
            "hmms": tokens,
            "length_aa": spec.get("length_aa"),
        }]
    return (cluster_protein or {}).get(seg_name) or []


def build_concatenated_sequences(
    complete_isolates: dict[str, list[Sequence]],
    segment_names: Optional[list[str]] = None,
    cluster_protein: Optional[dict[str, Any]] = None,
    require_protein: bool = False,
    report: Optional[QCReport] = None,
    *,
    segment_markers: Optional[dict[str, dict]] = None,
    hmm_active: bool = False,
    ga_cutoffs: Optional[dict[str, Optional[float]]] = None,
    hmm_cfg: Optional[dict[str, Any]] = None,
    policy=None,
) -> list[Sequence]:
    """Return one concatenated Sequence per complete isolate.

    When ``require_protein`` is True (protein-alphabet clustering), the
    marker protein for each segment is selected via
    ``select_marker_protein`` and concatenated in ``segment_names`` order
    onto the result's ``protein_sequence``.

    Per-segment marker specs are looked up in this priority order:
        1. ``segment_markers[seg_name]`` (new HMM-aware form: a dict
           ``{aliases?, hmms?}``) — wins when present.
        2. ``cluster_protein[seg_name]`` (legacy: list of alias strings,
           or with v0.13+ a mixed list including ``{name, aliases?,
           hmms?}`` dicts).
        3. Empty → longest-CDS fallback inside ``select_marker_protein``.

    Isolates whose any segment fails marker selection are dropped.
    HMM-tier failures (configured ``hmms`` had no passing CDS) are
    counted under ``report.removed_hmm_failed`` with per-marker
    breakdown in ``report.removed_hmm_by_marker``; non-HMM failures
    (no aliases matched / no proteins / no translation) keep the
    legacy ``removed_incomplete_isolates`` counter.
    """
    out: list[Sequence] = []
    cp = cluster_protein or {}
    sm = segment_markers or {}
    for isolate_id, segs in complete_isolates.items():
        protein_concat: Optional[str] = None
        marker_ids: Optional[list[str]] = None
        if require_protein:
            if segment_names is None:
                segment_names = [seg.segment or "" for seg in segs]
            seg_by_name = {s.segment: s for s in segs if s.segment}
            parts: list[str] = []
            ids: list[str] = []
            failed_segment: Optional[str] = None
            failure_info: Optional[MarkerFailure] = None
            for seg_name in segment_names:
                seg = seg_by_name.get(seg_name)
                if seg is None:
                    failed_segment = seg_name
                    failure_info = MarkerFailure(
                        reason="no_proteins", marker_name=seg_name
                    )
                    break
                # segment_markers wins over cluster_protein for the same segment.
                marker_specs = segment_marker_specs(seg_name, sm, cp)
                marker, failure = select_marker_protein(
                    seg.proteins,
                    marker_specs,
                    hmm_active=hmm_active,
                    ga_cutoffs=ga_cutoffs,
                    hmm_cfg=hmm_cfg,
                )
                if (
                    marker is None
                    and failure is not None
                    and failure.reason == "length_out_of_bounds"
                    and policy is not None
                    and policy.protects_any(segs, "marker_length")
                ):
                    # Protected against the marker-length stage: pick as if
                    # no bound existed, exactly like the output pickers do
                    # for this isolate. The QC stage already recorded the
                    # protection for isolates it saw; record it here only
                    # for the rare isolate that reached assembly without
                    # passing through that stage (no isolate_id upstream).
                    marker, _ = select_marker_protein(
                        seg.proteins, marker_specs,
                        hmm_active=hmm_active, ga_cutoffs=ga_cutoffs,
                        hmm_cfg=hmm_cfg, ignore_bounds=True,
                    )
                    if marker is not None:
                        if report is not None and not any(
                            pr["stage"] == "marker_length"
                            and pr["id"] in {s.id for s in segs}
                            for pr in report.protected
                        ):
                            report.add_protected(
                                seg.id, "marker_length",
                                _format_failure_reason(failure),
                            )
                        failure = None
                if marker is None:
                    # Force-keep whitelist: this isolate was protected
                    # against the `hmm` stage upstream, but the marker
                    # selector re-applies the HMM gate here and would re-drop
                    # it. Honour the protection by retrying WITHOUT the gate
                    # (longest-CDS / alias marker) so the isolate still gets
                    # a marker to concatenate. Only HMM-gate failures are
                    # rescued — a genuine "no usable CDS" still drops (it
                    # can't be concatenated, like an incomplete isolate).
                    if (
                        failure is not None
                        and failure.reason == "hmm_failed"
                        and policy is not None
                        and policy.protects_any(segs, "hmm")
                    ):
                        # The HMM-off pick still honours length_aa, unless
                        # the isolate is protected against that stage too.
                        marker, retry_failure = select_marker_protein(
                            seg.proteins, marker_specs, hmm_active=False,
                            ignore_bounds=policy.protects_any(
                                segs, "marker_length"
                            ),
                        )
                        if (
                            marker is None
                            and retry_failure is not None
                            and retry_failure.reason == "length_out_of_bounds"
                        ):
                            # HMM protection rescued the identity, but the
                            # protein is outside its range: a length drop.
                            failure = retry_failure
                        if marker is not None:
                            if report is not None:
                                report.add_protected(
                                    seg.accession or seg.id,
                                    "hmm",
                                    f"hmm_failed:{failure.marker_name or seg_name}",
                                )
                            parts.append(marker["sequence"])
                            if marker.get("protein_id"):
                                ids.append(marker["protein_id"])
                            continue
                    failed_segment = seg_name
                    failure_info = failure
                    break
                parts.append(marker["sequence"])
                if marker.get("protein_id"):
                    ids.append(marker["protein_id"])
            if failed_segment is not None:
                if report is not None:
                    is_hmm = (
                        failure_info is not None
                        and failure_info.reason == "hmm_failed"
                    )
                    is_length = (
                        failure_info is not None
                        and failure_info.reason == "length_out_of_bounds"
                    )
                    if is_length:
                        # Same reason shape as the marker-length QC stage;
                        # one isolate, first failing segment (this backstop
                        # only sees isolates the stage could not group).
                        detail = describe_length_miss(
                            failure_info.length_aa, failure_info.bounds
                        )
                        report.removed_marker_length += 1
                        report.removed_marker_length_by_marker[failed_segment] = (
                            report.removed_marker_length_by_marker.get(
                                failed_segment, 0) + 1
                        )
                        for seq in segs:
                            if seq.segment == failed_segment:
                                reason = f"marker_length:{failed_segment}:{detail}"
                            else:
                                reason = (
                                    f"marker_length_sibling:{failed_segment}:"
                                    f"{detail}"
                                )
                            seq.qc_passed = False
                            seq.qc_fail_reason = reason
                            report.add_removed(seq.id, reason)
                        continue
                    if is_hmm:
                        sub = _format_failure_reason(failure_info)
                        reason = f"incomplete_isolate:{sub}:{failed_segment}"
                        marker_key = (
                            failure_info.marker_name or failed_segment
                        )
                        # Count once per isolate (not per segment) so the
                        # per-marker breakdown matches user expectation
                        # ("12 isolates lost their L because RdRp_4 failed").
                        report.removed_hmm_failed += 1
                        report.removed_hmm_by_marker[marker_key] = (
                            report.removed_hmm_by_marker.get(marker_key, 0) + 1
                        )
                    else:
                        reason = (
                            f"incomplete_isolate:missing_marker_protein:"
                            f"{failed_segment}"
                        )
                    for seq in segs:
                        seq.qc_passed = False
                        seq.qc_fail_reason = reason
                        if not is_hmm:
                            report.removed_incomplete_isolates += 1
                        report.add_removed(seq.id, reason)
                continue
            protein_concat = "".join(parts)
            marker_ids = ids or None
        concat = concatenate_isolate(segs, isolate_id, protein_sequence=protein_concat)
        concat.marker_protein_ids = marker_ids
        out.append(concat)
    return out
