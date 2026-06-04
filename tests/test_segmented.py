"""Segmented virus completeness filter and concatenation."""
from __future__ import annotations

from repseq.models import QCReport
from repseq.segmented.completeness import (
    build_concatenated_sequences,
    concatenate_isolate,
    detect_strain_collisions,
    extract_isolate_id,
    filter_complete_isolates,
    identify_segment,
    segment_length_filter,
)


# ---------------------------------------------------------------------------
# Per-helper tests
# ---------------------------------------------------------------------------

ISOLATE_RE = r"(?P<isolate>[A-Z]/[^/]+/[^/]+/[^/]+/\d{4})"


def test_extract_isolate_id_named_group(make_seq):
    s = make_seq("a", "ACGT", header="A/duck/HongKong/1/1997 segment 4")
    assert extract_isolate_id(s, ISOLATE_RE) == "A/duck/HongKong/1/1997"


def test_extract_isolate_id_returns_none_on_miss(make_seq):
    s = make_seq("a", "ACGT", header="totally unrelated header")
    assert extract_isolate_id(s, ISOLATE_RE) is None


def test_identify_segment_numeric_to_name(make_seq):
    s = make_seq("a", "ACGT", segment="4")
    names = ["PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"]
    assert identify_segment(s, names) == "HA"


def test_identify_segment_by_header_name(make_seq):
    s = make_seq("a", "ACGT", header="A/duck/HK/1/97 hemagglutinin HA gene")
    names = ["PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"]
    assert identify_segment(s, names) == "HA"


def test_identify_segment_segment_n_pattern(make_seq):
    s = make_seq("a", "ACGT", header="A/duck/HK/1/97 segment 6 thing")
    names = ["PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"]
    assert identify_segment(s, names) == "NA"


def test_identify_segment_ignores_stray_single_char_name(make_seq):
    # A lone "M" in the header must NOT be read as the "M" segment by the
    # bare word-boundary search — single-char segment names are too
    # ambiguous and are excluded from that step.
    s = make_seq("a", "ACGT", header="some virus isolate M/3/2019 complete genome")
    names = ["PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"]
    assert identify_segment(s, names) is None


# ---------------------------------------------------------------------------
# segment_aliases
# ---------------------------------------------------------------------------

_BUNYA = ["L", "M", "S"]
_BUNYA_ALIASES = {
    "L": ["large", "large segment"],
    "M": ["medium", "medium segment", "glycoprotein"],
    "S": ["small", "small segment", "nucleoprotein"],
}


def test_alias_single_word_in_header(make_seq):
    s = make_seq("a", "ACGT", header="Some bunyavirus large gene something")
    assert identify_segment(s, _BUNYA, segment_aliases=_BUNYA_ALIASES) == "L"


def test_alias_multi_word_in_header(make_seq):
    s = make_seq("a", "ACGT", header="Some bunyavirus large segment polymerase")
    assert identify_segment(s, _BUNYA, segment_aliases=_BUNYA_ALIASES) == "L"


def test_alias_longer_term_wins_over_shorter(make_seq):
    """When both 'large' and 'large segment' could match, the longer alias
    should be tried first so the canonical resolution is consistent."""
    s = make_seq("a", "ACGT", header="virus large segment gene description")
    # Both terms map to "L" here, but the test still verifies the ordering
    # path: if we had aliases mapping "large" → "X" and "large segment" → "L",
    # the longer one would win.
    aliases = {"L": ["large segment"], "X": ["large"]}
    assert identify_segment(s, ["L", "X"], segment_aliases=aliases) == "L"


def test_alias_seq_segment_string_resolved_to_canonical(make_seq):
    """When seq.segment carries an alias string, identify_segment should
    return the canonical name, not the alias verbatim."""
    s = make_seq("a", "ACGT", segment="hemagglutinin")
    names = ["PB2", "HA", "NA"]
    aliases = {"HA": ["hemagglutinin", "haemagglutinin"]}
    assert identify_segment(s, names, segment_aliases=aliases) == "HA"


def test_alias_none_falls_back_to_existing_behavior(make_seq):
    """With no aliases, behavior should match the pre-alias implementation."""
    s = make_seq("a", "ACGT", header="A/duck/HK/1/97 segment 4 HA gene")
    names = ["PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"]
    assert identify_segment(s, names, segment_aliases=None) == "HA"


def test_alias_used_via_completeness_filter(make_seq):
    """End-to-end: aliases should let filter_complete_isolates group
    sequences that only mention protein names in their headers."""
    seqs = [
        make_seq("s1", "A" * 100, header="bunya/X/1/2020 large polymerase"),
        make_seq("s2", "A" * 100, header="bunya/X/1/2020 glycoprotein gene"),
        make_seq("s3", "A" * 100, header="bunya/X/1/2020 nucleoprotein gene"),
    ]
    virus_cfg = {
        "expected_segments": 3,
        "segments": _BUNYA,
        "isolate_regex": r"(?P<isolate>bunya/[^/]+/[^/]+/\d{4})",
        "segment_aliases": _BUNYA_ALIASES,
    }
    report = QCReport()
    kept, complete, _extras = filter_complete_isolates(seqs, virus_cfg, report)
    assert len(complete) == 1
    [(iso_id, ordered)] = complete.items()
    # Ordering follows cfg["segments"] = [L, M, S]
    assert [s.segment for s in ordered] == ["L", "M", "S"]
    assert {s.id for s in kept} == {"s1", "s2", "s3"}


# ---------------------------------------------------------------------------
# filter_complete_isolates
# ---------------------------------------------------------------------------

def _influenza_cfg(segments: list[str]) -> dict:
    return {
        "expected_segments": len(segments),
        "segments": segments,
        "isolate_regex": ISOLATE_RE,
    }


def test_filter_complete_isolates_keeps_complete_drops_incomplete(make_seq):
    segs = ["HA", "NA", "NP"]
    complete = [
        make_seq("c1", "A" * 100, header="A/duck/HK/1/1997 segment 4", segment="4"),
        make_seq("c2", "A" * 100, header="A/duck/HK/1/1997 segment 6", segment="6"),
        make_seq("c3", "A" * 100, header="A/duck/HK/1/1997 segment 5", segment="5"),
    ]
    incomplete = [
        make_seq("i1", "A" * 100, header="A/cat/US/9/2003 segment 4", segment="4"),
        # missing the other two
    ]
    report = QCReport()
    cfg = _influenza_cfg(segs)
    # filter expects only the segments listed in cfg["segments"]; reuse segment names as-is
    # (we set segment numerically; identify_segment maps 4->HA, 5->NP, 6->NA per the list above)
    # Override segments with the 8-name canonical list so 4/5/6 map to HA/NP/NA:
    cfg = {
        "expected_segments": 3,
        "segments": ["HA", "NP", "NA"],
        "isolate_regex": ISOLATE_RE,
    }
    # Re-map: segment "4" -> index 3 (out of range for 3-name list). So set explicit segment names instead.
    for s, name in zip(complete, ["HA", "NA", "NP"]):
        s.segment = name
    for s in incomplete:
        s.segment = "HA"

    kept, complete_iso, _extras = filter_complete_isolates(
        complete + incomplete, cfg, report
    )

    assert len(complete_iso) == 1
    [(iso_id, ordered)] = complete_iso.items()
    # Ordering follows cfg["segments"]: HA, NP, NA
    assert [s.segment for s in ordered] == ["HA", "NP", "NA"]
    assert {s.id for s in kept} == {"c1", "c2", "c3"}
    assert report.removed_incomplete_isolates == 1
    assert incomplete[0].qc_passed is False


def test_filter_complete_isolates_unresolved_marked_incomplete(make_seq):
    cfg = {
        "expected_segments": 2,
        "segments": ["HA", "NA"],
        "isolate_regex": ISOLATE_RE,
    }
    # header doesn't match isolate regex
    s = make_seq("u", "A" * 100, header="some unrelated header", segment="HA")
    report = QCReport()
    kept, complete_iso, _extras = filter_complete_isolates([s], cfg, report)
    assert kept == []
    assert complete_iso == {}
    assert report.removed_incomplete_isolates == 1
    assert s.qc_fail_reason and "could_not_identify" in s.qc_fail_reason


# ---------------------------------------------------------------------------
# extra_segments_action
# ---------------------------------------------------------------------------

def _extras_cfg() -> dict:
    return {
        "expected_segments": 3,
        "segments": ["L", "M", "S"],
        "isolate_regex": ISOLATE_RE,
    }


def _build_extras_seqs(make_seq):
    """L/M/S complete isolate plus a stray fourth segment named 'X'."""
    seqs = [
        make_seq("a_l", "A" * 100, header="A/duck/HK/1/1997 segment L"),
        make_seq("a_m", "A" * 100, header="A/duck/HK/1/1997 segment M"),
        make_seq("a_s", "A" * 100, header="A/duck/HK/1/1997 segment S"),
        make_seq("a_x", "A" * 100, header="A/duck/HK/1/1997 segment X"),
    ]
    for s, name in zip(seqs, ["L", "M", "S", "X"]):
        s.segment = name
    return seqs


def test_extra_segments_warn_keeps_isolate_with_expected_segments(make_seq):
    """warn (default): the isolate passes completeness with just its
    expected segments. The extra segment is silently pruned from the
    concat — same behaviour as before v0.14.1 — but the detection IS
    surfaced via the returned extras_by_isolate dict."""
    seqs = _build_extras_seqs(make_seq)
    report = QCReport()
    kept, complete, extras = filter_complete_isolates(
        seqs, _extras_cfg(), report, extra_segments_action="warn"
    )
    # Isolate kept, only expected segments in the ordered list.
    assert len(complete) == 1
    [(iso_id, ordered)] = complete.items()
    assert [s.segment for s in ordered] == ["L", "M", "S"]
    # Extras dict surfaces the detection.
    assert extras == {iso_id: ["X"]}
    # No counter increment, no record dropped.
    assert report.removed_extra_segments == 0
    assert report.removed_incomplete_isolates == 0


def test_extra_segments_drop_removes_whole_isolate(make_seq):
    """drop: the whole isolate is removed. Every segment (expected AND
    extra) lands in _qc_removed.tsv with reason
    extra_segments:<extras>, and the counter increments once per
    isolate (units: isolates)."""
    seqs = _build_extras_seqs(make_seq)
    report = QCReport()
    kept, complete, extras = filter_complete_isolates(
        seqs, _extras_cfg(), report, extra_segments_action="drop"
    )
    assert complete == {}
    assert kept == []
    assert len(extras) == 1
    assert list(extras.values()) == [["X"]]
    # Counted in isolates.
    assert report.removed_extra_segments == 1
    # No incomplete-isolate bookkeeping for this drop path.
    assert report.removed_incomplete_isolates == 0
    # Every segment dropped, reason carries the extras list.
    assert all(s.qc_passed is False for s in seqs)
    reasons = {s.qc_fail_reason for s in seqs}
    assert reasons == {"extra_segments:X"}
    # _qc_removed.tsv entries cover every segment.
    detail_reasons = {d["reason"] for d in report.details}
    assert detail_reasons == {"extra_segments:X"}


def test_extra_segments_default_is_warn_when_kwarg_omitted(make_seq):
    """The default action is 'warn' so back-compat is preserved for
    callers that don't pass the kwarg."""
    seqs = _build_extras_seqs(make_seq)
    report = QCReport()
    _, complete, extras = filter_complete_isolates(seqs, _extras_cfg(), report)
    assert len(complete) == 1
    assert extras and report.removed_extra_segments == 0


def test_extra_segments_returns_empty_dict_when_no_extras(make_seq):
    """An L/M/S-only isolate produces an empty extras dict — both
    actions are no-ops in this case."""
    seqs = _build_extras_seqs(make_seq)[:3]  # drop the X segment
    report = QCReport()
    _, complete, extras = filter_complete_isolates(
        seqs, _extras_cfg(), report, extra_segments_action="drop"
    )
    assert extras == {}
    assert len(complete) == 1
    assert report.removed_extra_segments == 0


def test_extra_segments_drop_multiple_extras_join_in_reason(make_seq):
    """When an isolate carries more than one extra segment, the
    _qc_removed.tsv reason joins the names in sorted order."""
    seqs = [
        make_seq("a_l", "A" * 100, header="A/duck/HK/1/1997 segment L"),
        make_seq("a_m", "A" * 100, header="A/duck/HK/1/1997 segment M"),
        make_seq("a_s", "A" * 100, header="A/duck/HK/1/1997 segment S"),
        make_seq("a_x", "A" * 100, header="A/duck/HK/1/1997 segment X"),
        make_seq("a_y", "A" * 100, header="A/duck/HK/1/1997 segment Y"),
    ]
    for s, name in zip(seqs, ["L", "M", "S", "X", "Y"]):
        s.segment = name
    report = QCReport()
    _, _, extras = filter_complete_isolates(
        seqs, _extras_cfg(), report, extra_segments_action="drop"
    )
    [(iso_id, names)] = extras.items()
    assert names == ["X", "Y"]
    reasons = {s.qc_fail_reason for s in seqs}
    assert reasons == {"extra_segments:X,Y"}


# ---------------------------------------------------------------------------
# Concatenation
# ---------------------------------------------------------------------------

def test_concatenate_isolate_joins_in_order(make_seq):
    a = make_seq("a", "AAA")
    b = make_seq("b", "CCC")
    c = make_seq("c", "GGG")
    out = concatenate_isolate([a, b, c], "ISO1")
    assert out.sequence == "AAACCCGGG"
    assert out.id == "CONCAT|ISO1"
    assert "ISO1" in out.header
    assert out.isolate_id == "ISO1"


def test_filter_complete_isolates_whitespace_in_isolate_id_is_underscored(make_seq):
    # Regression: viral isolate names like "yaba-7 virus strain yaba 7"
    # contain spaces. The normalised isolate_id (and therefore the
    # CONCAT seq.id derived from it) must be whitespace-free, otherwise
    # MMseqs2 truncates the FASTA header token and the cluster TSV
    # round-trip silently drops every cluster.
    cfg = {
        "expected_segments": 2,
        "segments": ["L", "S"],
        "isolate_regex": r"\|(?P<isolate>.+?)(?:\s+segment|\|)",
    }
    seqs = [
        make_seq("acc1", "A" * 100,
                 header="X|yaba-7 virus strain yaba 7|acc1 segment L",
                 segment="L"),
        make_seq("acc2", "C" * 100,
                 header="X|yaba-7 virus strain yaba 7|acc2 segment S",
                 segment="S"),
    ]
    report = QCReport()
    kept, complete, _extras = filter_complete_isolates(seqs, cfg, report)
    assert len(complete) == 1
    [(iso_id, _segs)] = complete.items()
    assert " " not in iso_id
    assert iso_id == "yaba-7_virus_strain_yaba_7"
    concats = build_concatenated_sequences(complete)
    assert len(concats) == 1
    assert " " not in concats[0].id


def test_concatenate_isolate_has_no_single_accession(make_seq):
    # A concatenated isolate spans several segments — it must not claim a
    # single segment's accession. Identity is carried by isolate_id; the
    # per-segment accessions remain recoverable from the header.
    a = make_seq("acc1", "AAA", accession="acc1")
    b = make_seq("acc2", "CCC", accession="acc2")
    out = concatenate_isolate([a, b], "ISO1")
    assert out.accession is None
    assert out.isolate_id == "ISO1"
    assert "acc1" in out.header and "acc2" in out.header


def test_concatenate_isolate_inherits_first_non_empty_metadata(make_seq):
    """When per-isolate metadata is blank on segment 0 but present on a
    later segment, the concat should pick up the later value rather than
    leave the field blank. The typical case is a single submission where
    all segments share metadata; this guards the long tail where only
    one segment was annotated."""
    a = make_seq("acc1", "AAA", accession="acc1",
                 host=None, country=None, organism=None)
    b = make_seq("acc2", "CCC", accession="acc2",
                 host="Apodemus agrarius", country="South Korea",
                 organism="Hantaan virus")
    out = concatenate_isolate([a, b], "ISO1")
    assert out.host == "Apodemus agrarius"
    assert out.country == "South Korea"
    assert out.organism == "Hantaan virus"


def test_concatenate_isolate_prefers_segment_0_when_all_set(make_seq):
    """When every segment has metadata, segment 0 still wins — first
    non-empty preserves the historical 'take the first segment' default."""
    a = make_seq("acc1", "AAA", accession="acc1", host="rodent A")
    b = make_seq("acc2", "CCC", accession="acc2", host="rodent B")
    out = concatenate_isolate([a, b], "ISO1")
    assert out.host == "rodent A"


def test_concatenate_isolate_stores_segments_for_downstream(make_seq):
    """The CONCAT record carries references to its source segments so
    the phyloXML writer can list every nuc accession + protein without
    re-fetching."""
    a = make_seq("acc1", "AAA", accession="acc1")
    b = make_seq("acc2", "CCC", accession="acc2")
    out = concatenate_isolate([a, b], "ISO1")
    assert out.concat_segments == [a, b]
    assert out.concat_segments[0] is a  # identity, not copy
    assert out.concat_segments[1] is b


def test_build_concat_protein_records_marker_protein_ids(make_seq):
    """For protein-alphabet runs, the concat record's
    marker_protein_ids lists the chosen marker's protein_id per
    segment in segment order."""
    L = make_seq("L_acc", "A" * 100, accession="L_acc", segment="L")
    L.proteins = [
        {"protein_id": "L_pol", "product": "polymerase",
         "length": 100, "sequence": "M" * 100},
    ]
    S = make_seq("S_acc", "G" * 100, accession="S_acc", segment="S")
    S.proteins = [
        {"protein_id": "S_N", "product": "nucleoprotein",
         "length": 40, "sequence": "M" * 40},
        {"protein_id": "S_NSs", "product": "NSs",
         "length": 20, "sequence": "M" * 20},
    ]
    complete = {"iso1": [L, S]}
    report = QCReport()
    out = build_concatenated_sequences(
        complete, segment_names=["L", "S"], require_protein=True, report=report,
    )
    assert len(out) == 1
    concat = out[0]
    # Markers: L's longest CDS, S's longest CDS.
    assert concat.marker_protein_ids == ["L_pol", "S_N"]


# ---------------------------------------------------------------------------
# segment_length_filter
# ---------------------------------------------------------------------------

def test_segment_length_filter_keeps_in_bounds(make_seq):
    seg_names = ["HA", "NA"]
    iso = {"iso1": [make_seq("h", "A" * 1700, segment="HA"),
                    make_seq("n", "A" * 1400, segment="NA")]}
    bounds = {"HA": {"min": 1600, "max": 1800}, "NA": {"min": 1300, "max": 1500}}
    report = QCReport()
    result = segment_length_filter(iso, seg_names, bounds, report)
    assert "iso1" in result
    assert report.removed_length == 0


def test_segment_length_filter_drops_too_short(make_seq):
    seg_names = ["HA", "NA"]
    iso = {"iso1": [make_seq("h", "A" * 500, segment="HA"),
                    make_seq("n", "A" * 1400, segment="NA")]}
    bounds = {"HA": {"min": 1600}}
    report = QCReport()
    result = segment_length_filter(iso, seg_names, bounds, report)
    assert result == {}
    assert report.removed_length == 2  # both seqs in the dropped isolate


def test_segment_length_filter_drops_too_long(make_seq):
    seg_names = ["HA", "NA"]
    iso = {"iso1": [make_seq("h", "A" * 2000, segment="HA"),
                    make_seq("n", "A" * 1400, segment="NA")]}
    bounds = {"HA": {"max": 1800}}
    report = QCReport()
    result = segment_length_filter(iso, seg_names, bounds, report)
    assert result == {}
    assert report.removed_length == 2


def test_segment_length_filter_partial_bounds_and_mixed_isolates(make_seq):
    seg_names = ["HA", "NA"]
    iso = {
        "good": [make_seq("g_ha", "A" * 1700, segment="HA"),
                 make_seq("g_na", "A" * 1400, segment="NA")],
        "bad":  [make_seq("b_ha", "A" * 500,  segment="HA"),
                 make_seq("b_na", "A" * 1400, segment="NA")],
    }
    bounds = {"HA": {"min": 1600, "max": 1800}}
    report = QCReport()
    result = segment_length_filter(iso, seg_names, bounds, report)
    assert set(result.keys()) == {"good"}
    assert report.removed_length == 2  # two seqs from the bad isolate


def test_segment_length_filter_per_segment_counter_too_short(make_seq):
    seg_names = ["HA", "NA"]
    iso = {
        "i1": [make_seq("a", "A" * 500,  segment="HA"),
               make_seq("b", "A" * 1400, segment="NA")],
        "i2": [make_seq("c", "A" * 400,  segment="HA"),
               make_seq("d", "A" * 1400, segment="NA")],
    }
    bounds = {"HA": {"min": 1600}}
    report = QCReport()
    segment_length_filter(iso, seg_names, bounds, report)
    # Counted in isolates, not segments: 2 isolates lost their HA.
    assert report.removed_length_by_segment == {
        "HA": {"too_short": 2, "too_long": 0}
    }


def test_segment_length_filter_per_segment_counter_too_long(make_seq):
    seg_names = ["HA", "NA"]
    iso = {
        "i1": [make_seq("a", "A" * 2000, segment="HA"),
               make_seq("b", "A" * 1400, segment="NA")],
    }
    bounds = {"HA": {"max": 1800}}
    report = QCReport()
    segment_length_filter(iso, seg_names, bounds, report)
    assert report.removed_length_by_segment == {
        "HA": {"too_short": 0, "too_long": 1}
    }


def test_segment_length_filter_per_segment_counter_mixed_segments(make_seq):
    seg_names = ["HA", "NA"]
    iso = {
        "i1": [make_seq("a", "A" * 500,  segment="HA"),
               make_seq("b", "A" * 1400, segment="NA")],
        "i2": [make_seq("c", "A" * 1700, segment="HA"),
               make_seq("d", "A" * 2000, segment="NA")],
    }
    bounds = {"HA": {"min": 1600}, "NA": {"max": 1500}}
    report = QCReport()
    segment_length_filter(iso, seg_names, bounds, report)
    assert report.removed_length_by_segment == {
        "HA": {"too_short": 1, "too_long": 0},
        "NA": {"too_short": 0, "too_long": 1},
    }


def test_segment_length_filter_first_failing_segment_wins(make_seq):
    # When an isolate has multiple bad segments, only the first one
    # encountered (in segment_names order) is counted — the iteration
    # short-circuits via ``break``, which keeps the qc_removed reason
    # and the per-segment counter consistent.
    seg_names = ["HA", "NA"]
    iso = {
        "i1": [make_seq("a", "A" * 500,  segment="HA"),
               make_seq("b", "A" * 100,  segment="NA")],
    }
    bounds = {"HA": {"min": 1600}, "NA": {"min": 1000}}
    report = QCReport()
    segment_length_filter(iso, seg_names, bounds, report)
    assert report.removed_length_by_segment == {
        "HA": {"too_short": 1, "too_long": 0}
    }


def test_segment_length_filter_no_bounds_keeps_counter_empty(make_seq):
    seg_names = ["HA", "NA"]
    iso = {"i1": [make_seq("a", "A" * 1700, segment="HA"),
                  make_seq("b", "A" * 1400, segment="NA")]}
    bounds = {}
    report = QCReport()
    result = segment_length_filter(iso, seg_names, bounds, report)
    assert "i1" in result
    assert report.removed_length_by_segment == {}


def test_qc_summary_per_segment_breakdown_replaces_skipped_line():
    # When the segmented length filter actually fired, the summary must
    # show the per-segment breakdown, not the misleading "skipped" line.
    report = QCReport(length_filter_skipped=True)
    report.removed_length_by_segment = {
        "L": {"too_short": 257, "too_long": 0},
        "M": {"too_short": 28, "too_long": 6},
    }
    out = report.summary()
    assert "skipped (segmented mode)" not in out
    assert "Removed (length)    : 291 isolate(s) (per-segment, isolate-level)" in out
    assert "L too short  : 257" in out
    assert "M too short  : 28" in out
    assert "M too long   : 6" in out


def test_qc_summary_keeps_skipped_when_no_segment_bounds_configured():
    # No per-segment counter populated -> fall back to the "skipped" line.
    report = QCReport(length_filter_skipped=True)
    out = report.summary()
    assert "Removed (length)    : skipped (segmented mode)" in out


def test_build_concatenated_sequences_one_per_isolate(make_seq):
    iso_map = {
        "ISO1": [make_seq("a1", "AAA"), make_seq("a2", "TTT")],
        "ISO2": [make_seq("b1", "GGG"), make_seq("b2", "CCC")],
    }
    out = build_concatenated_sequences(iso_map)
    assert len(out) == 2
    by_iso = {s.isolate_id: s for s in out}
    assert by_iso["ISO1"].sequence == "AAATTT"
    assert by_iso["ISO2"].sequence == "GGGCCC"
    # No protein concat unless explicitly requested.
    assert all(s.protein_sequence is None for s in out)


# ---------------------------------------------------------------------------
# Protein-alphabet concat (segmented mode)
# ---------------------------------------------------------------------------

def _cds(product: str, aa: str):
    return {"protein_id": "P", "product": product, "length": len(aa), "sequence": aa}


def test_build_concat_protein_uses_longest_cds_per_segment(make_seq):
    s_seg = make_seq("a_S", "AAA", segment="S")
    l_seg = make_seq("a_L", "TTT", segment="L")
    s_seg.proteins = [_cds("nucleoprotein", "N" * 400)]
    l_seg.proteins = [
        _cds("ns1 short", "M" * 50),
        _cds("polymerase", "M" * 2200),
    ]
    iso = {"iso1": [s_seg, l_seg]}
    out = build_concatenated_sequences(
        iso, segment_names=["S", "L"], require_protein=True,
    )
    assert len(out) == 1
    # AA concat is N*400 (S segment marker) + M*2200 (L segment, longest CDS).
    assert out[0].protein_sequence == "N" * 400 + "M" * 2200


def test_build_concat_protein_alias_per_segment(make_seq):
    s_seg = make_seq("a_S", "AAA", segment="S")
    l_seg = make_seq("a_L", "TTT", segment="L")
    s_seg.proteins = [
        _cds("nucleoprotein", "N" * 400),
        _cds("NSs", "X" * 100),
    ]
    l_seg.proteins = [_cds("RNA-dependent RNA polymerase", "M" * 2200)]
    iso = {"iso1": [s_seg, l_seg]}
    out = build_concatenated_sequences(
        iso, segment_names=["S", "L"], require_protein=True,
        cluster_protein={"S": ["NSs"], "L": ["polymerase"]},
    )
    # S marker is NSs (alias wins despite N being longer); L is polymerase.
    assert out[0].protein_sequence == "X" * 100 + "M" * 2200


def test_build_concat_protein_drops_isolate_when_marker_missing(make_seq):
    s_seg = make_seq("a_S", "AAA", segment="S")
    l_seg = make_seq("a_L", "TTT", segment="L")
    s_seg.proteins = [_cds("nucleoprotein", "N" * 400)]
    l_seg.proteins = []  # fetched, no CDS — no viable marker for L
    iso = {"iso1": [s_seg, l_seg]}
    report = QCReport()
    out = build_concatenated_sequences(
        iso, segment_names=["S", "L"], require_protein=True, report=report,
    )
    assert out == []
    assert report.removed_incomplete_isolates == 2
    assert all(
        "missing_marker_protein:L" in d["reason"] for d in report.details
    )


def test_build_concat_protein_no_aliases_match_falls_back_to_longest(make_seq):
    s_seg = make_seq("a_S", "AAA", segment="S")
    l_seg = make_seq("a_L", "TTT", segment="L")
    s_seg.proteins = [_cds("nucleoprotein", "N" * 400)]
    l_seg.proteins = [_cds("L protein", "M" * 2200)]
    iso = {"iso1": [s_seg, l_seg]}
    out = build_concatenated_sequences(
        iso, segment_names=["S", "L"], require_protein=True,
        # Alias for S doesn't match nucleoprotein — should fall back to longest.
        cluster_protein={"S": ["does_not_match"]},
    )
    assert out[0].protein_sequence == "N" * 400 + "M" * 2200


def _hmm_hit(target, passing):
    return {
        "target": target, "passing": passing, "dom_evalue": 1e-3,
        "dom_score": 50.0, "hmm_len": 300, "ali_span": 280,
        "ali_from": 1, "ali_to": 280,
    }


def _hmm_gated_isolate(make_seq):
    """An isolate whose L segment fails the HMM gate (hit present but not
    passing), so build_concatenated_sequences would normally re-drop it."""
    s_seg = make_seq("a_S", "AAA", segment="S")
    l_seg = make_seq("a_L", "TTT", segment="L")
    s_seg.accession = "ACC_S"
    l_seg.accession = "ACC_L"
    s_seg.isolate_id = l_seg.isolate_id = "iso1"
    s_seg.proteins = [
        {"protein_id": "P_N", "product": "nucleoprotein", "length": 400,
         "sequence": "N" * 400, "hmm_hits": [_hmm_hit("Bunya_N", True)]},
    ]
    l_seg.proteins = [
        {"protein_id": "P_L", "product": "polymerase", "length": 2200,
         "sequence": "M" * 2200, "hmm_hits": [_hmm_hit("RdRP_4", False)]},
    ]
    sm = {"S": {"hmms": ["Bunya_N"]}, "L": {"hmms": ["RdRP_4"]}}
    return {"iso1": [s_seg, l_seg]}, sm


def test_build_concat_hmm_gate_drops_unprotected_isolate(make_seq):
    """Baseline: with no override policy, an isolate failing the HMM gate
    at marker selection is dropped under removed_hmm_failed."""
    iso, sm = _hmm_gated_isolate(make_seq)
    report = QCReport()
    out = build_concatenated_sequences(
        iso, segment_names=["S", "L"], require_protein=True, report=report,
        segment_markers=sm, hmm_active=True,
    )
    assert out == []
    assert report.removed_hmm_failed == 1


def test_build_concat_hmm_protection_keeps_isolate_with_longest_cds(make_seq):
    """The force-keep whitelist must survive the marker-selection backstop:
    a protected isolate failing the HMM gate is rescued by retrying without
    the gate (longest CDS), kept, and recorded on report.protected."""
    from repseq.overrides import ProtectionPolicy, resolve_ids, resolve_stages

    iso, sm = _hmm_gated_isolate(make_seq)
    report = QCReport()
    policy = ProtectionPolicy(
        resolve_ids({"ids": ["iso1"]}), resolve_stages("all"), enabled=True,
    )
    out = build_concatenated_sequences(
        iso, segment_names=["S", "L"], require_protein=True, report=report,
        segment_markers=sm, hmm_active=True, policy=policy,
    )
    assert len(out) == 1
    # L marker is the longest CDS (M*2200), gate bypassed for the protected isolate.
    assert out[0].protein_sequence == "N" * 400 + "M" * 2200
    assert report.removed_hmm_failed == 0
    assert any(p["stage"] == "hmm" for p in report.protected)


def test_build_concat_hmm_protection_only_for_hmm_stage(make_seq):
    """A policy that protects a DIFFERENT stage must not rescue the HMM
    backstop — the isolate still drops."""
    from repseq.overrides import ProtectionPolicy, resolve_ids, resolve_stages

    iso, sm = _hmm_gated_isolate(make_seq)
    report = QCReport()
    policy = ProtectionPolicy(
        resolve_ids({"ids": ["iso1"]}), resolve_stages(["ambiguous"]),
        enabled=True,
    )
    out = build_concatenated_sequences(
        iso, segment_names=["S", "L"], require_protein=True, report=report,
        segment_markers=sm, hmm_active=True, policy=policy,
    )
    assert out == []
    assert report.removed_hmm_failed == 1


# ---------------------------------------------------------------------------
# isolate_id_source provenance + strain-collision detector
# ---------------------------------------------------------------------------

_ISO_RE_SIMPLE = r"isolate=(?P<isolate>[A-Za-z0-9_-]+)"


def test_filter_complete_isolates_tags_regex_source(make_seq):
    # Headers carry isolate=X / segment=Y — GenBank pre-pass did NOT run,
    # so the regex fires inside filter_complete_isolates and tags source.
    a_l = make_seq("a_l", "AAA", header="isolate=I1 segment=L")
    a_m = make_seq("a_m", "CCC", header="isolate=I1 segment=M")
    a_s = make_seq("a_s", "GGG", header="isolate=I1 segment=S")
    virus_cfg = {
        "segments": ["L", "M", "S"],
        "isolate_regex": _ISO_RE_SIMPLE,
        "segment_regex": r"segment=(?P<segment>[LMS])",
    }
    report = QCReport()
    _, complete, _extras = filter_complete_isolates(
        [a_l, a_m, a_s], virus_cfg, report
    )
    assert "i1" in complete
    for seg in complete["i1"]:
        assert seg.isolate_id_source == "regex"


def test_filter_complete_isolates_preserves_existing_source(make_seq):
    # If the GenBank pre-pass already populated isolate_id_source, the
    # regex stage must not overwrite it.
    a_l = make_seq("a_l", "AAA", header="isolate=I1 segment=L", isolate_id="I1")
    a_l.isolate_id_source = "strain"
    a_m = make_seq("a_m", "CCC", header="isolate=I1 segment=M", isolate_id="I1")
    a_m.isolate_id_source = "strain"
    a_s = make_seq("a_s", "GGG", header="isolate=I1 segment=S", isolate_id="I1")
    a_s.isolate_id_source = "strain"
    virus_cfg = {
        "segments": ["L", "M", "S"],
        "isolate_regex": _ISO_RE_SIMPLE,
        "segment_regex": r"segment=(?P<segment>[LMS])",
    }
    report = QCReport()
    _, complete, _extras = filter_complete_isolates(
        [a_l, a_m, a_s], virus_cfg, report
    )
    for seg in complete["i1"]:
        assert seg.isolate_id_source == "strain"


def test_concatenate_isolate_inherits_isolate_id_source(make_seq):
    segs = [
        make_seq("a_L", "A", segment="L", isolate_id="I1"),
        make_seq("a_M", "C", segment="M", isolate_id="I1"),
    ]
    for s in segs:
        s.isolate_id_source = "strain"
    concat = concatenate_isolate(segs, "I1")
    assert concat.isolate_id_source == "strain"


def test_concatenate_isolate_id_source_none_when_unset(make_seq):
    segs = [
        make_seq("a_L", "A", segment="L", isolate_id="I1"),
        make_seq("a_M", "C", segment="M", isolate_id="I1"),
    ]
    concat = concatenate_isolate(segs, "I1")
    assert concat.isolate_id_source is None


def test_detect_strain_collisions_flags_same_strain_same_segment(make_seq):
    # Two distinct accessions, both with /strain=L99, both on segment S.
    a = make_seq("acc1", "AAA", segment="S", isolate_id="L99", accession="ACC1")
    b = make_seq("acc2", "CCC", segment="S", isolate_id="L99", accession="ACC2")
    a.isolate_id_source = "strain"
    b.isolate_id_source = "strain"
    out = detect_strain_collisions([a, b])
    assert out == {("L99", "S"): ["ACC1", "ACC2"]}


def test_detect_strain_collisions_ignores_isolate_source(make_seq):
    # Two accessions with same id+segment but source=isolate are NOT
    # flagged — /isolate is submitter-asserted unique.
    a = make_seq("acc1", "AAA", segment="S", isolate_id="L99", accession="ACC1")
    b = make_seq("acc2", "CCC", segment="S", isolate_id="L99", accession="ACC2")
    a.isolate_id_source = "isolate"
    b.isolate_id_source = "isolate"
    assert detect_strain_collisions([a, b]) == {}


def test_detect_strain_collisions_ignores_different_segments(make_seq):
    # Two strain-sourced accessions sharing isolate_id but on different
    # segments — that's the normal segmented-virus pattern, not a
    # collision.
    a = make_seq("acc1", "AAA", segment="L", isolate_id="L99", accession="ACC1")
    b = make_seq("acc2", "CCC", segment="M", isolate_id="L99", accession="ACC2")
    a.isolate_id_source = "strain"
    b.isolate_id_source = "strain"
    assert detect_strain_collisions([a, b]) == {}


def test_detect_strain_collisions_skips_unresolved_records(make_seq):
    # Records missing isolate_id or segment are silently skipped — they
    # fail completeness for a different reason and shouldn't appear here.
    a = make_seq("acc1", "AAA", segment="S", isolate_id=None, accession="ACC1")
    b = make_seq("acc2", "CCC", segment=None, isolate_id="L99", accession="ACC2")
    a.isolate_id_source = "strain"
    b.isolate_id_source = "strain"
    assert detect_strain_collisions([a, b]) == {}


def test_detect_strain_collisions_multiple_groups_independent(make_seq):
    # Two separate collisions on different (isolate, segment) pairs;
    # detector returns both, sorted accession lists.
    a1 = make_seq("a1", "A", segment="S", isolate_id="L99", accession="A1")
    a2 = make_seq("a2", "A", segment="S", isolate_id="L99", accession="A2")
    b1 = make_seq("b1", "C", segment="L", isolate_id="M77", accession="B2")
    b2 = make_seq("b2", "C", segment="L", isolate_id="M77", accession="B1")
    for s in (a1, a2, b1, b2):
        s.isolate_id_source = "strain"
    out = detect_strain_collisions([a1, a2, b1, b2])
    assert out == {
        ("L99", "S"): ["A1", "A2"],
        ("M77", "L"): ["B1", "B2"],
    }


def test_detect_strain_collisions_dedups_repeated_accessions(make_seq):
    # Two records with the same accession (e.g. the same nuc record
    # cited twice in input) should not be reported as a collision —
    # the detector dedups by accession before counting.
    a = make_seq("a", "A", segment="S", isolate_id="L99", accession="DUP")
    b = make_seq("b", "C", segment="S", isolate_id="L99", accession="DUP")
    a.isolate_id_source = "strain"
    b.isolate_id_source = "strain"
    assert detect_strain_collisions([a, b]) == {}
