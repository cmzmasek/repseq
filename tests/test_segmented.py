"""Segmented virus completeness filter and concatenation."""
from __future__ import annotations

from repseq.models import QCReport
from repseq.segmented.completeness import (
    build_concatenated_sequences,
    concatenate_isolate,
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
    kept, complete = filter_complete_isolates(seqs, virus_cfg, report)
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

    kept, complete_iso = filter_complete_isolates(
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
    kept, complete_iso = filter_complete_isolates([s], cfg, report)
    assert kept == []
    assert complete_iso == {}
    assert report.removed_incomplete_isolates == 1
    assert s.qc_fail_reason and "could_not_identify" in s.qc_fail_reason


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
    kept, complete = filter_complete_isolates(seqs, cfg, report)
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
