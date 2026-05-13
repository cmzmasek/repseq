"""Segmented virus completeness filter and concatenation."""
from __future__ import annotations

from repseq.models import QCReport
from repseq.segmented.completeness import (
    build_concatenated_sequences,
    concatenate_isolate,
    extract_isolate_id,
    filter_complete_isolates,
    identify_segment,
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
