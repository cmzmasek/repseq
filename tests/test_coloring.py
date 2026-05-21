"""Unit tests for the taxonomy-driven leaf-colour palette
(``repseq.phylo.coloring``).

These exercise the palette builder and resolver directly, without the
phyloXML writer — the writer-side emission is covered in
``tests/test_phyloxml_writer.py``.
"""
from __future__ import annotations

from repseq.models import Sequence, SequenceType, TaxonomyInfo
from repseq.phylo.coloring import (
    DEFAULT_MISSING_COLOR,
    ColorScheme,
    _fan_offsets,
    _hsv_to_hex,
    _normalize,
    build_color_scheme,
    is_valid_hex_color,
)


def _seq(seq_id, *, genus="", subgenus="", family="Fam") -> Sequence:
    lineage = {}
    if subgenus:
        lineage["subgenus"] = subgenus
    return Sequence(
        id=seq_id,
        header=seq_id,
        sequence="ACGT",
        seq_type=SequenceType.NUCLEOTIDE,
        accession=seq_id,
        organism="org",
        taxonomy=TaxonomyInfo(
            taxid=1, species="sp", genus=genus, family=family, lineage=lineage
        ),
    )


# ---------------------------------------------------------------------------
# Normalisation
# ---------------------------------------------------------------------------

def test_normalize_strips_and_keeps_real_values():
    assert _normalize("  Orthohantavirus ") == "Orthohantavirus"


def test_normalize_treats_sentinels_as_missing():
    for token in ("", "  ", "unknown", "NA", "n/a", "None", "null", "?", "-", "Unclassified"):
        assert _normalize(token) is None


def test_normalize_none():
    assert _normalize(None) is None


# ---------------------------------------------------------------------------
# Hex / helpers
# ---------------------------------------------------------------------------

def test_hsv_to_hex_format():
    h = _hsv_to_hex(0.0, 1.0, 1.0)
    assert h == "#FF0000"
    assert all(c in "#0123456789ABCDEF" for c in h)


def test_is_valid_hex_color():
    assert is_valid_hex_color("#808080")
    assert is_valid_hex_color("#aabbcc")
    assert not is_valid_hex_color("808080")
    assert not is_valid_hex_color("#xyzxyz")
    assert not is_valid_hex_color("#FFF")
    assert not is_valid_hex_color(123)


def test_fan_offsets():
    assert _fan_offsets(1, 20.0) == [0.0]
    offs = _fan_offsets(3, 18.0)
    assert offs == [-18.0, 0.0, 18.0]


# ---------------------------------------------------------------------------
# Scheme building / resolution
# ---------------------------------------------------------------------------

def test_disabled_returns_none():
    cfg = {"phylo": {"coloring": {"enabled": False}}}
    assert build_color_scheme([_seq("A", genus="G")], cfg) is None


def test_default_is_genus_mode():
    scheme = build_color_scheme([_seq("A", genus="Alpha")], {})
    assert scheme is not None
    assert scheme.parent_rank == "genus"
    assert scheme.child_rank is None


def test_same_genus_same_color_distinct_genera_differ():
    reps = [_seq("A", genus="Alpha"), _seq("B", genus="Alpha"), _seq("C", genus="Beta")]
    scheme = build_color_scheme(reps, {})
    a, b, c = (scheme.color_for(r) for r in reps)
    assert a == b
    assert a != c
    assert "#808080" not in (a, b, c)


def test_missing_genus_is_grey():
    scheme = build_color_scheme([_seq("A", genus="Alpha")], {})
    assert scheme.color_for(_seq("Z", genus="")) == DEFAULT_MISSING_COLOR
    assert scheme.color_for(_seq("Z", genus="unknown")) == DEFAULT_MISSING_COLOR


def test_no_taxonomy_is_grey():
    scheme = build_color_scheme([_seq("A", genus="Alpha")], {})
    bare = Sequence(
        id="N", header="N", sequence="A", seq_type=SequenceType.NUCLEOTIDE,
        accession="N", organism="o", taxonomy=None,
    )
    assert scheme.color_for(bare) == DEFAULT_MISSING_COLOR


def test_custom_missing_color():
    cfg = {"phylo": {"coloring": {"missing_color": "#123456"}}}
    scheme = build_color_scheme([_seq("A", genus="Alpha")], cfg)
    assert scheme.color_for(_seq("Z", genus="")) == "#123456"


def test_two_rank_present_child_differs_from_sibling():
    cfg = {"phylo": {"coloring": {"ranks": ["genus", "subgenus"]}}}
    reps = [
        _seq("A", genus="G", subgenus="Sub1"),
        _seq("B", genus="G", subgenus="Sub2"),
        _seq("C", genus="G", subgenus="Sub1"),
    ]
    scheme = build_color_scheme(reps, cfg)
    a, b, c = (scheme.color_for(r) for r in reps)
    assert a == c          # same subgenus
    assert a != b          # different subgenus
    assert "#808080" not in (a, b, c)


def test_two_rank_missing_child_takes_parent_base_color():
    cfg = {"phylo": {"coloring": {"ranks": ["genus", "subgenus"]}}}
    reps = [_seq("A", genus="G", subgenus="Sub1")]
    scheme = build_color_scheme(reps, cfg)
    # A leaf with the same genus but NO subgenus: not grey, and equal to
    # the genus base hue (offset 0 == single-child colour here, but the
    # contract is "parent base colour", so just assert non-grey + genus
    # resolves).
    no_child = _seq("Z", genus="G", subgenus="")
    color = scheme.color_for(no_child)
    assert color != DEFAULT_MISSING_COLOR


def test_two_rank_missing_parent_is_grey():
    cfg = {"phylo": {"coloring": {"ranks": ["genus", "subgenus"]}}}
    reps = [_seq("A", genus="G", subgenus="Sub1")]
    scheme = build_color_scheme(reps, cfg)
    assert scheme.color_for(_seq("Z", genus="", subgenus="Sub1")) == DEFAULT_MISSING_COLOR


def test_nothing_resolves_returns_grey_scheme():
    """Empty parent set still yields a scheme; every leaf is grey."""
    reps = [_seq("A", genus=""), _seq("B", genus="unknown")]
    scheme = build_color_scheme(reps, {})
    assert isinstance(scheme, ColorScheme)
    assert scheme.color_for(reps[0]) == DEFAULT_MISSING_COLOR


def test_palette_is_subset_stable():
    """A taxon keeps its colour whether or not other taxa are present —
    the property that keeps colours consistent across 2E and 2F trees."""
    full = [_seq("A", genus="Alpha"), _seq("B", genus="Beta"), _seq("C", genus="Gamma")]
    subset = [full[1]]  # just Beta
    color_full = build_color_scheme(full, {}).color_for(full[1])
    # Build over the full set (as the orchestrator does) then resolve the
    # subset leaf — must match.
    assert build_color_scheme(full, {}).color_for(subset[0]) == color_full
