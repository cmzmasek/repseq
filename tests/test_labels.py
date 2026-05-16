"""Leaf-label formatter: placeholder substitution, fallback, separator drop.

These tests pin the contract that a single template can serve both
segmented and non-segmented runs — segmented gets the strain→isolate_id
fallback, non-segmented just resolves placeholders as-is, and empty
fields never produce ``||`` or trailing ``|`` in the rendered label.
"""
from __future__ import annotations

from repseq.models import Sequence, SequenceType, TaxonomyInfo
from repseq.phylo.labels import (
    _parse_year,
    format_leaf_label,
    labeling_options,
    pick_format_string,
)


def _make_seq(**overrides) -> Sequence:
    base = dict(
        id="ACC123",
        header="ACC123 fake header",
        sequence="ACGT",
        seq_type=SequenceType.NUCLEOTIDE,
        accession="ACC123",
        organism="Hantaan virus",
        host="Apodemus agrarius",
        strain="76-118",
        collection_date="1976-04",
        country="South Korea",
        isolate_id="76-118",
        taxonomy=TaxonomyInfo(
            taxid=1980456,
            species="Hantaan orthohantavirus",
            genus="Orthohantavirus",
            family="Hantaviridae",
        ),
    )
    base.update(overrides)
    return Sequence(**base)


def test_format_all_fields_present():
    seq = _make_seq()
    out = format_leaf_label(seq, "{species}|{id}|{host}")
    assert out == "Hantaan_orthohantavirus|ACC123|Apodemus_agrarius"


def test_format_replace_whitespace_off():
    seq = _make_seq()
    out = format_leaf_label(
        seq, "{species}|{host}", replace_whitespace=False,
    )
    assert out == "Hantaan orthohantavirus|Apodemus agrarius"


def test_empty_host_drops_separator():
    """The '|' immediately before an empty field is removed too."""
    seq = _make_seq(host=None)
    out = format_leaf_label(seq, "{species}|{id}|{host}")
    assert out == "Hantaan_orthohantavirus|ACC123"


def test_unknown_token_treated_as_empty():
    """NCBI metadata is full of 'unknown' / 'n/a' — they count as empty."""
    seq = _make_seq(host="unknown")
    out = format_leaf_label(seq, "{species}|{id}|{host}")
    assert out == "Hantaan_orthohantavirus|ACC123"


def test_strain_falls_back_to_isolate_id():
    """When {strain} is requested but absent, isolate_id substitutes."""
    seq = _make_seq(strain=None, isolate_id="ISO-99")
    out = format_leaf_label(seq, "{species}|{strain}|{host}")
    assert out == "Hantaan_orthohantavirus|ISO-99|Apodemus_agrarius"


def test_strain_absent_and_no_isolate_id_drops_separator():
    seq = _make_seq(strain=None, isolate_id=None)
    out = format_leaf_label(seq, "{species}|{strain}|{host}")
    assert out == "Hantaan_orthohantavirus|Apodemus_agrarius"


def test_keep_separator_on_empty():
    seq = _make_seq(host=None)
    out = format_leaf_label(
        seq, "{species}|{id}|{host}", keep_separator_on_empty=True,
    )
    assert out == "Hantaan_orthohantavirus|ACC123|"


def test_year_is_parsed_from_collection_date():
    seq = _make_seq(collection_date="04-Apr-1976")
    out = format_leaf_label(seq, "{species}|{year}")
    assert out == "Hantaan_orthohantavirus|1976"


def test_year_missing_drops_separator():
    seq = _make_seq(collection_date=None)
    out = format_leaf_label(seq, "{species}|{year}")
    assert out == "Hantaan_orthohantavirus"


def test_taxonomy_ranks_resolved_via_get_rank():
    seq = _make_seq()
    out = format_leaf_label(seq, "{family}|{genus}|{species}")
    assert out == "Hantaviridae|Orthohantavirus|Hantaan_orthohantavirus"


def test_no_taxonomy_drops_all_rank_placeholders():
    seq = _make_seq(taxonomy=None, organism="Hantaan virus")
    out = format_leaf_label(seq, "{species}|{id}")
    assert out == "ACC123"  # species dropped, leading separator gone


def test_unknown_placeholder_is_dropped_too():
    """An unknown field name resolves to empty, same drop policy."""
    seq = _make_seq()
    out = format_leaf_label(seq, "{species}|{nonexistent}|{id}")
    assert out == "Hantaan_orthohantavirus|ACC123"


def test_literal_text_is_preserved():
    seq = _make_seq()
    out = format_leaf_label(seq, "leaf_{id}_v1")
    assert out == "leaf_ACC123_v1"


def test_id_fallback_to_seq_id_when_no_accession():
    seq = _make_seq(accession=None)
    out = format_leaf_label(seq, "{accession}|{id}")
    # No accession → use seq.id for both
    assert out == "ACC123|ACC123"


def test_isolate_id_placeholder_in_segmented_label():
    seq = _make_seq(strain=None, isolate_id="ISO-42")
    out = format_leaf_label(seq, "{species}|{isolate_id}")
    assert out == "Hantaan_orthohantavirus|ISO-42"


def test_parse_year_finds_first_4digit_run():
    assert _parse_year("1976-04") == "1976"
    assert _parse_year("04/12/2020") == "2020"
    assert _parse_year("Apr-1976") == "1976"
    assert _parse_year("unknown") is None
    assert _parse_year(None) is None
    assert _parse_year("") is None


# ---------------------------------------------------------------------------
# pick_format_string / labeling_options helpers
# ---------------------------------------------------------------------------

def test_pick_format_string_uses_segmented_when_enabled():
    cfg = {
        "phylo": {
            "labeling": {
                "format": "{species}|{id}|{host}",
                "segmented_format": "{species}|{strain}|{host}",
            }
        }
    }
    assert pick_format_string(cfg, segmented=True) == "{species}|{strain}|{host}"
    assert pick_format_string(cfg, segmented=False) == "{species}|{id}|{host}"


def test_pick_format_string_segmented_falls_back_to_format():
    """segmented_format=None means: use the regular format."""
    cfg = {
        "phylo": {
            "labeling": {
                "format": "{species}|{id}|{host}",
                "segmented_format": None,
            }
        }
    }
    assert pick_format_string(cfg, segmented=True) == "{species}|{id}|{host}"


def test_pick_format_string_default_when_cfg_empty():
    assert pick_format_string(None, segmented=False) == "{species}|{id}|{host}"
    assert pick_format_string({}, segmented=True) == "{species}|{id}|{host}"


def test_labeling_options_defaults():
    opts = labeling_options(None)
    assert opts == {"replace_whitespace": True, "keep_separator_on_empty": False}


def test_labeling_options_overrides():
    cfg = {
        "phylo": {
            "labeling": {
                "replace_whitespace": False,
                "keep_separator_on_empty": True,
            }
        }
    }
    opts = labeling_options(cfg)
    assert opts == {"replace_whitespace": False, "keep_separator_on_empty": True}
