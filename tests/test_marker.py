"""Marker-protein selection + non-segmented protein_sequence population."""
from __future__ import annotations

from repseq.clustering.marker import (
    populate_protein_sequences,
    select_marker_protein,
)
from repseq.models import QCReport


def _cds(product: str, sequence: str, protein_id: str = "P_x"):
    return {
        "protein_id": protein_id,
        "product": product,
        "length": len(sequence),
        "sequence": sequence,
    }


def test_select_marker_returns_none_for_empty_proteins():
    assert select_marker_protein(None) is None
    assert select_marker_protein([]) is None


def test_select_marker_returns_longest_when_no_aliases():
    proteins = [
        _cds("nucleoprotein", "M" * 50),
        _cds("RNA-dependent RNA polymerase", "M" * 2200),
        _cds("glycoprotein", "M" * 500),
    ]
    chosen = select_marker_protein(proteins)
    assert chosen["product"] == "RNA-dependent RNA polymerase"


def test_select_marker_alias_order_encodes_preference():
    proteins = [
        _cds("nucleoprotein", "M" * 400),
        _cds("RNA-dependent RNA polymerase", "M" * 2200),
    ]
    # Even though polymerase is longer, alias order picks nucleoprotein first.
    chosen = select_marker_protein(proteins, ["nucleoprotein", "polymerase"])
    assert chosen["product"] == "nucleoprotein"


def test_select_marker_alias_case_insensitive_substring():
    proteins = [
        _cds("hypothetical protein", "M" * 100),
        _cds("Hemagglutinin precursor", "M" * 500),
    ]
    chosen = select_marker_protein(proteins, ["HEMAGGLUTININ"])
    assert chosen["product"] == "Hemagglutinin precursor"


def test_select_marker_falls_back_to_longest_when_no_alias_matches():
    proteins = [
        _cds("nucleoprotein", "M" * 400),
        _cds("polymerase", "M" * 2200),
    ]
    chosen = select_marker_protein(proteins, ["glycoprotein"])
    # No alias matched, fall back to longest.
    assert chosen["product"] == "polymerase"


def test_select_marker_skips_proteins_without_sequence():
    proteins = [
        {"protein_id": "p1", "product": "polymerase", "length": 2200, "sequence": None},
        _cds("nucleoprotein", "M" * 400),
    ]
    chosen = select_marker_protein(proteins)
    # The polymerase has no translation; the nucleoprotein is the only viable
    # marker, not the longest among "all CDSes".
    assert chosen["product"] == "nucleoprotein"


def test_select_marker_returns_none_if_no_cds_has_translation():
    proteins = [
        {"protein_id": "p1", "product": "polymerase", "length": 2200, "sequence": None},
    ]
    assert select_marker_protein(proteins) is None


def test_select_marker_breaks_alias_ties_by_length():
    """Two CDSes match the same alias — keep the longer one."""
    proteins = [
        _cds("polymerase short", "M" * 100),
        _cds("polymerase full", "M" * 2000),
    ]
    chosen = select_marker_protein(proteins, ["polymerase"])
    assert chosen["product"] == "polymerase full"


# ---------------------------------------------------------------------------
# populate_protein_sequences — non-segmented input
# ---------------------------------------------------------------------------

def test_populate_sets_protein_sequence_per_input(make_seq):
    a = make_seq("a", "ACGT")
    a.proteins = [_cds("polymerase", "MMMMM")]
    b = make_seq("b", "ACGT")
    b.proteins = [_cds("nucleoprotein", "NNNNN")]
    kept = populate_protein_sequences([a, b])
    assert kept == [a, b]
    assert a.protein_sequence == "MMMMM"
    assert b.protein_sequence == "NNNNN"


def test_populate_honours_alias_list(make_seq):
    a = make_seq("a", "ACGT")
    a.proteins = [
        _cds("polymerase", "M" * 2000),
        _cds("nucleoprotein", "N" * 400),
    ]
    populate_protein_sequences([a], aliases=["nucleoprotein"])
    assert a.protein_sequence == "N" * 400


def test_populate_drops_sequences_with_no_marker(make_seq):
    a = make_seq("a", "ACGT")
    a.proteins = [_cds("polymerase", "MMMMM")]
    b = make_seq("b", "ACGT")
    b.proteins = []  # fetched, none found
    c = make_seq("c", "ACGT")
    c.proteins = None  # never fetched
    report = QCReport()
    kept = populate_protein_sequences([a, b, c], report=report)
    assert kept == [a]
    assert report.removed_proteins == 2
    assert all("no_marker_protein_for_clustering" in d["reason"] for d in report.details)
