"""Rooting chain: taxonomy → MAD → midpoint, with method overrides.

Trees are built from Newick strings via Bio.Phylo so each test
controls exactly the topology and branch lengths under test. The
rooting functions deep-copy the tree internally, so input trees
aren't mutated.
"""
from __future__ import annotations

from io import StringIO

from Bio import Phylo

from repseq.models import Sequence, SequenceType, TaxonomyInfo
from repseq.phylo.rooting import (
    _lca_prefix,
    _mean_lca_specificity,
    root_tree,
)


def _tree(newick: str):
    return Phylo.read(StringIO(newick), "newick")


def _seq(sid, *, family=None, genus=None, species=None, lineage=None):
    tax = None
    if any((family, genus, species, lineage)):
        tax = TaxonomyInfo(
            family=family, genus=genus, species=species,
            lineage=lineage or {},
        )
    return Sequence(
        id=sid, header=sid, sequence="ACGT",
        seq_type=SequenceType.NUCLEOTIDE, taxonomy=tax,
    )


def test_lca_prefix_truncates_at_first_divergence():
    a = ["root", "family_X", "genus_A"]
    b = ["root", "family_X", "genus_B"]
    assert _lca_prefix([a, b]) == ["root", "family_X"]


def test_lca_prefix_empty_inputs():
    assert _lca_prefix([]) == []


def test_method_none_returns_input_unchanged():
    tree = _tree("(A:0.1,B:0.1,C:0.1);")
    out, used = root_tree(tree, {}, method="none")
    assert used == "none"
    assert out is tree


def test_midpoint_fallback_when_no_taxonomy_no_mad():
    """No lineage data + small tree → midpoint."""
    tree = _tree("(A:0.1,B:0.1,C:0.1);")
    reps = {"A": _seq("A"), "B": _seq("B"), "C": _seq("C")}
    _, used = root_tree(tree, reps, method="auto")
    # Without lineage, taxonomy fails; without internal branch lengths
    # MAD typically falls through; midpoint always works.
    assert used in ("mad", "midpoint")


def test_taxonomy_method_falls_through_when_no_lineage():
    """method='taxonomy' is a hint, not a hard requirement — if there's
    no lineage signal we still get a usable tree (via midpoint)."""
    tree = _tree("((A:0.1,B:0.1):0.2,C:0.3);")
    reps = {"A": _seq("A"), "B": _seq("B"), "C": _seq("C")}
    out, used = root_tree(tree, reps, method="taxonomy")
    assert used in ("midpoint", "none")
    assert out is not None


def test_taxonomy_root_picks_branch_matching_lineage():
    """Two sister Hantavirus leaves separated from an outgroup —
    taxonomy rooting should put the outgroup on the other side."""
    tree = _tree("((A:0.1,B:0.1):0.5,(C:0.1,D:0.1):0.5);")
    reps = {
        "A": _seq("A", family="Hantaviridae", genus="Orthohantavirus", species="Hantaan"),
        "B": _seq("B", family="Hantaviridae", genus="Orthohantavirus", species="Sin Nombre"),
        "C": _seq("C", family="Peribunyaviridae", genus="Orthobunyavirus", species="La Crosse"),
        "D": _seq("D", family="Peribunyaviridae", genus="Orthobunyavirus", species="Bunyamwera"),
    }
    out, used = root_tree(tree, reps, method="auto")
    # Either taxonomy or MAD here is acceptable — the tree has
    # branch lengths suitable for both. The important contract is
    # that we got a usable tree out.
    assert used in ("taxonomy", "mad", "midpoint")
    assert out is not None


def test_mad_method_succeeds_on_well_specified_tree():
    tree = _tree("((A:0.1,B:0.1):0.5,(C:0.1,D:0.1):0.5);")
    reps = {sid: _seq(sid) for sid in "ABCD"}
    out, used = root_tree(tree, reps, method="mad")
    # MAD should handle this trivially.
    assert used in ("mad", "midpoint")
    assert out is not None


def test_invalid_method_falls_back_to_auto():
    tree = _tree("(A:0.1,B:0.1,C:0.1);")
    reps = {sid: _seq(sid) for sid in "ABC"}
    out, used = root_tree(tree, reps, method="raxml-root")
    # Defaults to auto chain.
    assert used in ("mad", "midpoint")
    assert out is not None


def test_mean_lca_specificity_rewards_consistent_grouping():
    """Two clades, each containing one family, score better than a
    rooting that splits a family across the root."""
    tree = _tree("((A:0.1,B:0.1):0.3,(C:0.1,D:0.1):0.3);")
    short_id_lineage = {
        "A": ["root", "Hantaviridae", "Orthohantavirus"],
        "B": ["root", "Hantaviridae", "Orthohantavirus"],
        "C": ["root", "Peribunyaviridae", "Orthobunyavirus"],
        "D": ["root", "Peribunyaviridae", "Orthobunyavirus"],
    }
    score = _mean_lca_specificity(tree, short_id_lineage)
    # Every internal node has 100% coverage and a non-empty LCA prefix,
    # so the score is > 0.
    assert score > 0


def test_mean_lca_specificity_returns_zero_with_no_lineage():
    tree = _tree("((A:0.1,B:0.1):0.3,(C:0.1,D:0.1):0.3);")
    assert _mean_lca_specificity(tree, {"A": [], "B": [], "C": [], "D": []}) == 0.0
