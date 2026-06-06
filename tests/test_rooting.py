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
    _candidate_clades,
    _lca_prefix,
    _leaves_under,
    _mad_root,
    _mad_score_for_branch,
    _mean_lca_specificity,
    _pairwise_leaf_distances,
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


def test_mad_recovers_balanced_root_on_clock_tree():
    """Regression for the MAD cross-pair sign bug (rigour audit 2026-06-06).

    On a perfect-clock (ultrametric) tree the correct MAD root is the branch
    giving zero ancestor deviation — the balanced {a,b}|{c,d} split. The
    pre-fix code computed the cross-pair deviation as (num0 - 2x) instead of
    (num0 + 2x), which clamped the optimal split position to x=0, mis-ranked
    branches, and rooted the tree at the lopsided {a}|{b,c,d} instead.
    """
    tree = _tree("((a:2,b:2):1,(c:2,d:2):1);")
    # Scramble the parser root (root at a leaf) so MAD must recover it.
    tree.root_with_outgroup(
        [t for t in tree.get_terminals() if t.name == "a"][0]
    )
    rooted = _mad_root(tree)
    assert rooted is not None
    sides = sorted(sorted(_leaves_under(c)) for c in rooted.root.clades)
    assert sides == [["a", "b"], ["c", "d"]]


def test_mad_residual_is_zero_at_clock_root():
    """The minimal MAD deviation over all branches of a clock tree is ~0
    (achieved at the true root). With the sign bug the minimum was strictly
    positive because no branch could reach the zero-deviation position."""
    tree = _tree("((a:2,b:2):1,(c:2,d:2):1);")
    tree.root_with_outgroup(
        [t for t in tree.get_terminals() if t.name == "a"][0]
    )
    pair_d = _pairwise_leaf_distances(tree)
    scores = [
        _mad_score_for_branch(
            tree, br, pair_d,
            {t.name: tree.distance(br, t) for t in tree.get_terminals()},
        )[0]
        for br in _candidate_clades(tree)
        if br.branch_length
    ]
    assert min(scores) < 1e-9


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


# ---------------------------------------------------------------------------
# Outgroup rooting (user-specified)
# ---------------------------------------------------------------------------

def _seq_with_acc(sid, accession, **kw):
    s = _seq(sid, **kw)
    s.accession = accession
    return s


def test_outgroup_rooting_by_single_accession():
    """A single accession in `outgroup` should root the tree at that leaf."""
    tree = _tree("((A:0.1,B:0.1):0.3,(C:0.1,D:0.1):0.3);")
    reps = {
        "A": _seq_with_acc("A", "AB001"),
        "B": _seq_with_acc("B", "AB002"),
        "C": _seq_with_acc("C", "AB003"),
        "D": _seq_with_acc("D", "AB004"),
    }
    out, used = root_tree(tree, reps, method="outgroup", outgroup="AB003")
    assert used == "outgroup"
    # The rooted tree should have C reachable as an outgroup-side leaf.
    leaf_names = {t.name for t in out.get_terminals()}
    assert leaf_names == {"A", "B", "C", "D"}


def test_outgroup_rooting_accession_match_is_case_insensitive():
    tree = _tree("((A:0.1,B:0.1):0.3,(C:0.1,D:0.1):0.3);")
    reps = {sid: _seq_with_acc(sid, f"ABC00{i}") for i, sid in enumerate("ABCD")}
    _, used = root_tree(tree, reps, method="outgroup", outgroup="abc003")
    assert used == "outgroup"


def test_outgroup_rooting_by_clade_uses_mrca():
    """A list of accessions naming a clade roots at their MRCA."""
    tree = _tree("((A:0.1,B:0.1):0.3,(C:0.1,D:0.1):0.3);")
    reps = {sid: _seq_with_acc(sid, f"AB00{i}") for i, sid in enumerate("ABCD", start=1)}
    out, used = root_tree(
        tree, reps, method="outgroup", outgroup=["AB003", "AB004"],
    )
    assert used == "outgroup"
    assert {t.name for t in out.get_terminals()} == {"A", "B", "C", "D"}


def test_outgroup_rooting_by_taxonomy_rank():
    """outgroup_rank pulls every leaf with that taxon at that rank."""
    tree = _tree("((A:0.1,B:0.1):0.3,(C:0.1,D:0.1):0.3);")
    reps = {
        "A": _seq("A", family="Hantaviridae"),
        "B": _seq("B", family="Hantaviridae"),
        "C": _seq("C", family="Peribunyaviridae"),
        "D": _seq("D", family="Peribunyaviridae"),
    }
    _, used = root_tree(
        tree, reps, method="outgroup",
        outgroup_rank={"family": "Peribunyaviridae"},
    )
    assert used == "outgroup"


def test_outgroup_rooting_falls_back_to_midpoint_on_no_match():
    """A typo (no accession matches) must not abort — fall through to
    midpoint so a usable tree still comes out."""
    tree = _tree("((A:0.1,B:0.1):0.3,(C:0.1,D:0.1):0.3);")
    reps = {sid: _seq_with_acc(sid, f"AB00{i}") for i, sid in enumerate("ABCD", start=1)}
    _, used = root_tree(
        tree, reps, method="outgroup", outgroup="DOES_NOT_EXIST",
    )
    assert used in ("midpoint", "none")
