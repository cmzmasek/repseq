"""Internal-node LCA annotation, label cleanup, suffix-rank inference,
same-species suppression."""
from __future__ import annotations

from io import StringIO

from Bio import Phylo

from repseq.models import Sequence, SequenceType, TaxonomyInfo
from repseq.phylo.lca import (
    _infer_rank_from_name,
    _lca_prefix,
    _reaches_min_rank,
    annotate_internal_nodes,
    keep_deepest_labels,
    phyloxml_rank,
    suppress_same_species_pairs,
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


# ---------------------------------------------------------------------------
# Suffix-rank inference (ICTV convention)
# ---------------------------------------------------------------------------

def test_infer_rank_from_ictv_suffix():
    assert _infer_rank_from_name("Hantaviridae") == "family"
    assert _infer_rank_from_name("Bunyavirales") == "order"
    assert _infer_rank_from_name("Negarnaviricota") == "phylum"
    assert _infer_rank_from_name("Ellioviricetes") == "class"
    assert _infer_rank_from_name("Orthohantavirinae") == "subfamily"


def test_infer_rank_single_word_virus_is_genus():
    assert _infer_rank_from_name("Orthohantavirus") == "genus"


def test_infer_rank_multi_word_virus_no_inference():
    """Multi-word names like 'Hantaan virus' don't follow the
    single-word genus convention."""
    assert _infer_rank_from_name("Hantaan virus") == "no rank"


# ---------------------------------------------------------------------------
# min_rank gate
# ---------------------------------------------------------------------------

def test_reaches_min_rank_true_when_species_present():
    lineage = [
        ("family", "Hantaviridae"),
        ("genus", "Orthohantavirus"),
        ("species", "Hantaan virus"),
    ]
    assert _reaches_min_rank(lineage, "species")


def test_reaches_min_rank_false_when_too_coarse():
    lineage = [("family", "Hantaviridae")]
    assert not _reaches_min_rank(lineage, "genus")


def test_reaches_min_rank_none_always_true():
    assert _reaches_min_rank([], "none")
    assert _reaches_min_rank([("family", "X")], "none")


# ---------------------------------------------------------------------------
# LCA prefix
# ---------------------------------------------------------------------------

def test_lca_prefix_truncates_at_first_divergence():
    a = [("family", "X"), ("genus", "A")]
    b = [("family", "X"), ("genus", "B")]
    assert _lca_prefix([a, b]) == [("family", "X")]


# ---------------------------------------------------------------------------
# annotate_internal_nodes
# ---------------------------------------------------------------------------

def test_annotate_labels_internal_with_family_when_genera_differ():
    """Two sister leaves with different genera in the same family →
    internal node gets labelled with the family."""
    tree = _tree("((A:0.1,B:0.1):0.3,C:0.5);")
    reps = {
        "A": _seq("A", family="Hantaviridae", genus="Orthohantavirus",
                  species="Hantaan", lineage={
                      "family": "Hantaviridae",
                      "genus": "Orthohantavirus",
                      "species": "Hantaan",
                  }),
        "B": _seq("B", family="Hantaviridae", genus="Loanvirus",
                  species="Loanvirus brnoense", lineage={
                      "family": "Hantaviridae",
                      "genus": "Loanvirus",
                      "species": "Loanvirus brnoense",
                  }),
        "C": _seq("C", family="Peribunyaviridae", genus="Orthobunyavirus",
                  species="La Crosse", lineage={
                      "family": "Peribunyaviridae",
                      "genus": "Orthobunyavirus",
                      "species": "La Crosse",
                  }),
    }
    annotate_internal_nodes(tree, reps, min_rank="genus")
    # The A/B sister clade should be labelled "Hantaviridae".
    sister = [
        n for n in tree.get_nonterminals()
        if {t.name for t in n.get_terminals()} == {"A", "B"}
    ][0]
    assert getattr(sister, "_lca_name", None) == "Hantaviridae"
    assert getattr(sister, "_lca_rank", None) == "family"


def test_annotate_skips_when_coverage_below_threshold():
    """A clade where most leaves lack lineage shouldn't be labelled."""
    tree = _tree("((A:0.1,B:0.1):0.3,(C:0.1,D:0.1):0.3);")
    reps = {
        "A": _seq("A", family="Hantaviridae", genus="Orthohantavirus",
                  species="Hantaan", lineage={
                      "family": "Hantaviridae",
                      "genus": "Orthohantavirus",
                      "species": "Hantaan",
                  }),
        "B": _seq("B"),   # no taxonomy
        "C": _seq("C"),
        "D": _seq("D"),
    }
    annotate_internal_nodes(
        tree, reps, min_rank="genus", coverage_threshold=0.5,
    )
    # The (C,D) clade has 0% lineage coverage → no label.
    cd = [
        n for n in tree.get_nonterminals()
        if {t.name for t in n.get_terminals()} == {"C", "D"}
    ][0]
    assert getattr(cd, "_lca_name", None) is None
    # The root has 1/4 = 25% < 50% → no label either.
    assert getattr(tree.root, "_lca_name", None) is None


def test_annotate_filters_leaves_below_min_rank():
    """A leaf whose lineage stops at family is excluded from the vote
    when min_rank='genus' — but still stays in the tree."""
    tree = _tree("((A:0.1,B:0.1):0.3,C:0.5);")
    reps = {
        "A": _seq("A", family="Hantaviridae", genus="Orthohantavirus",
                  species="Hantaan", lineage={
                      "family": "Hantaviridae",
                      "genus": "Orthohantavirus",
                      "species": "Hantaan",
                  }),
        "B": _seq("B", family="Hantaviridae", genus="Orthohantavirus",
                  species="Sin Nombre", lineage={
                      "family": "Hantaviridae",
                      "genus": "Orthohantavirus",
                      "species": "Sin Nombre",
                  }),
        "C": _seq("C", family="Peribunyaviridae", lineage={
            "family": "Peribunyaviridae",
        }),  # no genus → excluded from vote with min_rank=genus
    }
    annotate_internal_nodes(tree, reps, min_rank="genus")
    # The A/B sister clade: both vote, both Orthohantavirus → labelled.
    ab = [
        n for n in tree.get_nonterminals()
        if {t.name for t in n.get_terminals()} == {"A", "B"}
    ][0]
    assert getattr(ab, "_lca_name", None) == "Orthohantavirus"
    # Tree still has all three leaves.
    assert {t.name for t in tree.get_terminals()} == {"A", "B", "C"}


# ---------------------------------------------------------------------------
# keep_deepest_labels
# ---------------------------------------------------------------------------

def test_keep_deepest_clears_nested_duplicate_labels():
    tree = _tree("(((A:0.1,B:0.1):0.1,C:0.2):0.1,D:0.3);")
    # Set the same LCA name on multiple nested internals.
    for n in tree.get_nonterminals():
        n._lca_name = "Hantaviridae"
        n._lca_rank = "family"
    keep_deepest_labels(tree)
    labeled = [
        n for n in tree.get_nonterminals()
        if getattr(n, "_lca_name", None) == "Hantaviridae"
    ]
    # Only the largest (root) clade should keep the label.
    assert len(labeled) == 1
    assert len(labeled[0].get_terminals()) == 4


def test_keep_deepest_preserves_distinct_labels():
    tree = _tree("(((A:0.1,B:0.1):0.1,C:0.2):0.1,D:0.3);")
    ints = tree.get_nonterminals()
    ints[0]._lca_name = "Hantaviridae"; ints[0]._lca_rank = "family"
    ints[1]._lca_name = "Orthohantavirus"; ints[1]._lca_rank = "genus"
    ints[2]._lca_name = "Hantaan virus"; ints[2]._lca_rank = "species"
    keep_deepest_labels(tree)
    # Three distinct names → all three survive.
    names = {
        getattr(n, "_lca_name", None) for n in tree.get_nonterminals()
    }
    assert names == {"Hantaviridae", "Orthohantavirus", "Hantaan virus"}


# ---------------------------------------------------------------------------
# Same-species suppression
# ---------------------------------------------------------------------------

def test_suppress_same_species_pair_on_two_leaf_internal():
    tree = _tree("((A:0.1,B:0.1):0.3,C:0.5);")
    reps = {
        "A": _seq("A", species="Hantaan virus"),
        "B": _seq("B", species="Hantaan virus"),
        "C": _seq("C", species="Other virus"),
    }
    # Pretend the LCA annotator labelled the (A,B) internal.
    ab = [
        n for n in tree.get_nonterminals()
        if {t.name for t in n.get_terminals()} == {"A", "B"}
    ][0]
    ab._lca_name = "Hantaan virus"
    ab._lca_rank = "species"
    suppress_same_species_pairs(tree, reps)
    assert getattr(ab, "_lca_name", None) is None


def test_does_not_suppress_when_species_differ():
    tree = _tree("((A:0.1,B:0.1):0.3,C:0.5);")
    reps = {
        "A": _seq("A", species="Hantaan virus"),
        "B": _seq("B", species="Sin Nombre virus"),
        "C": _seq("C", species="Other virus"),
    }
    ab = [
        n for n in tree.get_nonterminals()
        if {t.name for t in n.get_terminals()} == {"A", "B"}
    ][0]
    ab._lca_name = "Orthohantavirus"
    suppress_same_species_pairs(tree, reps)
    assert getattr(ab, "_lca_name", None) == "Orthohantavirus"


def test_does_not_suppress_when_more_than_two_children():
    tree = _tree("((A:0.1,B:0.1,X:0.1):0.3,C:0.5);")
    reps = {
        "A": _seq("A", species="Hantaan"),
        "B": _seq("B", species="Hantaan"),
        "X": _seq("X", species="Hantaan"),
        "C": _seq("C", species="Other"),
    }
    abx = [
        n for n in tree.get_nonterminals()
        if {t.name for t in n.get_terminals()} == {"A", "B", "X"}
    ][0]
    abx._lca_name = "Hantaan"
    suppress_same_species_pairs(tree, reps)
    # Three children, even if all same species → keep label (we only
    # suppress the 2-leaf case where the species name is unambiguously
    # already on both leaves).
    assert getattr(abx, "_lca_name", None) == "Hantaan"


# ---------------------------------------------------------------------------
# PhyloXML rank validation
# ---------------------------------------------------------------------------

def test_phyloxml_rank_accepts_standard():
    assert phyloxml_rank("family") == "family"
    assert phyloxml_rank("species") == "species"
    assert phyloxml_rank("genus") == "genus"


def test_phyloxml_rank_falls_back_to_other_for_unknown():
    # NCBI "no rank" / "clade" / weird custom ranks aren't in the
    # PhyloXML schema enumeration; rank must round-trip as "other".
    assert phyloxml_rank("no rank") == "other"
    assert phyloxml_rank("clade") == "other"
    assert phyloxml_rank("subrealm") == "other"  # not in PhyloXML enum
    assert phyloxml_rank(None) == "other"
    assert phyloxml_rank("") == "other"
