"""Tests for repseq.phylo.basis.describe_tree_basis — the plain-English
"what is this tree based on" descriptions that make every phyloXML
self-describing."""
from __future__ import annotations

from repseq.phylo.basis import ROLES, describe_tree_basis


def _props(role, **kw):
    _, props = describe_tree_basis(role, **kw)
    return props


def test_every_role_returns_sentence_and_five_properties():
    for role in ROLES:
        sentence, props = describe_tree_basis(
            role, alphabet="protein", segmented=False,
        )
        assert sentence and sentence.endswith(".")
        assert set(props) == {
            "tree_basis", "analysis_mode", "substrate", "alphabet", "leaf_unit",
        }
        assert props["tree_basis"] == sentence


def test_genome_segmented_protein_is_segment_marker_concat():
    sentence, props = describe_tree_basis(
        "genome", alphabet="protein", segmented=True,
        markers="L:RdRp, M:GPC, S:N",
    )
    assert props["substrate"] == "segment_marker_concat"
    assert props["analysis_mode"] == "segmented"
    assert props["leaf_unit"] == "isolate"
    assert "one marker protein per segment" in sentence
    assert "L:RdRp, M:GPC, S:N" in sentence
    assert "representative isolate" in sentence


def test_genome_segmented_nucleotide_is_segment_nt_concat():
    _, props = describe_tree_basis(
        "genome", alphabet="nucleotide", segmented=True,
    )
    assert props["substrate"] == "segment_nt_concat"
    assert props["alphabet"] == "nucleotide"


def test_genome_nonseg_single_marker_vs_concat():
    single = _props(
        "genome", alphabet="protein", segmented=False,
        markers="Spike", concat_markers=False,
    )
    assert single["substrate"] == "single_marker"
    concat = _props(
        "genome", alphabet="protein", segmented=False,
        markers="Spike, Nucleocapsid", concat_markers=True,
    )
    assert concat["substrate"] == "marker_protein_concat"


def test_genome_nonseg_nucleotide_is_genome_nt():
    sentence, props = describe_tree_basis(
        "genome", alphabet="nucleotide", segmented=False,
    )
    assert props["substrate"] == "genome_nt"
    assert "whole-genome nucleotide" in sentence
    assert props["leaf_unit"] == "sequence"


def test_partitioned_supermatrix_lists_families():
    sentence, props = describe_tree_basis(
        "genome_partitioned", alphabet="protein", segmented=False,
        families=["Spike", "Nucleocapsid"],
    )
    assert props["substrate"] == "supermatrix"
    assert "Spike, Nucleocapsid" in sentence
    assert "partitioned supermatrix" in sentence


def test_marker_tree_names_family_and_architecture():
    sentence, props = describe_tree_basis(
        "marker", alphabet="protein", segmented=False,
        family="Spike", architecture="CoV_S1--CoV_S2 OR bCoV_S1_N--CoV_S2",
    )
    assert props["substrate"] == "marker"
    assert "Spike marker protein only" in sentence
    assert "HMM architecture CoV_S1--CoV_S2 OR bCoV_S1_N--CoV_S2" in sentence


def test_segment_nt_tree_names_segment():
    sentence, props = describe_tree_basis(
        "segment_nt", alphabet="nucleotide", segmented=True, segment="L",
    )
    assert props["substrate"] == "segment_nt"
    assert "L segment nucleotide sequence only" in sentence


def test_extra_protein_tree_flags_accessory():
    sentence, props = describe_tree_basis(
        "extra_protein", alphabet="protein", segmented=False, family="ORF7",
    )
    assert props["substrate"] == "accessory_protein"
    assert "accessory protein ORF7" in sentence
    assert "do not drive clustering" in sentence


def test_peptide_tree_names_parent_polyprotein():
    sentence, props = describe_tree_basis(
        "peptide", alphabet="protein", segmented=False,
        family="NSP12", parent="ORF1ab",
    )
    assert props["substrate"] == "peptide"
    assert "NSP12 mature peptide" in sentence
    assert "sliced from the ORF1ab polyprotein" in sentence


def test_pre_cluster_overview_marks_reps():
    sentence, props = describe_tree_basis(
        "pre_cluster", alphabet="nucleotide", segmented=True,
    )
    assert props["substrate"] == "overview"
    assert props["leaf_unit"] == "isolate"
    assert "overview of every post-QC isolate" in sentence
    assert "[repr]" in sentence
