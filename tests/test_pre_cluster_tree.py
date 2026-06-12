"""Pre-cluster overview tree (2H) — helpers + end-to-end integration.

The end-to-end test mocks MAFFT and FastTree at their wrapper level so
we exercise the orchestration (id map + label-prefix dict + temp MSA
handling + Newick parse + midpoint root + phyloXML write) without
needing real binaries.
"""
from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from repseq.phylo.pre_cluster import (
    _build_id_map_pre,
    _use_protein_sequence,
    _write_id_map_with_rep_flag,
    run_pre_cluster_phylogeny,
)
from repseq.phylo.pipeline import PhyloError


# ---------------------------------------------------------------------------
# _use_protein_sequence
# ---------------------------------------------------------------------------

def test_use_protein_when_alphabet_protein_and_all_populated(make_seq):
    seqs = [make_seq("a", "ACGT"), make_seq("b", "ACGT")]
    for s in seqs:
        s.protein_sequence = "MKL"
    cfg = {"clustering": {"alphabet_for_clustering": "protein"}}
    assert _use_protein_sequence(seqs, cfg) is True


def test_use_protein_false_when_any_protein_sequence_missing(make_seq):
    """Half-protein / half-NT alignment is never acceptable — drop to
    NT if any sequence lacks the protein body."""
    seqs = [make_seq("a", "ACGT"), make_seq("b", "ACGT")]
    seqs[0].protein_sequence = "MKL"
    # seqs[1].protein_sequence stays None
    cfg = {"clustering": {"alphabet_for_clustering": "protein"}}
    assert _use_protein_sequence(seqs, cfg) is False


def test_use_protein_false_when_alphabet_nucleotide(make_seq):
    seqs = [make_seq("a", "ACGT")]
    seqs[0].protein_sequence = "MKL"
    cfg = {"clustering": {"alphabet_for_clustering": "nucleotide"}}
    assert _use_protein_sequence(seqs, cfg) is False


# ---------------------------------------------------------------------------
# _build_id_map_pre — order + count
# ---------------------------------------------------------------------------

def test_build_id_map_pre_preserves_input_order(make_seq):
    seqs = [make_seq(f"seq{i}", "ACGT") for i in range(4)]
    id_map = _build_id_map_pre(seqs)
    short_ids = list(id_map.keys())
    assert short_ids == ["S0001", "S0002", "S0003", "S0004"]
    assert [id_map[k] for k in short_ids] == ["seq0", "seq1", "seq2", "seq3"]


# ---------------------------------------------------------------------------
# _write_id_map_with_rep_flag — three-column TSV
# ---------------------------------------------------------------------------

def test_id_map_tsv_has_three_columns_and_is_rep_flag(tmp_path):
    id_map = {"S0001": "rep_a", "S0002": "non_rep_b", "S0003": "rep_c"}
    rep_ids = {"rep_a", "rep_c"}
    out = tmp_path / "id_map.tsv"
    _write_id_map_with_rep_flag(id_map, rep_ids, out)
    text = out.read_text()
    lines = text.strip().split("\n")
    assert lines[0] == "short_id\taccession\tis_rep"
    assert "S0001\trep_a\tTRUE" in text
    assert "S0002\tnon_rep_b\tFALSE" in text
    assert "S0003\trep_c\tTRUE" in text


# ---------------------------------------------------------------------------
# label_prefix_by_id on write_phyloxml — verify the additive parameter
# ---------------------------------------------------------------------------

def test_phyloxml_label_prefix_appears_on_marked_leaves(tmp_path, make_seq):
    """Direct test of the v0.32.0 label_prefix_by_id parameter on
    write_phyloxml: when a leaf's seq.id is in the prefix dict, the
    prefix gets prepended to its phyloXML <name>. Other leaves are
    unaffected."""
    pytest.importorskip("Bio")
    from repseq.phylo.phyloxml_writer import write_phyloxml

    # Tiny tree: 3 leaves, two are "reps" (carry the prefix), one is not.
    reps = [
        make_seq("rep1", "ACGT", organism="Virus A"),
        make_seq("seq2", "ACGT", organism="Virus B"),
        make_seq("rep3", "ACGT", organism="Virus C"),
    ]
    id_map = {"S0001": "rep1", "S0002": "seq2", "S0003": "rep3"}
    nwk = tmp_path / "tree.nwk"
    nwk.write_text("(S0001:0.1,S0002:0.1,S0003:0.1);\n")

    out = tmp_path / "tree.xml"
    write_phyloxml(
        nwk, out, reps, id_map,
        cfg={}, prefix="x", alphabet="nucleotide",
        msa_tool="MAFFT", msa_version="?",
        tree_tool="FastTree", tree_version="?",
        model="GTR", ufboot=None,
        label_prefix_by_id={"rep1": "[repr] ", "rep3": "[repr] "},
    )
    xml = out.read_text()
    # rep1 and rep3 leaves carry the prefix in their <name>; seq2 does not.
    assert "[repr] " in xml
    # Make sure there are at least 2 prefixed names (one per rep).
    assert xml.count("[repr] ") >= 2


def test_phyloxml_representative_ids_emits_boolean_property(tmp_path, make_seq):
    """write_phyloxml(representative_ids=...) emits a
    repseq:is_representative xsd:boolean <property> on EVERY leaf —
    'true' for ids in the set, 'false' otherwise."""
    import re

    pytest.importorskip("Bio")
    from repseq.phylo.phyloxml_writer import write_phyloxml

    reps = [
        make_seq("rep1", "ACGT", organism="Virus A"),
        make_seq("seq2", "ACGT", organism="Virus B"),
        make_seq("rep3", "ACGT", organism="Virus C"),
    ]
    id_map = {"S0001": "rep1", "S0002": "seq2", "S0003": "rep3"}
    nwk = tmp_path / "tree.nwk"
    nwk.write_text("(S0001:0.1,S0002:0.1,S0003:0.1);\n")

    out = tmp_path / "tree.xml"
    write_phyloxml(
        nwk, out, reps, id_map,
        cfg={}, prefix="x", alphabet="nucleotide",
        msa_tool="MAFFT", msa_version="?",
        tree_tool="FastTree", tree_version="?",
        model="GTR", ufboot=None,
        representative_ids={"rep1", "rep3"},
    )
    xml = out.read_text()
    assert 'datatype="xsd:boolean"' in xml
    # One property per leaf (3), values match membership.
    vals = re.findall(r'ref="repseq:is_representative"[^>]*>(\w+)<', xml)
    assert sorted(vals) == ["false", "true", "true"]


def test_phyloxml_no_representative_property_when_ids_none(tmp_path, make_seq):
    """The default (representative_ids=None) — every other tree — emits
    NO repseq:is_representative property at all."""
    pytest.importorskip("Bio")
    from repseq.phylo.phyloxml_writer import write_phyloxml

    reps = [
        make_seq("a", "ACGT", organism="Virus A"),
        make_seq("b", "ACGT", organism="Virus B"),
        make_seq("c", "ACGT", organism="Virus C"),
    ]
    id_map = {"S0001": "a", "S0002": "b", "S0003": "c"}
    nwk = tmp_path / "tree.nwk"
    nwk.write_text("(S0001:0.1,S0002:0.1,S0003:0.1);\n")

    out = tmp_path / "tree.xml"
    write_phyloxml(
        nwk, out, reps, id_map,
        cfg={}, prefix="x", alphabet="nucleotide",
        msa_tool="MAFFT", msa_version="?",
        tree_tool="FastTree", tree_version="?",
        model="GTR", ufboot=None,
    )
    assert "repseq:is_representative" not in out.read_text()


# ---------------------------------------------------------------------------
# run_pre_cluster_phylogeny — end-to-end (MAFFT/FastTree mocked)
# ---------------------------------------------------------------------------

def _stub_mafft(in_fa, out_fa, cfg, extra_args=None, use_auto=True):
    """Drop the input into the output as a stub MSA."""
    out_fa.parent.mkdir(parents=True, exist_ok=True)
    out_fa.write_text(in_fa.read_text())


def _stub_fasttree(msa_fa, out_nwk, cfg, is_protein):
    """Emit a simple star Newick whose leaf labels = the MSA's first-token ids.

    A real FastTree output is more nuanced, but the orchestrator just
    parses the Newick + roots it; a star topology is enough to verify
    the wiring end-to-end."""
    short_ids = []
    for line in msa_fa.read_text().splitlines():
        if line.startswith(">"):
            short_ids.append(line[1:].split()[0])
    leaves = ",".join(f"{sid}:0.1" for sid in short_ids)
    out_nwk.parent.mkdir(parents=True, exist_ok=True)
    out_nwk.write_text(f"({leaves});\n")


def test_run_pre_cluster_full_pipeline_writes_three_files(tmp_path, make_seq):
    pytest.importorskip("Bio")
    seqs = [
        make_seq("rep1", "ACGTACGT", organism="Virus A"),
        make_seq("seq2", "ACGTAAGT", organism="Virus A"),
        make_seq("rep3", "ACGTAAGT", organism="Virus B"),
        make_seq("seq4", "ACGTAAGT", organism="Virus B"),
    ]
    reps = [seqs[0], seqs[2]]
    cfg = {
        "clustering": {"alphabet_for_clustering": "nucleotide"},
        "segmented": {"enabled": False},
        "output": {"dir": str(tmp_path), "prefix": "test"},
        "temp_dir": str(tmp_path / "tmp"),
    }
    with patch("repseq.phylo.pre_cluster.run_mafft", side_effect=_stub_mafft), \
         patch("repseq.phylo.pre_cluster.run_fasttree", side_effect=_stub_fasttree):
        files = run_pre_cluster_phylogeny(seqs, reps, cfg, tmp_path, "test")
    written = {p.name for p in files}
    assert "test_pre_cluster_tree.nwk" in written
    assert "test_pre_cluster_tree.xml" in written
    assert "test_pre_cluster_tree_id_map.tsv" in written

    # id_map_tsv carries is_rep flags for the elected reps.
    id_map_txt = (tmp_path / "test_pre_cluster_tree_id_map.tsv").read_text()
    assert "rep1\tTRUE" in id_map_txt
    assert "seq2\tFALSE" in id_map_txt
    assert "rep3\tTRUE" in id_map_txt
    assert "seq4\tFALSE" in id_map_txt

    # phyloXML carries the [repr] prefix on the two rep leaves.
    xml = (tmp_path / "test_pre_cluster_tree.xml").read_text()
    assert "[repr] " in xml
    assert xml.count("[repr] ") >= 2

    # ...and a machine-readable repseq:is_representative boolean on every
    # leaf (4 leaves: 2 reps → true, 2 non-reps → false).
    import re
    vals = re.findall(r'ref="repseq:is_representative"[^>]*>(\w+)<', xml)
    assert sorted(vals) == ["false", "false", "true", "true"]


def test_run_pre_cluster_raises_under_3_sequences(tmp_path, make_seq):
    seqs = [make_seq("a", "ACGT"), make_seq("b", "ACGT")]  # only 2
    cfg = {"output": {"dir": str(tmp_path), "prefix": "test"}}
    with pytest.raises(PhyloError, match=">= 3 sequences"):
        run_pre_cluster_phylogeny(seqs, [], cfg, tmp_path, "test")


def test_run_pre_cluster_uses_protein_when_alphabet_protein(tmp_path, make_seq):
    """When clustering ran on protein and every sequence carries a
    protein_sequence, the pre-cluster tree is built on AA. Verified by
    inspecting the FASTA the wrapped MAFFT stub receives."""
    pytest.importorskip("Bio")
    seqs = [make_seq(f"s{i}", "ACGTACGT") for i in range(3)]
    for s in seqs:
        s.protein_sequence = "MKLPQE"
    cfg = {
        "clustering": {"alphabet_for_clustering": "protein"},
        "segmented": {"enabled": False},
        "output": {"dir": str(tmp_path), "prefix": "test"},
        "temp_dir": str(tmp_path / "tmp"),
    }
    captured = {}
    def _capture_mafft(in_fa, out_fa, cfg, extra_args=None, use_auto=True):
        captured["fasta"] = in_fa.read_text()
        out_fa.parent.mkdir(parents=True, exist_ok=True)
        out_fa.write_text(in_fa.read_text())
    with patch("repseq.phylo.pre_cluster.run_mafft", side_effect=_capture_mafft), \
         patch("repseq.phylo.pre_cluster.run_fasttree", side_effect=_stub_fasttree):
        run_pre_cluster_phylogeny(seqs, [], cfg, tmp_path, "test")
    # Body should be the protein, not the NT.
    assert "MKLPQE" in captured["fasta"]
    assert "ACGTACGT" not in captured["fasta"]


def test_run_pre_cluster_passes_retree1_to_mafft(tmp_path, make_seq):
    """MAFFT must always be called with --retree 1 + use_auto=False,
    regardless of phylo.mafft.* config."""
    pytest.importorskip("Bio")
    seqs = [make_seq(f"s{i}", "ACGTACGT") for i in range(3)]
    cfg = {
        "clustering": {"alphabet_for_clustering": "nucleotide"},
        "segmented": {"enabled": False},
        # User has set L-INS-i for the post-cluster tree — pre-cluster
        # must ignore this.
        "phylo": {"mafft": {"extra_args": ["--maxiterate", "1000", "--localpair"]}},
        "output": {"dir": str(tmp_path), "prefix": "test"},
        "temp_dir": str(tmp_path / "tmp"),
    }
    captured = {}
    def _capture_mafft(in_fa, out_fa, cfg, extra_args=None, use_auto=True):
        captured["extra_args"] = extra_args
        captured["use_auto"] = use_auto
        out_fa.parent.mkdir(parents=True, exist_ok=True)
        out_fa.write_text(in_fa.read_text())
    with patch("repseq.phylo.pre_cluster.run_mafft", side_effect=_capture_mafft), \
         patch("repseq.phylo.pre_cluster.run_fasttree", side_effect=_stub_fasttree):
        run_pre_cluster_phylogeny(seqs, [], cfg, tmp_path, "test")
    assert captured["extra_args"] == ["--retree", "1"]
    assert captured["use_auto"] is False


def _capture_mafft_factory(captured):
    def _capture(in_fa, out_fa, cfg, extra_args=None, use_auto=True):
        captured["extra_args"] = extra_args
        out_fa.parent.mkdir(parents=True, exist_ok=True)
        out_fa.write_text(in_fa.read_text())
    return _capture


def test_run_pre_cluster_switches_to_parttree_above_threshold(tmp_path, make_seq):
    """At/above phylo.pre_cluster_tree.parttree_threshold MAFFT must switch
    from --retree 1 to the PartTree guide (--retree 2 --parttree) so a huge
    pool builds instead of OOM-ing on the O(N^2) distance matrix."""
    pytest.importorskip("Bio")
    seqs = [make_seq(f"s{i}", "ACGTACGT") for i in range(4)]
    cfg = {
        "clustering": {"alphabet_for_clustering": "nucleotide"},
        "segmented": {"enabled": False},
        # 4 sequences >= threshold 3 -> PartTree path.
        "phylo": {"pre_cluster_tree": {"parttree_threshold": 3}},
        "output": {"dir": str(tmp_path), "prefix": "test"},
        "temp_dir": str(tmp_path / "tmp"),
    }
    captured = {}
    with patch(
        "repseq.phylo.pre_cluster.run_mafft",
        side_effect=_capture_mafft_factory(captured),
    ), patch(
        "repseq.phylo.pre_cluster.run_fasttree", side_effect=_stub_fasttree,
    ):
        run_pre_cluster_phylogeny(seqs, [], cfg, tmp_path, "test")
    assert captured["extra_args"] == ["--retree", "2", "--parttree"]
    # The phyloXML <description> records the strategy that actually ran.
    xml = (tmp_path / "test_pre_cluster_tree.xml").read_text()
    assert "--parttree" in xml


def test_run_pre_cluster_keeps_retree1_below_threshold(tmp_path, make_seq):
    """Below the threshold the standard --retree 1 pass is used (default
    threshold 10000 — a small pool never triggers PartTree)."""
    pytest.importorskip("Bio")
    seqs = [make_seq(f"s{i}", "ACGTACGT") for i in range(4)]
    cfg = {
        "clustering": {"alphabet_for_clustering": "nucleotide"},
        "segmented": {"enabled": False},
        "phylo": {"pre_cluster_tree": {"parttree_threshold": 10000}},
        "output": {"dir": str(tmp_path), "prefix": "test"},
        "temp_dir": str(tmp_path / "tmp"),
    }
    captured = {}
    with patch(
        "repseq.phylo.pre_cluster.run_mafft",
        side_effect=_capture_mafft_factory(captured),
    ), patch(
        "repseq.phylo.pre_cluster.run_fasttree", side_effect=_stub_fasttree,
    ):
        run_pre_cluster_phylogeny(seqs, [], cfg, tmp_path, "test")
    assert captured["extra_args"] == ["--retree", "1"]
