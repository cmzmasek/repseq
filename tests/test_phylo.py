"""Phylogeny step: id-remap round-trip, skip rules, Newick→phyloXML.

The MAFFT and FastTree subprocesses are mocked: each test stubs out
``run_mafft`` and ``run_fasttree`` to write deterministic intermediate
files, so the orchestrator's id assignment, name restoration, and skip
behaviour can be locked without touching real binaries.
"""
from __future__ import annotations

from pathlib import Path
from unittest.mock import patch
from xml.etree import ElementTree as ET

import pytest

from repseq.models import Sequence, SequenceType
from repseq.phylo.fasttree import FastTreeError
from repseq.phylo.mafft import MafftError
from repseq.phylo.pipeline import (
    PhyloError,
    _build_id_map,
    _write_id_map,
    _write_short_id_fasta,
    run_phylogeny,
)


def _seq(sid: str, seq: str, seq_type: SequenceType = SequenceType.PROTEIN) -> Sequence:
    return Sequence(id=sid, header=sid, sequence=seq, seq_type=seq_type, accession=sid)


# ---------------------------------------------------------------------------
# Short-id assignment
# ---------------------------------------------------------------------------

def test_build_id_map_deterministic_order():
    reps = [_seq("alpha", "MK"), _seq("CONCAT|iso-1", "MK"), _seq("β-name", "MK")]
    m = _build_id_map(reps)
    assert list(m.items()) == [
        ("S0001", "alpha"),
        ("S0002", "CONCAT|iso-1"),
        ("S0003", "β-name"),
    ]


def test_write_short_id_fasta_uses_short_id_as_sole_header(tmp_path):
    reps = [_seq("CONCAT|iso 1 with space", "MKLPQE"),
            _seq("very long descriptor that would choke FastTree", "MMMMMM")]
    id_map = _build_id_map(reps)
    out = tmp_path / "input.fasta"
    _write_short_id_fasta(reps, id_map, out)

    header_tokens = [
        line[1:].strip()
        for line in out.read_text().splitlines()
        if line.startswith(">")
    ]
    assert header_tokens == ["S0001", "S0002"]


def test_write_id_map_writes_round_trip_tsv(tmp_path):
    id_map = {"S0001": "alpha", "S0002": "CONCAT|iso-1"}
    out = tmp_path / "id_map.tsv"
    _write_id_map(id_map, out)

    rows = out.read_text().splitlines()
    assert rows[0] == "short_id\toriginal_id"
    assert set(rows[1:]) == {"S0001\talpha", "S0002\tCONCAT|iso-1"}


# ---------------------------------------------------------------------------
# Orchestrator: skip rules + happy path with mocked binaries
# ---------------------------------------------------------------------------

def test_run_phylogeny_skips_with_fewer_than_three_reps(tmp_path):
    reps = [_seq("a", "MK"), _seq("b", "MK")]
    with pytest.raises(PhyloError, match="need >= 3"):
        run_phylogeny(reps, {}, tmp_path, "test")


def _stub_mafft_writes_alignment(input_fasta: Path, output_fasta: Path, cfg):
    # Echo the input back as a "trivial alignment" — same headers, same
    # sequences. The orchestrator never inspects the alignment content,
    # only its path, so an identity copy is enough.
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    output_fasta.write_text(input_fasta.read_text())


def _stub_fasttree_writes_newick(short_ids: list[str]):
    """Build a stub that writes a balanced Newick over the given short ids."""

    def _run(msa_fasta: Path, output_newick: Path, cfg, is_protein):
        # Ladder topology — every leaf is named with its short id so the
        # orchestrator's rename step has something to match on.
        body = short_ids[0]
        for sid in short_ids[1:]:
            body = f"({body}:0.1,{sid}:0.1)"
        output_newick.parent.mkdir(parents=True, exist_ok=True)
        output_newick.write_text(body + ";\n")

    return _run


def test_run_phylogeny_happy_path_writes_all_outputs(tmp_path):
    reps = [
        _seq("alpha", "MKLPQEFIL"),
        _seq("beta",  "MKLPQEFIA"),
        _seq("gamma", "MKLPQEFIY"),
    ]
    short_ids = ["S0001", "S0002", "S0003"]
    # Pin FastTree so the test stays alphabet-agnostic; the IQ-TREE
    # dispatch path has its own dedicated test.
    cfg = {"phylo": {"tool": "fasttree"}}

    with patch("repseq.phylo.pipeline.run_mafft", side_effect=_stub_mafft_writes_alignment), \
         patch("repseq.phylo.pipeline.run_fasttree",
               side_effect=_stub_fasttree_writes_newick(short_ids)):
        files = run_phylogeny(reps, cfg, tmp_path, "test")

    names = [f.name for f in files]
    assert names == [
        "test_msa.fasta",
        "test_tree.nwk",
        "test_tree.xml",
        "test_tree_id_map.tsv",
    ]
    # phyloXML: terminal names restored.
    xml = (tmp_path / "test_tree.xml").read_text()
    for original in ("alpha", "beta", "gamma"):
        assert original in xml
    for short in short_ids:
        assert short not in xml
    # Newick keeps the short ids — by design, since it is what FastTree
    # produced and the id_map.tsv lets downstream tools decode it.
    nwk = (tmp_path / "test_tree.nwk").read_text()
    for short in short_ids:
        assert short in nwk
    # The temp input FASTA should have been cleaned up.
    assert not (tmp_path / "test_phylo_input.fasta").exists()


def test_run_phylogeny_picks_nucleotide_model_for_nt_reps(tmp_path):
    reps = [
        _seq(f"n{i}", "ACGTACGT", seq_type=SequenceType.NUCLEOTIDE) for i in range(3)
    ]
    seen_is_protein: list[bool] = []

    def _fasttree(msa_fasta, output_newick, cfg, is_protein):
        seen_is_protein.append(is_protein)
        # Write a valid 3-leaf tree so the orchestrator can proceed.
        output_newick.write_text("(S0001:0.1,(S0002:0.1,S0003:0.1):0.1);\n")

    # auto picks FastTree for NT input, so no need to pin the tool.
    with patch("repseq.phylo.pipeline.run_mafft", side_effect=_stub_mafft_writes_alignment), \
         patch("repseq.phylo.pipeline.run_fasttree", side_effect=_fasttree):
        run_phylogeny(reps, {}, tmp_path, "test")

    assert seen_is_protein == [False]


def test_run_phylogeny_uses_protein_sequence_when_alphabet_protein(tmp_path):
    """alphabet=protein on NT-typed concat reps: MSA input is AA, model is protein."""
    reps = [
        _seq(f"CONCAT|iso{i}", "ACGTACGTACGT", seq_type=SequenceType.NUCLEOTIDE)
        for i in range(3)
    ]
    for i, rep in enumerate(reps):
        rep.protein_sequence = "M" + ("K" * (10 + i))
    seen_is_protein: list[bool] = []
    seen_input_bodies: list[str] = []

    def _mafft(input_fasta, output_fasta, cfg):
        # Capture what's being aligned.
        body = "".join(
            line for line in input_fasta.read_text().splitlines()
            if not line.startswith(">")
        )
        seen_input_bodies.append(body)
        _stub_mafft_writes_alignment(input_fasta, output_fasta, cfg)

    def _fasttree(msa_fasta, output_newick, cfg, is_protein):
        seen_is_protein.append(is_protein)
        output_newick.write_text("(S0001:0.1,(S0002:0.1,S0003:0.1):0.1);\n")

    with patch("repseq.phylo.pipeline.run_mafft", side_effect=_mafft), \
         patch("repseq.phylo.pipeline.run_fasttree", side_effect=_fasttree):
        # Pin FastTree to keep the assertions on this code path; the
        # IQ-TREE-on-protein dispatch is covered separately.
        run_phylogeny(
            reps,
            {"clustering": {"alphabet": "protein"}, "phylo": {"tool": "fasttree"}},
            tmp_path, "test",
        )

    assert seen_is_protein == [True]
    # MSA input must come from protein_sequence, not seq.sequence (the NT one).
    assert all("ACGT" not in body for body in seen_input_bodies)
    assert all("M" in body and "K" in body for body in seen_input_bodies)


def test_run_phylogeny_wraps_mafft_error_as_phyloerror(tmp_path):
    reps = [_seq(f"p{i}", "MK") for i in range(3)]

    def _boom(input_fasta, output_fasta, cfg):
        raise MafftError("mafft segfaulted")

    with patch("repseq.phylo.pipeline.run_mafft", side_effect=_boom):
        with pytest.raises(PhyloError, match="mafft segfaulted"):
            run_phylogeny(reps, {}, tmp_path, "test")


def test_run_phylogeny_wraps_fasttree_error_as_phyloerror(tmp_path):
    reps = [_seq(f"p{i}", "MK") for i in range(3)]

    def _boom(msa_fasta, output_newick, cfg, is_protein):
        raise FastTreeError("FastTree exit 137")

    with patch("repseq.phylo.pipeline.run_mafft", side_effect=_stub_mafft_writes_alignment), \
         patch("repseq.phylo.pipeline.run_fasttree", side_effect=_boom):
        with pytest.raises(PhyloError, match="FastTree exit 137"):
            run_phylogeny(reps, {"phylo": {"tool": "fasttree"}}, tmp_path, "test")


# ---------------------------------------------------------------------------
# IQ-TREE dispatch and run
# ---------------------------------------------------------------------------

def _stub_iqtree_writes_newick(short_ids: list[str], summary_text: str = ""):
    """Stub for run_iqtree: mirrors the FastTree stub but also writes a
    minimal summary file when summary_path is provided."""

    def _run(msa_fasta, output_newick, cfg, is_protein, summary_path=None):
        body = short_ids[0]
        for sid in short_ids[1:]:
            body = f"({body}:0.1,{sid}:0.1)"
        output_newick.parent.mkdir(parents=True, exist_ok=True)
        output_newick.write_text(body + ";\n")
        if summary_path is not None and summary_text:
            summary_path.parent.mkdir(parents=True, exist_ok=True)
            summary_path.write_text(summary_text)

    return _run


def test_run_phylogeny_picks_iqtree_for_protein_by_default(tmp_path):
    """`phylo.tool=auto` (the default) picks IQ-TREE for protein reps."""
    reps = [_seq(f"p{i}", "MK") for i in range(3)]
    short_ids = ["S0001", "S0002", "S0003"]
    called: list[str] = []

    def _iqtree(msa_fasta, output_newick, cfg, is_protein, summary_path=None):
        called.append("iqtree")
        _stub_iqtree_writes_newick(short_ids)(
            msa_fasta, output_newick, cfg, is_protein, summary_path=summary_path,
        )

    def _fasttree(*a, **kw):
        called.append("fasttree")

    with patch("repseq.phylo.pipeline.run_mafft", side_effect=_stub_mafft_writes_alignment), \
         patch("repseq.phylo.pipeline.run_iqtree", side_effect=_iqtree), \
         patch("repseq.phylo.pipeline.run_fasttree", side_effect=_fasttree):
        run_phylogeny(reps, {}, tmp_path, "test")

    assert called == ["iqtree"]


def test_run_phylogeny_picks_fasttree_for_nucleotide_by_default(tmp_path):
    """`phylo.tool=auto` picks FastTree for nucleotide reps."""
    reps = [
        _seq(f"n{i}", "ACGTACGT", seq_type=SequenceType.NUCLEOTIDE) for i in range(3)
    ]
    short_ids = ["S0001", "S0002", "S0003"]
    called: list[str] = []

    def _iqtree(*a, **kw):
        called.append("iqtree")

    def _fasttree(msa_fasta, output_newick, cfg, is_protein):
        called.append("fasttree")
        _stub_fasttree_writes_newick(short_ids)(msa_fasta, output_newick, cfg, is_protein)

    with patch("repseq.phylo.pipeline.run_mafft", side_effect=_stub_mafft_writes_alignment), \
         patch("repseq.phylo.pipeline.run_iqtree", side_effect=_iqtree), \
         patch("repseq.phylo.pipeline.run_fasttree", side_effect=_fasttree):
        run_phylogeny(reps, {}, tmp_path, "test")

    assert called == ["fasttree"]


def test_run_phylogeny_iqtree_summary_file_appended_to_output_list(tmp_path):
    reps = [_seq(f"p{i}", "MK") for i in range(3)]
    short_ids = ["S0001", "S0002", "S0003"]

    with patch("repseq.phylo.pipeline.run_mafft", side_effect=_stub_mafft_writes_alignment), \
         patch(
             "repseq.phylo.pipeline.run_iqtree",
             side_effect=_stub_iqtree_writes_newick(short_ids, "BIC model: LG+G4\n"),
         ):
        files = run_phylogeny(reps, {}, tmp_path, "test")

    names = [f.name for f in files]
    assert "test_iqtree_summary.txt" in names
    assert (tmp_path / "test_iqtree_summary.txt").read_text() == "BIC model: LG+G4\n"


def test_run_phylogeny_iqtree_no_summary_means_no_extra_file(tmp_path):
    """The `.iqtree` summary is optional — if the wrapper didn't write one,
    the output list does not include the path."""
    reps = [_seq(f"p{i}", "MK") for i in range(3)]
    short_ids = ["S0001", "S0002", "S0003"]

    with patch("repseq.phylo.pipeline.run_mafft", side_effect=_stub_mafft_writes_alignment), \
         patch(
             "repseq.phylo.pipeline.run_iqtree",
             side_effect=_stub_iqtree_writes_newick(short_ids),  # no summary text
         ):
        files = run_phylogeny(reps, {}, tmp_path, "test")

    assert all(f.name != "test_iqtree_summary.txt" for f in files)


def test_run_phylogeny_explicit_tool_iqtree_overrides_alphabet(tmp_path):
    """`phylo.tool=iqtree` picks IQ-TREE even for nucleotide reps."""
    reps = [
        _seq(f"n{i}", "ACGTACGT", seq_type=SequenceType.NUCLEOTIDE) for i in range(3)
    ]
    short_ids = ["S0001", "S0002", "S0003"]
    called: list[str] = []

    def _iqtree(msa_fasta, output_newick, cfg, is_protein, summary_path=None):
        called.append("iqtree")
        _stub_iqtree_writes_newick(short_ids)(
            msa_fasta, output_newick, cfg, is_protein, summary_path=summary_path,
        )

    with patch("repseq.phylo.pipeline.run_mafft", side_effect=_stub_mafft_writes_alignment), \
         patch("repseq.phylo.pipeline.run_iqtree", side_effect=_iqtree):
        run_phylogeny(reps, {"phylo": {"tool": "iqtree"}}, tmp_path, "test")

    assert called == ["iqtree"]


def test_run_phylogeny_wraps_iqtree_error_as_phyloerror(tmp_path):
    from repseq.phylo.iqtree import IQTreeError
    reps = [_seq(f"p{i}", "MK") for i in range(3)]

    def _boom(msa_fasta, output_newick, cfg, is_protein, summary_path=None):
        raise IQTreeError("IQ-TREE segfaulted")

    with patch("repseq.phylo.pipeline.run_mafft", side_effect=_stub_mafft_writes_alignment), \
         patch("repseq.phylo.pipeline.run_iqtree", side_effect=_boom):
        with pytest.raises(PhyloError, match="IQ-TREE segfaulted"):
            run_phylogeny(reps, {}, tmp_path, "test")


# ---------------------------------------------------------------------------
# CLI integration: --phylo skip is surfaced to stderr without aborting
# ---------------------------------------------------------------------------

def test_write_output_phylo_skip_does_not_abort(tmp_path, capsys):
    # When the phylo step raises PhyloError, _write_output should swallow
    # it and emit "[phylo skipped] ..." on stderr — matching the existing
    # --plot behaviour. The rest of _write_output must still complete.
    from repseq.cli import _write_output
    from repseq.models import Cluster, RunResult

    rep = _seq("only_one", "MK")
    result = RunResult(
        mode="global:count",
        representatives=[rep],
        clusters=[Cluster(cluster_id="c1", representative=rep)],
        group_stats=[],
        config_snapshot={},
    )
    cfg = {
        "output": {"dir": str(tmp_path), "prefix": "t"},
        "qc": {"remove_duplicates": False},
        "seed": 42, "threads": 1,
    }

    # write_results would itself need a full config; stub it so we only
    # exercise the phylo branch.
    with patch("repseq.cli.write_results", return_value=[]), \
         patch("repseq.cli.write_all_reports"):
        _write_output(
            result, qc_report=None, cfg=cfg, input_paths=[],
            complete_isolates=None, segment_names=None,
            phylo=True,
        )

    err = capsys.readouterr().err
    assert "[phylo skipped]" in err
    assert "need >= 3" in err
