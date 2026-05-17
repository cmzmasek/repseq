"""Bench-scientist progress messages from the MAFFT and FastTree wrappers.

These lock the ``[phylo] starting … (…)`` / ``[phylo] X finished (Ys)``
format. The user runs the pipeline interactively, often on large
inputs where each step can take minutes; a silent terminal would let
the user assume the pipeline froze. The IQ-TREE equivalent lives in
``test_iqtree.py`` next to its other wrapper tests.
"""
from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

from repseq.phylo.fasttree import run_fasttree
from repseq.phylo.mafft import run_mafft


def _stub_subprocess_run_with_output(payload: str = "(A,B,C);\n"):
    """A subprocess.run stub that writes ``payload`` to the requested
    stdout file handle (the wrappers redirect stdout into the output
    file)."""
    def _run(cmd, check, stdout, stderr, text):
        stdout.write(payload)

        class _R:
            returncode = 0
        return _R()

    return _run


def test_mafft_prints_start_and_finish(tmp_path, capsys):
    input_fasta = tmp_path / "in.fasta"
    input_fasta.write_text(">a\nMKL\n>b\nMKL\n")
    out = tmp_path / "msa.fasta"

    with patch("repseq.phylo.mafft._check_mafft", return_value="/fake/mafft"), \
         patch(
             "repseq.phylo.mafft.subprocess.run",
             side_effect=_stub_subprocess_run_with_output(">a\nMKL\n>b\nMKL\n"),
         ):
        run_mafft(input_fasta, out, {"threads": 4})

    err = capsys.readouterr().err
    assert "[phylo] starting MAFFT (" in err
    assert "[phylo] MAFFT finished (" in err
    # Binary path and input filename are stripped.
    assert "/fake/mafft" not in err
    assert "in.fasta" not in err
    # But the user-meaningful args are shown.
    assert "--auto" in err
    assert "--thread 4" in err


def test_fasttree_prints_start_and_finish_nt(tmp_path, capsys):
    msa = tmp_path / "msa.fasta"
    msa.write_text(">a\nACGT\n>b\nACGT\n")
    out = tmp_path / "tree.nwk"

    with patch("repseq.phylo.fasttree._check_fasttree", return_value="/fake/FastTree"), \
         patch(
             "repseq.phylo.fasttree.subprocess.run",
             side_effect=_stub_subprocess_run_with_output("(a,b);\n"),
         ):
        run_fasttree(msa, out, {}, is_protein=False)

    err = capsys.readouterr().err
    assert "[phylo] starting FastTree (" in err
    assert "[phylo] FastTree finished (" in err
    assert "/fake/FastTree" not in err
    assert "msa.fasta" not in err
    # NT model flags are visible.
    assert "-nt" in err
    assert "-gtr" in err


def test_fasttree_protein_omits_nt_flags(tmp_path, capsys):
    """Protein runs don't pass -nt/-gtr; the start message reflects that."""
    msa = tmp_path / "msa.fasta"
    msa.write_text(">a\nMKL\n>b\nMKL\n")
    out = tmp_path / "tree.nwk"

    with patch("repseq.phylo.fasttree._check_fasttree", return_value="/fake/FastTree"), \
         patch(
             "repseq.phylo.fasttree.subprocess.run",
             side_effect=_stub_subprocess_run_with_output("(a,b);\n"),
         ):
        run_fasttree(msa, out, {}, is_protein=True)

    err = capsys.readouterr().err
    assert "[phylo] starting FastTree (" in err
    # No NT flags when running on protein.
    assert "-nt" not in err
    assert "-gtr" not in err
