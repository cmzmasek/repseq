"""Output writer: protein-alphabet representatives FASTA."""
from __future__ import annotations

from repseq.models import RunResult
from repseq.output.writer import write_results


def test_write_results_emits_protein_fasta_when_reps_have_protein_sequence(
    tmp_path, make_seq,
):
    out_dir = tmp_path / "out"
    a = make_seq("a", "ACGTACGT")
    a.protein_sequence = "MEEPMEEP"
    b = make_seq("b", "TTTTGGGG")
    b.protein_sequence = "QQQQRRRR"
    result = RunResult(mode="global", representatives=[a, b])
    cfg = {"output": {"dir": str(out_dir), "prefix": "p"}}
    written = write_results(result, cfg)
    names = {p.name for p in written}
    assert "p_representatives.fasta" in names
    assert "p_representatives_protein.fasta" in names
    aa = (out_dir / "p_representatives_protein.fasta").read_text()
    assert ">a" in aa and "MEEPMEEP" in aa
    assert ">b" in aa and "QQQQRRRR" in aa


def test_write_results_skips_protein_fasta_when_no_protein_sequence(
    tmp_path, make_seq,
):
    """Pure nucleotide run — protein FASTA must not be emitted."""
    out_dir = tmp_path / "out"
    a = make_seq("a", "ACGTACGT")
    result = RunResult(mode="global", representatives=[a])
    cfg = {"output": {"dir": str(out_dir), "prefix": "p"}}
    written = write_results(result, cfg)
    names = {p.name for p in written}
    assert "p_representatives.fasta" in names
    assert "p_representatives_protein.fasta" not in names
    assert not (out_dir / "p_representatives_protein.fasta").exists()
