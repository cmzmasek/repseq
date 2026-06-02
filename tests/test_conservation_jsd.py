"""Tests for the per-MSA JSD conservation scorer (phylo/conservation.py)."""

from __future__ import annotations

import random

import pytest

from repseq.phylo.conservation import (
    _classify,
    _untrimmed_sibling,
    detect_alphabet,
    henikoff_weights,
    score_rows,
    write_msa_conservation_report,
)


def test_identical_protein_scores_near_ceiling():
    """Invariant protein columns approach (but never reach) 1."""
    m = score_rows(["MKVLLACGTW"] * 8)
    assert m["alphabet"] == "protein"
    assert 0.80 <= m["mean_conservation"] <= 0.97
    assert m["mean_conservation"] < 1.0


def test_random_protein_scores_low():
    """Random sequences land well below a real family's conservation."""
    random.seed(1)
    aa = "ACDEFGHIKLMNPQRSTVWY"
    rows = ["".join(random.choice(aa) for _ in range(200)) for _ in range(40)]
    assert score_rows(rows)["mean_conservation"] < 0.30


def test_identical_nt_has_lower_ceiling_than_protein():
    """A perfectly conserved NT column is ~0.55 (4-symbol uniform bg)."""
    m = score_rows(["ACGTACGTACGT"] * 8)
    assert m["alphabet"] == "nucleotide"
    assert 0.50 <= m["mean_conservation"] <= 0.60


def test_random_nt_scores_near_zero():
    random.seed(2)
    rows = ["".join(random.choice("ACGT") for _ in range(300)) for _ in range(40)]
    assert score_rows(rows)["mean_conservation"] < 0.10


def test_henikoff_weighting_downweights_redundancy():
    """7 identical + 1 outlier must score below all-identical."""
    pure = score_rows(["MKVLLACGTW"] * 8)["mean_conservation"]
    biased = score_rows(["MKVLLACGTW"] * 7 + ["QQQQQQQQQQ"])["mean_conservation"]
    assert biased < pure


def test_henikoff_weights_sum_to_one():
    w = henikoff_weights(["MKVL", "MKVD", "QKVL"])
    assert pytest.approx(sum(w), abs=1e-9) == 1.0
    assert all(x > 0 for x in w)


def test_henikoff_single_sequence():
    assert henikoff_weights(["MKVL"]) == [1.0]


def test_gap_penalty_lowers_overall_but_core_stays_high():
    """A mostly-gap column drags mean_conservation below the core mean."""
    rows = ["MKVLLA----", "MKVLLA----", "MKVLDA-CGT", "MKVLLA----"]
    m = score_rows(rows)
    assert m["mean_conservation_core"] is not None
    assert m["mean_conservation"] < m["mean_conservation_core"]


def test_score_rows_empty_and_zero_width():
    assert score_rows([]) is None
    assert score_rows(["", ""]) is None


def test_detect_alphabet():
    assert detect_alphabet(["ACGTACGT", "ACGTACGN"]) == "nucleotide"
    assert detect_alphabet(["MKVLLWYF", "MKILLWYF"]) == "protein"
    # RNA folds to nucleotide too.
    assert detect_alphabet(["ACGUACGU"]) == "nucleotide"


def test_rna_u_folds_onto_t():
    dna = score_rows(["ACGTACGT"] * 5)["mean_conservation"]
    rna = score_rows(["ACGUACGU"] * 5)["mean_conservation"]
    assert pytest.approx(dna, abs=1e-9) == rna


# Third tuple element is is_untrimmed (the `_untrimmed` filename infix),
# NOT whether the alignment was trimmed — that's decided by sibling
# existence in the sweep (see test_trimmed_column_*).
@pytest.mark.parametrize(
    "rel, expect",
    [
        ("run_msa.fasta", ("genome", "", False)),
        ("run_msa_untrimmed.fasta", ("genome", "", True)),
        ("run_msa_Spike.fasta", ("partition_family", "Spike", False)),
        ("run_msa_Spike_untrimmed.fasta", ("partition_family", "Spike", True)),
        ("run_per_protein/Spike_msa.fasta", ("marker", "Spike", False)),
        ("run_per_protein/run_Spike_msa.fasta", ("marker", "Spike", False)),
        ("run_extra_protein/ORF7_msa.fasta", ("extra_protein", "ORF7", False)),
        ("run_polyprotein/ORF1ab_NSP3_msa.fasta", ("peptide", "ORF1ab_NSP3", False)),
        ("run_per_segment/L_msa.fasta", ("segment_nt", "L", False)),
    ],
)
def test_classify(rel, expect):
    from pathlib import Path

    assert _classify(Path(rel), "run") == expect


def test_untrimmed_sibling_paths():
    from pathlib import Path

    assert _untrimmed_sibling(Path("run_msa.fasta")) == Path("run_msa_untrimmed.fasta")
    assert _untrimmed_sibling(Path("run_msa_Spike.fasta")) == Path(
        "run_msa_Spike_untrimmed.fasta"
    )
    assert _untrimmed_sibling(Path("p/Spike_msa.fasta")) == Path(
        "p/Spike_msa_untrimmed.fasta"
    )


def _trimmed_col(tmp_path, prefix="run"):
    """Run the sweep and return {msa_path: trimmed_value} from the TSV."""
    path = write_msa_conservation_report(tmp_path, prefix, cfg=None)
    assert path is not None
    lines = path.read_text().strip().splitlines()
    header = lines[0].split("\t")
    i_msa, i_trim = header.index("msa"), header.index("trimmed")
    out = {}
    for line in lines[1:]:
        cells = line.split("\t")
        out[cells[i_msa]] = cells[i_trim]
    return out


def test_trimmed_column_false_when_trimal_off(tmp_path):
    """Default run (no `_untrimmed` companion): the tree-input `_msa.fasta`
    was never trimmed, so `trimmed` must be FALSE — the bug fix."""
    (tmp_path / "run_msa.fasta").write_text(">a\nMKVL\n>b\nMKVD\n>c\nMKVL\n")
    cols = _trimmed_col(tmp_path)
    assert cols["run_msa.fasta"] == "FALSE"


def test_trimmed_column_true_only_with_untrimmed_sibling(tmp_path):
    """With trimAl on, the trimmed tree input has an `_untrimmed` companion:
    the `_msa.fasta` row is TRUE, the `_untrimmed` row is FALSE."""
    (tmp_path / "run_msa.fasta").write_text(">a\nMKVL\n>b\nMKVD\n>c\nMKVL\n")
    (tmp_path / "run_msa_untrimmed.fasta").write_text(
        ">a\nMK--VL\n>b\nMK--VD\n>c\nMK--VL\n"
    )
    cols = _trimmed_col(tmp_path)
    assert cols["run_msa.fasta"] == "TRUE"
    assert cols["run_msa_untrimmed.fasta"] == "FALSE"


def test_write_report_sweeps_all_msas(tmp_path):
    prefix = "run"
    (tmp_path / f"{prefix}_per_protein").mkdir()
    # genome MSA
    (tmp_path / f"{prefix}_msa.fasta").write_text(
        ">a\nMKVLLACGTW\n>b\nMKVLLACGTW\n>c\nMKVLLACGTW\n"
    )
    # a marker MSA in the subdir
    (tmp_path / f"{prefix}_per_protein" / "Spike_msa.fasta").write_text(
        ">a\nMKVLDACGTW\n>b\nMKVLLACGTW\n>c\nMKVLLACGTW\n"
    )
    # a non-MSA file that must be ignored
    (tmp_path / f"{prefix}_tree.nwk").write_text("(a,b,c);\n")

    path = write_msa_conservation_report(tmp_path, prefix, cfg=None)
    assert path is not None
    lines = path.read_text().strip().splitlines()
    header = lines[0].split("\t")
    assert header[:4] == ["msa", "role", "label", "alphabet"]
    assert "mean_conservation" in header
    body = lines[1:]
    assert len(body) == 2  # genome + marker, tree.nwk ignored
    # Genome row sorts first.
    assert body[0].split("\t")[1] == "genome"
    assert body[1].split("\t")[1] == "marker"


def test_write_report_no_msas_returns_none(tmp_path):
    assert write_msa_conservation_report(tmp_path, "run", cfg=None) is None


def test_write_report_respects_disabled_config(tmp_path):
    (tmp_path / "run_msa.fasta").write_text(">a\nMKVL\n>b\nMKVL\n")
    cfg = {"phylo": {"conservation": {"enabled": False}}}
    assert write_msa_conservation_report(tmp_path, "run", cfg=cfg) is None
