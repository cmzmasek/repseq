"""Writer tests for :func:`repseq.output.report.write_polyprotein_outputs`.

Exercises the cross-rep aggregation, FASTA header tags, audit TSV
contents, and segmented vs non-segmented dispatch. As with
``test_polyprotein.py``, synthetic HMM hits mirror the producer schema
in ``hmm/hmmscan.py:_parse_domtblout`` (``target`` for the profile name).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from repseq.models import RunResult, Sequence, SequenceType
from repseq.output.report import write_polyprotein_outputs


def _hit(target: str, ali_from: int, ali_to: int, evalue: float = 1e-30) -> dict:
    return {
        "target": target,
        "ali_from": ali_from,
        "ali_to": ali_to,
        "dom_evalue": evalue,
        "evalue": evalue,
        "ali_span": ali_to - ali_from + 1,
    }


def _rep_with_polyprotein(seq_id: str = "AB12345") -> Sequence:
    poly_seq = "M" + "A" * 49 + "LQ" + "B" * 50 + "LQ" + "C" * 46
    protein = {
        "protein_id": f"YP_{seq_id}.1",
        "product": "polyprotein",
        "length": len(poly_seq),
        "sequence": poly_seq,
        "hmm_hits": [
            _hit("P_VP4", 1, 50),
            _hit("P_VP2", 53, 103),
            _hit("P_VP3", 106, 150),
        ],
    }
    return Sequence(
        id=seq_id,
        header=seq_id,
        sequence="ACGT" * len(poly_seq),  # any NT placeholder
        seq_type=SequenceType.NUCLEOTIDE,
        accession=seq_id,
        organism="Test virus",
        proteins=[protein],
    )


def _cfg_with_polyprotein() -> dict:
    return {
        "_hmm_runtime": {"active": True},
        "output": {"dir": "/dev/null", "prefix": "test"},
        "clustering": {
            "polyprotein": [
                {
                    "name": "P1",
                    "peptides": [
                        {"name": "VP4", "hmm": "P_VP4"},
                        {"name": "VP2", "hmm": "P_VP2", "cleavage_motif": "LQ"},
                        {"name": "VP3", "hmm": "P_VP3", "cleavage_motif": "LQ"},
                    ],
                    "cut_strategy": "motif",
                    "motif_window_aa": 10,
                    "min_peptides_hit": 2,
                }
            ]
        },
        "segmented": {"enabled": False},
        "hmm": {"enabled": True},
    }


def test_writes_one_fasta_per_peptide(tmp_path):
    reps = [_rep_with_polyprotein("AB000001"), _rep_with_polyprotein("CD000002")]
    result = RunResult(mode="global", representatives=reps, clusters=[])
    cfg = _cfg_with_polyprotein()
    written = write_polyprotein_outputs(result, cfg, tmp_path, "test")

    # Expect 3 peptide FASTAs + 1 audit TSV.
    fasta_paths = [p for p in written if p.suffix == ".fasta"]
    tsv_paths = [p for p in written if p.suffix == ".tsv"]
    assert len(fasta_paths) == 3
    assert len(tsv_paths) == 1
    sub = tmp_path / "test_polyprotein"
    assert sub.exists()
    assert (sub / "test_P1_VP4.fasta").exists()
    assert (sub / "test_P1_VP2.fasta").exists()
    assert (sub / "test_P1_VP3.fasta").exists()
    assert (sub / "test_P1_peptides.tsv").exists()


def test_fasta_header_carries_expected_bracket_tags(tmp_path):
    reps = [_rep_with_polyprotein("AB000001")]
    result = RunResult(mode="global", representatives=reps, clusters=[])
    cfg = _cfg_with_polyprotein()
    write_polyprotein_outputs(result, cfg, tmp_path, "test")

    body = (tmp_path / "test_polyprotein" / "test_P1_VP2.fasta").read_text()
    assert body.startswith(">YP_AB000001.1:VP2")
    # Polyprotein-specific tags must be present.
    assert "[polyprotein=P1]" in body
    assert "[peptide=VP2]" in body
    assert "[peptide_range_aa=" in body
    assert "[cut_method=motif:LQ]" in body
    # Shared bracket-tag set from _write_protein_fasta_record family.
    assert "[organism=Test virus]" in body
    assert "[length=" in body
    assert "[parent=AB000001]" in body


def test_audit_tsv_records_status_per_rep_per_peptide(tmp_path):
    # One rep with a clean slice + one rep with no parent CDS (no HMM hits).
    clean = _rep_with_polyprotein("CLEAN001")
    empty = Sequence(
        id="EMPTY001",
        header="EMPTY001",
        sequence="ACGT" * 50,
        seq_type=SequenceType.NUCLEOTIDE,
        accession="EMPTY001",
        organism="Test virus",
        proteins=[{
            "protein_id": "YP_EMPTY001.1",
            "product": "polyprotein",
            "length": 100,
            "sequence": "X" * 100,
            "hmm_hits": [],
        }],
    )
    reps = [clean, empty]
    result = RunResult(mode="global", representatives=reps, clusters=[])
    cfg = _cfg_with_polyprotein()
    write_polyprotein_outputs(result, cfg, tmp_path, "test")

    tsv = (tmp_path / "test_polyprotein" / "test_P1_peptides.tsv").read_text()
    lines = tsv.strip().split("\n")
    header = lines[0].split("\t")
    assert header[:5] == [
        "isolate_id", "peptide_name", "parent_accession",
        "parent_protein_id", "range_aa_from",
    ]
    # 2 reps × 3 peptides = 6 audit rows.
    rows = [line.split("\t") for line in lines[1:]]
    assert len(rows) == 6
    # The clean rep's rows should be all ok; the empty rep's all no_parent_cds.
    clean_rows = [r for r in rows if r[0] == "CLEAN001"]
    empty_rows = [r for r in rows if r[0] == "EMPTY001"]
    assert all(r[header.index("status")] == "ok" for r in clean_rows)
    assert all(r[header.index("status")] == "no_parent_cds" for r in empty_rows)


def test_soft_fails_when_hmm_tier_inactive(tmp_path, capsys):
    reps = [_rep_with_polyprotein("AB000001")]
    result = RunResult(mode="global", representatives=reps, clusters=[])
    cfg = _cfg_with_polyprotein()
    cfg["_hmm_runtime"]["active"] = False
    # Strip HMM hits to ensure _hmm_tier_ran returns False.
    for rep in reps:
        for prot in rep.proteins or []:
            prot["hmm_hits"] = []

    written = write_polyprotein_outputs(result, cfg, tmp_path, "test")
    assert written == []
    err = capsys.readouterr().err
    assert "polyprotein cutting skipped" in err
    # No outputs directory should have been created.
    assert not (tmp_path / "test_polyprotein").exists()


def test_no_specs_declared_is_a_silent_noop(tmp_path):
    reps = [_rep_with_polyprotein("AB000001")]
    result = RunResult(mode="global", representatives=reps, clusters=[])
    cfg = _cfg_with_polyprotein()
    cfg["clustering"]["polyprotein"] = []
    assert write_polyprotein_outputs(result, cfg, tmp_path, "test") == []
    assert not (tmp_path / "test_polyprotein").exists()
