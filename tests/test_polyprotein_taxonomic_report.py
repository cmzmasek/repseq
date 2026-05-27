"""Writer tests for
:func:`repseq.output.report.write_polyprotein_taxonomic_report`.

Mirrors the synthetic-HMM-hit shape used by ``test_polyprotein.py`` and
``test_polyprotein_output.py`` (the ``target`` key — what
``hmm/hmmscan.py:_parse_domtblout`` actually emits — drives the slicer).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from repseq.models import Sequence, SequenceType, TaxonomyInfo
from repseq.output.report import write_polyprotein_taxonomic_report


def _hit(target: str, ali_from: int, ali_to: int, evalue: float = 1e-30) -> dict:
    return {
        "target": target,
        "ali_from": ali_from,
        "ali_to": ali_to,
        "dom_evalue": evalue,
        "evalue": evalue,
        "ali_span": ali_to - ali_from + 1,
    }


def _rep(seq_id: str, genus: str, *, hits: list[dict] | None = None) -> Sequence:
    poly_seq = "M" + "A" * 49 + "LQ" + "B" * 50 + "LQ" + "C" * 46
    protein = {
        "protein_id": f"YP_{seq_id}.1",
        "product": "polyprotein",
        "length": len(poly_seq),
        "sequence": poly_seq,
        "hmm_hits": hits if hits is not None else [
            _hit("P_VP4", 1, 50),
            _hit("P_VP2", 53, 103),
            _hit("P_VP3", 106, 150),
        ],
    }
    return Sequence(
        id=seq_id,
        header=seq_id,
        sequence="ACGT" * len(poly_seq),
        seq_type=SequenceType.NUCLEOTIDE,
        accession=seq_id,
        organism="Test virus",
        proteins=[protein],
        taxonomy=TaxonomyInfo(genus=genus),
    )


def _cfg() -> dict:
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


def test_soft_fails_when_no_polyprotein_spec(tmp_path):
    cfg = _cfg()
    cfg["clustering"]["polyprotein"] = []
    out = tmp_path / "report.txt"
    assert write_polyprotein_taxonomic_report(
        [_rep("A001", "Alphavirus")], [_rep("A001", "Alphavirus")],
        cfg, segmented=False, path=out,
    ) is False
    assert not out.exists()


def test_soft_fails_when_hmm_tier_off(tmp_path):
    """When the HMM tier didn't run AND no rep carries hmm_hits, the
    writer mirrors :func:`write_polyprotein_outputs` and emits nothing.
    The fallback in ``_hmm_tier_ran`` scans for cached hits, so we
    also strip them from the test reps."""
    cfg = _cfg()
    cfg["_hmm_runtime"] = {"active": False}
    rep = _rep("A001", "Alphavirus", hits=[])
    out = tmp_path / "report.txt"
    assert write_polyprotein_taxonomic_report(
        [rep], [rep], cfg, segmented=False, path=out,
    ) is False
    assert not out.exists()


def test_writes_spec_section_with_all_peptide_columns(tmp_path):
    reps = [
        _rep("A001", "Alphavirus"),
        _rep("A002", "Alphavirus"),
        _rep("B001", "Betavirus"),
    ]
    out = tmp_path / "report.txt"
    ok = write_polyprotein_taxonomic_report(
        reps, reps, _cfg(), segmented=False, path=out,
    )
    assert ok is True
    body = out.read_text()
    # H2 banner per spec
    assert "========== P1 ==========" in body
    # All three declared peptides as columns (even if some had no coverage).
    assert "VP4" in body
    assert "VP2" in body
    assert "VP3" in body
    # Two sub-table titles at genus rank (the only rank populated here).
    assert "coverage (post-QC pool" in body
    assert "coverage (representatives" in body
    assert "peptide length statistics" in body
    # Trailing architectures block.
    assert "== Peptide architectures ==" in body
    assert "VP2: P_VP2  [cleavage_motif=LQ]" in body


def test_missing_peptide_keeps_column_with_zero_coverage(tmp_path):
    # Rep carries VP4 and VP2 but not VP3 (no P_VP3 hit).
    partial = _rep("A001", "Alphavirus", hits=[
        _hit("P_VP4", 1, 50),
        _hit("P_VP2", 53, 103),
    ])
    out = tmp_path / "report.txt"
    write_polyprotein_taxonomic_report(
        [partial], [partial], _cfg(), segmented=False, path=out,
    )
    body = out.read_text()
    # VP3 must still appear as a header column ...
    assert "VP3" in body
    # ... but its coverage cell for the only genus row is 0.
    lines = [ln for ln in body.splitlines() if ln.lstrip().startswith("Alphavirus")]
    assert lines, "expected at least one Alphavirus row"
    # The first such row is the post-QC coverage row; its last column is VP3.
    cov_row = lines[0]
    assert cov_row.rstrip().endswith("0 0%") or " 0 0%" in cov_row


def test_no_parent_cds_excluded_from_coverage(tmp_path):
    """A rep whose polyprotein CDS has no peptide HMM hits at all gets a
    `no_parent_cds` audit row per peptide — it must NOT count toward
    coverage (since no FASTA record is written for it either)."""
    no_hits = _rep("A001", "Alphavirus", hits=[])
    out = tmp_path / "report.txt"
    write_polyprotein_taxonomic_report(
        [no_hits], [no_hits], _cfg(), segmented=False, path=out,
    )
    body = out.read_text()
    # The genus row exists (1 item at Alphavirus) but VP4/VP2/VP3 all 0.
    rows = [ln for ln in body.splitlines() if ln.lstrip().startswith("Alphavirus")]
    assert rows, "expected an Alphavirus row in the coverage table"
    cov_row = rows[0]
    # All three peptide cells should be `0 0%`.
    assert cov_row.count("0 0%") == 3


def test_returns_true_and_writes_file(tmp_path):
    reps = [_rep("A001", "Alphavirus")]
    out = tmp_path / "report.txt"
    assert write_polyprotein_taxonomic_report(
        reps, reps, _cfg(), segmented=False, path=out,
    ) is True
    assert out.exists()
    assert out.stat().st_size > 0


def test_per_spec_h2_for_multiple_specs(tmp_path):
    cfg = _cfg()
    cfg["clustering"]["polyprotein"].append({
        "name": "ORF2",
        "peptides": [
            {"name": "X1", "hmm": "P_VP4"},
            {"name": "X2", "hmm": "P_VP2"},
        ],
        "cut_strategy": "bisect",
        "min_peptides_hit": 1,
    })
    reps = [_rep("A001", "Alphavirus")]
    out = tmp_path / "report.txt"
    write_polyprotein_taxonomic_report(reps, reps, cfg, segmented=False, path=out)
    body = out.read_text()
    assert "========== P1 ==========" in body
    assert "========== ORF2 ==========" in body
    # Each spec has its own architecture block.
    assert body.count("== Peptide architectures ==") == 2
