"""Tests for the tidy long-format TSV companions to the four
``_*_taxonomic_report.txt`` files.

Schema invariant under test: every TSV emits the same 8-column header
in the same order, one observation per row, parsable by a vanilla CSV
reader. The data-gathering helpers are shared with the `.txt` writers,
so these tests focus on the TSV shape, not the underlying counts.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from repseq.models import Sequence, SequenceType, TaxonomyInfo
from repseq.output.report import (
    _TIDY_TSV_COLUMNS,
    write_nucleotide_taxonomic_report_tsv,
    write_polyprotein_taxonomic_report_tsv,
    write_protein_taxonomic_report_tsv,
    write_taxonomic_report_tsv,
)


# -----------------------------------------------------------
# Fixtures
# -----------------------------------------------------------

def _hit(target: str, ali_from: int, ali_to: int, evalue: float = 1e-30) -> dict:
    return {
        "target": target,
        "ali_from": ali_from,
        "ali_to": ali_to,
        "dom_evalue": evalue,
        "evalue": evalue,
        "ali_span": ali_to - ali_from + 1,
    }


def _rep(seq_id: str, genus: str, *, hits=None, length: int = 800) -> Sequence:
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
        sequence="A" * length,
        seq_type=SequenceType.NUCLEOTIDE,
        accession=seq_id,
        organism="Test virus",
        proteins=[protein],
        taxonomy=TaxonomyInfo(genus=genus),
    )


def _polyprotein_cfg() -> dict:
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
                },
            ]
        },
        "segmented": {"enabled": False},
        "hmm": {"enabled": True},
    }


# -----------------------------------------------------------
# Diversity TSV
# -----------------------------------------------------------

def test_diversity_tsv_header_and_distinct_rows(tmp_path):
    reps = [
        _rep("A001", "Alphavirus"),
        _rep("A002", "Alphavirus"),
        _rep("B001", "Betavirus"),
    ]
    out = tmp_path / "tax.tsv"
    assert write_taxonomic_report_tsv(reps, reps[:2], path=out) is True
    with out.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    # Header order must match the canonical schema.
    assert list(csv.reader(out.open(), delimiter="\t").__next__()) == list(_TIDY_TSV_COLUMNS)
    # Every row is "diversity" type.
    assert all(r["report"] == "diversity" for r in rows)
    # The pool column only carries our two pools.
    assert set(r["pool"] for r in rows) == {"post_qc", "reps"}
    # *ALL* distinct_taxa rows: one per (rank, pool); 9 ranks × 2 pools = 18.
    all_rows = [r for r in rows if r["taxon"] == "*ALL*"]
    assert len(all_rows) == 18
    # The genus distinct row should match what we set up.
    genus_post_qc = next(
        r for r in all_rows
        if r["rank"] == "genus" and r["pool"] == "post_qc"
    )
    assert genus_post_qc["metric"] == "distinct_taxa"
    assert genus_post_qc["value"] == "2"


def test_diversity_tsv_member_count_present_in_both_pools(tmp_path):
    reps = [_rep("A001", "Alphavirus"), _rep("B001", "Betavirus")]
    out = tmp_path / "tax.tsv"
    write_taxonomic_report_tsv(reps, [reps[0]], path=out)
    with out.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    # Betavirus exists in post_qc=1 but reps=0; both rows should be present.
    beta_post = [r for r in rows if r["taxon"] == "Betavirus" and r["pool"] == "post_qc"]
    beta_reps = [r for r in rows if r["taxon"] == "Betavirus" and r["pool"] == "reps"]
    assert beta_post and beta_reps
    assert any(r["metric"] == "member_count" and r["value"] == "1" for r in beta_post)
    assert any(r["metric"] == "member_count" and r["value"] == "0" for r in beta_reps)


# -----------------------------------------------------------
# Nucleotide TSV
# -----------------------------------------------------------

def test_nucleotide_tsv_emits_length_metrics_only(tmp_path):
    reps = [
        _rep("A001", "Alphavirus", length=800),
        _rep("A002", "Alphavirus", length=820),
    ]
    out = tmp_path / "nt.tsv"
    cfg = {"output": {"dir": "/dev/null", "prefix": "test"}}
    assert write_nucleotide_taxonomic_report_tsv(
        reps, reps, cfg, segmented=False, path=out,
    ) is True
    with out.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert all(r["report"] == "nucleotide" for r in rows)
    assert all(r["spec"] == "genome" for r in rows)
    # No coverage_* metrics in the NT report (segments always 100%).
    assert not any(r["metric"].startswith("coverage_") for r in rows)
    # We should see at least one length_min row at genus.
    assert any(
        r["rank"] == "genus" and r["metric"] == "length_min" and r["value"] == "800"
        for r in rows
    )


# -----------------------------------------------------------
# Polyprotein TSV
# -----------------------------------------------------------

def test_polyprotein_tsv_soft_fails_when_no_spec(tmp_path):
    cfg = _polyprotein_cfg()
    cfg["clustering"]["polyprotein"] = []
    out = tmp_path / "pp.tsv"
    assert write_polyprotein_taxonomic_report_tsv(
        [_rep("A001", "Alphavirus")], [_rep("A001", "Alphavirus")],
        cfg, path=out,
    ) is False
    assert not out.exists()


def test_polyprotein_tsv_soft_fails_when_hmm_off(tmp_path):
    cfg = _polyprotein_cfg()
    cfg["_hmm_runtime"] = {"active": False}
    out = tmp_path / "pp.tsv"
    rep = _rep("A001", "Alphavirus", hits=[])
    assert write_polyprotein_taxonomic_report_tsv(
        [rep], [rep], cfg, path=out,
    ) is False


def test_polyprotein_tsv_uses_composite_spec_column(tmp_path):
    reps = [_rep("A001", "Alphavirus")]
    out = tmp_path / "pp.tsv"
    write_polyprotein_taxonomic_report_tsv(
        reps, reps, _polyprotein_cfg(), path=out,
    )
    with out.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert all(r["report"] == "polyprotein" for r in rows)
    # Composite spec: <poly_name>:<peptide_name>.
    specs = {r["spec"] for r in rows}
    assert "P1:VP4" in specs and "P1:VP2" in specs and "P1:VP3" in specs


def test_polyprotein_tsv_skips_length_when_coverage_zero(tmp_path):
    """A peptide with no successful slice should emit coverage_count=0,
    coverage_pct=0, length_n=0 — and NO length_min/max/median/iqr rows
    (those would be fake values that skew downstream aggregations)."""
    partial = _rep("A001", "Alphavirus", hits=[
        _hit("P_VP4", 1, 50),
        _hit("P_VP2", 53, 103),
        # VP3 missing.
    ])
    out = tmp_path / "pp.tsv"
    write_polyprotein_taxonomic_report_tsv(
        [partial], [partial], _polyprotein_cfg(), path=out,
    )
    with out.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    vp3_rows = [r for r in rows if r["spec"] == "P1:VP3"]
    assert vp3_rows, "expected at least one VP3 row"
    metrics_emitted = {r["metric"] for r in vp3_rows}
    # The 0-coverage spec should carry coverage_count, coverage_pct, length_n.
    assert "coverage_count" in metrics_emitted
    assert "coverage_pct" in metrics_emitted
    assert "length_n" in metrics_emitted
    # And NOT the other length stats.
    assert "length_min" not in metrics_emitted
    assert "length_max" not in metrics_emitted
    assert "length_median" not in metrics_emitted
    assert "length_iqr" not in metrics_emitted


# -----------------------------------------------------------
# Protein TSV — covers the dispatch shared between cluster_protein
# / extra_protein. Uses a minimal alias-only spec (no HMM gate) so we
# don't need to wire HMM hits matching the marker.
# -----------------------------------------------------------

def test_protein_tsv_emits_per_spec_rows(tmp_path):
    rep = _rep("A001", "Alphavirus", hits=[])
    # Add a satisfying CDS for an alias-only marker.
    rep.proteins.append({
        "protein_id": "YP_marker.1",
        "product": "spike glycoprotein",
        "length": 1200,
        "sequence": "S" * 1200,
        "hmm_hits": [],
    })
    cfg = {
        "_hmm_runtime": {"active": False},
        "output": {"dir": "/dev/null", "prefix": "test"},
        "clustering": {
            "cluster_protein": [
                {"name": "Spike", "aliases": ["spike"]},
            ],
        },
        "segmented": {"enabled": False},
        "hmm": {"enabled": False},
    }
    out = tmp_path / "pr.tsv"
    assert write_protein_taxonomic_report_tsv(
        [rep], [rep], cfg, path=out,
    ) is True
    with out.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert all(r["report"] == "protein" for r in rows)
    # Cluster-driving markers carry the * suffix (parity with the .txt).
    assert any(r["spec"] == "Spike*" for r in rows)
    # Should see coverage_pct=100 since the only seq carries the marker.
    cov_pct_rows = [r for r in rows if r["metric"] == "coverage_pct"]
    assert cov_pct_rows
    assert any(r["value"] == "100" for r in cov_pct_rows)


def test_protein_tsv_returns_false_when_no_specs(tmp_path):
    out = tmp_path / "pr.tsv"
    rep = _rep("A001", "Alphavirus", hits=[])
    cfg = {
        "_hmm_runtime": {"active": False},
        "output": {"dir": "/dev/null", "prefix": "test"},
        "clustering": {},
        "segmented": {"enabled": False},
    }
    assert write_protein_taxonomic_report_tsv(
        [rep], [rep], cfg, path=out,
    ) is False
    assert not out.exists()
