"""Tests for the v0.13 ``hmmscan`` TSV column + matching FASTA tag.

Covers the format (``Name(E=val,cov=val);...``), ordering (best E-value
first), the passing-only filter, the FASTA-header tag, and the
back-compat behaviour when no protein carries ``hmm_hits``.
"""
from __future__ import annotations

from pathlib import Path

from repseq.models import Cluster, RunResult, TaxonomyInfo
from repseq.output.report import (
    _format_hmmscan_cell,
    write_isolate_proteins_tsv,
    write_proteins_fasta,
)


def _hit(target, *, dom_evalue, ali_span, hmm_len=300, passing=True):
    # Coverage is now measured on the HMM model span (hmm_to - hmm_from + 1),
    # so express the covered span via hmm coords. ali_span is kept (it's a
    # real field used for domain ordering) but no longer drives coverage.
    return {
        "target": target,
        "dom_evalue": dom_evalue,
        "dom_score": 200.0,
        "ali_span": ali_span,
        "hmm_from": 1,
        "hmm_to": ali_span,
        "hmm_len": hmm_len,
        "passing": passing,
    }


def test_format_hmmscan_cell_empty_when_no_hits():
    assert _format_hmmscan_cell({"product": "x"}) == ""
    assert _format_hmmscan_cell({"hmm_hits": []}) == ""


def test_format_hmmscan_cell_excludes_failing_hits():
    """v0.13 design: only show passing hits — failing hits are internal
    diagnostics, not part of the user-facing summary."""
    prot = {"hmm_hits": [
        _hit("RdRP_4", dom_evalue=0.5, ali_span=20, passing=False),
        _hit("Helicase", dom_evalue=1e-30, ali_span=280, passing=True),
    ]}
    out = _format_hmmscan_cell(prot)
    assert "RdRP_4" not in out
    assert out.startswith("Helicase(")


def test_format_hmmscan_cell_orders_by_best_evalue_first():
    prot = {"hmm_hits": [
        _hit("Z_weak", dom_evalue=1e-10, ali_span=200),
        _hit("A_strong", dom_evalue=1e-50, ali_span=290),
        _hit("M_mid", dom_evalue=1e-30, ali_span=250),
    ]}
    out = _format_hmmscan_cell(prot)
    assert out.startswith("A_strong(")
    # Semicolon-separated, three entries
    assert out.count(";") == 2
    # Middle entry second, weakest last
    parts = out.split(";")
    assert parts[1].startswith("M_mid(")
    assert parts[2].startswith("Z_weak(")


def test_format_hmmscan_cell_renders_coverage_two_decimal():
    prot = {"hmm_hits": [_hit("X", dom_evalue=1e-30, ali_span=150, hmm_len=300)]}
    out = _format_hmmscan_cell(prot)
    # 150/300 = 0.50
    assert "cov=0.50" in out
    assert out.startswith("X(E=1e-30,")


def test_format_hmmscan_cell_format_string_matches_spec():
    """Pinned to the user-confirmed format:
    ``Bunya_Gn(E=1e-30,cov=0.85);Bunya_Gc(E=3e-25,cov=0.78)``"""
    prot = {"hmm_hits": [
        _hit("Bunya_Gn", dom_evalue=1e-30, ali_span=170, hmm_len=200),
        _hit("Bunya_Gc", dom_evalue=3e-25, ali_span=156, hmm_len=200),
    ]}
    out = _format_hmmscan_cell(prot)
    assert out == "Bunya_Gn(E=1e-30,cov=0.85);Bunya_Gc(E=3e-25,cov=0.78)"


def test_isolate_proteins_tsv_emits_hmmscan_values(tmp_path: Path, make_seq):
    """The column carries hit data when proteins have passing hits."""
    s = make_seq("s1", "ACGT", segment="L", accession="ACC.1")
    s.proteins = [{
        "protein_id": "P_pol",
        "product": "polymerase",
        "length": 2200,
        "hmm_hits": [_hit("RdRP_4", dom_evalue=1e-50, ali_span=280, hmm_len=300)],
    }]
    path = tmp_path / "iso.tsv"
    assert write_isolate_proteins_tsv({"ISO1": [s]}, path) is True
    lines = path.read_text().splitlines()
    assert "\thmmscan\t" in lines[0]
    row = lines[1].split("\t")
    # hmmscan is the column right of representative (index 9 — see TSV
    # schema in repseq/output/report.py).
    assert row[9].startswith("RdRP_4(E=1e-50,cov=")


def test_isolate_proteins_tsv_blank_hmmscan_when_no_hits(tmp_path: Path, make_seq):
    s = make_seq("s1", "ACGT", segment="HA", accession="ACC.1")
    s.proteins = [{"protein_id": "P", "product": "x", "length": 100}]
    path = tmp_path / "iso.tsv"
    write_isolate_proteins_tsv({"ISO1": [s]}, path)
    row = path.read_text().splitlines()[1].split("\t")
    assert row[9] == ""  # hmmscan column blank — back-compat


def test_proteins_fasta_emits_hmmscan_tag(tmp_path: Path, make_seq):
    """FASTA header tag matches the TSV column verbatim."""
    s = make_seq("s1", "ACGT", segment="L", accession="ACC.1",
                 organism="Test virus")
    s.proteins = [{
        "protein_id": "P_pol",
        "product": "polymerase",
        "sequence": "M" * 50,
        "length": 50,
        "hmm_hits": [_hit("RdRP_4", dom_evalue=1e-50, ali_span=280, hmm_len=300)],
    }]
    result = RunResult(
        mode="x", representatives=[s],
        clusters=[Cluster(cluster_id="c1", representative=s)],
    )
    out = tmp_path / "proteins.fasta"
    write_proteins_fasta(result, complete_isolates=None, path=out)
    header = out.read_text().splitlines()[0]
    assert "[hmmscan=RdRP_4(E=1e-50,cov=" in header


def test_proteins_fasta_omits_hmmscan_tag_when_no_passing_hits(
    tmp_path: Path, make_seq
):
    s = make_seq("s1", "ACGT", segment="L", accession="ACC.1")
    s.proteins = [{
        "protein_id": "P_pol",
        "product": "polymerase",
        "sequence": "M" * 50,
        "length": 50,
        "hmm_hits": [_hit("RdRP_4", dom_evalue=1.0, ali_span=20, passing=False)],
    }]
    result = RunResult(
        mode="x", representatives=[s],
        clusters=[Cluster(cluster_id="c1", representative=s)],
    )
    out = tmp_path / "proteins.fasta"
    write_proteins_fasta(result, complete_isolates=None, path=out)
    header = out.read_text().splitlines()[0]
    assert "[hmmscan=" not in header
