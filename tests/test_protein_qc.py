"""Protein-annotation QC: batched GenBank fetch, caching, filtering, TSV output."""
from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from repseq.models import QCReport, SequenceSource
from repseq.output.report import write_isolate_proteins_tsv
from repseq.qc.protein_qc import (
    attach_proteins,
    filter_by_protein_count,
    run_protein_qc,
)
from repseq.taxonomy.cache import TaxonomyCache
from repseq.taxonomy.ncbi import NCBITaxonomy


# ---------------------------------------------------------------------------
# Minimal GenBank record (single CDS) for parser tests
# ---------------------------------------------------------------------------

_GB_TWO_CDS = """\
LOCUS       MW626064                1234 bp    cRNA    linear   VRL 01-JAN-2020
DEFINITION  Influenza A virus segment 7.
ACCESSION   MW626064
VERSION     MW626064.1
KEYWORDS    .
SOURCE      Influenza A virus
FEATURES             Location/Qualifiers
     source          1..1234
                     /organism="Influenza A virus"
     CDS             1..759
                     /protein_id="QQQ12345.1"
                     /product="matrix protein 1"
                     /translation="MSLLTEVETYVLSIVPSGPLKAEIAQRLEDVFAGKNTDLEALMEW"
     CDS             1..294
                     /protein_id="QQQ12346.1"
                     /product="matrix protein 2"
                     /translation="MSLLTEVETPIRNEWGCRCNDSSDPLVVAASIIGILHLILWILDR"
ORIGIN
        1 atgagccttc taaccgaggt cgaaacgtac gttctctcta tcgtcccgtc aggccccctc
//
"""

_GB_NO_CDS = """\
LOCUS       XX000001                 100 bp    DNA     linear   VRL 01-JAN-2020
DEFINITION  Fragment with no CDS.
ACCESSION   XX000001
VERSION     XX000001.1
FEATURES             Location/Qualifiers
     source          1..100
                     /organism="Test virus"
ORIGIN
        1 atgcatgcat gcatgcatgc atgcatgcat gcatgcatgc atgcatgcat gcatgcatgc
//
"""


# ---------------------------------------------------------------------------
# fetch_proteins_batch
# ---------------------------------------------------------------------------

def test_fetch_proteins_batch_parses_cds_features(tmp_cache_dir, monkeypatch):
    """Parser test: feed a fake efetch response, verify CDS features extracted."""
    cache = TaxonomyCache(tmp_cache_dir)
    ncbi = NCBITaxonomy(cache, email="t@e.com")

    fake_resp = MagicMock()
    fake_resp.text = _GB_TWO_CDS
    fake_resp.raise_for_status = lambda: None

    with patch("repseq.taxonomy.ncbi.requests.get", return_value=fake_resp) as mock_get:
        out = ncbi.fetch_proteins_batch(["MW626064.1"])

    mock_get.assert_called_once()
    # Inspect the request actually made
    _, kwargs = mock_get.call_args
    assert kwargs["params"]["db"] == "nuccore"
    assert kwargs["params"]["rettype"] == "gb"
    assert kwargs["params"]["id"] == "MW626064.1"

    proteins = out["MW626064.1"]
    assert len(proteins) == 2
    assert proteins[0]["protein_id"] == "QQQ12345.1"
    assert proteins[0]["product"] == "matrix protein 1"
    assert proteins[0]["length"] == 45  # length of the /translation= string
    assert proteins[1]["product"] == "matrix protein 2"


def test_fetch_proteins_batch_uses_cache_on_second_call(tmp_cache_dir):
    cache = TaxonomyCache(tmp_cache_dir)
    ncbi = NCBITaxonomy(cache, email="t@e.com")

    fake_resp = MagicMock()
    fake_resp.text = _GB_TWO_CDS
    fake_resp.raise_for_status = lambda: None

    with patch("repseq.taxonomy.ncbi.requests.get", return_value=fake_resp) as mock_get:
        ncbi.fetch_proteins_batch(["MW626064.1"])
        # Second call should hit the cache, not the network
        out2 = ncbi.fetch_proteins_batch(["MW626064.1"])

    assert mock_get.call_count == 1
    assert len(out2["MW626064.1"]) == 2


def test_fetch_proteins_batch_chunks_large_inputs(tmp_cache_dir):
    """Batched: accessions are split into chunks of batch_size."""
    cache = TaxonomyCache(tmp_cache_dir)
    ncbi = NCBITaxonomy(cache, email="t@e.com")

    accessions = [f"FAKE{i:03d}.1" for i in range(450)]
    fake_resp = MagicMock()
    fake_resp.text = ""  # empty response → all accessions get []
    fake_resp.raise_for_status = lambda: None

    with patch("repseq.taxonomy.ncbi.requests.get", return_value=fake_resp) as mock_get:
        out = ncbi.fetch_proteins_batch(accessions, batch_size=200)

    # 450 with batch_size=200 → 3 batches (200, 200, 50)
    assert mock_get.call_count == 3
    assert len(out) == 450
    assert all(out[a] == [] for a in accessions)


def test_fetch_proteins_batch_records_with_no_cds(tmp_cache_dir):
    cache = TaxonomyCache(tmp_cache_dir)
    ncbi = NCBITaxonomy(cache, email="t@e.com")
    fake_resp = MagicMock()
    fake_resp.text = _GB_NO_CDS
    fake_resp.raise_for_status = lambda: None

    with patch("repseq.taxonomy.ncbi.requests.get", return_value=fake_resp):
        out = ncbi.fetch_proteins_batch(["XX000001.1"])
    assert out["XX000001.1"] == []


def test_fetch_proteins_batch_network_failure_caches_empty(tmp_cache_dir):
    """A failed fetch shouldn't crash; missing accessions cache as []."""
    cache = TaxonomyCache(tmp_cache_dir)
    ncbi = NCBITaxonomy(cache, email="t@e.com")

    with patch(
        "repseq.taxonomy.ncbi.requests.get",
        side_effect=RuntimeError("network down"),
    ):
        out = ncbi.fetch_proteins_batch(["A.1", "B.1"])

    assert out == {"A.1": [], "B.1": []}


# ---------------------------------------------------------------------------
# attach_proteins
# ---------------------------------------------------------------------------

def test_attach_proteins_skips_uniprot_sequences(make_seq):
    """UniProt records are proteins themselves — they have no GenBank CDS table."""
    ncbi = MagicMock()
    ncbi.fetch_proteins_batch.return_value = {"NC_001.1": [{"protein_id": "X"}]}

    uniprot_seq = make_seq("p1", "MEEP", source=SequenceSource.UNIPROT, accession="P12345")
    ncbi_seq = make_seq("n1", "ACGT", source=SequenceSource.NCBI, accession="NC_001.1")
    attach_proteins([uniprot_seq, ncbi_seq], ncbi)

    assert uniprot_seq.proteins is None  # untouched
    assert ncbi_seq.proteins == [{"protein_id": "X"}]
    # NCBI fetcher only received the ncbi accession
    ncbi.fetch_proteins_batch.assert_called_once_with(["NC_001.1"])


def test_attach_proteins_noop_when_no_accessions(make_seq):
    ncbi = MagicMock()
    s = make_seq("p1", "MEEP", source=SequenceSource.UNIPROT)
    attach_proteins([s], ncbi)
    ncbi.fetch_proteins_batch.assert_not_called()


# ---------------------------------------------------------------------------
# filter_by_protein_count
# ---------------------------------------------------------------------------

def test_filter_min_proteins_drops_under_annotated(make_seq):
    a = make_seq("a", "ACGT")
    a.proteins = []  # no CDS
    b = make_seq("b", "ACGT")
    b.proteins = [{"protein_id": "X"}]
    report = QCReport()

    kept = filter_by_protein_count(
        [a, b],
        qc_cfg={"protein_annotation": {"enabled": True, "min_proteins": 1}},
        virus_cfg=None,
        report=report,
    )
    assert [s.id for s in kept] == ["b"]
    assert report.removed_proteins == 1
    assert a.qc_passed is False
    assert "protein_count_below_min" in a.qc_fail_reason


def test_filter_per_segment_expected_count(make_seq):
    """Segment-specific expected count: drop segments that don't match."""
    # HA segment expects exactly 1 protein
    ha_good = make_seq("ha1", "ACGT", segment="HA")
    ha_good.proteins = [{"protein_id": "HA1"}]
    ha_bad = make_seq("ha2", "ACGT", segment="HA")
    ha_bad.proteins = []  # 0 != 1 expected → fail
    # M segment expects exactly 2 proteins
    m_good = make_seq("m1", "ACGT", segment="M")
    m_good.proteins = [{"protein_id": "M1"}, {"protein_id": "M2"}]
    m_bad = make_seq("m2", "ACGT", segment="M")
    m_bad.proteins = [{"protein_id": "M1_only"}]  # 1 != 2 expected → fail

    virus_cfg = {
        "segments": ["HA", "M"],
        "expected_proteins_per_segment": {"HA": 1, "M": 2},
    }
    report = QCReport()
    kept = filter_by_protein_count(
        [ha_good, ha_bad, m_good, m_bad],
        qc_cfg={},
        virus_cfg=virus_cfg,
        report=report,
    )
    assert {s.id for s in kept} == {"ha1", "m1"}
    assert report.removed_proteins == 2


def test_filter_passes_through_when_proteins_none(make_seq):
    """Sequences without protein data (e.g. UniProt) shouldn't be filtered."""
    s = make_seq("a", "MEEP")
    s.proteins = None
    report = QCReport()
    kept = filter_by_protein_count(
        [s],
        qc_cfg={"protein_annotation": {"enabled": True, "min_proteins": 1}},
        virus_cfg=None,
        report=report,
    )
    assert kept == [s]
    assert report.removed_proteins == 0


def test_filter_noop_when_nothing_configured(make_seq):
    s = make_seq("a", "ACGT")
    s.proteins = []
    report = QCReport()
    kept = filter_by_protein_count([s], qc_cfg={}, virus_cfg=None, report=report)
    assert kept == [s]
    assert report.removed_proteins == 0


# ---------------------------------------------------------------------------
# run_protein_qc (integration of attach + filter)
# ---------------------------------------------------------------------------

def test_run_protein_qc_disabled_is_noop(make_seq):
    ncbi = MagicMock()
    s = make_seq("a", "ACGT", source=SequenceSource.NCBI)
    report = QCReport()
    kept = run_protein_qc([s], ncbi, cfg={"qc": {}}, virus_cfg=None, report=report)
    assert kept == [s]
    ncbi.fetch_proteins_batch.assert_not_called()


def test_run_protein_qc_end_to_end_min_proteins(make_seq):
    ncbi = MagicMock()
    ncbi.fetch_proteins_batch.return_value = {
        "A.1": [],
        "B.1": [{"protein_id": "P1"}],
    }
    a = make_seq("a", "ACGT", source=SequenceSource.NCBI, accession="A.1")
    b = make_seq("b", "ACGT", source=SequenceSource.NCBI, accession="B.1")
    report = QCReport()
    cfg = {"qc": {"protein_annotation": {"enabled": True, "min_proteins": 1}}}
    kept = run_protein_qc([a, b], ncbi, cfg=cfg, virus_cfg=None, report=report)
    assert [s.id for s in kept] == ["b"]
    assert report.removed_proteins == 1


# ---------------------------------------------------------------------------
# Output TSV
# ---------------------------------------------------------------------------

def test_write_isolate_proteins_tsv(tmp_path: Path, make_seq):
    s1 = make_seq("s1", "ACGT", segment="HA", accession="NC_001.1")
    s1.proteins = [{"protein_id": "HA_P1", "product": "hemagglutinin", "length": 566}]
    s2 = make_seq("s2", "ACGT", segment="NA", accession="NC_002.1")
    s2.proteins = [{"protein_id": "NA_P1", "product": "neuraminidase", "length": 469}]
    complete_isolates = {"A/duck/HK/1/97": [s1, s2]}

    path = tmp_path / "iso_proteins.tsv"
    wrote = write_isolate_proteins_tsv(complete_isolates, path)
    assert wrote is True

    lines = path.read_text().strip().splitlines()
    assert lines[0] == "isolate_id\tsegment\taccession\tprotein_id\tproduct\tlength"
    assert len(lines) == 3  # header + 2 protein rows
    assert "HA_P1\themagglutinin\t566" in lines[1]
    assert "NA_P1\tneuraminidase\t469" in lines[2]


def test_write_isolate_proteins_tsv_skips_when_no_proteins(tmp_path: Path, make_seq):
    s = make_seq("s1", "ACGT", segment="HA")
    s.proteins = None
    path = tmp_path / "iso_proteins.tsv"
    wrote = write_isolate_proteins_tsv({"iso1": [s]}, path)
    assert wrote is False
    assert not path.exists()
