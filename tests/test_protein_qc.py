"""Protein-annotation QC: batched GenBank fetch, caching, filtering, TSV output."""
from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from repseq.models import Cluster, QCReport, RunResult, SequenceSource, TaxonomyInfo
from repseq.output.report import write_isolate_proteins_tsv, write_proteins_fasta
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

_GB_WITH_SOURCE_QUALIFIERS = """\
LOCUS       MW626064                1234 bp    cRNA    linear   VRL 01-JAN-2020
DEFINITION  Hantavirus L segment.
ACCESSION   MW626064
VERSION     MW626064.1
FEATURES             Location/Qualifiers
     source          1..1234
                     /organism="Sin Nombre virus"
                     /isolate="SNV-NM-H10"
                     /strain="Convict Creek 107"
                     /segment="L"
     CDS             1..759
                     /protein_id="QQQ12345.1"
                     /product="RNA-dependent RNA polymerase"
                     /translation="MSLLTEVETYVLSIVPSGPLKAEIAQRLEDVFAGKNTDLEALMEW"
ORIGIN
        1 atgagccttc taaccgaggt cgaaacgtac gttctctcta tcgtcccgtc aggccccctc
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
    assert proteins[0]["sequence"].startswith("MSLLTEVET")  # translation captured
    assert proteins[1]["product"] == "matrix protein 2"
    assert proteins[1]["sequence"].startswith("MSLLTEVETPIRNEW")


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
# fetch_source_metadata_batch — isolate/strain/segment from /source feature
# ---------------------------------------------------------------------------

def test_fetch_source_metadata_batch_extracts_qualifiers(tmp_cache_dir):
    """Source feature carries /isolate, /strain, /segment — all captured."""
    cache = TaxonomyCache(tmp_cache_dir)
    ncbi = NCBITaxonomy(cache, email="t@e.com")

    fake_resp = MagicMock()
    fake_resp.text = _GB_WITH_SOURCE_QUALIFIERS
    fake_resp.raise_for_status = lambda: None

    with patch("repseq.taxonomy.ncbi.requests.get", return_value=fake_resp):
        out = ncbi.fetch_source_metadata_batch(["MW626064.1"])

    meta = out["MW626064.1"]
    assert meta["isolate"] == "SNV-NM-H10"
    assert meta["strain"] == "Convict Creek 107"
    assert meta["segment"] == "L"


def test_fetch_source_metadata_batch_handles_missing_qualifiers(tmp_cache_dir):
    """A source feature with only /organism returns all-None for our keys."""
    cache = TaxonomyCache(tmp_cache_dir)
    ncbi = NCBITaxonomy(cache, email="t@e.com")

    fake_resp = MagicMock()
    fake_resp.text = _GB_NO_CDS
    fake_resp.raise_for_status = lambda: None

    with patch("repseq.taxonomy.ncbi.requests.get", return_value=fake_resp):
        out = ncbi.fetch_source_metadata_batch(["XX000001.1"])

    assert out["XX000001.1"] == {"isolate": None, "strain": None, "segment": None}


def test_fetch_source_metadata_batch_shares_cache_with_proteins(tmp_cache_dir):
    """One efetch populates both protein list and source metadata."""
    cache = TaxonomyCache(tmp_cache_dir)
    ncbi = NCBITaxonomy(cache, email="t@e.com")

    fake_resp = MagicMock()
    fake_resp.text = _GB_WITH_SOURCE_QUALIFIERS
    fake_resp.raise_for_status = lambda: None

    with patch(
        "repseq.taxonomy.ncbi.requests.get", return_value=fake_resp,
    ) as mock_get:
        # First call: hits network.
        proteins = ncbi.fetch_proteins_batch(["MW626064.1"])
        # Second call for source metadata on the same accession: cache hit.
        meta = ncbi.fetch_source_metadata_batch(["MW626064.1"])

    assert mock_get.call_count == 1
    assert len(proteins["MW626064.1"]) == 1
    assert meta["MW626064.1"]["segment"] == "L"


def test_fetch_source_metadata_batch_legacy_cache_returns_none(tmp_cache_dir):
    """Old cache entries (proteins only, no 'source' key) yield all-None."""
    from repseq.taxonomy.cache import TaxonomyCache as _TC
    cache = _TC(tmp_cache_dir)
    cache.set("ncbi_proteins", "ZZ000001.1", {"proteins": []})  # v0.5.9 shape
    ncbi = NCBITaxonomy(cache, email="t@e.com")

    # No network call should be made — entry is "cached", just incomplete.
    with patch(
        "repseq.taxonomy.ncbi.requests.get",
        side_effect=AssertionError("should not hit network"),
    ):
        out = ncbi.fetch_source_metadata_batch(["ZZ000001.1"])

    assert out["ZZ000001.1"] == {"isolate": None, "strain": None, "segment": None}


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
    # NCBI fetcher only received the ncbi accession.
    ncbi.fetch_proteins_batch.assert_called_once()
    args, _ = ncbi.fetch_proteins_batch.call_args
    assert args[0] == ["NC_001.1"]


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


def test_filter_per_segment_accepts_list_of_counts(make_seq):
    """List-valued expected count: any value in the list is acceptable."""
    one_protein = make_seq("a", "ACGT", segment="PB1")
    one_protein.proteins = [{"protein_id": "PB1"}]
    two_proteins = make_seq("b", "ACGT", segment="PB1")
    two_proteins.proteins = [{"protein_id": "PB1"}, {"protein_id": "PB1-F2"}]
    three_proteins = make_seq("c", "ACGT", segment="PB1")
    three_proteins.proteins = [{}, {}, {}]  # 3 not in [1, 2]

    virus_cfg = {
        "segments": ["PB1"],
        "expected_proteins_per_segment": {"PB1": [1, 2]},
    }
    report = QCReport()
    kept = filter_by_protein_count(
        [one_protein, two_proteins, three_proteins],
        qc_cfg={}, virus_cfg=virus_cfg, report=report,
    )
    assert {s.id for s in kept} == {"a", "b"}
    assert report.removed_proteins == 1
    assert "expected_one_of=[1, 2]" in three_proteins.qc_fail_reason


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
    tax = TaxonomyInfo(
        species="Influenza A virus",
        genus="Alphainfluenzavirus",
        family="Orthomyxoviridae",
        order="Articulavirales",
        class_="Insthoviricetes",
        lineage={
            "species": "Influenza A virus",
            "genus": "Alphainfluenzavirus",
            "family": "Orthomyxoviridae",
            "order": "Articulavirales",
            "class": "Insthoviricetes",
        },
    )
    s1 = make_seq("s1", "ACGT" * 400, segment="HA", accession="NC_001.1", taxonomy=tax)
    s1.proteins = [{"protein_id": "HA_P1", "product": "hemagglutinin", "length": 566}]
    s2 = make_seq("s2", "ACGT" * 350, segment="NA", accession="NC_002.1", taxonomy=tax)
    s2.proteins = [{"protein_id": "NA_P1", "product": "neuraminidase", "length": 469}]
    complete_isolates = {"A/duck/HK/1/97": [s1, s2]}

    path = tmp_path / "iso_proteins.tsv"
    wrote = write_isolate_proteins_tsv(complete_isolates, path)
    assert wrote is True

    lines = path.read_text().strip().splitlines()
    assert lines[0] == (
        "protein_id\tproduct\tlength\tisolate_id\tsegment\tsegment_length\t"
        "accession\tspecies\tsubgenus\tgenus\tsubfamily\tfamily\tsuborder\t"
        "order\tsubclass\tclass"
    )
    assert len(lines) == 3  # header + 2 protein rows

    row1 = lines[1].split("\t")
    assert row1[0] == "HA_P1"
    assert row1[1] == "hemagglutinin"
    assert row1[2] == "566"
    assert row1[3] == "A/duck/HK/1/97"
    assert row1[4] == "HA"
    assert row1[5] == "1600"           # len("ACGT" * 400)
    assert row1[6] == "NC_001.1"
    assert row1[7] == "Influenza A virus"
    assert row1[8] == ""                # subgenus (absent)
    assert row1[9] == "Alphainfluenzavirus"
    assert row1[10] == ""               # subfamily (absent)
    assert row1[11] == "Orthomyxoviridae"
    assert row1[12] == ""               # suborder (absent)
    assert row1[13] == "Articulavirales"
    assert row1[14] == ""               # subclass (absent)
    assert row1[15] == "Insthoviricetes"

    row2 = lines[2].split("\t")
    assert row2[0] == "NA_P1"
    assert row2[4] == "NA"
    assert row2[5] == "1400"            # len("ACGT" * 350)


def test_write_isolate_proteins_tsv_emits_sub_ranks_from_lineage(
    tmp_path: Path, make_seq
):
    """Sub-ranks (subgenus/subfamily/suborder/subclass) come only via the
    lineage map — they have no standard TaxonomyInfo field. Confirm they
    round-trip into the output."""
    tax = TaxonomyInfo(
        species="Schmallenberg virus",
        genus="Orthobunyavirus",
        family="Peribunyaviridae",
        order="Elliovirales",
        class_="Ellioviricetes",
        lineage={
            "species": "Schmallenberg virus",
            "subgenus": "Simbu serogroup",
            "genus": "Orthobunyavirus",
            "subfamily": "Bunyavirinae",
            "family": "Peribunyaviridae",
            "suborder": "Bunyavirales-suborder",
            "order": "Elliovirales",
            "subclass": "Some-subclass",
            "class": "Ellioviricetes",
        },
    )
    s = make_seq("s1", "ACGT", segment="L", accession="ACC.1", taxonomy=tax)
    s.proteins = [{"protein_id": "P1", "product": "L protein", "length": 2200}]

    path = tmp_path / "iso_proteins.tsv"
    assert write_isolate_proteins_tsv({"ISO1": [s]}, path) is True
    row = path.read_text().strip().splitlines()[1].split("\t")
    assert row[8] == "Simbu serogroup"           # subgenus
    assert row[10] == "Bunyavirinae"             # subfamily
    assert row[12] == "Bunyavirales-suborder"    # suborder
    assert row[14] == "Some-subclass"            # subclass


def test_write_isolate_proteins_tsv_no_taxonomy_leaves_rank_cells_blank(
    tmp_path: Path, make_seq
):
    s = make_seq("s1", "ACGT", segment="HA", accession="ACC.1", taxonomy=None)
    s.proteins = [{"protein_id": "P1", "product": "x", "length": 10}]
    path = tmp_path / "iso_proteins.tsv"
    assert write_isolate_proteins_tsv({"ISO1": [s]}, path) is True
    # Use splitlines() without strip() — trailing empty TSV cells matter
    # and would otherwise be eaten as trailing whitespace on the last row.
    rows = path.read_text().splitlines()
    row = rows[1].split("\t")
    assert row[0] == "P1"
    assert row[5] == "4"   # segment_length
    # All 9 taxonomy cells (indices 7..15) are blank
    assert row[7:] == [""] * 9


def test_write_proteins_fasta_segmented(tmp_path: Path, make_seq):
    """Segmented path: emits proteins from segments of represented isolates."""
    # Build segments with translations
    s1 = make_seq("s1", "ACGT", segment="HA", accession="NC_001.1", isolate_id="ISO1")
    s1.proteins = [{
        "protein_id": "YP_HA1",
        "product": "hemagglutinin",
        "length": 5,
        "sequence": "MAKLM",
    }]
    s2 = make_seq("s2", "ACGT", segment="M", accession="NC_002.1", isolate_id="ISO1")
    s2.proteins = [
        {"protein_id": "YP_M1", "product": "matrix protein 1",
         "length": 4, "sequence": "MAAA"},
        {"protein_id": "YP_M2", "product": "matrix protein 2",
         "length": 4, "sequence": "MBBB"},
    ]
    # Another isolate that was NOT selected — its proteins must be excluded
    s3 = make_seq("s3", "ACGT", segment="HA", accession="NC_003.1", isolate_id="ISO2")
    s3.proteins = [{
        "protein_id": "YP_HA_other", "product": "hemagglutinin",
        "length": 5, "sequence": "MSKIP",
    }]
    complete_isolates = {"ISO1": [s1, s2], "ISO2": [s3]}

    # The selected representative is the CONCAT for ISO1 only
    concat = make_seq("concat", "ACGTACGT", accession="ISO1", isolate_id="ISO1")
    concat.id = "CONCAT|ISO1"
    result = RunResult(
        mode="global:threshold",
        representatives=[concat],
        clusters=[Cluster(cluster_id="c1", representative=concat)],
    )

    out = tmp_path / "proteins.fasta"
    wrote = write_proteins_fasta(result, complete_isolates, out)
    assert wrote is True

    body = out.read_text()
    # ISO1 proteins present
    assert ">YP_HA1 hemagglutinin" in body
    assert ">YP_M1 matrix protein 1" in body
    assert ">YP_M2 matrix protein 2" in body
    # Tags carry isolate, segment, parent accession
    assert "[isolate=ISO1]" in body
    assert "[segment=HA]" in body
    assert "[parent=NC_001.1]" in body
    # ISO2 (not selected) must NOT appear
    assert "YP_HA_other" not in body


def test_write_proteins_fasta_non_segmented(tmp_path: Path, make_seq):
    """Non-segmented path: proteins from result.representatives directly."""
    rep = make_seq("r1", "MEEP", accession="P12345")
    rep.proteins = [{
        "protein_id": "P12345", "product": "test protein",
        "length": 4, "sequence": "MEEP",
    }]
    result = RunResult(
        mode="global:count",
        representatives=[rep],
        clusters=[Cluster(cluster_id="c1", representative=rep)],
    )

    out = tmp_path / "proteins.fasta"
    wrote = write_proteins_fasta(result, complete_isolates=None, path=out)
    assert wrote is True
    body = out.read_text()
    assert ">P12345 test protein" in body
    # No isolate tag in non-segmented mode
    assert "[isolate=" not in body
    assert "[parent=P12345]" in body
    assert "MEEP" in body


def test_write_proteins_fasta_skips_when_no_translations(tmp_path: Path, make_seq):
    """If proteins lack 'sequence' (older cache), the file is not written."""
    rep = make_seq("r1", "ACGT", accession="NC_001.1")
    rep.proteins = [{"protein_id": "X", "product": "y", "length": 4}]  # no sequence
    result = RunResult(mode="x", representatives=[rep], clusters=[])
    out = tmp_path / "proteins.fasta"
    wrote = write_proteins_fasta(result, complete_isolates=None, path=out)
    assert wrote is False
    assert not out.exists()


def test_write_proteins_fasta_wraps_at_70_chars(tmp_path: Path, make_seq):
    rep = make_seq("r1", "X", accession="P1")
    long_seq = "M" + ("A" * 200)
    rep.proteins = [{"protein_id": "P1", "product": "x", "length": 201, "sequence": long_seq}]
    result = RunResult(mode="x", representatives=[rep], clusters=[])
    out = tmp_path / "p.fasta"
    write_proteins_fasta(result, complete_isolates=None, path=out)
    body_lines = [ln for ln in out.read_text().splitlines() if not ln.startswith(">")]
    assert all(len(ln) <= 70 for ln in body_lines)


def test_write_isolate_proteins_tsv_skips_when_no_proteins(tmp_path: Path, make_seq):
    s = make_seq("s1", "ACGT", segment="HA")
    s.proteins = None
    path = tmp_path / "iso_proteins.tsv"
    wrote = write_isolate_proteins_tsv({"iso1": [s]}, path)
    assert wrote is False
    assert not path.exists()
