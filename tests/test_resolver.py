"""MetadataResolver: strain-label fallback and failure tracking."""
from __future__ import annotations

from unittest.mock import MagicMock

from repseq.models import SequenceSource, TaxonomyInfo
from repseq.taxonomy.cache import TaxonomyCache
from repseq.taxonomy.resolver import MetadataResolver, _parse_strain_label


# ---------------------------------------------------------------------------
# Strain-label parsing
# ---------------------------------------------------------------------------

def test_parse_strain_label_influenza_extracts_host_country_year(make_seq):
    s = make_seq(
        "a", "ACGT",
        header="A/duck/Hong_Kong/1/1997(H5N1) hemagglutinin",
    )
    s.description = s.header
    out = _parse_strain_label(s)
    assert out.get("host") == "duck"
    assert out.get("country") == "Hong Kong"
    assert out.get("collection_date") == "1997"


def test_parse_strain_label_generic_year_only(make_seq):
    s = make_seq("a", "ACGT")
    s.description = "Some virus isolate from 2014 outbreak"
    out = _parse_strain_label(s)
    assert out == {"collection_date": "2014"}


def test_parse_strain_label_human_normalised_to_homo_sapiens(make_seq):
    s = make_seq("a", "ACGT", header="A/human/USA/1/2020(H1N1)")
    s.description = s.header
    out = _parse_strain_label(s)
    assert out.get("host") == "Homo sapiens"


def test_parse_strain_label_does_not_overwrite_existing(make_seq):
    s = make_seq("a", "ACGT", header="A/duck/Hong_Kong/1/1997(H5N1)",
                 host="Existing host", country="Existing country")
    s.description = s.header
    out = _parse_strain_label(s)
    # Resolver only adds keys that aren't set on the seq
    assert "host" not in out
    assert "country" not in out


# ---------------------------------------------------------------------------
# Resolve falls back to strain label when DB has nothing
# ---------------------------------------------------------------------------

def test_resolve_uses_strain_label_when_db_returns_none(tmp_cache_dir, make_seq):
    cache = TaxonomyCache(tmp_cache_dir)
    ncbi = MagicMock()
    ncbi.fetch_accession_metadata.return_value = None
    uniprot = MagicMock()
    uniprot.fetch_entry.return_value = None

    resolver = MetadataResolver(cache, ncbi, uniprot, threads=1)
    s = make_seq(
        "x", "ACGT",
        header="A/chicken/Vietnam/1/2005(H5N1)",
        source=SequenceSource.NCBI,
    )
    s.description = s.header
    resolver.resolve(s)

    assert s.host == "chicken"
    assert s.country == "Vietnam"
    assert s.collection_date == "2005"


def test_resolve_populates_taxonomy_from_db(tmp_cache_dir, make_seq):
    cache = TaxonomyCache(tmp_cache_dir)
    ncbi = MagicMock()
    ncbi.fetch_accession_metadata.return_value = {
        "organism": "Influenza A virus",
        "taxid": 11320,
        "lineage": {"genus": "Alphainfluenzavirus", "family": "Orthomyxoviridae"},
    }
    uniprot = MagicMock()

    resolver = MetadataResolver(cache, ncbi, uniprot, threads=1)
    s = make_seq("y", "ACGT", source=SequenceSource.NCBI)
    resolver.resolve(s)

    assert s.organism == "Influenza A virus"
    assert isinstance(s.taxonomy, TaxonomyInfo)
    assert s.taxonomy.genus == "Alphainfluenzavirus"
    assert s.taxonomy.family == "Orthomyxoviridae"


# ---------------------------------------------------------------------------
# Failure tracking — the bug we just fixed
# ---------------------------------------------------------------------------

def test_resolve_batch_records_failures_instead_of_swallowing(tmp_cache_dir, make_seq, caplog):
    cache = TaxonomyCache(tmp_cache_dir)
    ncbi = MagicMock()

    def explode(_acc):
        raise RuntimeError("simulated network failure")

    ncbi.fetch_accession_metadata.side_effect = explode
    uniprot = MagicMock()

    resolver = MetadataResolver(cache, ncbi, uniprot, threads=2)
    seqs = [
        make_seq("s1", "ACGT", source=SequenceSource.NCBI),
        make_seq("s2", "ACGT", source=SequenceSource.NCBI),
    ]
    with caplog.at_level("WARNING", logger="repseq.taxonomy.resolver"):
        resolver.resolve_batch(seqs)

    # Failures were recorded (not silently swallowed)
    assert len(resolver.failures) == 2
    failed_ids = {sid for sid, _ in resolver.failures}
    assert failed_ids == {"s1", "s2"}
    assert all("simulated network failure" in msg for _, msg in resolver.failures)
    # And surfaced via logging
    assert any("Resolution failed" in r.message for r in caplog.records)


def test_resolve_batch_success_does_not_record_failures(tmp_cache_dir, make_seq):
    cache = TaxonomyCache(tmp_cache_dir)
    ncbi = MagicMock()
    ncbi.fetch_accession_metadata.return_value = {"organism": "Test"}
    uniprot = MagicMock()

    resolver = MetadataResolver(cache, ncbi, uniprot, threads=2)
    seqs = [make_seq(f"s{i}", "ACGT", source=SequenceSource.NCBI) for i in range(3)]
    resolver.resolve_batch(seqs)
    assert resolver.failures == []
