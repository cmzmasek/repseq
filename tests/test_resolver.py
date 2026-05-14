"""MetadataResolver: strain-label fallback and failure tracking."""
from __future__ import annotations

from unittest.mock import MagicMock

from repseq.models import SequenceSource, TaxonomyInfo
from repseq.taxonomy.cache import TaxonomyCache
from repseq.taxonomy.ncbi import _looks_like_protein_acc
from repseq.taxonomy.resolver import MetadataResolver, _parse_strain_label


# ---------------------------------------------------------------------------
# Accession-type detection (picks the Entrez db: nuccore vs protein)
# ---------------------------------------------------------------------------

def test_looks_like_protein_acc_refseq_protein_prefixes():
    # NP_/XP_/YP_/WP_/AP_ are RefSeq *protein* accessions.
    for acc in ("NP_040980.1", "XP_001234.2", "YP_009725.1", "WP_000123.1"):
        assert _looks_like_protein_acc(acc) is True


def test_looks_like_protein_acc_refseq_nucleotide_not_protein():
    # Regression: RefSeq *nucleotide* prefixes must NOT be treated as
    # protein, or fetch_accession_metadata queries the wrong Entrez db
    # and resolves no metadata for curated RefSeq genomes.
    for acc in ("NC_026433.1", "NM_000546.6", "NG_007114.1",
                "NR_003286.4", "NW_000001.1", "NZ_CP000001.1", "XM_011.1"):
        assert _looks_like_protein_acc(acc) is False


def test_looks_like_protein_acc_genbank():
    assert _looks_like_protein_acc("AAA12345") is True       # GenBank protein
    assert _looks_like_protein_acc("MW626064.1") is False    # GenBank nucleotide


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


def test_resolve_db_value_overrides_header_parsed_field(tmp_cache_dir, make_seq):
    # H2 regression: the authoritative database value must win over the
    # heuristic header parse. Previously the resolver only filled EMPTY
    # fields, so a wrong header guess blocked the correct DB value.
    cache = TaxonomyCache(tmp_cache_dir)
    ncbi = MagicMock()
    ncbi.fetch_accession_metadata.return_value = {
        "organism": "Influenza A virus",
        "host": "Gallus gallus",
        "country": "Viet Nam",
    }
    uniprot = MagicMock()
    resolver = MetadataResolver(cache, ncbi, uniprot, threads=1)
    # seq carries values a fragile header parse might have produced
    s = make_seq("s1", "ACGT", source=SequenceSource.NCBI,
                 host="H5N1", country="USA")
    resolver.resolve(s)
    assert s.host == "Gallus gallus"
    assert s.country == "Viet Nam"


def test_resolve_keeps_header_value_when_db_field_absent(tmp_cache_dir, make_seq):
    # The DB only overrides fields it actually provides; a header-derived
    # value is retained when the database has nothing for that field.
    cache = TaxonomyCache(tmp_cache_dir)
    ncbi = MagicMock()
    ncbi.fetch_accession_metadata.return_value = {"organism": "Foo virus"}
    uniprot = MagicMock()
    resolver = MetadataResolver(cache, ncbi, uniprot, threads=1)
    s = make_seq("s1", "ACGT", source=SequenceSource.NCBI, host="duck")
    resolver.resolve(s)
    assert s.host == "duck"


# Taxonomy efetch XML — the *only* Entrez endpoint that returns the ranked
# lineage. (Taxonomy esummary carries none: genus/species are blank for
# viruses, which silently grouped every viral sequence as "Unknown".)
_TAXONOMY_XML = """<?xml version="1.0"?>
<TaxaSet><Taxon>
  <TaxId>111</TaxId>
  <ScientificName>Influenza A virus</ScientificName>
  <Rank>species</Rank>
  <LineageEx>
    <Taxon><ScientificName>Orthomyxoviridae</ScientificName><Rank>family</Rank></Taxon>
    <Taxon><ScientificName>Alphainfluenzavirus</ScientificName><Rank>genus</Rank></Taxon>
  </LineageEx>
</Taxon></TaxaSet>"""


def test_fetch_accession_metadata_parses_source_qualifiers(tmp_cache_dir):
    # H1 regression: host/country/collection_date/strain must be harvested
    # from the esummary subtype/subname fields, not left to header parsing.
    from repseq.taxonomy.ncbi import NCBITaxonomy

    cache = TaxonomyCache(tmp_cache_dir)
    ncbi = NCBITaxonomy(cache)

    def fake_get(endpoint, params):
        if endpoint == "esearch.fcgi":
            return {"esearchresult": {"idlist": ["999"]}}
        if endpoint == "esummary.fcgi":
            return {"result": {"999": {
                "organism": "Influenza A virus",
                "taxid": 111,
                "title": "Influenza A virus segment 4",
                "subtype": "strain|host|country|collection_date",
                "subname": "A/duck/Vietnam/1/2005|Gallus gallus|Viet Nam|2005-01",
            }}}
        raise AssertionError(f"unexpected _get({endpoint!r}, {params!r})")

    ncbi._get = fake_get  # type: ignore[assignment]
    # Lineage comes from taxonomy efetch (XML), not esummary.
    ncbi._get_text = lambda endpoint, params: _TAXONOMY_XML  # type: ignore[assignment]
    meta = ncbi.fetch_accession_metadata("MW000001.1")

    assert meta["organism"] == "Influenza A virus"
    assert meta["host"] == "Gallus gallus"
    assert meta["country"] == "Viet Nam"
    assert meta["collection_date"] == "2005-01"
    assert meta["strain"] == "A/duck/Vietnam/1/2005"
    assert meta["lineage"]["genus"] == "Alphainfluenzavirus"
    assert meta["lineage"]["family"] == "Orthomyxoviridae"


def test_parse_taxonomy_xml_extracts_ranked_lineage():
    # Regression: virus lineage must be read from efetch XML <LineageEx>.
    # The taxonomy esummary endpoint returns no lineage at all, so the old
    # esummary-based parse gave every viral sequence an empty lineage and
    # taxonomic modes grouped everything under "Unknown".
    from repseq.taxonomy.ncbi import _parse_taxonomy_xml

    xml_text = """<?xml version="1.0"?>
    <TaxaSet><Taxon>
      <TaxId>159137</TaxId>
      <ScientificName>Yaba-7 virus</ScientificName>
      <Rank>no rank</Rank>
      <LineageEx>
        <Taxon><ScientificName>Viruses</ScientificName><Rank>acellular root</Rank></Taxon>
        <Taxon><ScientificName>Bunyaviricetes</ScientificName><Rank>class</Rank></Taxon>
        <Taxon><ScientificName>Elliovirales</ScientificName><Rank>order</Rank></Taxon>
        <Taxon><ScientificName>Peribunyaviridae</ScientificName><Rank>family</Rank></Taxon>
        <Taxon><ScientificName>Orthobunyavirus</ScientificName><Rank>genus</Rank></Taxon>
        <Taxon><ScientificName>Orthobunyavirus heptayabaense</ScientificName><Rank>species</Rank></Taxon>
      </LineageEx>
    </Taxon></TaxaSet>"""

    ranks = _parse_taxonomy_xml(xml_text)
    assert ranks["genus"] == "Orthobunyavirus"
    assert ranks["family"] == "Peribunyaviridae"
    assert ranks["order"] == "Elliovirales"
    assert ranks["class"] == "Bunyaviricetes"
    assert ranks["species"] == "Orthobunyavirus heptayabaense"
    # 'no rank' entries (here the queried taxon itself) are skipped.
    assert "no rank" not in ranks


def test_parse_taxonomy_xml_handles_garbage():
    from repseq.taxonomy.ncbi import _parse_taxonomy_xml

    assert _parse_taxonomy_xml("not xml at all") == {}
    assert _parse_taxonomy_xml("<TaxaSet></TaxaSet>") == {}


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
