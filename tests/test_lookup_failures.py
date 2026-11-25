"""Database lookups that fail for a network reason: retried, then reported.

Found on the Orthoherpesviridae pass 1 (2026-11-12): a transient Entrez error
was swallowed by ``except Exception: return None`` — no retry, no count — so
9 of 5,425 records silently lost their taxonomy and formed an ``Unknown``
group. The GenBank batch fetch was worse: a failed batch CACHED its (up to
200) accessions as "record with no proteins".
"""
from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest
import requests

from repseq.models import SequenceSource
from repseq.output.report import write_lookup_failures_tsv
from repseq.taxonomy import http
from repseq.taxonomy.cache import TaxonomyCache
from repseq.taxonomy.ncbi import NCBITaxonomy
from repseq.taxonomy.resolver import MetadataResolver
from repseq.taxonomy.uniprot import UniProtAPI


@pytest.fixture(autouse=True)
def _fresh_log():
    http.reset_failures()
    with patch.object(http, "_sleep") as fake_sleep:
        yield fake_sleep
    http.reset_failures()


def _http_error(status, headers=None):
    resp = MagicMock()
    resp.status_code = status
    resp.headers = headers or {}
    return requests.HTTPError(f"{status}", response=resp)


def _ok(json_body=None, text=""):
    resp = MagicMock()
    resp.raise_for_status.return_value = None
    resp.json.return_value = json_body or {}
    resp.text = text
    return resp


def _failing(status):
    resp = MagicMock()
    resp.raise_for_status.side_effect = _http_error(status)
    return resp


# --- get_with_retry --------------------------------------------------------

def test_transient_error_is_retried_and_then_succeeds(_fresh_log):
    good = _ok({"x": 1})
    with patch("requests.get", side_effect=[_failing(429), _failing(503), good]) as g:
        resp = http.get_with_retry("http://x", params={})
    assert resp is good
    assert g.call_count == 3
    assert _fresh_log.call_count == 2          # waited between attempts


def test_gives_up_after_max_tries():
    with patch("requests.get", return_value=_failing(503)) as g:
        with pytest.raises(requests.HTTPError):
            http.get_with_retry("http://x", params={})
    assert g.call_count == http.MAX_TRIES


def test_not_found_is_not_retried():
    with patch("requests.get", return_value=_failing(404)) as g:
        with pytest.raises(requests.HTTPError):
            http.get_with_retry("http://x", params={})
    assert g.call_count == 1


def test_timeout_and_connection_errors_are_retried():
    good = _ok()
    with patch(
        "requests.get",
        side_effect=[requests.Timeout("slow"), requests.ConnectionError("x"), good],
    ) as g:
        assert http.get_with_retry("http://x", params={}) is good
    assert g.call_count == 3


def test_retry_after_header_sets_the_wait(_fresh_log):
    bad = MagicMock()
    bad.raise_for_status.side_effect = _http_error(429, {"Retry-After": "7"})
    with patch("requests.get", side_effect=[bad, _ok()]):
        http.get_with_retry("http://x", params={})
    _fresh_log.assert_called_once_with(7.0)


def test_throttle_runs_before_every_attempt():
    throttle = MagicMock()
    with patch("requests.get", side_effect=[_failing(500), _ok()]):
        http.get_with_retry("http://x", params={}, throttle=throttle)
    assert throttle.call_count == 2


# --- NCBI client -----------------------------------------------------------

def test_metadata_failure_is_recorded_and_not_cached(tmp_cache_dir):
    cache = TaxonomyCache(tmp_cache_dir)
    ncbi = NCBITaxonomy(cache, email="t@e.com")
    with patch("requests.get", return_value=_failing(429)):
        assert ncbi.fetch_accession_metadata("AB000001.1") is None
    assert http.failure_for(http.STAGE_METADATA, "AB000001.1") == (
        "HTTP 429 (rate limit reached)"
    )
    # Nothing cached, so a re-run retries it.
    assert cache.get("ncbi_nuccore", "AB000001.1") is None


def test_genuine_not_found_is_cached_and_not_a_failure(tmp_cache_dir):
    cache = TaxonomyCache(tmp_cache_dir)
    ncbi = NCBITaxonomy(cache, email="t@e.com")
    empty = _ok({"esearchresult": {"idlist": []}})
    with patch("requests.get", return_value=empty):
        assert ncbi.fetch_accession_metadata("ZZ999999.1") is None
    assert http.failures() == []
    assert cache.get("ncbi_nuccore", "ZZ999999.1") == {"_not_found": True}


def test_lineage_failure_is_recorded(tmp_cache_dir):
    ncbi = NCBITaxonomy(TaxonomyCache(tmp_cache_dir), email="t@e.com")
    with patch("requests.get", side_effect=requests.Timeout("slow")):
        assert ncbi.fetch_lineage(10298) is None
    # Keyed by taxid under an internal stage; not itself a report row.
    assert http.failure_for(http.STAGE_LINEAGE_TAXID, 10298) == (
        "no answer within the time limit"
    )
    assert http.failures() == []


def test_failed_genbank_batch_is_not_cached_as_empty_records(tmp_cache_dir):
    """The cache-poisoning half: a failed batch must not make its
    accessions look like records that annotate no proteins."""
    cache = TaxonomyCache(tmp_cache_dir)
    ncbi = NCBITaxonomy(cache, email="t@e.com")
    with patch("requests.get", return_value=_failing(502)):
        out = ncbi.fetch_proteins_batch(["A.1", "B.1"])
    assert out == {"A.1": [], "B.1": []}          # this run: nothing usable
    assert cache.get("ncbi_proteins", "A.1") is None
    assert cache.get("ncbi_proteins", "B.1") is None
    stages = {(f["stage"], f["id"]) for f in http.failures()}
    assert stages == {("genbank_record", "A.1"), ("genbank_record", "B.1")}


def test_failed_nucleotide_batch_is_not_cached_as_not_found(tmp_cache_dir):
    cache = TaxonomyCache(tmp_cache_dir)
    ncbi = NCBITaxonomy(cache, email="t@e.com")
    with patch("requests.get", return_value=_failing(503)):
        out = ncbi.fetch_nucleotide_batch(["A.1"])
    assert out == {"A.1": None}
    assert cache.get("ncbi_nuc_seq", "A.1") is None
    assert http.failure_for(http.STAGE_SEQUENCE, "A.1")


# --- resolver --------------------------------------------------------------

def test_resolver_counts_a_swallowed_lookup_failure(tmp_cache_dir, make_seq):
    cache = TaxonomyCache(tmp_cache_dir)
    ncbi = NCBITaxonomy(cache, email="t@e.com")
    resolver = MetadataResolver(cache, ncbi, UniProtAPI(cache), threads=1)
    seq = make_seq("AB000001.1", "ACGT" * 10)
    seq.accession = "AB000001.1"
    seq.source = SequenceSource.NCBI
    with patch("requests.get", return_value=_failing(429)):
        resolver.resolve_batch([seq])
    assert len(resolver.failures) == 1
    sid, msg = resolver.failures[0]
    assert sid == seq.id and "rate limit" in msg


def test_resolver_does_not_count_a_genuine_miss(tmp_cache_dir, make_seq):
    cache = TaxonomyCache(tmp_cache_dir)
    ncbi = NCBITaxonomy(cache, email="t@e.com")
    resolver = MetadataResolver(cache, ncbi, UniProtAPI(cache), threads=1)
    seq = make_seq("ZZ999999.1", "ACGT" * 10)
    seq.accession = "ZZ999999.1"
    seq.source = SequenceSource.NCBI
    with patch("requests.get", return_value=_ok({"esearchresult": {"idlist": []}})):
        resolver.resolve_batch([seq])
    assert resolver.failures == []


# --- audit file ------------------------------------------------------------

def test_lookup_failures_tsv_written_only_when_there_are_failures(tmp_path):
    path = tmp_path / "x_lookup_failures.tsv"
    assert write_lookup_failures_tsv({}, path) is False
    assert not path.exists()

    http.record_failure(http.STAGE_METADATA, "AB000001.1", "HTTP 429 (rate limit reached)")
    http.record_failure(http.STAGE_RECORD, "CD000002.1", "could not connect")
    cfg = {"_lookup_failures_runtime": {"failures": http.failures()}}
    assert write_lookup_failures_tsv(cfg, path) is True
    rows = path.read_text().splitlines()
    assert rows[0] == "stage\tid\terror"
    assert rows[1:] == [
        "genbank_record\tCD000002.1\tcould not connect",
        "metadata\tAB000001.1\tHTTP 429 (rate limit reached)",
    ]


# --- review fixes (2026-11-25) ---------------------------------------------

def test_retries_stop_once_the_network_is_clearly_gone(_fresh_log):
    """Offline: after a few unreachable requests in a row, fail fast instead
    of sleeping 6 s on every one of 120,000 lookups."""
    with patch("requests.get", side_effect=requests.ConnectionError("down")) as g:
        for _ in range(10):
            with pytest.raises(requests.ConnectionError):
                http.get_with_retry("http://x", params={})
    # First request: 3 attempts. By the 5th consecutive failure retries are
    # suspended, so the remaining requests cost one attempt each.
    assert g.call_count < 10 * http.MAX_TRIES
    assert g.call_count <= http._OFFLINE_AFTER + 2 + 9
    sleeps_while_offline = _fresh_log.call_count
    with patch("requests.get", side_effect=requests.ConnectionError("down")):
        with pytest.raises(requests.ConnectionError):
            http.get_with_retry("http://x", params={})
    assert _fresh_log.call_count == sleeps_while_offline      # no more waiting


def test_retries_resume_after_a_request_gets_through():
    with patch("requests.get", side_effect=requests.ConnectionError("down")):
        for _ in range(6):
            with pytest.raises(requests.ConnectionError):
                http.get_with_retry("http://x", params={})
    with patch("requests.get", return_value=_ok()):
        http.get_with_retry("http://x", params={})
    with patch("requests.get", side_effect=[requests.Timeout("slow"), _ok()]) as g:
        http.get_with_retry("http://x", params={})
    assert g.call_count == 2


def test_server_errors_do_not_trip_the_offline_switch():
    """A 503 is an ANSWER: the network is there, keep retrying."""
    with patch("requests.get", return_value=_failing(503)) as g:
        for _ in range(8):
            with pytest.raises(requests.HTTPError):
                http.get_with_retry("http://x", params={})
    assert g.call_count == 8 * http.MAX_TRIES


def test_uniprot_bad_request_is_not_a_network_failure(tmp_cache_dir):
    """UniProt answering 400 to an id it does not recognise is a fact about
    the record, not a failed lookup."""
    uni = UniProtAPI(TaxonomyCache(tmp_cache_dir))
    with patch("requests.get", return_value=_failing(400)):
        assert uni.fetch_entry("XYZ_001") is None
    assert http.failures() == []
    with patch("requests.get", return_value=_failing(503)):
        assert uni.fetch_entry("P12345") is None
    assert [f["id"] for f in http.failures()] == ["P12345"]


def test_later_successful_genbank_fetch_clears_the_failure(tmp_cache_dir):
    """The GenBank fetch runs at several pipeline steps; a record that got
    through on a later one must not stay on the failure report."""
    ncbi = NCBITaxonomy(TaxonomyCache(tmp_cache_dir), email="t@e.com")
    with patch("requests.get", return_value=_failing(502)):
        ncbi.fetch_proteins_batch(["A.1"])
    assert http.failure_for(http.STAGE_RECORD, "A.1")
    with patch("requests.get", return_value=_ok(text="")):
        ncbi.fetch_proteins_batch(["A.1"])
    assert http.failure_for(http.STAGE_RECORD, "A.1") is None


def test_lineage_failure_stays_pinned_to_the_record(tmp_cache_dir, make_seq):
    """Another record's later success for the same taxid clears the shared
    taxid entry, but the record already resolved without its lineage must
    stay on the report."""
    cache = TaxonomyCache(tmp_cache_dir)
    ncbi = NCBITaxonomy(cache, email="t@e.com")
    resolver = MetadataResolver(cache, ncbi, UniProtAPI(cache), threads=1)
    seq = make_seq("AB000001.1", "ACGT" * 10)
    seq.accession = "AB000001.1"
    seq.source = SequenceSource.NCBI
    esearch = _ok({"esearchresult": {"idlist": ["1"]}})
    esummary = _ok({"result": {"1": {"taxid": 10298, "organism": "HSV-1"}}})
    lineage_down = [_failing(503)] * (2 * http.MAX_TRIES)
    with patch("requests.get", side_effect=[esearch, esummary] + lineage_down):
        resolver.resolve_batch([seq])
    assert len(resolver.failures) == 1
    rows = http.failures()
    assert [(r["stage"], r["id"]) for r in rows] == [("taxonomy", "AB000001.1")]
    # Another record's fetch of the same taxid now succeeds ...
    http.clear_failure(http.STAGE_LINEAGE_TAXID, 10298)
    # ... and this record is still reported.
    assert [(r["stage"], r["id"]) for r in http.failures()] == [
        ("taxonomy", "AB000001.1")
    ]


def test_unknown_source_record_recovers_its_lineage_on_a_rerun(tmp_cache_dir, make_seq):
    """Metadata is cached even when the lineage fetch failed; the next run
    must ask for the lineage again on the unknown-source path too."""
    cache = TaxonomyCache(tmp_cache_dir)
    cache.set("ncbi_nuccore", "AB000001.1",
              {"accession": "AB000001.1", "taxid": 10298, "organism": "HSV-1"})
    ncbi = NCBITaxonomy(cache, email="t@e.com")
    resolver = MetadataResolver(cache, ncbi, UniProtAPI(cache), threads=1)
    seq = make_seq("AB000001.1", "ACGT" * 10)
    seq.accession = "AB000001.1"
    seq.source = SequenceSource.UNKNOWN
    with patch.object(ncbi, "fetch_lineage", return_value={
        "taxid": 10298, "species": "Simplexvirus humanalpha1", "lineage": {},
    }) as fl:
        resolver.resolve_batch([seq])
    fl.assert_called_once_with(10298)
    assert resolver.failures == []
