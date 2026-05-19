"""Tests for repseq.hmm.runner (cache-aware batch scan + cutoff logic)."""
from __future__ import annotations

from unittest.mock import patch

import pytest

from repseq.hmm.errors import HMMScanError
from repseq.hmm.runner import (
    CACHE_SOURCE,
    _cache_key,
    coverage_of,
    passes_cutoffs,
    scan_proteins,
)


def _hit(target="RdRp", dom_score=200.0, dom_evalue=1e-50, hmm_len=300, ali_span=280):
    # Coverage is measured on the HMM model span (hmm_to - hmm_from + 1).
    # These tests express the covered span via ali_span for readability, so
    # mirror it onto hmm coords (hmm_from=1, hmm_to=ali_span). ali_span stays
    # as a real field (used for domain ordering elsewhere).
    return {
        "target": target,
        "dom_score": dom_score,
        "dom_evalue": dom_evalue,
        "hmm_len": hmm_len,
        "ali_span": ali_span,
        "hmm_from": 1,
        "hmm_to": ali_span,
    }


def test_passes_cutoffs_ga_path():
    """When the HMM has a GA cutoff and use_ga is True, the dom_score is
    compared against the GA — E-value is ignored for the similarity gate."""
    h = _hit(dom_score=30.0, dom_evalue=10.0)  # huge E but passes GA
    assert passes_cutoffs(
        h, ga_cutoffs={"RdRp": 25.0},
        default_evalue=1e-5, relative_length_cutoff=0.5, use_ga_when_available=True,
    ) is True


def test_passes_cutoffs_ga_fails_when_score_below_threshold():
    h = _hit(dom_score=20.0)
    assert passes_cutoffs(
        h, ga_cutoffs={"RdRp": 25.0},
        default_evalue=1e-5, relative_length_cutoff=0.5, use_ga_when_available=True,
    ) is False


def test_passes_cutoffs_evalue_fallback_when_no_ga():
    """Profile has no GA → fall back to default_evalue (the E-value gate)."""
    h = _hit(dom_evalue=1e-30)
    assert passes_cutoffs(
        h, ga_cutoffs={"RdRp": None},
        default_evalue=1e-5, relative_length_cutoff=0.5, use_ga_when_available=True,
    ) is True
    h_bad = _hit(dom_evalue=1.0)
    assert passes_cutoffs(
        h_bad, ga_cutoffs={"RdRp": None},
        default_evalue=1e-5, relative_length_cutoff=0.5, use_ga_when_available=True,
    ) is False


def test_passes_cutoffs_ignores_ga_when_disabled():
    """use_ga_when_available=False forces the E-value gate even when GA is set."""
    h = _hit(dom_score=100.0, dom_evalue=1.0)
    assert passes_cutoffs(
        h, ga_cutoffs={"RdRp": 25.0},
        default_evalue=1e-5, relative_length_cutoff=0.5, use_ga_when_available=False,
    ) is False


def test_passes_cutoffs_length_gate_independent_of_similarity():
    """Even with a great E-value, a short alignment over a long HMM fails."""
    h = _hit(hmm_len=1000, ali_span=200)  # 20% coverage
    assert passes_cutoffs(
        h, ga_cutoffs={},
        default_evalue=1e-5, relative_length_cutoff=0.5, use_ga_when_available=True,
    ) is False


def test_passes_cutoffs_both_gates_required():
    """Both similarity AND length must pass."""
    h_short = _hit(hmm_len=500, ali_span=100)  # similarity OK, length fail
    h_weak = _hit(dom_evalue=1.0)               # length OK, similarity fail
    cfg = dict(default_evalue=1e-5, relative_length_cutoff=0.5, use_ga_when_available=True)
    assert passes_cutoffs(h_short, ga_cutoffs={}, **cfg) is False
    assert passes_cutoffs(h_weak, ga_cutoffs={}, **cfg) is False


def test_coverage_of_uses_hmm_model_span():
    """Coverage is (hmm_to - hmm_from + 1) / hmm_len, NOT the protein
    alignment span. A hit with a big insertion (large ali_span) but a
    small HMM-model span must report the small model coverage."""
    # Model span 150 of a 300-long model = 0.50, regardless of a large
    # protein alignment span (insertions don't inflate coverage).
    assert coverage_of(
        {"hmm_from": 1, "hmm_to": 150, "hmm_len": 300, "ali_span": 290}
    ) == 0.5
    # Pathological: model span > hmm_len would exceed 1.0 — clip.
    assert coverage_of({"hmm_from": 1, "hmm_to": 200, "hmm_len": 100}) == 1.0
    # hmm_len=0 (defensive): max(hmm_len, 1) → 1, so span/1 clipped to 1.0.
    assert coverage_of({"hmm_from": 1, "hmm_to": 50, "hmm_len": 0}) == 1.0


def test_cache_key_changes_with_db_signature():
    """Different DB signature must yield different cache keys (so a new
    DB invalidates old hits)."""
    seq = "MMMMM"
    k1 = _cache_key(seq, "sig_a")
    k2 = _cache_key(seq, "sig_b")
    assert k1 != k2


def test_scan_proteins_empty_input_returns_empty():
    assert scan_proteins({}, cfg={"hmm": {}}) == {}


def test_scan_proteins_raises_when_hmmscan_missing():
    with patch("repseq.hmm.runner.is_available", return_value=False):
        with pytest.raises(HMMScanError, match="not on PATH"):
            scan_proteins({"q": "M"}, cfg={"hmm": {}})


def test_scan_proteins_uses_cache_to_skip_already_seen(tmp_path):
    """Sequences with a cached entry must not be re-scanned. Mocks
    scan() to confirm only the uncached sequence is passed through."""
    from repseq.taxonomy.cache import TaxonomyCache
    cache = TaxonomyCache(tmp_path, ttl_days=0)
    # Pre-seed a hit for one sequence.
    seq_cached = "M" * 50
    seq_new = "K" * 50
    db_path = tmp_path / "fake.hmm"
    db_path.write_text("HMMER3/f\nNAME Foo\n//\n")
    # The bundled lookup in scan_proteins will resolve_database_path -> the
    # bundled .hmm (the user-supplied null path). Patch to a known temp path
    # so db_signature is deterministic.
    with patch("repseq.hmm.runner.resolve_database_path", return_value=db_path), \
         patch("repseq.hmm.runner.ensure_pressed"), \
         patch("repseq.hmm.runner.db_signature", return_value="sig123"), \
         patch("repseq.hmm.runner.is_available", return_value=True):
        # Cache the hit for the cached sequence under the right key.
        cache.set(
            CACHE_SOURCE, _cache_key(seq_cached, "sig123"),
            {"hits": [{"target": "CachedHit", "dom_evalue": 1e-30}]},
        )

        def fake_scan(db, queries, threads):
            # Only seq_new should be batched (the cached one is filtered out).
            assert len(queries) == 1
            # Return one synthesised hit for the uncached sequence.
            return {next(iter(queries)): [{"target": "FreshHit", "dom_evalue": 1e-25}]}

        with patch("repseq.hmm.runner.scan", side_effect=fake_scan):
            results = scan_proteins(
                {"q1": seq_cached, "q2": seq_new},
                cfg={"hmm": {}}, cache=cache,
            )
    assert results["q1"][0]["target"] == "CachedHit"
    assert results["q2"][0]["target"] == "FreshHit"


def test_scan_proteins_dedupes_identical_sequences(tmp_path):
    """Two qids pointing at the same AA sequence must result in a single
    scan call carrying one query (key collapse)."""
    db_path = tmp_path / "fake.hmm"
    db_path.write_text("HMMER3/f\nNAME Foo\n//\n")
    seq = "M" * 100
    captured_queries: list[dict] = []

    def fake_scan(db, queries, threads):
        captured_queries.append(queries)
        return {next(iter(queries)): [{"target": "DedupHit", "dom_evalue": 1e-30}]}

    with patch("repseq.hmm.runner.resolve_database_path", return_value=db_path), \
         patch("repseq.hmm.runner.ensure_pressed"), \
         patch("repseq.hmm.runner.db_signature", return_value="sig"), \
         patch("repseq.hmm.runner.is_available", return_value=True), \
         patch("repseq.hmm.runner.scan", side_effect=fake_scan):
        results = scan_proteins(
            {"a": seq, "b": seq, "c": seq},
            cfg={"hmm": {}},
        )
    assert len(captured_queries) == 1
    assert len(captured_queries[0]) == 1  # one unique seq → one scan query
    # All three qids end up with the same hit list.
    assert results["a"] == results["b"] == results["c"]
    assert results["a"][0]["target"] == "DedupHit"
