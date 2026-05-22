"""Tests for the HMM-gated branch of select_marker_protein.

The legacy alias-only behaviour is covered in test_marker.py; this file
focuses on the v0.13 strict HMM gate added for dict-form marker specs
with ``hmms: [...]``.
"""
from __future__ import annotations

from repseq.clustering.marker import (
    MarkerFailure,
    _format_failure_reason,
    populate_protein_sequences,
    select_marker_protein,
)
from repseq.models import QCReport


def _cds(product, sequence, *, protein_id="P_x", hmm_hits=None):
    p = {
        "protein_id": protein_id,
        "product": product,
        "length": len(sequence),
        "sequence": sequence,
    }
    if hmm_hits is not None:
        p["hmm_hits"] = hmm_hits
    return p


def _hit(target, passing=True, dom_evalue=1e-30, hmm_len=300, ali_span=280,
         ali_from=1, ali_to=None):
    if ali_to is None:
        ali_to = ali_from + ali_span - 1
    return {
        "target": target,
        "passing": passing,
        "dom_evalue": dom_evalue,
        "dom_score": 200.0,
        "hmm_len": hmm_len,
        "ali_span": ali_span,
        "ali_from": ali_from,
        "ali_to": ali_to,
    }


# ---------------------------------------------------------------------------
# Strict HMM gate
# ---------------------------------------------------------------------------

def test_hmm_gate_picks_only_cds_with_passing_hit():
    """The CDS labelled 'polymerase' but missing the HMM hit must NOT
    be chosen, even though its /product matches the alias the spec
    also lists. The whole point of the strict gate."""
    proteins = [
        _cds("RNA polymerase", "M" * 2200, protein_id="P_pol_misnamed"),
        _cds("hypothetical protein", "M" * 1000, protein_id="P_real_pol",
             hmm_hits=[_hit("RdRP_4", passing=True)]),
    ]
    spec = [{"name": "L", "aliases": ["polymerase"], "hmms": ["RdRP_4"]}]
    marker, failure = select_marker_protein(proteins, spec, hmm_active=True)
    assert failure is None
    assert marker["protein_id"] == "P_real_pol"


def test_hmm_gate_rejects_when_no_cds_passes():
    """The polymerase CDS hits the HMM but BELOW cutoffs (passing=False).
    The strict gate drops the sequence rather than falling back to the
    alias match — that's the v0.13 design."""
    proteins = [
        _cds("RNA polymerase", "M" * 2200,
             hmm_hits=[_hit("RdRP_4", passing=False, dom_evalue=0.01)]),
    ]
    spec = [{"name": "L", "aliases": ["polymerase"], "hmms": ["RdRP_4"]}]
    marker, failure = select_marker_protein(proteins, spec, hmm_active=True)
    assert marker is None
    assert failure.reason == "hmm_failed"
    assert failure.marker_name == "L"
    assert "RdRP_4" in failure.failed_tokens


def test_multidomain_token_requires_all_hmms_on_same_cds():
    """v0.14.0 multidomain token semantic: 'A--B' requires a SINGLE CDS
    with passing hits to both HMMs in N-to-C order (A N-terminal to B).
    The partial CDS satisfying only one domain must NOT be picked even
    though it is longer; only the full polyprotein qualifies."""
    proteins = [
        # Only one of two required HMMs hits — does NOT satisfy "A--B".
        _cds("glycoprotein", "M" * 1000, protein_id="gly_partial",
             hmm_hits=[_hit("Bunya_G1", passing=True, ali_from=500, ali_to=800)]),
        # Both HMMs hit in correct N-to-C order — does satisfy "G2--G1".
        _cds("glycoprotein", "M" * 800, protein_id="gly_full",
             hmm_hits=[
                 # Bunya_G2 N-terminal (positions 1-300).
                 _hit("Bunya_G2", passing=True, ali_from=1, ali_to=300),
                 # Bunya_G1 C-terminal (positions 500-800) → strictly after G2.
                 _hit("Bunya_G1", passing=True, ali_from=500, ali_to=800),
             ]),
    ]
    spec = [{"name": "M", "aliases": ["glycoprotein"], "hmms": ["Bunya_G2--Bunya_G1"]}]
    marker, failure = select_marker_protein(proteins, spec, hmm_active=True)
    assert failure is None
    assert marker["protein_id"] == "gly_full"


def test_separate_hmm_tokens_are_or_semantic_in_marker_selection():
    """v0.14.0 hard cutover: ``hmms: ['A', 'B']`` is TWO separate tokens,
    not the v0.13.0 list-AND. For marker selection, the longest CDS
    satisfying ANY token is chosen — so a CDS hitting only A is still
    a satisfying CDS for the marker (different intent from the
    multidomain token form 'A--B')."""
    proteins = [
        _cds("glycoprotein", "M" * 1500, protein_id="gly_A_only",
             hmm_hits=[_hit("Bunya_G1", passing=True)]),
        _cds("glycoprotein", "M" * 800, protein_id="gly_AB",
             hmm_hits=[
                 _hit("Bunya_G1", passing=True),
                 _hit("Bunya_G2", passing=True),
             ]),
    ]
    spec = [{"name": "M", "hmms": ["Bunya_G1", "Bunya_G2"]}]
    marker, failure = select_marker_protein(proteins, spec, hmm_active=True)
    assert failure is None
    # gly_A_only satisfies the 'Bunya_G1' token. gly_AB satisfies both.
    # Longest satisfying CDS wins → gly_A_only (1500 > 800).
    assert marker["protein_id"] == "gly_A_only"


def test_hmm_gate_inactive_falls_back_to_alias():
    """When hmm_active=False (tool missing / disabled), HMM specs
    behave like alias-only specs, preserving alias matching."""
    proteins = [
        _cds("RNA polymerase", "M" * 2200, protein_id="P_pol"),
        _cds("nucleoprotein", "M" * 400),
    ]
    spec = [{"name": "L", "aliases": ["polymerase"], "hmms": ["RdRP_4"]}]
    marker, failure = select_marker_protein(proteins, spec, hmm_active=False)
    assert failure is None
    assert marker["protein_id"] == "P_pol"


def test_hmm_gate_picks_longest_among_multiple_passing():
    """Tiebreaker among HMM-passing CDSes is length (consistent with
    the alias-tier tiebreaker)."""
    proteins = [
        _cds("polymerase", "M" * 100,
             hmm_hits=[_hit("RdRP_4", passing=True)], protein_id="short_pol"),
        _cds("polymerase", "M" * 2000,
             hmm_hits=[_hit("RdRP_4", passing=True)], protein_id="long_pol"),
    ]
    spec = [{"name": "L", "aliases": ["polymerase"], "hmms": ["RdRP_4"]}]
    marker, _ = select_marker_protein(proteins, spec, hmm_active=True)
    assert marker["protein_id"] == "long_pol"


def test_legacy_alias_fall_through_unchanged_for_alias_only_specs():
    """v0.12-compatible: alias-only specs that don't match fall through
    to longest CDS (this is the behaviour test_marker.py also pins;
    duplicated here to guard against accidental coupling with the new
    HMM gate logic)."""
    proteins = [
        _cds("polymerase short", "M" * 100),
        _cds("polymerase full", "M" * 2000),
    ]
    marker, _ = select_marker_protein(proteins, ["glycoprotein"])
    # No alias matched → fall through to longest.
    assert marker["product"] == "polymerase full"


def test_mixed_specs_hmm_fails_alias_only_succeeds():
    """When the first (HMM) spec fails but a later alias-only spec
    matches, the alias spec wins. Confirms specs are tried in order."""
    proteins = [
        _cds("nucleocapsid", "M" * 400, protein_id="N_pid",
             hmm_hits=[_hit("RdRP_4", passing=False)]),
    ]
    specs = [
        {"name": "L", "aliases": ["polymerase"], "hmms": ["RdRP_4"]},
        {"name": "N", "aliases": ["nucleocapsid"], "hmms": []},
    ]
    marker, failure = select_marker_protein(proteins, specs, hmm_active=True)
    assert failure is None
    assert marker["protein_id"] == "N_pid"


def test_format_failure_reason_renders_hmm_details():
    f = MarkerFailure(
        reason="hmm_failed", marker_name="L",
        failed_tokens=["RdRP_4"], best_evalue=0.01,
    )
    out = _format_failure_reason(f)
    assert out.startswith("hmm_failed:L:RdRP_4")
    assert "E=" in out


def test_format_failure_reason_renders_multidomain_token():
    """Multidomain tokens are rendered intact (the '--' separator is
    preserved) so the _qc_removed.tsv reason shows the exact token
    string the user configured."""
    f = MarkerFailure(
        reason="hmm_failed", marker_name="M",
        failed_tokens=["Bunya_G1--Bunya_G2"], best_evalue=None,
    )
    out = _format_failure_reason(f)
    assert out == "hmm_failed:M:Bunya_G1--Bunya_G2"


def test_format_failure_reason_non_hmm_uses_legacy_string():
    """Non-HMM failures keep the legacy string for back-compat with
    downstream parsers (e.g. anything that greps _qc_removed.tsv for
    'no_marker_protein_for_clustering')."""
    f = MarkerFailure(reason="no_proteins")
    assert _format_failure_reason(f) == "no_marker_protein_for_clustering"
    f = MarkerFailure(reason="no_alias_match")
    assert _format_failure_reason(f) == "no_marker_protein_for_clustering"


# ---------------------------------------------------------------------------
# populate_protein_sequences with HMM tier
# ---------------------------------------------------------------------------

def test_populate_increments_hmm_counters_on_drop(make_seq):
    """A sequence whose HMM-gated spec fails must be counted under
    removed_hmm_failed + removed_hmm_by_marker — NOT under the generic
    removed_proteins (which would blur it with the missing-marker case)."""
    s = make_seq("s1", "ACGT")
    s.proteins = [_cds("polymerase", "M" * 2200,
                       hmm_hits=[_hit("RdRP_4", passing=False)])]
    spec = [{"name": "L", "aliases": ["polymerase"], "hmms": ["RdRP_4"]}]
    report = QCReport()
    kept = populate_protein_sequences([s], spec, report, hmm_active=True)
    assert kept == []
    assert report.removed_hmm_failed == 1
    assert report.removed_hmm_by_marker.get("L") == 1
    # The legacy counter must NOT be incremented for HMM drops.
    assert report.removed_proteins == 0
    # The reason in _qc_removed.tsv carries the HMM detail.
    assert any("hmm_failed:L" in d["reason"] for d in report.details)


def test_populate_legacy_counter_still_used_for_non_hmm_drops(make_seq):
    """No-translation / no-proteins drops continue to land under
    removed_proteins so existing behaviour is preserved."""
    a = make_seq("a", "ACGT")
    a.proteins = None
    report = QCReport()
    kept = populate_protein_sequences([a], None, report)
    assert kept == []
    assert report.removed_proteins == 1
    assert report.removed_hmm_failed == 0
