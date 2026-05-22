"""Tests for the v0.14.0 HMM-QC pipeline step (``cli._run_hmm_qc`` and
its segmented / non-segmented branches).

The HMM scan itself is mocked: we pre-populate ``seq.proteins[*]
['hmm_hits']`` with pass/fail data and call the drop helpers directly,
so these tests run without hmmscan installed and without network.
"""
from __future__ import annotations

from repseq.cli import (
    _run_hmm_qc_non_segmented,
    _run_hmm_qc_segmented,
    _resolve_segment_hmms,
    _segment_fails_hmm_gate,
)
from repseq.models import QCReport


def _cds(product, sequence, *, protein_id="P_x", hmm_hits=None):
    return {
        "protein_id": protein_id,
        "product": product,
        "length": len(sequence),
        "sequence": sequence,
        "hmm_hits": hmm_hits or [],
    }


def _hit(target, ali_from, ali_to, *, passing=True, dom_evalue=1e-30):
    return {
        "target": target,
        "passing": passing,
        "dom_evalue": dom_evalue,
        "ali_from": ali_from,
        "ali_to": ali_to,
        "ali_span": ali_to - ali_from + 1,
    }


def _attach_proteins(seq, proteins):
    seq.proteins = proteins


# ---------------------------------------------------------------------------
# _resolve_segment_hmms — lookup priority
# ---------------------------------------------------------------------------

def test_resolve_segment_hmms_prefers_segment_markers_over_cluster_protein():
    sm = {"S": {"hmms": ["Bunya_nucleocap"]}}
    cp = {"S": [{"name": "S_legacy", "hmms": ["LegacyHMM"]}]}
    name, tokens = _resolve_segment_hmms("S", sm, cp)
    assert name == "S"
    assert tokens == ["Bunya_nucleocap"]


def test_resolve_segment_hmms_falls_back_to_cluster_protein():
    cp = {"L": [{"name": "RdRp", "hmms": ["RdRP_4"]}]}
    name, tokens = _resolve_segment_hmms("L", {}, cp)
    assert name == "RdRp"
    assert tokens == ["RdRP_4"]


def test_resolve_segment_hmms_returns_empty_when_no_spec():
    name, tokens = _resolve_segment_hmms("M", {}, {})
    assert tokens == []


def test_resolve_segment_hmms_returns_empty_when_spec_has_no_hmms():
    """Alias-only specs don't trigger the HMM QC gate."""
    sm = {"S": {"aliases": ["nucleocapsid"]}}
    name, tokens = _resolve_segment_hmms("S", sm, {})
    assert tokens == []


# ---------------------------------------------------------------------------
# _segment_fails_hmm_gate — per-segment token check
# ---------------------------------------------------------------------------

def test_segment_passes_when_every_token_satisfied(make_seq):
    seq = make_seq("seg1", "ACGT", segment="S")
    _attach_proteins(seq, [
        _cds("nucleocapsid", "M" * 300, hmm_hits=[_hit("Bunya_nucleocap", 1, 200)]),
    ])
    assert _segment_fails_hmm_gate(seq, ["Bunya_nucleocap"]) is None


def test_segment_fails_when_token_unsatisfied(make_seq):
    seq = make_seq("seg1", "ACGT", segment="S")
    _attach_proteins(seq, [
        _cds("nucleocapsid", "M" * 300, hmm_hits=[_hit("WrongHMM", 1, 200)]),
    ])
    assert _segment_fails_hmm_gate(seq, ["Bunya_nucleocap"]) == "Bunya_nucleocap"


def test_segment_passes_multidomain_token(make_seq):
    seq = make_seq("seg1", "ACGT", segment="M")
    _attach_proteins(seq, [
        _cds("glycoprotein polyprotein", "M" * 1100, hmm_hits=[
            _hit("Bunya_G2", 1, 300),       # N-terminal
            _hit("Bunya_G1", 500, 1000),    # C-terminal
        ]),
    ])
    # Token 'Bunya_G2--Bunya_G1' → G2 N-terminal, G1 C-terminal.
    assert _segment_fails_hmm_gate(seq, ["Bunya_G2--Bunya_G1"]) is None


def test_segment_fails_multidomain_when_only_one_domain_present(make_seq):
    seq = make_seq("seg1", "ACGT", segment="M")
    _attach_proteins(seq, [
        _cds("partial gly", "M" * 500, hmm_hits=[_hit("Bunya_G2", 1, 300)]),
    ])
    failed = _segment_fails_hmm_gate(seq, ["Bunya_G2--Bunya_G1"])
    assert failed == "Bunya_G2--Bunya_G1"


def test_segment_passes_when_any_alternative_token_satisfied(make_seq):
    """Tokens in one list are ALTERNATIVES (OR): satisfying either passes.

    Here both alternatives happen to be present (on different CDSes), which
    of course passes — the key OR cases (only one present, neither present)
    are the two tests below.
    """
    seq = make_seq("seg1", "ACGT", segment="M")
    _attach_proteins(seq, [
        _cds("Gn", "M" * 400, protein_id="gn", hmm_hits=[_hit("Bunya_G1", 1, 300)]),
        _cds("Gc", "M" * 400, protein_id="gc", hmm_hits=[_hit("Bunya_G2", 1, 300)]),
    ])
    assert _segment_fails_hmm_gate(seq, ["Bunya_G1", "Bunya_G2"]) is None


def test_segment_passes_when_only_one_alternative_present(make_seq):
    """OR: a segment carrying just one of the alternative architectures
    passes — the whole point of listing alternatives (e.g. a Spike that is
    CoV_S1--CoV_S2 OR bCoV_S1_N--bCoV_S1_RBD--CoV_S2)."""
    seq = make_seq("seg1", "ACGT", segment="M")
    _attach_proteins(seq, [
        _cds("Gn", "M" * 400, hmm_hits=[_hit("Bunya_G1", 1, 300)]),
    ])
    assert _segment_fails_hmm_gate(seq, ["Bunya_G1", "Bunya_G2"]) is None


def test_segment_fails_when_no_alternative_present(make_seq):
    """OR: only when NONE of the alternatives match does the segment fail;
    the reason names the alternatives joined with '|'."""
    seq = make_seq("seg1", "ACGT", segment="M")
    _attach_proteins(seq, [
        _cds("other", "M" * 400, hmm_hits=[_hit("WrongHMM", 1, 300)]),
    ])
    failed = _segment_fails_hmm_gate(seq, ["Bunya_G1", "Bunya_G2"])
    assert failed == "Bunya_G1|Bunya_G2"


# ---------------------------------------------------------------------------
# _run_hmm_qc_segmented — full integration
# ---------------------------------------------------------------------------

def _segmented_cfg():
    return {
        "segmented": {
            "enabled": True,
            "virus": "test_virus",
            "viruses": {
                "test_virus": {
                    "expected_segments": 3,
                    "segments": ["S", "M", "L"],
                    "segment_markers": {
                        "S": {"hmms": ["Bunya_nucleocap"]},
                        "M": {"hmms": ["Bunya_G2--Bunya_G1"]},
                        # L deliberately omitted — phased rollout case.
                    },
                },
            },
        },
        "_hmm_runtime": {"active": True},
    }


def _isolate_segs(make_seq, iso_id, *,
                  s_hits=None, m_hits=None, l_hits=None):
    s = make_seq(f"{iso_id}_S", "ACGT", segment="S", isolate_id=iso_id)
    m = make_seq(f"{iso_id}_M", "ACGT", segment="M", isolate_id=iso_id)
    l = make_seq(f"{iso_id}_L", "ACGT", segment="L", isolate_id=iso_id)
    _attach_proteins(s, [_cds("N", "M" * 300, hmm_hits=s_hits or [])])
    _attach_proteins(m, [_cds("GPC", "M" * 1100, hmm_hits=m_hits or [])])
    _attach_proteins(l, [_cds("RdRp", "M" * 2200, hmm_hits=l_hits or [])])
    return [s, m, l]


def test_segmented_qc_keeps_isolate_with_all_segments_passing(make_seq):
    cfg = _segmented_cfg()
    iso = _isolate_segs(
        make_seq, "iso1",
        s_hits=[_hit("Bunya_nucleocap", 1, 200)],
        m_hits=[
            _hit("Bunya_G2", 1, 300),
            _hit("Bunya_G1", 500, 1000),
        ],
        # L not gated; whatever hits or none.
        l_hits=[],
    )
    qc = QCReport()
    kept = _run_hmm_qc_segmented(iso, cfg, qc)
    assert len(kept) == 3
    assert qc.removed_hmm_failed == 0


def test_segmented_qc_drops_isolate_when_one_segment_fails(make_seq):
    cfg = _segmented_cfg()
    iso = _isolate_segs(
        make_seq, "iso2",
        s_hits=[_hit("Bunya_nucleocap", 1, 200)],
        # M is missing Bunya_G1 → token unsatisfied.
        m_hits=[_hit("Bunya_G2", 1, 300)],
        l_hits=[],
    )
    qc = QCReport()
    kept = _run_hmm_qc_segmented(iso, cfg, qc)
    assert kept == []
    assert qc.removed_hmm_failed == 1
    # Per-marker breakdown: M segment's multidomain token failed.
    assert "M:Bunya_G2--Bunya_G1" in qc.removed_hmm_by_marker
    assert qc.removed_hmm_by_marker["M:Bunya_G2--Bunya_G1"] == 1


def test_segmented_qc_records_per_segment_reasons_in_qc_removed(make_seq):
    cfg = _segmented_cfg()
    iso = _isolate_segs(
        make_seq, "iso3",
        # S fails its single-HMM token.
        s_hits=[_hit("WrongHMM", 1, 200)],
        m_hits=[
            _hit("Bunya_G2", 1, 300),
            _hit("Bunya_G1", 500, 1000),
        ],
        l_hits=[],
    )
    qc = QCReport()
    _run_hmm_qc_segmented(iso, cfg, qc)
    # The S segment carries the primary failure reason.
    s_entry = next(d for d in qc.details if d["id"] == "iso3_S")
    assert s_entry["reason"] == "hmm_failed:S:Bunya_nucleocap"
    # M and L carry sibling-drop reasons referencing S.
    m_entry = next(d for d in qc.details if d["id"] == "iso3_M")
    l_entry = next(d for d in qc.details if d["id"] == "iso3_L")
    assert "hmm_failed_sibling:S:Bunya_nucleocap" in m_entry["reason"]
    assert "hmm_failed_sibling:S:Bunya_nucleocap" in l_entry["reason"]


def test_segmented_qc_keeps_isolate_when_only_ungated_segment_lacks_hits(make_seq):
    """L has no spec → no gate. Even with zero hits, isolate must survive."""
    cfg = _segmented_cfg()
    iso = _isolate_segs(
        make_seq, "iso4",
        s_hits=[_hit("Bunya_nucleocap", 1, 200)],
        m_hits=[
            _hit("Bunya_G2", 1, 300),
            _hit("Bunya_G1", 500, 1000),
        ],
        l_hits=[],  # no hits, but L isn't gated
    )
    qc = QCReport()
    kept = _run_hmm_qc_segmented(iso, cfg, qc)
    assert len(kept) == 3
    assert qc.removed_hmm_failed == 0


def test_segmented_qc_mixed_isolates(make_seq):
    """Three isolates; only one passes."""
    cfg = _segmented_cfg()
    iso_good = _isolate_segs(
        make_seq, "good",
        s_hits=[_hit("Bunya_nucleocap", 1, 200)],
        m_hits=[
            _hit("Bunya_G2", 1, 300),
            _hit("Bunya_G1", 500, 1000),
        ],
    )
    iso_bad_s = _isolate_segs(
        make_seq, "bad_s",
        s_hits=[],  # no N protein hit
        m_hits=[
            _hit("Bunya_G2", 1, 300),
            _hit("Bunya_G1", 500, 1000),
        ],
    )
    iso_bad_m = _isolate_segs(
        make_seq, "bad_m",
        s_hits=[_hit("Bunya_nucleocap", 1, 200)],
        # M domains in wrong order for 'G2--G1' → multidomain token fails.
        m_hits=[
            _hit("Bunya_G1", 1, 300),     # G1 N-terminal (should be C-term)
            _hit("Bunya_G2", 500, 1000),  # G2 C-terminal (should be N-term)
        ],
    )
    all_seqs = iso_good + iso_bad_s + iso_bad_m
    qc = QCReport()
    kept = _run_hmm_qc_segmented(all_seqs, cfg, qc)
    kept_ids = {s.id for s in kept}
    # All 3 good isolate seqs survive; both bad isolates drop all 3 seqs.
    assert kept_ids == {"good_S", "good_M", "good_L"}
    assert qc.removed_hmm_failed == 2  # two isolates dropped


def test_segmented_qc_no_op_when_no_isolate_id(make_seq):
    """Sequences without isolate_id are not gated here — the regex
    fallback in _handle_segmented can still group them later, and
    build_concatenated_sequences applies its own backstop gate."""
    cfg = _segmented_cfg()
    seq = make_seq("orphan", "ACGT", segment="S")  # no isolate_id
    _attach_proteins(seq, [_cds("N", "M" * 300, hmm_hits=[])])
    qc = QCReport()
    kept = _run_hmm_qc_segmented([seq], cfg, qc)
    assert kept == [seq]
    assert qc.removed_hmm_failed == 0


# ---------------------------------------------------------------------------
# _run_hmm_qc_non_segmented
# ---------------------------------------------------------------------------

def _non_segmented_cfg(specs):
    return {
        "segmented": {"enabled": False},
        "clustering": {"cluster_protein": specs},
        "_hmm_runtime": {"active": True},
    }


def test_non_segmented_keeps_sequence_when_spec_satisfied(make_seq):
    cfg = _non_segmented_cfg([
        {"name": "RdRp", "hmms": ["RdRP_4"]},
    ])
    s = make_seq("s1", "ACGT")
    _attach_proteins(s, [
        _cds("polymerase", "M" * 2200, hmm_hits=[_hit("RdRP_4", 1, 280)]),
    ])
    qc = QCReport()
    kept = _run_hmm_qc_non_segmented([s], cfg, qc)
    assert kept == [s]
    assert qc.removed_hmm_failed == 0


def test_non_segmented_drops_sequence_when_spec_unsatisfied(make_seq):
    cfg = _non_segmented_cfg([
        {"name": "RdRp", "hmms": ["RdRP_4"]},
    ])
    s = make_seq("s1", "ACGT")
    _attach_proteins(s, [_cds("polymerase", "M" * 2200, hmm_hits=[])])
    qc = QCReport()
    kept = _run_hmm_qc_non_segmented([s], cfg, qc)
    assert kept == []
    assert qc.removed_hmm_failed == 1
    assert qc.removed_hmm_by_marker == {"RdRp": 1}


def test_non_segmented_requires_all_specs_to_pass(make_seq):
    """Multiple HMM-defining specs: a sequence must satisfy ALL of them
    (each is an independent marker the user wants verified)."""
    cfg = _non_segmented_cfg([
        {"name": "RdRp", "hmms": ["RdRP_4"]},
        {"name": "Nucleocap", "hmms": ["Bunya_nucleocap"]},
    ])
    # Sequence has RdRp but not Nucleocap.
    s = make_seq("s1", "ACGT")
    _attach_proteins(s, [
        _cds("polymerase", "M" * 2200, hmm_hits=[_hit("RdRP_4", 1, 280)]),
    ])
    qc = QCReport()
    kept = _run_hmm_qc_non_segmented([s], cfg, qc)
    assert kept == []
    assert qc.removed_hmm_by_marker == {"Nucleocap": 1}


def test_non_segmented_alias_only_specs_do_not_gate(make_seq):
    """Specs without ``hmms`` are advisory selection hints, not QC gates."""
    cfg = _non_segmented_cfg([
        {"name": "RdRp", "aliases": ["polymerase"]},
    ])
    s = make_seq("s1", "ACGT")
    # No HMM hits at all — alias-only spec shouldn't drop the sequence.
    _attach_proteins(s, [_cds("polymerase", "M" * 2200, hmm_hits=[])])
    qc = QCReport()
    kept = _run_hmm_qc_non_segmented([s], cfg, qc)
    assert kept == [s]


def test_non_segmented_no_hmm_specs_is_noop(make_seq):
    """When no spec carries hmms, the step is a no-op."""
    cfg = _non_segmented_cfg([])
    s = make_seq("s1", "ACGT")
    _attach_proteins(s, [_cds("anything", "M" * 100, hmm_hits=[])])
    qc = QCReport()
    kept = _run_hmm_qc_non_segmented([s], cfg, qc)
    assert kept == [s]
    assert qc.removed_hmm_failed == 0
