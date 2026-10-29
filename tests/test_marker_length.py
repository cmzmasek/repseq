"""Marker ``length_aa`` bounds (v0.76.0): a cluster marker's protein must lie
inside its range, enforced in EVERY picker.

The motivating data is real (global Alphainfluenzavirus 0.95 run,
2026-10-07): 197 H5N8 isolates from one submission batch carry a 702 aa PA
and a 481 aa NP — partial records truncated at the sequence edge — while
real NA length polymorphism (the 449 aa avian N1 stalk deletion) must
survive. The fixtures reproduce that shape.

Design points pinned here (agreed with the designer 2026-10-28):

* the bound acts at clustering ("step 1"), so the per-protein FASTA, the
  trees, the supermatrix and the coverage report all show exactly the CDS
  clustering used;
* a separate QC stage after the HMM gate drops the record, so a partial NP
  never reads as "marker fails HMM identity";
* ``removed_marker_length`` counts each isolate once, while the per-marker
  breakdown counts EVERY failing marker (shows which bound bites);
* a bounded alias-only marker never falls back to "the longest CDS";
* ``overrides.protect_qc`` (stage ``marker_length``) keeps the isolate and
  every picker then uses its real, out-of-range protein.

Protein dicts mirror the production schema (``protein_id``, ``product``,
``length``, ``sequence``, ``hmm_hits`` with ``target`` / ``dom_evalue`` /
``ali_from`` / ``ali_to`` / ``passing``).
"""
from __future__ import annotations

import csv
from pathlib import Path
from unittest.mock import patch

import pytest

from repseq.cli import _run_marker_length_qc
from repseq.clustering.marker import (
    _format_failure_reason,
    describe_length_miss,
    populate_protein_sequences,
    select_concatenated_markers,
    select_marker_protein,
)
from repseq.config import load_config, validate_config
from repseq.models import QCReport, RunResult, Sequence, SequenceType, TaxonomyInfo
from repseq.output.report import (
    _coverage_data_per_taxon,
    write_per_protein_fastas,
)
from repseq.output.summary import _render_qc
from repseq.phylo import per_protein
from repseq.phylo.per_protein import (
    build_extra_pick_context,
    build_marker_bounds_context,
    collect_marker_length_bounds,
    collect_marker_specs,
    pick_marker_cds,
)
from repseq.segmented.completeness import (
    build_concatenated_sequences,
    group_isolate_records,
    segment_marker_specs,
)


# ---------------------------------------------------------------------------
# Builders (production schema)
# ---------------------------------------------------------------------------

def _hit(target, e=1e-20, ali_from=1, ali_to=100, passing=True):
    return {
        "target": target, "dom_evalue": e,
        "ali_from": ali_from, "ali_to": ali_to, "passing": passing,
    }


def _prot(pid, product, n_aa, hits=()):
    return {
        "protein_id": pid, "product": product,
        "length": n_aa, "sequence": "M" * n_aa, "hmm_hits": list(hits),
    }


def _pa(iso, n=716):
    return _prot(f"{iso}_PA", "polymerase PA", n, [_hit("Flu_PA")])


def _np(iso, n=498, suffix=""):
    return _prot(f"{iso}_NP{suffix}", "nucleocapsid protein", n, [_hit("Flu_NP")])


def _ha(iso, n=566):
    return _prot(f"{iso}_HA", "hemagglutinin", n, [_hit("Hemagglutinin")])


def _segrec(iso, segment, proteins):
    return Sequence(
        id=f"{iso}_{segment}", header=f"{iso}_{segment}", sequence="ACGT" * 30,
        seq_type=SequenceType.NUCLEOTIDE, accession=f"{iso}_{segment}",
        segment=segment, isolate_id=iso, proteins=proteins,
    )


def _isolate(iso, *, pa=716, np=498, ha=566, extra_np=None):
    nps = [_np(iso, np)] + ([_np(iso, extra_np, "b")] if extra_np else [])
    return [
        _segrec(iso, "3", [_pa(iso, pa)]),
        _segrec(iso, "4", [_ha(iso, ha)]),
        _segrec(iso, "5", nps),
    ]


def _concat_rep(iso, **kw):
    segs = _isolate(iso, **kw)
    return Sequence(
        id=f"CONCAT|{iso}", header=iso, sequence="ACGT" * 90,
        seq_type=SequenceType.NUCLEOTIDE, isolate_id=iso, concat_segments=segs,
        taxonomy=TaxonomyInfo(genus="Alphainfluenzavirus", family="Orthomyxoviridae"),
    )


PA_BOUNDS = [705, 740]
NP_BOUNDS = [485, 510]


def _seg_cfg(*, protect=None, hmm_active=True):
    cfg = {
        "phylo": {"tool": "fasttree", "per_protein": {"min_taxa": 3}},
        "segmented": {
            "enabled": True, "virus": "FluA",
            "viruses": {"FluA": {
                "expected_segments": 3, "segments": ["3", "4", "5"],
                "isolate_regex": r"(?P<isolate>.+)",
                "segment_markers": {
                    "3": {"hmms": ["Flu_PA"], "length_aa": PA_BOUNDS},
                    "4": {"hmms": ["Hemagglutinin"]},  # deliberately unbounded
                    "5": {"hmms": ["Flu_NP"], "length_aa": NP_BOUNDS},
                },
            }},
        },
        "hmm": {"multidomain_overlap_tolerance": 30},
    }
    if hmm_active:
        cfg["_hmm_runtime"] = {"active": True}
    if protect:
        cfg["overrides"] = {
            "ids": list(protect), "protect_qc": True,
            "protect_stages": ["marker_length"],
        }
    return cfg


NP_SPEC = {"name": "NP", "hmms": ["Flu_NP"], "length_aa": NP_BOUNDS}


# ---------------------------------------------------------------------------
# The clustering picker (select_marker_protein)
# ---------------------------------------------------------------------------

def test_in_range_cds_wins_over_a_longer_out_of_range_one():
    proteins = [_np("x", 498), _np("x", 560, "fused")]
    marker, failure = select_marker_protein(proteins, [NP_SPEC], hmm_active=True)
    assert failure is None
    assert marker["protein_id"] == "x_NP"


def test_partial_marker_is_a_length_failure_not_an_hmm_failure():
    marker, failure = select_marker_protein([_np("x", 481)], [NP_SPEC], hmm_active=True)
    assert marker is None
    assert failure.reason == "length_out_of_bounds"
    assert failure.length_aa == 481
    assert failure.bounds == (485, 510)
    assert _format_failure_reason(failure) == "marker_length:NP:481aa_below_485"


def test_no_hmm_hit_is_still_an_hmm_failure():
    bare = _prot("x_NP", "nucleocapsid protein", 481)
    _marker, failure = select_marker_protein([bare], [NP_SPEC], hmm_active=True)
    assert failure.reason == "hmm_failed"


def test_ignore_bounds_picks_the_out_of_range_protein():
    marker, failure = select_marker_protein(
        [_np("x", 481)], [NP_SPEC], hmm_active=True, ignore_bounds=True,
    )
    assert failure is None and len(marker["sequence"]) == 481


def test_bounded_alias_only_marker_never_falls_back_to_the_longest_cds():
    spec = {"name": "G", "aliases": ["attachment glycoprotein"], "length_aa": [250, 350]}
    proteins = [
        _prot("L", "large protein", 2165),
        _prot("F", "fusion glycoprotein", 574),
        _prot("G", "attachment glycoprotein", 420),
    ]
    marker, failure = select_marker_protein(proteins, [spec], hmm_active=False)
    assert marker is None
    assert failure.reason == "length_out_of_bounds"
    assert _format_failure_reason(failure) == "marker_length:G:420aa_above_350"
    # No alias match at all: still no fallback, and it is reported as a
    # missing marker, not as a length failure (no protein was identified).
    marker, failure = select_marker_protein(proteins[:2], [spec], hmm_active=False)
    assert marker is None
    assert failure.reason == "no_alias_match"
    assert _format_failure_reason(failure) == "no_marker_protein_for_clustering"


def test_unmatched_bounded_alias_is_not_counted_by_the_length_stage():
    cfg = load_config(None)
    cfg["clustering"]["cluster_protein"] = [
        {"name": "NP", "aliases": ["nucleoprotein"], "length_aa": [485, 510]},
    ]
    seq = Sequence(
        id="x", header="x", sequence="ACGT", seq_type=SequenceType.NUCLEOTIDE,
        accession="x", proteins=[_prot("x_NP", "nucleocapsid protein", 498)],
    )
    qc = QCReport()
    assert _run_marker_length_qc([seq], cfg, qc, None) == [seq]
    assert qc.removed_marker_length == 0


def test_unbounded_alias_only_marker_keeps_the_legacy_fallback():
    spec = {"name": "G", "aliases": ["attachment glycoprotein"]}
    proteins = [_prot("L", "large protein", 2165), _prot("F", "fusion glycoprotein", 574)]
    marker, failure = select_marker_protein(proteins, [spec], hmm_active=False)
    assert failure is None and marker["protein_id"] == "L"


def test_hmm_only_spec_without_the_hmm_tier_falls_back_inside_the_range():
    """HMM tier soft-failed: legacy longest-CDS fallback, but in range."""
    proteins = [_prot("PB1", "polymerase PB1", 757), _prot("F2", "PB1-F2 protein", 90)]
    spec = {"name": "2", "hmms": ["Flu_PB1"], "length_aa": [80, 110]}
    marker, failure = select_marker_protein(proteins, [spec], hmm_active=False)
    assert failure is None and marker["protein_id"] == "F2"
    spec = {"name": "2", "hmms": ["Flu_PB1"], "length_aa": [800, 900]}
    marker, failure = select_marker_protein(proteins, [spec], hmm_active=False)
    assert marker is None and failure.reason == "length_out_of_bounds"


def test_concatenated_markers_keep_their_bounds():
    """select_concatenated_markers re-feeds normalised specs — the bound
    must survive that round trip."""
    proteins = [_np("x", 481), _prot("S", "spike", 1273, [_hit("CoV_S")])]
    specs = [{"name": "S", "hmms": ["CoV_S"]}, NP_SPEC]
    markers, failure = select_concatenated_markers(proteins, specs, hmm_active=True)
    assert markers is None and failure.reason == "length_out_of_bounds"


@pytest.mark.parametrize("length,bounds,expected", [
    (481, (485, 510), "481aa_below_485"),
    (820, (None, 790), "820aa_above_790"),
    (None, (485, 510), "no_candidate"),
])
def test_describe_length_miss(length, bounds, expected):
    assert describe_length_miss(length, bounds) == expected


# ---------------------------------------------------------------------------
# Clustering and output pickers agree
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("lengths", [
    [498], [481], [560], [498, 560], [481, 498, 505], [], [481, 560],
])
@pytest.mark.parametrize("kind", ["hmm", "alias"])
def test_clustering_and_output_pickers_pick_the_same_cds(lengths, kind):
    proteins = [_np("x", n, f"_{i}") for i, n in enumerate(lengths)]
    if kind == "hmm":
        spec = NP_SPEC
        tokens, aliases = ["Flu_NP"], []
    else:
        spec = {"name": "NP", "aliases": ["nucleocapsid"], "length_aa": NP_BOUNDS}
        tokens, aliases = [], ["nucleocapsid"]
    clustered, _f = select_marker_protein(proteins, [spec], hmm_active=True)
    shown = pick_marker_cds(
        proteins, tokens, aliases, True, overlap_tolerance=30,
        bounds=tuple(NP_BOUNDS),
    )
    assert clustered is shown


def test_exempt_output_pick_matches_the_protected_clustering_pick():
    proteins = [_np("x", 481)]
    clustered, _ = select_marker_protein(
        proteins, [NP_SPEC], hmm_active=True, ignore_bounds=True,
    )
    shown = pick_marker_cds(
        proteins, ["Flu_NP"], [], True, bounds=tuple(NP_BOUNDS), exempt=True,
    )
    assert clustered is shown is proteins[0]
    assert pick_marker_cds(proteins, ["Flu_NP"], [], True, bounds=tuple(NP_BOUNDS)) is None


# ---------------------------------------------------------------------------
# The marker-length QC stage — segmented
# ---------------------------------------------------------------------------

def _details(qc):
    return {d["id"]: d["reason"] for d in qc.details}


def test_partial_pa_and_np_drop_the_isolate_once_and_both_markers_count():
    seqs = _isolate("good") + _isolate("h5n8", pa=702, np=481)
    qc = QCReport()
    kept = _run_marker_length_qc(seqs, _seg_cfg(), qc, None)
    assert {s.isolate_id for s in kept} == {"good"}
    assert qc.removed_marker_length == 1                     # isolates, once
    assert qc.removed_marker_length_by_marker == {"3": 1, "5": 1}  # every marker
    assert qc.removed_hmm_failed == 0
    assert _details(qc) == {
        "h5n8_3": "marker_length:3:702aa_below_705",
        "h5n8_4": "marker_length_sibling:3:702aa_below_705",
        "h5n8_5": "marker_length:5:481aa_below_485",
    }


def test_unbounded_segment_is_not_judged():
    seqs = _isolate("short_ha", ha=467)  # HA has no length_aa
    qc = QCReport()
    assert _run_marker_length_qc(seqs, _seg_cfg(), qc, None) == seqs
    assert qc.removed_marker_length == 0


def test_missing_marker_is_left_to_the_hmm_gate():
    seqs = _isolate("no_np")
    seqs[2].proteins = [_prot("no_np_X", "hypothetical protein", 120)]
    qc = QCReport()
    assert len(_run_marker_length_qc(seqs, _seg_cfg(), qc, None)) == 3
    assert qc.removed_marker_length == 0


def test_records_without_isolate_id_are_grouped_by_the_regex_like_assembly():
    """No GenBank isolate id: the gate groups by isolate_regex exactly as
    assembly will, so the isolate is judged here in every alphabet."""
    seqs = _isolate("floater", np=481)
    for s in seqs:
        s.isolate_id = None
    cfg = _seg_cfg()
    cfg["segmented"]["viruses"]["FluA"]["isolate_regex"] = r"^(?P<isolate>[^_]+)_"
    qc = QCReport()
    assert _run_marker_length_qc(seqs, cfg, qc, None) == []
    assert qc.removed_marker_length == 1
    assert _details(qc)["floater_5"] == "marker_length:5:481aa_below_485"


def test_records_whose_isolate_cannot_be_resolved_are_left_alone():
    seqs = _isolate("floater", np=481)
    for s in seqs:
        s.isolate_id = None
    cfg = _seg_cfg()
    cfg["segmented"]["viruses"]["FluA"]["isolate_regex"] = r"(?P<isolate>NOMATCH)"
    assert len(_run_marker_length_qc(seqs, cfg, QCReport(), None)) == 3


def test_no_bounds_configured_is_a_no_op():
    cfg = _seg_cfg()
    for spec in cfg["segmented"]["viruses"]["FluA"]["segment_markers"].values():
        spec.pop("length_aa", None)
    seqs = _isolate("h5n8", pa=702, np=481)
    assert _run_marker_length_qc(seqs, cfg, QCReport(), None) is seqs


def test_protected_isolate_is_kept_and_recorded():
    seqs = _isolate("ref", np=481)
    qc = QCReport()
    kept = _run_marker_length_qc(seqs, _seg_cfg(protect=["ref"]), qc, None)
    assert len(kept) == 3
    assert qc.removed_marker_length == 0
    stages = {(p["id"], p["stage"], p["reason"]) for p in qc.protected}
    assert ("ref_5", "marker_length", "marker_length:5:481aa_below_485") in stages
    assert ("ref_3", "marker_length", "marker_length_sibling:5:481aa_below_485") in stages


def test_skipped_with_a_warning_under_no_resolve(capsys):
    seqs = _isolate("h5n8", np=481)
    seqs[0].proteins = None
    assert _run_marker_length_qc(seqs, _seg_cfg(), QCReport(), None) is seqs
    assert "--no-resolve was set" in capsys.readouterr().err


# ---------------------------------------------------------------------------
# The marker-length QC stage — non-segmented
# ---------------------------------------------------------------------------

def _nonseg_cfg():
    return {
        "clustering": {"cluster_protein": [
            {"name": "S", "hmms": ["CoV_S"], "length_aa": [1100, 1400]},
            {"name": "N", "hmms": ["CoV_N"], "length_aa": [350, 450]},
            {"name": "E", "hmms": ["CoV_E"]},  # unbounded: never judged here
        ]},
        "_hmm_runtime": {"active": True},
    }


def _cov(sid, s_len, n_len):
    return Sequence(
        id=sid, header=sid, sequence="ACGT" * 30,
        seq_type=SequenceType.NUCLEOTIDE, accession=sid,
        proteins=[
            _prot(f"{sid}_S", "spike", s_len, [_hit("CoV_S")]),
            _prot(f"{sid}_N", "nucleocapsid", n_len, [_hit("CoV_N")]),
        ],
    )


def test_non_segmented_requires_every_bounded_marker():
    seqs = [_cov("ok", 1273, 419), _cov("short_n", 1273, 300), _cov("both", 900, 300)]
    qc = QCReport()
    kept = _run_marker_length_qc(seqs, _nonseg_cfg(), qc, None)
    assert [s.id for s in kept] == ["ok"]
    assert qc.removed_marker_length == 2                      # sequences
    assert qc.removed_marker_length_by_marker == {"S": 1, "N": 2}
    assert _details(qc) == {
        "short_n": "marker_length:N:300aa_below_350",
        "both": "marker_length:S:900aa_below_1100",
    }


def test_non_segmented_protection():
    cfg = _nonseg_cfg()
    cfg["overrides"] = {"ids": ["short_n"], "protect_qc": True, "protect_stages": "all"}
    qc = QCReport()
    kept = _run_marker_length_qc([_cov("short_n", 1273, 300)], cfg, qc, None)
    assert len(kept) == 1 and qc.removed_marker_length == 0
    assert qc.protected[0]["stage"] == "marker_length"


# ---------------------------------------------------------------------------
# Clustering-side consumers
# ---------------------------------------------------------------------------

def test_populate_routes_length_failures_and_honours_protection():
    from repseq.overrides import ProtectionPolicy

    specs = [{"name": "N", "hmms": ["CoV_N"], "length_aa": [350, 450]}]
    qc = QCReport()
    kept = populate_protein_sequences(
        [_cov("short_n", 1273, 300)], specs, qc, hmm_active=True,
    )
    assert kept == [] and qc.removed_marker_length == 1
    assert qc.removed_proteins == 0
    assert _details(qc) == {"short_n": "marker_length:N:300aa_below_350"}

    policy = ProtectionPolicy.from_cfg(
        {"overrides": {"ids": ["short_n"], "protect_qc": True}}
    )
    kept = populate_protein_sequences(
        [_cov("short_n", 1273, 300)], specs, QCReport(), hmm_active=True,
        policy=policy,
    )
    assert len(kept) == 1 and len(kept[0].protein_sequence) == 300


def test_assembly_backstop_drops_and_counts_by_length():
    vcfg = _seg_cfg()["segmented"]["viruses"]["FluA"]
    qc = QCReport()
    out = build_concatenated_sequences(
        {"h5n8": _isolate("h5n8", np=481), "good": _isolate("good")},
        segment_names=["3", "4", "5"], require_protein=True, report=qc,
        segment_markers=vcfg["segment_markers"], hmm_active=True,
    )
    assert [s.isolate_id for s in out] == ["good"]
    assert qc.removed_marker_length == 1
    assert qc.removed_marker_length_by_marker == {"5": 1}
    assert qc.removed_incomplete_isolates == 0 and qc.removed_hmm_failed == 0
    assert _details(qc)["h5n8_5"] == "marker_length:5:481aa_below_485"
    assert _details(qc)["h5n8_3"] == "marker_length_sibling:5:481aa_below_485"


def test_assembly_backstop_uses_the_real_protein_for_a_protected_isolate():
    from repseq.overrides import ProtectionPolicy

    cfg = _seg_cfg(protect=["ref"])
    vcfg = cfg["segmented"]["viruses"]["FluA"]
    qc = QCReport()
    out = build_concatenated_sequences(
        {"ref": _isolate("ref", np=481)},
        segment_names=["3", "4", "5"], require_protein=True, report=qc,
        segment_markers=vcfg["segment_markers"], hmm_active=True,
        policy=ProtectionPolicy.from_cfg(cfg),
    )
    assert len(out) == 1
    assert len(out[0].protein_sequence) == 716 + 566 + 481
    assert any(p["stage"] == "marker_length" for p in qc.protected)


def test_assembly_takes_the_in_range_cds_when_a_longer_one_exists():
    vcfg = _seg_cfg()["segmented"]["viruses"]["FluA"]
    out = build_concatenated_sequences(
        {"dup": _isolate("dup", extra_np=560)},
        segment_names=["3", "4", "5"], require_protein=True, report=QCReport(),
        segment_markers=vcfg["segment_markers"], hmm_active=True,
    )
    assert out[0].marker_protein_ids == ["dup_PA", "dup_HA", "dup_NP"]


# ---------------------------------------------------------------------------
# Output-side consumers show exactly what clustering used
# ---------------------------------------------------------------------------

def test_collect_marker_length_bounds_uses_the_family_labels():
    bounds = collect_marker_length_bounds(_seg_cfg())
    labels = {(s[0], s[3]) for s in collect_marker_specs(_seg_cfg())}
    assert bounds == {("3_Flu_PA", "3"): (705, 740), ("5_Flu_NP", "5"): (485, 510)}
    assert set(bounds) <= labels
    assert collect_marker_length_bounds(_nonseg_cfg()) == {
        ("S", None): (1100, 1400), ("N", None): (350, 450),
    }


def _reps():
    return [
        _concat_rep("dup", extra_np=560),  # in-range 498 beside a 560 fusion
        _concat_rep("ok1"),
        _concat_rep("ok2", np=497),
    ]


def _fasta_ids(path: Path) -> list[str]:
    return [ln[1:].split()[0] for ln in path.read_text().splitlines() if ln.startswith(">")]


def test_per_protein_fasta_shows_the_clustering_cds(tmp_path):
    result = RunResult(mode="global", representatives=_reps())
    written = write_per_protein_fastas(result, _seg_cfg(), None, tmp_path, "pfx")
    by_name = {p.name: p for p in written}
    assert _fasta_ids(by_name["pfx_5_Flu_NP.fasta"]) == ["dup_NP", "ok1_NP", "ok2_NP"]


def test_per_protein_fasta_shows_a_protected_isolates_real_protein(tmp_path):
    reps = [_concat_rep("ref", np=481), _concat_rep("plain", np=481), _concat_rep("ok")]
    result = RunResult(mode="global", representatives=reps)
    written = write_per_protein_fastas(
        result, _seg_cfg(protect=["ref"]), None, tmp_path, "pfx",
    )
    by_name = {p.name: p for p in written}
    # "plain" could never reach here in a real run (QC drops it); the
    # output picker alone simply finds no in-range NP for it.
    assert _fasta_ids(by_name["pfx_5_Flu_NP.fasta"]) == ["ref_NP", "ok_NP"]


def test_per_protein_trees_use_the_same_picks(tmp_path):
    calls: dict[str, dict] = {}

    def _fake_build_tree(leaf_reps, bodies, **kw):
        calls[kw["file_prefix"]] = {k: len(v) for k, v in bodies.items()}
        return []

    with patch.object(per_protein, "_build_tree", _fake_build_tree):
        per_protein.run_per_protein_phylogeny(_reps(), _seg_cfg(), tmp_path, "pfx")
    assert calls["5_Flu_NP"] == {"CONCAT|dup": 498, "CONCAT|ok1": 498, "CONCAT|ok2": 497}


def test_coverage_report_counts_the_clustering_cds():
    cfg = _seg_cfg()
    specs = [(l, t, a, s, True) for l, t, a, s in collect_marker_specs(cfg)]
    ctx = build_extra_pick_context(cfg, True)
    _totals, lengths = _coverage_data_per_taxon(
        _reps(), specs, True, "genus", 30, ctx,
    )
    np_idx = [s[0] for s in specs].index("5_Flu_NP")
    assert sorted(lengths["Alphainfluenzavirus"][np_idx]) == [497, 498, 498]


def test_markers_claim_their_in_range_pick():
    cfg = _seg_cfg()
    ctx = build_extra_pick_context(cfg, True)
    proteins = [_np("dup", 498), _np("dup", 560, "b")]
    claimed = ctx.claimed(proteins, "5")
    assert set(claimed) == {id(proteins[0])}


def test_marker_bounds_context_exemption_follows_the_protection_list():
    ctx = build_marker_bounds_context(_seg_cfg(protect=["ref"]))
    assert ctx.exempt(_concat_rep("ref"))
    assert not ctx.exempt(_concat_rep("other"))
    # Protection off for this stage ⇒ never exempt.
    cfg = _seg_cfg(protect=["ref"])
    cfg["overrides"]["protect_stages"] = ["hmm"]
    assert not build_marker_bounds_context(cfg).exempt(_concat_rep("ref"))


# ---------------------------------------------------------------------------
# Reporting surfaces
# ---------------------------------------------------------------------------

def test_console_tally_shows_the_stage():
    qc = QCReport(total_input=12)
    qc.removed_marker_length = 1
    qc.final_survivors = 2
    qc.final_survivors_unit = "isolates"
    qc.display = {"segmented": True, "length_isolates_pre_assembly": 1}
    assert "marker protein outside length range" in qc.summary()
    assert "−1 isolates" in qc.summary()


def test_methods_summary_names_the_ranges_and_the_count():
    qc = QCReport(total_input=6)
    qc.removed_marker_length = 1
    qc.removed_marker_length_by_marker = {"3": 1, "5": 1}
    text = _render_qc(qc, _seg_cfg())
    assert (
        "configured length ranges (amino acids, inclusive): "
        "segment `3` 705–740 aa; segment `5` 485–510 aa" in text
    )
    assert "**1** isolates were removed this way" in text
    assert "failing markers: `3` (1), `5` (1)" in text
    assert "The same ranges decided which coding sequence" in text


def test_methods_summary_mentions_bounds_even_when_nothing_dropped():
    text = _render_qc(QCReport(total_input=6), _seg_cfg())
    assert "no isolates were removed this way" in text


# ---------------------------------------------------------------------------
# Config validation
# ---------------------------------------------------------------------------

def _validation_cfg(sm_np):
    cfg = load_config(None)
    cfg["segmented"] = {
        "enabled": True, "virus": "FluA",
        "viruses": {"FluA": {
            "expected_segments": 1, "segments": ["5"],
            "isolate_regex": r"(?P<isolate>X)",
            "segment_markers": {"5": sm_np},
        }},
    }
    return cfg


def test_validate_accepts_length_aa_on_markers():
    for bounds in ([485, 510], [485, None], [None, 510], None):
        sm = {"hmms": ["Flu_NP"], "length_aa": bounds}
        assert validate_config(_validation_cfg(sm)) == [], bounds
    cfg = load_config(None)
    cfg["clustering"]["cluster_protein"] = [
        {"name": "N", "hmms": ["CoV_N"], "length_aa": [350, 450]},
    ]
    assert validate_config(cfg) == []


def test_validate_rejects_malformed_marker_length_aa():
    errs = validate_config(_validation_cfg({"hmms": ["Flu_NP"], "length_aa": [510, 485]}))
    assert any("min (510) is larger than max (485)" in e for e in errs), errs
    cfg = load_config(None)
    cfg["clustering"]["cluster_protein"] = [
        {"name": "N", "hmms": ["CoV_N"], "length_aa": "350-450"},
    ]
    assert any("two-item list" in e for e in validate_config(cfg))


def test_validate_accepts_the_marker_length_protect_stage():
    cfg = load_config(None)
    cfg["overrides"]["protect_stages"] = ["hmm", "marker_length"]
    assert validate_config(cfg) == []


# ---------------------------------------------------------------------------
# "How many removed isolates had every segment" — both isolate-level gates
# ---------------------------------------------------------------------------
#
# Both gates run before isolate assembly, so part of what they remove was
# already missing a segment. Measured on the 2026-10-28 H5N1 rerun: 72 % of
# the failing markers belonged to isolates that were incomplete anyway (HA:
# 423 vs 14), so the headline count alone overstated the cost.

from repseq.cli import _count_isolates_with_all_segments, _run_hmm_qc_segmented


def _incomplete(iso, **kw):
    """The isolate minus its HA record (an HA-less submission)."""
    return [s for s in _isolate(iso, **kw) if s.segment != "4"]


def test_count_isolates_with_all_segments():
    vcfg = _seg_cfg()["segmented"]["viruses"]["FluA"]
    groups = group_isolate_records(_isolate("a") + _incomplete("b"), vcfg)
    assert _count_isolates_with_all_segments(groups, ["a", "b"], ["3", "4", "5"]) == 1
    assert _count_isolates_with_all_segments(groups, ["a", "b"], []) is None


def test_complete_count_uses_canonical_segments_and_normalised_ids():
    """A /segment alias and an isolate name differing only in case or
    whitespace still count as one complete isolate, as at assembly."""
    cfg = _seg_cfg()
    vcfg = cfg["segmented"]["viruses"]["FluA"]
    vcfg["segment_aliases"] = {"4": ["hemagglutinin segment"]}
    seqs = _isolate("A/Duck/1", np=481)
    seqs[1].segment = "hemagglutinin segment"
    seqs[2].isolate_id = "a/duck/1"
    groups = group_isolate_records(seqs, vcfg)
    assert list(groups) == ["a/duck/1"]
    assert _count_isolates_with_all_segments(groups, groups, ["3", "4", "5"]) == 1


def test_length_gate_reports_how_many_removed_isolates_were_complete(capsys):
    seqs = _isolate("full", np=481) + _incomplete("part", np=481) + _isolate("ok")
    qc = QCReport()
    _run_marker_length_qc(seqs, _seg_cfg(), qc, None)
    assert qc.removed_marker_length == 2
    assert qc.display["length_complete_isolates"] == 1
    assert "1 with all 3 segments" in capsys.readouterr().out


def test_hmm_gate_reports_how_many_removed_isolates_were_complete(capsys):
    full, part = _isolate("full"), _incomplete("part")
    for seg in (full[2], part[1]):  # strip the NP's HMM hit on both
        seg.proteins = [dict(p, hmm_hits=[]) for p in seg.proteins]
    qc = QCReport()
    _run_hmm_qc_segmented(full + part + _isolate("ok"), _seg_cfg(), qc)
    assert qc.removed_hmm_failed == 2
    assert qc.display["hmm_complete_isolates"] == 1
    assert "1 with all 3 segments" in capsys.readouterr().out


def test_console_tally_shows_the_complete_share_for_both_gates():
    qc = QCReport(total_input=30)
    qc.removed_hmm_failed = 4
    qc.removed_marker_length = 5
    qc.final_survivors = 2
    qc.final_survivors_unit = "isolates"
    qc.display = {
        "segmented": True, "n_segments": 3,
        "hmm_isolates_pre_assembly": 4, "hmm_complete_isolates": 1,
        "length_isolates_pre_assembly": 5, "length_complete_isolates": 2,
    }
    text = qc.summary()
    assert "−4 isolates (1 with all 3 segments)" in text
    assert "−5 isolates (2 with all 3 segments)" in text


def test_methods_summary_states_the_complete_share():
    qc = QCReport(total_input=30)
    qc.removed_marker_length = 3          # 2 at the gate + 1 at the backstop
    qc.removed_marker_length_by_marker = {"5": 3}
    qc.removed_hmm_failed = 4
    qc.removed_hmm_by_marker = {"5:Flu_NP": 4}
    qc.display = {
        "length_isolates_pre_assembly": 2, "length_complete_isolates": 1,
        "hmm_isolates_pre_assembly": 4, "hmm_complete_isolates": 1,
    }
    text = _render_qc(qc, _seg_cfg())
    # Backstop removals are complete by construction, so 1 + 1 = 2.
    assert (
        "Of these 3 isolates, **2** carried all 3 expected segments; the "
        "others were already missing a segment" in text
    )
    assert "Of these 4 isolates, **1** carried all 3 expected segments" in text


def test_methods_summary_omits_the_share_without_the_figures():
    qc = QCReport(total_input=30)
    qc.removed_hmm_failed = 4
    qc.removed_hmm_by_marker = {"5:Flu_NP": 4}
    assert "carried all" not in _render_qc(qc, _seg_cfg())


# ---------------------------------------------------------------------------
# Review fixes (2026-10-28): the gates judge what assembly keeps
# ---------------------------------------------------------------------------

def _dup_np(iso, n_aa, nt_len, suffix):
    """A second segment-5 record for ``iso`` (a resubmission)."""
    return Sequence(
        id=f"{iso}_5{suffix}", header=f"{iso}_5{suffix}", sequence="A" * nt_len,
        seq_type=SequenceType.NUCLEOTIDE, accession=f"{iso}_5{suffix}",
        segment="5", isolate_id=iso, proteins=[_np(iso, n_aa, suffix)],
    )


def _with_np_records(iso, *nps):
    """``_isolate`` with its NP record replaced by ``nps`` (in that order)."""
    return [s for s in _isolate(iso) if s.segment != "5"] + list(nps)


def test_partial_duplicate_does_not_remove_an_isolate_with_a_complete_record():
    # Complete NP record (1,500 nt) plus an older, shorter partial one.
    seqs = _with_np_records("dup", _dup_np("dup", 498, 1500, "a"), _dup_np("dup", 481, 1445, "b"))
    qc = QCReport()
    assert _run_marker_length_qc(seqs, _seg_cfg(), qc, None) == seqs
    assert qc.removed_marker_length == 0


def test_the_record_assembly_keeps_decides_when_it_is_the_partial_one():
    seqs = _with_np_records("dup", _dup_np("dup", 498, 1400, "a"), _dup_np("dup", 481, 1445, "b"))
    qc = QCReport()
    assert _run_marker_length_qc(seqs, _seg_cfg(), qc, None) == []
    assert qc.removed_marker_length == 1
    assert qc.removed_marker_length_by_marker == {"5": 1}
    d = _details(qc)
    assert d["dup_5b"] == "marker_length:5:481aa_below_485"
    # The discarded duplicate is a sibling, not a second failure.
    assert d["dup_5a"] == "marker_length_sibling:5:481aa_below_485"


def test_hmm_gate_judges_only_the_record_assembly_keeps():
    good = _dup_np("dup", 498, 1500, "a")
    bad = _dup_np("dup", 498, 1445, "b")
    bad.proteins = [dict(p, hmm_hits=[]) for p in bad.proteins]  # no NP hit
    seqs = _with_np_records("dup", good, bad)
    qc = QCReport()
    assert _run_hmm_qc_segmented(seqs, _seg_cfg(), qc) == seqs
    assert qc.removed_hmm_failed == 0
    # Swap the lengths: now the failing record is the one assembly keeps.
    good.sequence, bad.sequence = "A" * 1445, "A" * 1500
    qc = QCReport()
    assert _run_hmm_qc_segmented(seqs, _seg_cfg(), qc) == []
    assert qc.removed_hmm_failed == 1


def test_gates_resolve_named_segments_like_assembly():
    """GenBank /segment '5' on a config naming its segments: the bound on
    segment_markers['NP'] must still apply."""
    cfg = _seg_cfg()
    vcfg = cfg["segmented"]["viruses"]["FluA"]
    vcfg["segments"] = ["PB2", "PB1", "PA", "HA", "NP"]
    vcfg["segment_markers"] = {"NP": {"hmms": ["Flu_NP"], "length_aa": NP_BOUNDS}}
    seqs = _isolate("h5n8", np=481)
    qc = QCReport()
    assert _run_marker_length_qc(seqs, cfg, qc, None) == []
    assert _details(qc)["h5n8_5"] == "marker_length:NP:481aa_below_485"


def test_exact_segment_name_wins_over_the_position_lookup():
    from repseq.segmented.completeness import identify_segment

    seq = _isolate("x")[0]  # GenBank /segment "3"
    assert identify_segment(seq, ["3", "4", "5"]) == "3"
    assert identify_segment(seq, ["PB2", "PB1", "PA"]) == "PA"


def test_marker_length_stage_does_not_refetch_proteins_that_carry_hmm_hits():
    """One record lacking proteins must not trigger a re-fetch of the others
    (the cache hands back fresh dicts without hmm_hits)."""
    from repseq.cli import _attach_missing_proteins

    seqs = _isolate("ok")
    orphan = Sequence(
        id="orphan", header="orphan", sequence="ACGT",
        seq_type=SequenceType.NUCLEOTIDE, accession="ORPH1", segment="5",
    )
    fetched = []
    with patch("repseq.cli.attach_proteins", side_effect=lambda s, n: fetched.extend(s)):
        _attach_missing_proteins(seqs + [orphan], object())
    assert fetched == [orphan]
    # A record that cannot be fetched at all (no accession) never triggers
    # the --no-resolve skip.
    orphan.accession = None
    qc = QCReport()
    assert len(_run_marker_length_qc(seqs + [orphan], _seg_cfg(), qc, None)) == 4
    assert "marker_length_skipped" not in qc.display


def test_protected_isolate_with_hmm_and_length_failures_is_kept_at_assembly():
    from repseq.overrides import ProtectionPolicy

    cfg = _seg_cfg(protect=["ref"])
    cfg["overrides"]["protect_stages"] = "all"
    vcfg = cfg["segmented"]["viruses"]["FluA"]
    segs = _isolate("ref", np=481)
    segs[2].proteins = [dict(p, hmm_hits=[]) for p in segs[2].proteins]
    qc = QCReport()
    out = build_concatenated_sequences(
        {"ref": segs}, segment_names=["3", "4", "5"], require_protein=True,
        report=qc, segment_markers=vcfg["segment_markers"], hmm_active=True,
        policy=ProtectionPolicy.from_cfg(cfg),
    )
    assert len(out) == 1
    assert len(out[0].protein_sequence) == 716 + 566 + 481


def test_isolate_protected_only_against_hmm_drops_by_length_at_assembly():
    from repseq.overrides import ProtectionPolicy

    cfg = _seg_cfg(protect=["ref"])
    cfg["overrides"]["protect_stages"] = ["hmm"]
    vcfg = cfg["segmented"]["viruses"]["FluA"]
    segs = _isolate("ref", np=481)
    segs[2].proteins = [dict(p, hmm_hits=[]) for p in segs[2].proteins]
    qc = QCReport()
    out = build_concatenated_sequences(
        {"ref": segs}, segment_names=["3", "4", "5"], require_protein=True,
        report=qc, segment_markers=vcfg["segment_markers"], hmm_active=True,
        policy=ProtectionPolicy.from_cfg(cfg),
    )
    assert out == []
    assert qc.removed_marker_length == 1 and qc.removed_hmm_failed == 0


def test_split_config_merges_cluster_protein_hmms_into_the_bounded_spec():
    specs = segment_marker_specs(
        "S",
        {"S": {"aliases": ["nucleoprotein"], "length_aa": [200, 260]}},
        {"S": [{"name": "N", "hmms": ["Bunya_nucleocap"]}]},
    )
    assert specs == [{
        "name": "S", "aliases": ["nucleoprotein"],
        "hmms": ["Bunya_nucleocap"], "length_aa": [200, 260],
    }]


def test_output_exemption_is_only_looked_up_when_the_bounded_pick_misses():
    calls = []

    def exempt():
        calls.append(1)
        return True

    np_ok = [_np("x", 498)]
    assert pick_marker_cds(np_ok, ["Flu_NP"], [], True, bounds=(485, 510), exempt=exempt)
    assert calls == []
    np_short = [_np("x", 481)]
    assert pick_marker_cds(np_short, ["Flu_NP"], [], True, bounds=(485, 510), exempt=exempt)
    assert calls == [1]


def test_methods_summary_says_when_the_ranges_were_not_applied():
    qc = QCReport(total_input=6)
    qc.display["marker_length_skipped"] = True
    text = _render_qc(qc, _seg_cfg())
    assert "but were **not applied** in this run" in text
    assert "removed this way" not in text
    assert "The same ranges decided" not in text


def test_methods_summary_is_alphabet_aware():
    cfg = _seg_cfg()
    cfg.setdefault("clustering", {})["alphabet_for_clustering"] = "nucleotide"
    text = _render_qc(QCReport(total_input=6), cfg)
    assert "Clustering itself used nucleotide sequences" in text
    assert "represented each marker in clustering" not in text
