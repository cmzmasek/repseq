"""Marker 5'-completeness (``qc.marker_completeness.require_met_start``,
v0.77.0): a cluster marker's CDS must translate from an initial methionine,
enforced in EVERY picker.

Why the test exists at all. On the global Alphainfluenzavirus 0.95 run
(2026-10-28, v0.76.0) the ``length_aa`` ranges more than halved the
representatives carrying a 5'-partial marker protein — 17.97 % → 8.10 % —
but could not close the gap, because a short N-terminal truncation stays
*inside* the range. A probe over nine runs then showed the residual is
truncation, not biology: pooling 609 non-Met marker proteins across five
clades, the first residue spans all 20 amino acids with a maximum share of
13.3 % (Leu), i.e. a random truncation point rather than the concentration a
non-AUG initiation mechanism would produce. Four clades (Coronaviridae,
Betacoronaviridae, Filoviridae, Betainfluenzavirus) sit at exactly 0 % over
4,127 proteins.

Design points pinned here (agreed with the designer 2026-10-29, one by one):

* ONE global switch, default off, on in the influenza configs — not a
  per-spec flag, which would let a config gate seven of eight segments and
  look entirely normal;
* it shares the marker-quality QC stage with ``length_aa`` but keeps its own
  reason and counter, with LENGTH taking precedence so the counters sum to
  exactly what the stage dropped;
* the protect token is ``marker_quality`` and covers BOTH checks (so
  ``select_marker_protein`` takes one ``ignore_restrictions`` flag rather
  than letting a caller half-protect a reference strain);
* it applies to ``extra_protein`` specs by parity, where the consequence is
  "the protein is absent", not a drop;
* polyprotein peptides are exempt — an internal cleavage product never
  starts with Met — which holds by construction, since the slicer has its
  own parent-CDS picker.

Protein dicts mirror the production schema (``protein_id`` / ``product`` /
``length`` / ``sequence`` / ``hmm_hits`` with ``target`` / ``dom_evalue`` /
``ali_from`` / ``ali_to`` / ``passing``).
"""
from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

from repseq.cli import _run_marker_quality_qc
from repseq.clustering.marker import (
    MARKER_QUALITY_RUNTIME,
    _format_failure_reason,
    require_met_start_active,
    populate_protein_sequences,
    require_met_start_from_cfg,
    select_marker_protein,
    starts_with_met,
)
from repseq.config import DEFAULTS, load_config, validate_config
from repseq.models import QCReport, RunResult, Sequence, SequenceType, TaxonomyInfo
from repseq.output.report import write_per_protein_fastas
from repseq.output.summary import _render_qc
from repseq.phylo import per_protein
from repseq.phylo.per_protein import (
    build_extra_pick_context,
    build_marker_bounds_context,
    pick_extra_cds,
    pick_marker_cds,
)
from repseq.segmented.completeness import build_concatenated_sequences


# ---------------------------------------------------------------------------
# Builders (production schema)
# ---------------------------------------------------------------------------

def _hit(target, e=1e-20, ali_from=1, ali_to=100, passing=True):
    return {
        "target": target, "dom_evalue": e,
        "ali_from": ali_from, "ali_to": ali_to, "passing": passing,
    }


def _prot(pid, product, n_aa, hits=(), *, first="M"):
    """A CDS whose translation starts with ``first``.

    ``first="K"`` is the 5'-partial case: GenBank gives a /translation that
    begins mid-protein because the CDS feature is partial at the 5' end.
    """
    return {
        "protein_id": pid, "product": product, "length": n_aa,
        "sequence": first + "A" * (n_aa - 1), "hmm_hits": list(hits),
    }


def _pa(iso, n=716, *, first="M"):
    return _prot(f"{iso}_PA", "polymerase PA", n, [_hit("Flu_PA")], first=first)


def _np(iso, n=498, suffix="", *, first="M"):
    return _prot(
        f"{iso}_NP{suffix}", "nucleocapsid protein", n, [_hit("Flu_NP")],
        first=first,
    )


def _ha(iso, n=566, *, first="M"):
    return _prot(f"{iso}_HA", "hemagglutinin", n, [_hit("Hemagglutinin")], first=first)


def _segrec(iso, segment, proteins, nt=1500):
    return Sequence(
        id=f"{iso}_{segment}", header=f"{iso}_{segment}", sequence="A" * nt,
        seq_type=SequenceType.NUCLEOTIDE, accession=f"{iso}_{segment}",
        segment=segment, isolate_id=iso, proteins=proteins,
    )


def _isolate(iso, *, pa=716, np=498, ha=566, pa_first="M", np_first="M",
             ha_first="M", extra_np=None):
    """Three segments: 3 (PA, bounded), 4 (HA, UNBOUNDED), 5 (NP, bounded)."""
    nps = [_np(iso, np, first=np_first)]
    if extra_np:
        nps.append(_np(iso, extra_np, "b"))
    return [
        _segrec(iso, "3", [_pa(iso, pa, first=pa_first)]),
        _segrec(iso, "4", [_ha(iso, ha, first=ha_first)]),
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


def _seg_cfg(*, protect=None, met=True, hmm_active=True, extra=None):
    cfg = {
        "qc": {"marker_completeness": {"require_met_start": met}},
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
    if extra is not None:
        cfg["segmented"]["viruses"]["FluA"]["extra_protein"] = extra
    if hmm_active:
        cfg["_hmm_runtime"] = {"active": True}
    if protect:
        cfg["overrides"] = {
            "ids": list(protect), "protect_qc": True,
            "protect_stages": ["marker_quality"],
        }
    return cfg


NP_SPEC = {"name": "NP", "hmms": ["Flu_NP"], "length_aa": NP_BOUNDS}


def _details(qc: QCReport) -> dict:
    return {d["id"]: d["reason"] for d in qc.details}


# ---------------------------------------------------------------------------
# The predicate
# ---------------------------------------------------------------------------

def test_starts_with_met_reads_the_first_residue():
    assert starts_with_met(_np("x"))
    assert not starts_with_met(_np("x", first="K"))
    # Lower-case translations occur in hand-built records.
    assert starts_with_met({"sequence": "mdsntv"})
    # No translation is NOT judged here — callers filter those out first, and
    # a second "no translation" drop would misattribute the cause.
    assert starts_with_met({"sequence": None})
    assert starts_with_met({})


def test_the_switch_is_off_by_default():
    assert require_met_start_from_cfg(DEFAULTS) is False
    assert require_met_start_from_cfg({}) is False
    assert require_met_start_from_cfg(_seg_cfg()) is True


# ---------------------------------------------------------------------------
# The clustering picker (select_marker_protein)
# ---------------------------------------------------------------------------

def test_a_5prime_partial_marker_is_rejected():
    marker, failure = select_marker_protein(
        [_np("x", 497, first="K")], [NP_SPEC],
        hmm_active=True, require_met_start=True,
    )
    assert marker is None
    assert failure.reason == "no_met_start"
    assert failure.marker_name == "NP"
    assert failure.length_aa == 497           # the CDS it would have taken
    assert _format_failure_reason(failure) == "marker_completeness:NP:no_met_start"


def test_the_switch_off_keeps_the_partial_marker():
    marker, failure = select_marker_protein(
        [_np("x", 497, first="K")], [NP_SPEC], hmm_active=True,
    )
    assert failure is None and marker["protein_id"] == "x_NP"


def test_a_complete_sibling_wins_over_a_partial_longer_one():
    """The filter runs BEFORE "longest wins", so completeness is not traded
    away for length."""
    proteins = [_np("x", 497), _np("x", 505, "b", first="K")]
    marker, failure = select_marker_protein(
        proteins, [NP_SPEC], hmm_active=True, require_met_start=True,
    )
    assert failure is None and marker["protein_id"] == "x_NP"


def test_length_takes_precedence_over_completeness():
    """A candidate that fails both is reported as a length miss, so the two
    counters sum to exactly the records the stage dropped."""
    _marker, failure = select_marker_protein(
        [_np("x", 481, first="K")], [NP_SPEC],
        hmm_active=True, require_met_start=True,
    )
    assert failure.reason == "length_out_of_bounds"


def test_one_protect_flag_waives_both_restrictions():
    marker, failure = select_marker_protein(
        [_np("x", 481, first="K")], [NP_SPEC],
        hmm_active=True, require_met_start=True, ignore_restrictions=True,
    )
    assert failure is None and len(marker["sequence"]) == 481


def test_an_alias_only_spec_is_strict_about_completeness():
    """It must not fall through to "the longest CDS", which could substitute
    an unrelated protein that happens to be complete."""
    spec = {"name": "NP", "aliases": ["nucleocapsid"]}
    proteins = [
        _np("x", 497, first="K"),
        _prot("x_other", "hypothetical protein", 600),   # complete, unrelated
    ]
    marker, failure = select_marker_protein(
        proteins, [spec], hmm_active=False, require_met_start=True,
    )
    assert marker is None and failure.reason == "no_met_start"
    # With the switch off this is the legacy alias pick.
    marker, failure = select_marker_protein(proteins, [spec], hmm_active=False)
    assert failure is None and marker["protein_id"] == "x_NP"


def test_an_unmatched_alias_still_falls_back_to_a_complete_cds():
    """An alias that matched NOTHING stays advisory (legacy behaviour); the
    fallback is just restricted to a complete CDS."""
    spec = {"name": "NP", "aliases": ["no-such-product"]}
    proteins = [_prot("x_a", "hypothetical", 300, first="K"), _prot("x_b", "orf", 200)]
    marker, failure = select_marker_protein(
        proteins, [spec], hmm_active=False, require_met_start=True,
    )
    assert failure is None and marker["protein_id"] == "x_b"
    # Nothing complete at all → a completeness failure, not a silent pick.
    marker, failure = select_marker_protein(
        [proteins[0]], [spec], hmm_active=False, require_met_start=True,
    )
    assert marker is None and failure.reason == "no_met_start"


# ---------------------------------------------------------------------------
# The QC stage (segmented)
# ---------------------------------------------------------------------------

def test_stage_drops_the_isolate_and_counts_it():
    seqs = _isolate("part", np_first="K")
    qc = QCReport()
    assert _run_marker_quality_qc(seqs, _seg_cfg(), qc, None) == []
    assert qc.removed_marker_completeness == 1          # isolates, once
    assert qc.removed_marker_completeness_by_marker == {"5": 1}
    assert qc.removed_marker_length == 0
    d = _details(qc)
    assert d["part_5"] == "marker_completeness:5:no_met_start"
    assert d["part_3"] == "marker_completeness_sibling:5:no_met_start"
    assert d["part_4"] == "marker_completeness_sibling:5:no_met_start"


def test_an_unbounded_segment_is_judged_too():
    """The requirement is global, so segment 4 — which declares no
    ``length_aa`` and is therefore invisible to the length check — is still
    gated. This is the whole point of a global switch."""
    qc = QCReport()
    assert _run_marker_quality_qc(_isolate("part", ha_first="K"), _seg_cfg(), qc, None) == []
    assert qc.removed_marker_completeness_by_marker == {"4": 1}

    # With the switch off, that same isolate passes untouched.
    seqs = _isolate("part", ha_first="K")
    qc = QCReport()
    assert _run_marker_quality_qc(seqs, _seg_cfg(met=False), qc, None) == seqs
    assert qc.removed_marker_completeness == 0


def test_a_complete_isolate_passes():
    seqs = _isolate("ok")
    qc = QCReport()
    assert _run_marker_quality_qc(seqs, _seg_cfg(), qc, None) == seqs
    assert qc.removed_marker_completeness == 0 and qc.removed_marker_length == 0


def test_an_isolate_failing_both_counts_once_under_length():
    """Precedence has to hold at the isolate level too, or the two counters
    would double-count the same isolate."""
    seqs = _isolate("both", np=481, ha_first="K")
    qc = QCReport()
    assert _run_marker_quality_qc(seqs, _seg_cfg(), qc, None) == []
    assert qc.removed_marker_length == 1
    assert qc.removed_marker_completeness == 0
    # The per-marker breakdowns still show BOTH defects, each under its own
    # stage — that is what the breakdowns are for.
    assert qc.removed_marker_length_by_marker == {"5": 1}
    assert qc.removed_marker_completeness_by_marker == {"4": 1}
    # Every record of the isolate names the stage the COUNTERS used, not
    # whichever defect happened to come first in input order — otherwise
    # _qc_removed.tsv and the tally disagree about why the isolate went.
    assert _details(qc) == {
        "both_3": "marker_length_sibling:5:481aa_below_485",
        "both_4": "marker_completeness:4:no_met_start",
        "both_5": "marker_length:5:481aa_below_485",
    }
    # ... and the isolate stays visible as a completeness failure too.
    assert qc.display["marker_completeness_also_length"] == 1


def test_two_isolates_one_per_stage_sum_to_the_drop():
    seqs = _isolate("a", np=481) + _isolate("b", np_first="K")
    qc = QCReport()
    assert _run_marker_quality_qc(seqs, _seg_cfg(), qc, None) == []
    assert qc.removed_marker_length == 1
    assert qc.removed_marker_completeness == 1
    assert qc.display["length_complete_isolates"] == 1
    assert qc.display["completeness_complete_isolates"] == 1


def test_protection_keeps_the_isolate_and_records_the_stage():
    cfg = _seg_cfg(protect=["ref"])
    seqs = _isolate("ref", np_first="K")
    qc = QCReport()
    assert _run_marker_quality_qc(seqs, cfg, qc, None) == seqs
    assert qc.removed_marker_completeness == 0
    stages = {(p["id"], p["stage"], p["reason"]) for p in qc.protected}
    assert ("ref_5", "marker_quality", "marker_completeness:5:no_met_start") in stages


def test_only_the_record_assembly_keeps_is_judged():
    """A partial duplicate must not remove an isolate that also has a
    complete record for that segment (the 2026-10-28 review's finding)."""
    keep = _segrec("dup", "5", [_np("dup", 498, "a")], nt=1500)
    discard = _segrec("dup", "5", [_np("dup", 498, "b", first="K")], nt=1400)
    seqs = [s for s in _isolate("dup") if s.segment != "5"] + [keep, discard]
    qc = QCReport()
    assert _run_marker_quality_qc(seqs, _seg_cfg(), qc, None) == seqs
    assert qc.removed_marker_completeness == 0

    # Swap the NT lengths: now the partial record is the one assembly keeps.
    keep.sequence, discard.sequence = "A" * 1400, "A" * 1500
    qc = QCReport()
    assert _run_marker_quality_qc(seqs, _seg_cfg(), qc, None) == []
    assert qc.removed_marker_completeness == 1
    d = _details(qc)
    assert d["dup_5"] == "marker_completeness:5:no_met_start"


# ---------------------------------------------------------------------------
# The QC stage (non-segmented) — parity
# ---------------------------------------------------------------------------

def _cov(sid, s_len=1273, n_len=419, *, s_first="M", n_first="M"):
    return Sequence(
        id=sid, header=sid, sequence="ACGT" * 100,
        seq_type=SequenceType.NUCLEOTIDE, accession=sid,
        proteins=[
            _prot(f"{sid}_S", "spike", s_len, [_hit("CoV_S")], first=s_first),
            _prot(f"{sid}_N", "nucleocapsid", n_len, [_hit("CoV_N")], first=n_first),
        ],
    )


def _nonseg_cfg(*, met=True):
    return {
        "qc": {"marker_completeness": {"require_met_start": met}},
        "clustering": {"cluster_protein": [
            {"name": "S", "hmms": ["CoV_S"], "length_aa": [1100, 1400]},
            {"name": "N", "hmms": ["CoV_N"]},      # unbounded
        ]},
        "hmm": {"multidomain_overlap_tolerance": 30},
        "_hmm_runtime": {"active": True},
    }


def test_non_segmented_requires_every_marker_complete():
    seqs = [_cov("ok"), _cov("part_n", n_first="K"), _cov("part_s", s_first="K")]
    qc = QCReport()
    kept = _run_marker_quality_qc(seqs, _nonseg_cfg(), qc, None)
    assert [s.id for s in kept] == ["ok"]
    assert qc.removed_marker_completeness == 2            # sequences
    assert qc.removed_marker_completeness_by_marker == {"N": 1, "S": 1}
    assert _details(qc) == {
        "part_n": "marker_completeness:N:no_met_start",
        "part_s": "marker_completeness:S:no_met_start",
    }


def test_non_segmented_judges_only_already_mandatory_specs():
    """The global switch must not promote an advisory spec into a drop
    condition (designer's call 2026-10-29).

    Non-segmented clustering uses the FIRST matching spec, and the HMM gate
    already requires every ``hmms:`` spec, so judging those (and any spec
    that opted in with ``length_aa``) only strengthens a check the record had
    to pass. An unbounded alias-only spec was never gated by anything.
    """
    # Legacy bare-string form — alias-only, and it cannot even carry a range.
    cfg = _nonseg_cfg()
    cfg["clustering"]["cluster_protein"] = ["spike"]
    qc = QCReport()
    kept = _run_marker_quality_qc([_cov("part_s", s_first="K")], cfg, qc, None)
    assert [x.id for x in kept] == ["part_s"]
    assert qc.removed_marker_completeness == 0
    assert qc.removed_marker_completeness_by_marker == {}

    # Dict form, still alias-only and unbounded: same.
    cfg = _nonseg_cfg()
    cfg["clustering"]["cluster_protein"] = [{"name": "S", "aliases": ["spike"]}]
    qc = QCReport()
    assert len(_run_marker_quality_qc([_cov("part_s", s_first="K")], cfg, qc, None)) == 1
    assert qc.removed_marker_completeness == 0

    # Give that same spec a range and it has opted in — now it is judged.
    cfg = _nonseg_cfg()
    cfg["clustering"]["cluster_protein"] = [
        {"name": "S", "aliases": ["spike"], "length_aa": [1100, 1400]}
    ]
    qc = QCReport()
    assert _run_marker_quality_qc([_cov("part_s", s_first="K")], cfg, qc, None) == []
    assert qc.removed_marker_completeness_by_marker == {"S": 1}


def test_non_segmented_protection():
    cfg = _nonseg_cfg()
    cfg["overrides"] = {"ids": ["part_n"], "protect_qc": True, "protect_stages": "all"}
    qc = QCReport()
    kept = _run_marker_quality_qc([_cov("part_n", n_first="K")], cfg, qc, None)
    assert len(kept) == 1 and qc.removed_marker_completeness == 0
    assert qc.protected[0]["stage"] == "marker_quality"


def test_populate_routes_completeness_failures():
    specs = [{"name": "N", "hmms": ["CoV_N"]}]
    qc = QCReport()
    kept = populate_protein_sequences(
        [_cov("part_n", n_first="K")], specs, qc,
        hmm_active=True, require_met_start=True,
    )
    assert kept == []
    assert qc.removed_marker_completeness == 1
    assert qc.removed_marker_completeness_by_marker == {"N": 1}
    assert qc.removed_proteins == 0 and qc.removed_marker_length == 0
    assert _details(qc) == {"part_n": "marker_completeness:N:no_met_start"}


# ---------------------------------------------------------------------------
# The assembly backstop
# ---------------------------------------------------------------------------

def test_assembly_backstop_drops_and_counts_a_partial_marker():
    vcfg = _seg_cfg()["segmented"]["viruses"]["FluA"]
    qc = QCReport()
    out = build_concatenated_sequences(
        {"part": _isolate("part", np_first="K"), "good": _isolate("good")},
        segment_names=["3", "4", "5"], require_protein=True, report=qc,
        segment_markers=vcfg["segment_markers"], hmm_active=True,
        require_met_start=True,
    )
    assert [s.isolate_id for s in out] == ["good"]
    assert qc.removed_marker_completeness == 1
    assert qc.removed_marker_completeness_by_marker == {"5": 1}
    assert qc.removed_incomplete_isolates == 0 and qc.removed_hmm_failed == 0
    d = _details(qc)
    assert d["part_5"] == "marker_completeness:5:no_met_start"
    assert d["part_3"] == "marker_completeness_sibling:5:no_met_start"


def test_assembly_backstop_uses_the_real_protein_for_a_protected_isolate():
    from repseq.overrides import ProtectionPolicy

    cfg = _seg_cfg(protect=["ref"])
    vcfg = cfg["segmented"]["viruses"]["FluA"]
    qc = QCReport()
    out = build_concatenated_sequences(
        {"ref": _isolate("ref", np_first="K")},
        segment_names=["3", "4", "5"], require_protein=True, report=qc,
        segment_markers=vcfg["segment_markers"], hmm_active=True,
        policy=ProtectionPolicy.from_cfg(cfg), require_met_start=True,
    )
    assert len(out) == 1
    assert len(out[0].protein_sequence) == 716 + 566 + 498
    assert any(p["stage"] == "marker_quality" for p in qc.protected)


# ---------------------------------------------------------------------------
# Output-side pickers show exactly what clustering used
# ---------------------------------------------------------------------------

def _fasta_ids(path: Path) -> list[str]:
    return [
        ln[1:].split()[0] for ln in path.read_text().splitlines()
        if ln.startswith(">")
    ]


def test_pick_marker_cds_applies_the_requirement():
    proteins = [_np("x", 497, first="K")]
    assert pick_marker_cds(proteins, ["Flu_NP"], [], True) is not None
    assert pick_marker_cds(
        proteins, ["Flu_NP"], [], True, require_met_start=True,
    ) is None
    # exempt (a protected rep) falls back to the unrestricted pick.
    assert pick_marker_cds(
        proteins, ["Flu_NP"], [], True, require_met_start=True, exempt=True,
    ) is not None


def test_marker_bounds_context_carries_the_requirement():
    ctx = build_marker_bounds_context(_seg_cfg())
    assert ctx.require_met_start is True
    rep = _concat_rep("part", np_first="K")
    assert ctx.pick(
        rep, rep.concat_segments[2].proteins, "5_Flu_NP", ["Flu_NP"], [],
        "5", True, 30,
    ) is None
    assert build_marker_bounds_context(_seg_cfg(met=False)).require_met_start is False


def test_per_protein_fasta_omits_a_partial_cds(tmp_path):
    reps = [_concat_rep("part", np_first="K"), _concat_rep("ok1"), _concat_rep("ok2")]
    written = write_per_protein_fastas(
        RunResult(mode="global", representatives=reps), _seg_cfg(), None,
        tmp_path, "pfx",
    )
    by_name = {p.name: p for p in written}
    assert _fasta_ids(by_name["pfx_5_Flu_NP.fasta"]) == ["ok1_NP", "ok2_NP"]
    # The partial rep keeps its complete markers on the other segments.
    assert _fasta_ids(by_name["pfx_3_Flu_PA.fasta"]) == ["part_PA", "ok1_PA", "ok2_PA"]


def test_per_protein_fasta_shows_a_protected_reps_real_protein(tmp_path):
    reps = [_concat_rep("ref", np_first="K"), _concat_rep("ok1"), _concat_rep("ok2")]
    written = write_per_protein_fastas(
        RunResult(mode="global", representatives=reps),
        _seg_cfg(protect=["ref"]), None, tmp_path, "pfx",
    )
    by_name = {p.name: p for p in written}
    assert _fasta_ids(by_name["pfx_5_Flu_NP.fasta"]) == ["ref_NP", "ok1_NP", "ok2_NP"]


def test_per_protein_trees_use_the_same_picks(tmp_path):
    calls: dict[str, dict] = {}

    def _fake_build_tree(leaf_reps, bodies, **kw):
        calls[kw["file_prefix"]] = sorted(bodies)
        return []

    reps = [
        _concat_rep("part", np_first="K"),
        _concat_rep("ok1"), _concat_rep("ok2"), _concat_rep("ok3"),
    ]
    with patch.object(per_protein, "_build_tree", _fake_build_tree):
        per_protein.run_per_protein_phylogeny(reps, _seg_cfg(), tmp_path, "pfx")
    # The partial rep contributes to every marker tree but the NP one.
    assert calls["5_Flu_NP"] == ["CONCAT|ok1", "CONCAT|ok2", "CONCAT|ok3"]
    assert calls["3_Flu_PA"] == [
        "CONCAT|ok1", "CONCAT|ok2", "CONCAT|ok3", "CONCAT|part",
    ]


# ---------------------------------------------------------------------------
# extra_protein parity: absent, not dropped
# ---------------------------------------------------------------------------

_M2 = {"name": "M2", "aliases": ["matrix protein 2"]}


def _extra_cfg(**kw):
    return _seg_cfg(extra={"5": [_M2]}, **kw)


def test_a_partial_accessory_protein_is_absent_with_a_reason():
    ctx = build_extra_pick_context(_extra_cfg(), True)
    partial = _prot("x_M2", "matrix protein 2", 97, first="K")
    cds, set_aside = pick_extra_cds(
        [_np("x"), partial], "M2", [], ["matrix protein 2"], "5", ctx,
    )
    assert cds is None
    assert set_aside == [
        (partial, "5'-partial CDS (translation does not start with Met)")
    ]
    # Nothing was dropped: the isolate's cluster markers are untouched.
    qc = QCReport()
    seqs = _isolate("x")
    seqs[2].proteins = seqs[2].proteins + [partial]
    assert _run_marker_quality_qc(seqs, _extra_cfg(), qc, None) == seqs


def test_a_complete_accessory_protein_is_still_picked():
    ctx = build_extra_pick_context(_extra_cfg(), True)
    complete = _prot("x_M2", "matrix protein 2", 97)
    cds, set_aside = pick_extra_cds(
        [_np("x"), complete], "M2", [], ["matrix protein 2"], "5", ctx,
    )
    assert cds is complete and set_aside == []


def test_the_switch_off_leaves_accessory_picking_alone():
    ctx = build_extra_pick_context(_extra_cfg(met=False), True)
    partial = _prot("x_M2", "matrix protein 2", 97, first="K")
    cds, set_aside = pick_extra_cds(
        [_np("x"), partial], "M2", [], ["matrix protein 2"], "5", ctx,
    )
    assert cds is partial and set_aside == []


# ---------------------------------------------------------------------------
# Polyprotein peptides are exempt by construction
# ---------------------------------------------------------------------------

def test_polyprotein_slicing_is_untouched_by_the_requirement():
    """An internal cleavage product never starts with Met, so the slicer must
    not inherit the requirement. It cannot: ``slice_polyprotein`` has its own
    parent-CDS picker and no path to the flag. This pins that — a future
    refactor routing the slicer through the marker pickers would empty every
    per-peptide FASTA."""
    import inspect

    from repseq.polyprotein import collect_polyprotein_specs, slice_polyprotein

    assert "require_met_start" not in inspect.signature(slice_polyprotein).parameters

    cfg = {"clustering": {"polyprotein": [{
        "name": "P1",
        "peptides": [
            {"name": "VP4", "hmms": ["P_VP4"]},
            {"name": "VP2", "hmms": ["P_VP2"]},
        ],
    }]}}
    spec = collect_polyprotein_specs(cfg)[0]
    # A 5'-partial polyprotein parent: still elected, still sliced.
    parent = _prot(
        "x_pp", "polyprotein", 400,
        [_hit("P_VP4", ali_from=1, ali_to=100),
         _hit("P_VP2", ali_from=110, ali_to=300)],
        first="K",
    )
    cds, peptides = slice_polyprotein([parent], spec)
    assert cds is parent
    assert [p.status for p in peptides] == ["ok", "ok"]


# ---------------------------------------------------------------------------
# Config surface
# ---------------------------------------------------------------------------

def test_the_shipped_reference_config_documents_the_key():
    text = (
        Path(__file__).resolve().parents[1]
        / "repseq" / "data" / "default_config.yaml"
    ).read_text()
    assert "marker_completeness" in text
    assert "require_met_start" in text
    cfg = load_config(
        Path(__file__).resolve().parents[1]
        / "repseq" / "data" / "default_config.yaml"
    )
    assert validate_config(cfg) == []
    assert require_met_start_from_cfg(cfg) is False


def test_validate_rejects_a_non_boolean_switch():
    cfg = load_config(None)
    cfg["qc"]["marker_completeness"]["require_met_start"] = "yes"
    errors = validate_config(cfg)
    assert any("require_met_start" in e for e in errors)


def test_the_renamed_protect_token_hard_errors_with_a_pointer():
    cfg = load_config(None)
    cfg["overrides"]["protect_stages"] = ["hmm", "marker_length"]
    errors = validate_config(cfg)
    assert len(errors) == 1
    assert "marker_quality" in errors[0] and "0.77.0" in errors[0]

    cfg["overrides"]["protect_stages"] = ["hmm", "marker_quality"]
    assert validate_config(cfg) == []


# ---------------------------------------------------------------------------
# The Methods prose must not drift (output/summary.py)
# ---------------------------------------------------------------------------

def test_summary_reports_the_requirement_and_its_cost():
    qc = QCReport(total_input=10)
    qc.removed_marker_completeness = 2
    qc.removed_marker_completeness_by_marker = {"5": 2}
    text = _render_qc(qc, _seg_cfg())
    assert "5'-complete" in text
    assert "initial methionine" in text
    assert "**2** isolates were removed this way" in text
    assert "`5` (2)" in text


def test_summary_says_when_the_check_was_skipped():
    qc = QCReport(total_input=10)
    qc.display["marker_completeness_skipped"] = True
    text = _render_qc(qc, _seg_cfg())
    assert "not applied" in text and "--no-resolve" in text


def test_summary_is_silent_when_the_switch_is_off():
    text = _render_qc(QCReport(total_input=10), _seg_cfg(met=False))
    assert "5'-complete" not in text


def test_no_resolve_skips_the_check_and_flags_it():
    """The stage needs GenBank translations; without them it must skip AND
    say so, or the Methods prose would silently overclaim."""
    seqs = _isolate("part", np_first="K")
    for s in seqs:
        s.proteins = None                      # nothing fetched yet
    qc = QCReport()
    assert _run_marker_quality_qc(seqs, _seg_cfg(), qc, None) == seqs
    assert qc.display["marker_completeness_skipped"] is True
    assert qc.removed_marker_completeness == 0


# ---------------------------------------------------------------------------
# Follow-ups from the 2026-10-29 code review of this feature
# ---------------------------------------------------------------------------

def test_the_first_declared_spec_names_a_completeness_failure():
    """One attribution rule for both picker branches.

    The reported ``marker_name`` reaches ``_qc_removed.tsv`` and the
    per-marker breakdown, so it must not depend on whether the failing specs
    went through the HMM branch (which used to keep the LAST name) or the
    alias branch (which kept the first).
    """
    prots = [
        _prot("x_S", "spike", 1273, [_hit("CoV_S")], first="K"),
        _prot("x_N", "nucleocapsid", 419, [_hit("CoV_N")], first="K"),
    ]
    hmm_specs = [
        {"name": "S", "hmms": ["CoV_S"]}, {"name": "N", "hmms": ["CoV_N"]},
    ]
    marker, failure = select_marker_protein(
        prots, hmm_specs, hmm_active=True, require_met_start=True,
    )
    assert marker is None and failure.reason == "no_met_start"
    assert failure.marker_name == "S"

    alias_specs = [
        {"name": "S", "aliases": ["spike"]},
        {"name": "N", "aliases": ["nucleocapsid"]},
    ]
    marker, failure = select_marker_protein(
        prots, alias_specs, hmm_active=False, require_met_start=True,
    )
    assert marker is None and failure.reason == "no_met_start"
    assert failure.marker_name == "S"


def test_an_accessory_hmm_claimant_claims_the_cds_it_will_pick():
    """A claim must land on the CDS the claimant's own picker takes.

    Claiming exists so two specs cannot write the same CDS into two files.
    If an HMM ``extra_protein`` claimant claimed without the Met filter its
    claim would land on a 5'-partial CDS it then withholds, leaving its real
    pick unclaimed and free for an alias spec to take as well.
    """
    partial = _prot("x_M2a", "matrix protein 2", 97, [_hit("Flu_M2")], first="K")
    complete = _prot("x_M2b", "M2 ion channel", 94, [_hit("Flu_M2")])
    cfg = _seg_cfg(extra={"5": [
        {"name": "M2", "hmms": ["Flu_M2"]},
        {"name": "M2b", "aliases": ["ion channel"]},
    ]})
    ctx = build_extra_pick_context(cfg, True)
    prots = [_np("x"), partial, complete]

    pick, _ = pick_extra_cds(prots, "M2", ["Flu_M2"], [], "5", ctx)
    assert pick is complete                      # the partial is withheld
    # ... so the alias spec must find it already claimed.
    claimed = ctx.claimed(prots, "5")
    assert claimed.get(id(complete)) == "M2"
    pick_b, set_aside = pick_extra_cds(
        prots, "M2b", [], ["ion channel"], "5", ctx,
    )
    assert pick_b is None
    assert [r for _p, r in set_aside] == ["already selected by HMM spec M2"]


def test_a_skipped_check_is_not_enforced_by_the_output_pickers():
    """The Methods prose and the artifacts must agree.

    When the stage skips (no GenBank translations under ``--no-resolve``) the
    summary says the requirement was NOT applied, so no picker may go on
    omitting 5'-partial CDSes from the FASTAs, trees and coverage tables.
    """
    cfg = _seg_cfg()
    assert require_met_start_active(cfg) is True
    assert build_marker_bounds_context(cfg).require_met_start is True

    seqs = _isolate("part", np_first="K")
    for s in seqs:
        s.proteins = None
    qc = QCReport()
    assert _run_marker_quality_qc(seqs, cfg, qc, None) == seqs

    assert cfg[MARKER_QUALITY_RUNTIME]["require_met_start_skipped"] is True
    assert require_met_start_active(cfg) is False
    assert build_marker_bounds_context(cfg).require_met_start is False
    assert build_extra_pick_context(cfg, True).require_met_start is False
    # The CONFIGURED value is unchanged, which is what lets the prose say
    # "configured but not applied".
    assert require_met_start_from_cfg(cfg) is True


def test_the_prose_never_reads_the_check_as_a_no_op():
    """With length taking precedence, a run whose truncated records all break
    a range too would otherwise report "no isolates were removed this way"
    while _qc_removed.tsv carries completeness rows."""
    qc = QCReport(total_input=10)
    # The realistic shape: every 5'-partial marker also broke a range, so
    # precedence put all three isolates in the length counter.
    qc.removed_marker_length = 3
    qc.removed_marker_length_by_marker = {"5": 3}
    qc.removed_marker_completeness = 0
    qc.removed_marker_completeness_by_marker = {"4": 3}
    qc.display["marker_completeness_also_length"] = 3
    text = _render_qc(qc, _seg_cfg())
    assert "removed by this check alone" in text
    assert "**3** that carried a 5'-partial marker" in text

    # And when it did remove some, the co-counted ones are still named.
    qc2 = QCReport(total_input=10)
    qc2.removed_marker_completeness = 2
    qc2.removed_marker_completeness_by_marker = {"4": 2}
    qc2.display["marker_completeness_also_length"] = 3
    text2 = _render_qc(qc2, _seg_cfg())
    assert "**2** isolates were removed this way" in text2
    assert "A further **3** isolates" in text2


def test_the_console_tally_keeps_co_counted_isolates_visible():
    qc = QCReport(total_input=9)
    qc.display.update({
        "segmented": True, "n_segments": 3,
        "length_isolates_pre_assembly": 4,
        "marker_completeness_also_length": 1,
    })
    qc.removed_marker_length = 4
    assert "incl. 1 also 5'-partial" in qc.summary()


def test_the_concatenation_note_accounts_for_every_cause():
    """The parenthetical annotates a total, so it has to add up to it."""
    qc = QCReport(total_input=9)
    qc.display.update({
        "segmented": True, "n_segments": 3,
        "concat_marker_isolates": 10, "concat_hmm_isolates": 3,
        "concat_length_isolates": 2, "concat_completeness_isolates": 5,
    })
    text = qc.summary()
    assert "3 via HMM" in text and "2 via length" in text
    assert "5 via 5'-partial marker" in text
