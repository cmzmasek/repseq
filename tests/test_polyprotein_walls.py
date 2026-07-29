"""Polyprotein peptide-coverage "wall of zeros" alarm.

A clade that clusters fine (passes the marker gate) but whose polyprotein
peptides are almost entirely uncovered by its best-matching spec — the
Orthoflavivirus-only-slicing-against-Hepacivirus/Pegivirus/Pestivirus case.
Covers the shared-core-detector, the TSV-reading flags path, and the config
knobs.
"""
from __future__ import annotations

import repseq.output.report as report
from repseq.output.report import compute_polyprotein_walls, find_polyprotein_coverage_walls
from repseq.output.flags import _polyprotein_wall_flags, collect_flags, write_flags_report
from repseq.config import load_config, validate_config

FLAVI = ["C", "prM", "E", "NS1", "NS2A", "NS2B", "NS3", "NS4A", "NS4B", "NS5"]
HEP = ["Core", "E1", "E2", "NS2", "NS3", "NS4A", "NS4B", "NS5A", "NS5B"]
PES = ["Npro", "C", "E2", "NS2", "NS3", "NS5B"]


# ---------------------------------------------------------------------------
# compute_polyprotein_walls — the shared core
# ---------------------------------------------------------------------------

def test_old_single_spec_fires_wall_for_non_flavi():
    """The pre-fix state: one Orthoflavivirus spec, RdRP_3 in NS5, so the
    non-flavi genera realise only NS5 (1/10) → a wall."""
    cov = {
        ("flavivirus", "Orthoflavivirus"): {p: 148 for p in FLAVI},
        ("flavivirus", "Orthohepacivirus"): {**{p: 0 for p in FLAVI}, "NS5": 149},
        ("flavivirus", "Orthopestivirus"): {**{p: 0 for p in FLAVI}, "NS5": 145},
    }
    totals = {"Orthoflavivirus": 148, "Orthohepacivirus": 149, "Orthopestivirus": 145}
    walls = compute_polyprotein_walls(
        cov, totals, {"flavivirus": FLAVI}, rank="genus",
    )
    flagged = {w["taxon"] for w in walls}
    assert flagged == {"Orthohepacivirus", "Orthopestivirus"}
    for w in walls:
        assert w["kind"] == "mistuned_spec"
        assert w["spec"] == "flavivirus"
        assert (w["covered_peptides"], w["total_peptides"]) == (1, 10)
    # Orthoflavivirus (full coverage) is NOT flagged.


def test_new_per_clade_specs_no_false_wall_despite_shared_marker():
    """The fix: three per-clade specs. RdRP_3 is shared as NS5B in BOTH the
    hepacivirus and pestivirus specs, so a pestivirus is nominally "claimed"
    by the hepacivirus spec via that one peptide — the argmax-home rule must
    ignore it and raise NO wall."""
    cov = {
        ("flavivirus", "Orthoflavivirus"): {p: 148 for p in FLAVI},
        ("flavivirus", "Orthohepacivirus"): {p: 0 for p in FLAVI},
        ("flavivirus", "Orthopestivirus"): {p: 0 for p in FLAVI},
        ("hepacivirus", "Orthohepacivirus"): {p: 149 for p in HEP},
        # pestivirus hits the hepacivirus spec ONLY at the shared NS5B (RdRP_3):
        ("hepacivirus", "Orthopestivirus"): {**{p: 0 for p in HEP}, "NS5B": 145},
        # pegivirus legitimately lacks Core/E1 (2/9 zero = 22% < 60%):
        ("hepacivirus", "Pegivirus"): {**{p: 144 for p in HEP}, "Core": 0, "E1": 0},
        ("pestivirus", "Orthopestivirus"): {p: 145 for p in PES},
        # hepacivirus hits the pestivirus spec ONLY at the shared NS5B:
        ("pestivirus", "Orthohepacivirus"): {**{p: 0 for p in PES}, "NS5B": 149},
        ("pestivirus", "Pegivirus"): {**{p: 0 for p in PES}, "NS5B": 144},
    }
    totals = {"Orthoflavivirus": 148, "Orthohepacivirus": 149,
              "Orthopestivirus": 145, "Pegivirus": 144}
    specs = {"flavivirus": FLAVI, "hepacivirus": HEP, "pestivirus": PES}
    walls = compute_polyprotein_walls(cov, totals, specs, rank="genus")
    assert walls == []


def test_unsliced_taxon_when_no_spec_slices_anything():
    """A clade with reps but zero peptides sliced by any spec → the more
    severe unsliced_taxon wall."""
    cov = {
        ("flavivirus", "Orthoflavivirus"): {p: 50 for p in FLAVI},
        ("flavivirus", "Weirdvirus"): {p: 0 for p in FLAVI},
    }
    totals = {"Orthoflavivirus": 50, "Weirdvirus": 8}
    walls = compute_polyprotein_walls(cov, totals, {"flavivirus": FLAVI}, rank="genus")
    assert len(walls) == 1
    assert walls[0]["kind"] == "unsliced_taxon"
    assert walls[0]["taxon"] == "Weirdvirus"
    assert walls[0]["spec"] is None
    assert walls[0]["n_reps"] == 8


def test_min_reps_guards_tiny_taxa():
    """A 2-rep clade with an all-zero wall is below the default min_reps=3
    and is NOT flagged (noise guard)."""
    cov = {
        ("flavivirus", "Orthoflavivirus"): {p: 50 for p in FLAVI},
        ("flavivirus", "Tinygenus"): {**{p: 0 for p in FLAVI}, "NS5": 2},
    }
    totals = {"Orthoflavivirus": 50, "Tinygenus": 2}
    walls = compute_polyprotein_walls(cov, totals, {"flavivirus": FLAVI},
                                      rank="genus", min_reps=3)
    assert [w["taxon"] for w in walls] == []
    # With min_reps=1 the same clade IS flagged.
    walls = compute_polyprotein_walls(cov, totals, {"flavivirus": FLAVI},
                                      rank="genus", min_reps=1)
    assert [w["taxon"] for w in walls] == ["Tinygenus"]


def test_wall_fraction_boundary():
    """A clade at exactly wall_fraction zeros fires; just below does not."""
    # 6/10 zero = 0.6 → fires at threshold 0.6
    cov = {("flavivirus", "G"): {**{p: 0 for p in FLAVI[:6]},
                                 **{p: 20 for p in FLAVI[6:]}}}
    totals = {"G": 20}
    assert compute_polyprotein_walls(cov, totals, {"flavivirus": FLAVI},
                                     rank="genus", wall_fraction=0.6)
    # 5/10 zero = 0.5 < 0.6 → no wall
    cov = {("flavivirus", "G"): {**{p: 0 for p in FLAVI[:5]},
                                 **{p: 20 for p in FLAVI[5:]}}}
    assert compute_polyprotein_walls(cov, totals, {"flavivirus": FLAVI},
                                     rank="genus", wall_fraction=0.6) == []


def test_sort_order_unsliced_before_mistuned():
    cov = {
        ("flavivirus", "Orthoflavivirus"): {p: 50 for p in FLAVI},
        ("flavivirus", "Mistuned"): {**{p: 0 for p in FLAVI}, "NS5": 40},
        ("flavivirus", "Unsliced"): {p: 0 for p in FLAVI},
    }
    totals = {"Orthoflavivirus": 50, "Mistuned": 40, "Unsliced": 10}
    walls = compute_polyprotein_walls(cov, totals, {"flavivirus": FLAVI}, rank="genus")
    assert [w["kind"] for w in walls] == ["unsliced_taxon", "mistuned_spec"]


# ---------------------------------------------------------------------------
# find_polyprotein_coverage_walls — the in-memory wrapper (drives the console
# WARNING). Guards + coverage-dict construction + delegation. The real slicing
# helper is monkeypatched so we exercise the wrapper's plumbing without brittle
# HMM-hit fixtures.
# ---------------------------------------------------------------------------

def _cfg_with_flavi_spec():
    cfg = load_config(None)
    cfg["clustering"]["polyprotein"] = [{
        "name": "flavivirus", "cut_strategy": "bisect", "min_peptides_hit": 2,
        "peptides": [{"name": p, "hmm": f"Flavi_{p}"} for p in FLAVI],
    }]
    cfg["_hmm_runtime"] = {"active": True}  # make _hmm_tier_ran short-circuit True
    return cfg


def _fake_cov(totals_map, ns5_only=None):
    """Return a stand-in for _polyprotein_coverage_data_per_taxon: every taxon
    fully covered on all peptides, except taxa in ``ns5_only`` which get only
    the last peptide (NS5) — the mistuned-clade shape."""
    def _inner(seqs, spec, rank, overlap_tolerance):
        n = len(spec.peptides)
        totals = dict(totals_map)
        lengths = {}
        for taxon, count in totals_map.items():
            if ns5_only and taxon in ns5_only:
                lengths[taxon] = [[] for _ in range(n)]
                lengths[taxon][-1] = [800] * count   # NS5 only
            else:
                lengths[taxon] = [[500] * count for _ in range(n)]
        return totals, lengths
    return _inner


def test_wrapper_builds_coverage_and_delegates(monkeypatch):
    monkeypatch.setattr(
        report, "_polyprotein_coverage_data_per_taxon",
        _fake_cov({"Orthoflavivirus": 50, "Orthohepacivirus": 40},
                  ns5_only={"Orthohepacivirus"}),
    )
    walls = find_polyprotein_coverage_walls([], _cfg_with_flavi_spec())
    assert [(w["taxon"], w["kind"], w["covered_peptides"], w["total_peptides"])
            for w in walls] == [("Orthohepacivirus", "mistuned_spec", 1, 10)]


def test_wrapper_returns_empty_when_disabled(monkeypatch):
    monkeypatch.setattr(
        report, "_polyprotein_coverage_data_per_taxon",
        _fake_cov({"Orthohepacivirus": 40}, ns5_only={"Orthohepacivirus"}),
    )
    cfg = _cfg_with_flavi_spec()
    cfg["output"]["polyprotein_report"]["wall_warning"]["enabled"] = False
    assert find_polyprotein_coverage_walls([], cfg) == []


def test_wrapper_returns_empty_without_polyprotein_specs():
    cfg = load_config(None)
    cfg["_hmm_runtime"] = {"active": True}
    # no clustering.polyprotein declared
    assert find_polyprotein_coverage_walls([], cfg) == []


def test_wrapper_returns_empty_when_hmm_tier_did_not_run():
    cfg = _cfg_with_flavi_spec()
    cfg.pop("_hmm_runtime")  # and reps=[] carry no hmm_hits → tier didn't run
    assert find_polyprotein_coverage_walls([], cfg) == []


def test_wrapper_honours_cfg_thresholds(monkeypatch):
    # 40-rep clade with only NS5 → 9/10 zero. min_reps=50 excludes it.
    monkeypatch.setattr(
        report, "_polyprotein_coverage_data_per_taxon",
        _fake_cov({"Orthohepacivirus": 40}, ns5_only={"Orthohepacivirus"}),
    )
    cfg = _cfg_with_flavi_spec()
    cfg["output"]["polyprotein_report"]["wall_warning"]["min_reps"] = 50
    assert find_polyprotein_coverage_walls([], cfg) == []


# ---------------------------------------------------------------------------
# console path — _final_summary emits the WARNING from qc_report.polyprotein_walls
# ---------------------------------------------------------------------------

def test_console_warning_mistuned_and_unsliced(capsys):
    from repseq.cli import _final_summary
    from repseq.models import QCReport, RunResult, Sequence, SequenceType
    qc = QCReport(total_input=5, passed=5)
    qc.polyprotein_walls = [
        {"kind": "unsliced_taxon", "rank": "genus", "taxon": "Weirdvirus",
         "spec": None, "n_reps": 8, "covered_peptides": 0,
         "zero_peptides": 10, "total_peptides": 10},
        {"kind": "mistuned_spec", "rank": "genus", "taxon": "Orthohepacivirus",
         "spec": "flavivirus", "n_reps": 149, "covered_peptides": 1,
         "zero_peptides": 9, "total_peptides": 10},
    ]
    rep = Sequence(id="r", header="r", sequence="A" * 9, seq_type=SequenceType.NUCLEOTIDE)
    result = RunResult(mode="global", representatives=[rep], clusters=[])
    _final_summary(result, qc, {"segmented": {"enabled": False}})
    err = capsys.readouterr().err
    assert "peptide-coverage WALL" in err
    assert "Orthohepacivirus (genus) [flavivirus: 1/10 peptides]" in err
    assert "Weirdvirus (genus) [no spec slices it]" in err


def test_console_no_warning_when_no_walls(capsys):
    from repseq.cli import _final_summary
    from repseq.models import QCReport, RunResult, Sequence, SequenceType
    qc = QCReport(total_input=5, passed=5)  # polyprotein_walls defaults to []
    rep = Sequence(id="r", header="r", sequence="A" * 9, seq_type=SequenceType.NUCLEOTIDE)
    result = RunResult(mode="global", representatives=[rep], clusters=[])
    _final_summary(result, qc, {"segmented": {"enabled": False}})
    assert "peptide-coverage WALL" not in capsys.readouterr().err


# ---------------------------------------------------------------------------
# flags path — reads the polyprotein tidy TSV
# ---------------------------------------------------------------------------

def _write_poly_tsv(path, coverage, totals, spec_peptides):
    lines = ["\t".join(("report", "rank", "pool", "taxon", "taxon_count",
                        "spec", "metric", "value"))]
    for (spec, taxon), peps in coverage.items():
        for pep in spec_peptides[spec]:
            cnt = peps.get(pep, 0)
            lines.append("\t".join((
                "polyprotein", "genus", "reps", taxon, str(totals[taxon]),
                f"{spec}:{pep}", "coverage_count", str(cnt),
            )))
    path.write_text("\n".join(lines) + "\n")


def test_flags_reads_tsv_and_flags_wall(tmp_path):
    _write_poly_tsv(
        tmp_path / "x_polyprotein_taxonomic_report.tsv",
        {
            ("flavivirus", "Orthoflavivirus"): {p: 148 for p in FLAVI},
            ("flavivirus", "Orthohepacivirus"): {**{p: 0 for p in FLAVI}, "NS5": 149},
        },
        {"Orthoflavivirus": 148, "Orthohepacivirus": 149},
        {"flavivirus": FLAVI},
    )
    flags = _polyprotein_wall_flags(tmp_path, "x", cfg=None)
    assert len(flags) == 1
    assert flags[0].category == "polyprotein_wall"
    assert "Orthohepacivirus" in flags[0].message
    assert flags[0].severity == "warn"


def test_flags_disabled_toggle(tmp_path):
    _write_poly_tsv(
        tmp_path / "x_polyprotein_taxonomic_report.tsv",
        {("flavivirus", "Orthohepacivirus"): {**{p: 0 for p in FLAVI}, "NS5": 149}},
        {"Orthohepacivirus": 149},
        {"flavivirus": FLAVI},
    )
    cfg = {"output": {"polyprotein_report": {"wall_warning": {"enabled": False}}}}
    assert _polyprotein_wall_flags(tmp_path, "x", cfg=cfg) == []


def test_write_flags_report_fires_on_wall_only(tmp_path):
    """The flags file is written even with no conflict tables when a wall is
    detected (parity with the QC-elimination alarm)."""
    _write_poly_tsv(
        tmp_path / "x_polyprotein_taxonomic_report.tsv",
        {("flavivirus", "Orthohepacivirus"): {**{p: 0 for p in FLAVI}, "NS5": 149}},
        {"Orthohepacivirus": 149},
        {"flavivirus": FLAVI},
    )
    path = write_flags_report(tmp_path, "x", cfg=None)
    assert path is not None
    text = path.read_text()
    assert "Polyprotein peptide-coverage walls" in text
    assert "Orthohepacivirus" in text


# ---------------------------------------------------------------------------
# config validation
# ---------------------------------------------------------------------------

def test_default_config_has_wall_warning():
    cfg = load_config(None)
    ww = cfg["output"]["polyprotein_report"]["wall_warning"]
    assert ww == {"enabled": True, "rank": "genus", "wall_fraction": 0.6, "min_reps": 3}
    assert validate_config(cfg) == []


def test_validate_rejects_bad_wall_warning():
    cfg = load_config(None)
    cfg["output"]["polyprotein_report"]["wall_warning"]["wall_fraction"] = 1.5
    assert any("wall_fraction" in e for e in validate_config(cfg))
    cfg = load_config(None)
    cfg["output"]["polyprotein_report"]["wall_warning"]["min_reps"] = 0
    assert any("min_reps" in e for e in validate_config(cfg))
    cfg = load_config(None)
    cfg["output"]["polyprotein_report"]["wall_warning"]["rank"] = "nonsense"
    assert any("rank" in e for e in validate_config(cfg))
