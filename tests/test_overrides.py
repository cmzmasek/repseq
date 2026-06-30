"""Tests for the QC force-keep override (overrides.protect_qc)."""

from __future__ import annotations

import copy

import pytest

from repseq.config import DEFAULTS, validate_config
from repseq.models import (
    QCReport,
    Sequence,
    SequenceSource,
    SequenceType,
)
from repseq.overrides import (
    QC_PROTECT_STAGES,
    ProtectionPolicy,
    _norm_id,
    resolve_ids,
    resolve_stages,
)
from repseq.segmented.completeness import _normalise_isolate_id
from repseq.models import Cluster, RunResult
from repseq.overrides import apply_force_select
from repseq.output.report import write_force_selected_tsv, write_overrides_tsv
from repseq.qc.pipeline import ambiguous_filter, run_qc
from repseq.segmented.taxonomy_consistency import (
    filter_taxonomy_consistent_isolates,
)


def _seq(sid, sequence="ACGT", *, accession=None, isolate_id=None):
    return Sequence(
        id=sid,
        sequence=sequence,
        seq_type=SequenceType.NUCLEOTIDE,
        source=SequenceSource.NCBI,
        accession=accession,
        isolate_id=isolate_id,
        description="",
        header=sid,
    )


# ---------------------------------------------------------------------------
# id resolution + matching
# ---------------------------------------------------------------------------

def test_resolve_ids_normalises_and_version_augments():
    ids = resolve_ids({"ids": ["NC_045512.2", "  MixedCase  "]})
    # full + version-stripped forms are stored, lower-cased + trimmed
    assert "nc_045512.2" in ids
    assert "nc_045512" in ids
    assert "mixedcase" in ids


def test_policy_matches_version_insensitive_both_directions():
    # configured WITH version → matches a seq WITHOUT version
    p = ProtectionPolicy(resolve_ids({"ids": ["NC_045512.2"]}),
                         resolve_stages("all"), enabled=True)
    assert p.protects(_seq("s", accession="NC_045512"), "ambiguous")
    # configured WITHOUT version → matches a seq WITH version
    p2 = ProtectionPolicy(resolve_ids({"ids": ["NC_045512"]}),
                          resolve_stages("all"), enabled=True)
    assert p2.protects(_seq("s", accession="NC_045512.3"), "ambiguous")


def test_policy_matches_id_and_isolate_id():
    p = ProtectionPolicy(resolve_ids({"ids": ["ISO-1"]}),
                         resolve_stages("all"), enabled=True)
    assert p.protects(_seq("x", isolate_id="iso-1"), "hmm")
    assert not p.protects(_seq("x", isolate_id="iso-2"), "hmm")


def test_norm_id_canonical_outputs():
    # Pin the canonical normalisation to CONCRETE outputs (not a self-
    # referential comparison, which would be tautological now that
    # _normalise_isolate_id delegates to _norm_id): lower-case, every
    # whitespace run -> "_", every pipe -> "_".
    assert _norm_id("A/Louisiana/12/2024 OR_IR") == "a/louisiana/12/2024_or_ir"
    assert _norm_id("  Mixed   WS\tName  ") == "mixed_ws_name"
    assert _norm_id("a|b|c") == "a_b_c"
    assert _norm_id("NC_045512.2") == "nc_045512.2"  # internal "_" preserved
    assert _norm_id("") is None
    assert _norm_id("   ") is None


def test_normalise_isolate_id_concrete_outputs():
    # _normalise_isolate_id delegates to _norm_id but must keep producing
    # exactly these strings (grouping dict key + CONCAT seq.id). Pinning the
    # concrete output catches a future inlined re-implementation that drifts
    # (the regression class the delegation removed) — the previous test
    # compared the two functions, which can never fail while one calls the
    # other.
    assert _normalise_isolate_id("A/Foo Bar/1") == "a/foo_bar/1"
    assert _normalise_isolate_id("a|b|c") == "a_b_c"
    assert _normalise_isolate_id("") == ""  # empty -> "" (not None)


def test_override_matches_pipe_bearing_accession():
    # Regression (review finding 1/2): an UNKNOWN-source FASTA header can
    # yield a seq.accession carrying a leading/trailing pipe (io/fasta.py's
    # first-token fallback). _norm_id maps the pipe to "_"; version-
    # insensitive matching must still bind the bare / versioned override id —
    # the trailing "_" must not defeat _strip_version.
    p = ProtectionPolicy(resolve_ids({"ids": ["AB123456"]}),
                         resolve_stages("all"), enabled=True)
    assert p.protects(_seq("s", accession="AB123456.1|"), "ambiguous")
    assert p.protects(_seq("s", accession="|AB123456.1"), "ambiguous")
    assert p.protects(_seq("s", accession="AB123456.1"), "ambiguous")  # control
    assert not p.protects(_seq("s", accession="XY999999.1|"), "ambiguous")


def test_exclude_conflict_detected_across_pipe_and_version():
    # Regression (review finding 2): the exclude-vs-keep hard error must fire
    # even when the two spellings differ only by a trailing pipe + version,
    # i.e. their version-stripped match keys coincide.
    cfg = copy.deepcopy(DEFAULTS)
    cfg["overrides"]["protect_qc"] = True
    cfg["overrides"]["ids"] = ["AB123456.1"]
    cfg["overrides"]["exclude"]["enabled"] = True
    cfg["overrides"]["exclude"]["ids"] = ["AB123456.1|"]
    assert any("cannot be both" in e for e in validate_config(cfg))


def test_force_select_binds_segmented_isolate_by_space_form_name():
    # The user lists the natural strain name (with spaces); it must match
    # BOTH a raw per-segment isolate_id and the underscored CONCAT id that
    # _normalise_isolate_id produces.
    p = ProtectionPolicy(resolve_ids({"ids": ["A/Louisiana/12/2024 OR_IR"]}),
                         resolve_stages("all"), enabled=True, force_select=True)
    concat_id = _normalise_isolate_id("A/Louisiana/12/2024 OR_IR")
    assert concat_id == "a/louisiana/12/2024_or_ir"
    # CONCAT pool sequence (isolate_id already normalised to underscores)
    assert p.pins(_seq("CONCAT|" + concat_id, isolate_id=concat_id))
    # raw per-segment sequence (isolate_id is the original space-form strain)
    assert p.protects(
        _seq("seg8", isolate_id="A/Louisiana/12/2024 OR_IR"), "hmm")
    # a different isolate is not bound
    assert not p.pins(_seq("CONCAT|other", isolate_id="a/other/1"))


def test_policy_inactive_without_ids():
    p = ProtectionPolicy(frozenset(), resolve_stages("all"), enabled=True)
    assert p.enabled is False
    assert not p.protects(_seq("x", accession="A"), "ambiguous")


def test_policy_respects_stage_scoping():
    p = ProtectionPolicy(resolve_ids({"ids": ["A"]}),
                         resolve_stages(["ambiguous"]), enabled=True)
    s = _seq("x", accession="A")
    assert p.protects(s, "ambiguous")
    assert not p.protects(s, "annotation")  # not in protect_stages


def test_from_cfg_reads_runtime_cache():
    cfg = {"_overrides_runtime": {
        "ids": resolve_ids({"ids": ["A"]}),
        "stages": resolve_stages("all"),
        "protect_qc": True,
    }}
    p = ProtectionPolicy.from_cfg(cfg)
    assert p.enabled and p.protects(_seq("x", accession="A"), "hmm")


def test_from_cfg_resolves_inline_when_no_cache():
    cfg = {"overrides": {"ids": ["A"], "protect_qc": True}}
    p = ProtectionPolicy.from_cfg(cfg)
    assert p.enabled and p.protects(_seq("x", accession="A"), "duplicates")


# ---------------------------------------------------------------------------
# stage integration: ambiguous filter + full run_qc
# ---------------------------------------------------------------------------

def test_ambiguous_filter_protects_listed_sequence():
    report = QCReport()
    policy = ProtectionPolicy(resolve_ids({"ids": ["B"]}),
                              resolve_stages("all"), enabled=True)
    clean = _seq("A", "ACGT" * 25, accession="A")
    dirty = _seq("B", "NNNN" * 25, accession="B")
    kept = ambiguous_filter([clean, dirty], 0.05, report, policy=policy)
    assert sorted(s.id for s in kept) == ["A", "B"]
    assert report.removed_ambiguous == 0
    assert report.protected == [
        {"id": "B", "stage": "ambiguous", "reason": "ambiguous_fraction:1.000>0.05"}
    ]


def test_ambiguous_filter_drops_when_not_protected():
    report = QCReport()
    policy = ProtectionPolicy(resolve_ids({"ids": ["Z"]}),
                              resolve_stages("all"), enabled=True)
    dirty = _seq("B", "NNNN" * 25, accession="B")
    kept = ambiguous_filter([dirty], 0.05, report, policy=policy)
    assert kept == []
    assert report.removed_ambiguous == 1
    assert report.protected == []


def test_run_qc_protects_only_named_stage():
    # Protect against ambiguous but NOT annotation: a dirty+annotated seq
    # is rescued from ambiguous yet still dropped by annotation.
    clean = _seq("A", "ACGT" * 25, accession="A")
    dirty = _seq("B", "NNNN" * 25, accession="B")
    dirty.description = "synthetic construct"
    cfg = {
        "qc": {
            "ambiguous_threshold": 0.05,
            "annotation_filter": {"enabled": True, "keywords": ["synthetic"]},
            "remove_duplicates": True,
        },
        "segmented": {"enabled": False},
        "overrides": {"ids": ["B"], "protect_qc": True,
                      "protect_stages": ["ambiguous"]},
    }
    kept, report = run_qc([clean, dirty], cfg)
    assert sorted(s.id for s in kept) == ["A"]  # B dropped by annotation
    assert report.removed_ambiguous == 0       # but NOT by ambiguous
    assert report.removed_annotation == 1
    assert any(p["stage"] == "ambiguous" for p in report.protected)


def test_run_qc_no_overrides_is_unchanged():
    dirty = _seq("B", "NNNN" * 25, accession="B")
    cfg = {
        "qc": {"ambiguous_threshold": 0.05,
               "annotation_filter": {"enabled": False},
               "remove_duplicates": True},
        "segmented": {"enabled": False},
    }
    kept, report = run_qc([dirty], cfg)
    assert kept == []
    assert report.removed_ambiguous == 1
    assert report.protected == []


# ---------------------------------------------------------------------------
# isolate-level stage: taxonomy consistency
# ---------------------------------------------------------------------------

class _Tax:
    def __init__(self, species):
        self._species = species

    def get_rank(self, rank):
        return self._species if rank == "species" else None


def test_taxonomy_consistency_protects_whole_isolate():
    a = _seq("a", isolate_id="iso1", accession="a")
    b = _seq("b", isolate_id="iso1", accession="b")
    a.taxonomy = _Tax("Virus A")
    b.taxonomy = _Tax("Virus B")  # mismatch → isolate would drop
    policy = ProtectionPolicy(resolve_ids({"ids": ["iso1"]}),
                              resolve_stages("all"), enabled=True)
    protected_out: list = []
    kept, removed = filter_taxonomy_consistent_isolates(
        [a, b], rank="species", policy=policy, protected_out=protected_out,
    )
    assert sorted(s.id for s in kept) == ["a", "b"]
    assert removed == []
    assert {p[0] for p in protected_out} == {"a", "b"}
    assert all(p[1] == "taxonomy_consistency" for p in protected_out)


def test_taxonomy_consistency_drops_when_not_protected():
    a = _seq("a", isolate_id="iso1", accession="a")
    b = _seq("b", isolate_id="iso1", accession="b")
    a.taxonomy = _Tax("Virus A")
    b.taxonomy = _Tax("Virus B")
    kept, removed = filter_taxonomy_consistent_isolates([a, b], rank="species")
    assert kept == []
    assert len(removed) == 2


# ---------------------------------------------------------------------------
# config validation
# ---------------------------------------------------------------------------

def test_validate_rejects_unknown_stage_token():
    cfg = copy.deepcopy(DEFAULTS)
    cfg["overrides"] = {"ids": ["x"], "protect_qc": True,
                        "protect_stages": ["bogus"]}
    errs = [e for e in validate_config(cfg) if "protect_stages" in e]
    assert errs and "bogus" in errs[0]


def test_validate_rejects_bad_types():
    cfg = copy.deepcopy(DEFAULTS)
    cfg["overrides"] = {"ids": "notalist", "protect_qc": "yes"}
    errs = validate_config(cfg)
    assert any("overrides.ids must be a list" in e for e in errs)
    assert any("overrides.protect_qc must be a boolean" in e for e in errs)


def test_validate_accepts_all_and_subset():
    for stages in ("all", list(QC_PROTECT_STAGES), ["hmm", "ambiguous"]):
        cfg = copy.deepcopy(DEFAULTS)
        cfg["overrides"] = {"ids": ["x"], "protect_qc": True,
                            "protect_stages": stages}
        assert [e for e in validate_config(cfg) if "overrides" in e] == []


# ---------------------------------------------------------------------------
# report writer
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# force_select (apply_force_select)
# ---------------------------------------------------------------------------

def _S(sid, n=100, *, accession=None, isolate_id=None, refseq=False):
    s = _seq(sid, "A" * n, accession=accession or sid, isolate_id=isolate_id)
    s.is_refseq = refseq
    return s


def _fs_cfg(ids, on=True):
    return {
        "overrides": {"ids": ids, "force_select": on},
        "representative": {"priority": ["refseq", "reviewed_uniprot", "longest"]},
    }


def test_force_select_elects_pinned_member_over_representative():
    r1 = _S("R1", refseq=True)        # would normally win on refseq
    p1 = _S("P1", 120)                # pinned, not refseq
    m2 = _S("M2")
    c = Cluster("c1", representative=r1, members=[p1, m2])
    result = RunResult(mode="m", representatives=[r1], clusters=[c])
    apply_force_select(result, [r1, p1, m2], _fs_cfg(["P1"]))
    assert c.representative is p1                 # pin won election
    assert r1 in c.members                        # old rep demoted
    assert p1 in result.representatives and r1 not in result.representatives
    assert [e["action"] for e in result.force_selected] == ["elected_representative"]


def test_force_select_splits_colliding_pins_into_singletons():
    r1 = _S("R1", refseq=True)
    p1 = _S("P1", 120)
    p2 = _S("P2", 90)
    c = Cluster("c1", representative=r1, members=[p1, p2])
    result = RunResult(mode="m", representatives=[r1], clusters=[c])
    apply_force_select(result, [r1, p1, p2], _fs_cfg(["P1", "P2"]))
    # Longest pin (P1) wins; P2 split into its own singleton cluster.
    assert c.representative is p1
    assert any(cl.representative is p2 and cl.members == [] for cl in result.clusters)
    assert {"P1", "P2"} <= {s.id for s in result.representatives}
    actions = sorted(e["action"] for e in result.force_selected)
    assert actions == ["elected_representative", "split_singleton"]


def test_force_select_adds_diversity_deselected_orphan():
    # global -n style: selected reps are singleton clusters; D1 is in the
    # pool but in no cluster.
    r1 = _S("R1")
    d1 = _S("D1")
    c = Cluster("c1", representative=r1)
    result = RunResult(mode="global:count", representatives=[r1], clusters=[c])
    apply_force_select(result, [r1, d1], _fs_cfg(["D1"]))
    assert d1 in result.representatives
    assert any(cl.representative is d1 for cl in result.clusters)
    assert result.force_selected[0]["action"] == "added_representative"


def test_force_select_already_representative_is_noop_action():
    r1 = _S("R1")
    c = Cluster("c1", representative=r1)
    result = RunResult(mode="m", representatives=[r1], clusters=[c])
    apply_force_select(result, [r1], _fs_cfg(["R1"]))
    assert result.representatives == [r1]
    assert result.force_selected == [
        {"id": "R1", "action": "already_representative", "detail": ""}
    ]


def test_force_select_reports_unavailable_pin():
    r1 = _S("R1")
    c = Cluster("c1", representative=r1)
    result = RunResult(mode="m", representatives=[r1], clusters=[c])
    apply_force_select(result, [r1], _fs_cfg(["GHOST"]))
    assert [e["action"] for e in result.force_selected] == ["unavailable"]
    assert result.force_selected[0]["id"] == "GHOST"


def test_force_select_noop_when_disabled():
    r1 = _S("R1")
    p1 = _S("P1")
    c = Cluster("c1", representative=r1, members=[p1])
    result = RunResult(mode="m", representatives=[r1], clusters=[c])
    apply_force_select(result, [r1, p1], _fs_cfg(["P1"], on=False))
    assert result.representatives == [r1]
    assert result.force_selected == []


def test_force_select_matches_segmented_isolate_id():
    # CONCAT rep carries isolate_id; pin by isolate_id.
    r1 = _S("CONCAT|isoA", isolate_id="isoA")
    p1 = _S("CONCAT|isoB", isolate_id="isoB")
    c = Cluster("c1", representative=r1, members=[p1])
    result = RunResult(mode="m", representatives=[r1], clusters=[c])
    apply_force_select(result, [r1, p1], _fs_cfg(["isoB"]))
    assert c.representative is p1
    assert result.force_selected[0]["action"] == "elected_representative"


def test_write_force_selected_tsv(tmp_path):
    result = RunResult(mode="m")
    result.force_selected = [
        {"id": "P1", "action": "elected_representative", "detail": "cluster=c1"},
        {"id": "GHOST", "action": "unavailable", "detail": "no surviving match"},
    ]
    path = tmp_path / "x_force_selected.tsv"
    assert write_force_selected_tsv(result, path) is True
    lines = path.read_text().splitlines()
    assert lines[0] == "id\taction\tdetail"
    assert lines[1] == "P1\telected_representative\tcluster=c1"


def test_write_force_selected_tsv_skips_when_empty(tmp_path):
    result = RunResult(mode="m")
    path = tmp_path / "x_force_selected.tsv"
    assert write_force_selected_tsv(result, path) is False
    assert not path.exists()


def test_validate_rejects_bad_force_select_type():
    cfg = copy.deepcopy(DEFAULTS)
    cfg["overrides"] = {"ids": ["x"], "force_select": "yes"}
    errs = validate_config(cfg)
    assert any("overrides.force_select must be a boolean" in e for e in errs)


def test_pin_ids_flag_enables_force_select(tmp_path):
    from repseq.cli import _load_and_validate

    pins = tmp_path / "pins.txt"
    pins.write_text("NC_1.1\nNC_2\n")
    out = tmp_path / "out"
    cfg = _load_and_validate(
        config_path=None, output_dir=str(out), prefix="t",
        threads=None, seed=None, pin_ids=str(pins),
    )
    assert cfg["overrides"]["force_select"] is True
    rt = cfg["_overrides_runtime"]
    assert rt["force_select"] is True
    assert "nc_1.1" in rt["ids"] and "nc_2" in rt["ids"]
    assert "NC_1.1" in rt["raw_ids"]


def test_write_overrides_tsv_skips_when_empty(tmp_path):
    report = QCReport()
    path = tmp_path / "x_overrides.tsv"
    assert write_overrides_tsv(report, path) is False
    assert not path.exists()


# ---------------------------------------------------------------------------
# --protect-ids CLI flag (via _load_and_validate)
# ---------------------------------------------------------------------------

def test_protect_ids_flag_merges_file_and_enables_protect_qc(tmp_path):
    from repseq.cli import _load_and_validate

    vip = tmp_path / "vip.txt"
    vip.write_text("# comment\nNC_045512.2\n\nMN908947\n")
    out = tmp_path / "out"
    cfg = _load_and_validate(
        config_path=None, output_dir=str(out), prefix="t",
        threads=None, seed=None, protect_ids=str(vip),
    )
    # config default has protect_qc off + empty ids; the flag flips it on
    # and merges the file ids.
    assert cfg["overrides"]["protect_qc"] is True
    assert "NC_045512.2" in cfg["overrides"]["ids"]
    assert "MN908947" in cfg["overrides"]["ids"]
    rt = cfg["_overrides_runtime"]
    assert rt["protect_qc"] is True
    # version-augmented + normalised in the resolved runtime set
    assert "nc_045512.2" in rt["ids"] and "nc_045512" in rt["ids"]


def test_protect_ids_flag_unions_with_config_ids(tmp_path):
    from repseq.cli import _load_and_validate

    yaml_cfg = tmp_path / "c.yaml"
    yaml_cfg.write_text("overrides:\n  ids: [AB000001]\n")
    vip = tmp_path / "vip.txt"
    vip.write_text("CD000002\n")
    out = tmp_path / "out"
    cfg = _load_and_validate(
        config_path=str(yaml_cfg), output_dir=str(out), prefix="t",
        threads=None, seed=None, protect_ids=str(vip),
    )
    assert set(cfg["overrides"]["ids"]) == {"AB000001", "CD000002"}


def test_protect_ids_missing_file_exits(tmp_path):
    from repseq.cli import _load_and_validate

    out = tmp_path / "out"
    with pytest.raises(SystemExit):
        _load_and_validate(
            config_path=None, output_dir=str(out), prefix="t",
            threads=None, seed=None, protect_ids=str(tmp_path / "nope.txt"),
        )


def test_write_overrides_tsv_writes_rows(tmp_path):
    report = QCReport()
    report.add_protected("B", "ambiguous", "ambiguous_fraction:1.000>0.05")
    report.add_protected("B", "hmm", "hmm_failed:S:Foo")
    path = tmp_path / "x_overrides.tsv"
    assert write_overrides_tsv(report, path) is True
    lines = path.read_text().splitlines()
    assert lines[0] == "id\tstage\twould_be_reason"
    assert lines[1] == "B\tambiguous\tambiguous_fraction:1.000>0.05"
    assert lines[2] == "B\thmm\thmm_failed:S:Foo"


# ---------------------------------------------------------------------------
# exclude (input blocklist: apply_exclusions / validation / CLI flag)
# ---------------------------------------------------------------------------

from repseq.overrides import apply_exclusions, resolve_raw_ids
from repseq.output.report import write_excluded_tsv


def _excl_cfg(ids):
    """A cfg with the exclude runtime populated, as _resolve_overrides would."""
    ov = {"exclude": {"enabled": True, "ids": ids}}
    return {
        "_overrides_runtime": {
            "exclude_ids": resolve_ids(ov["exclude"]),
            "exclude_raw_ids": resolve_raw_ids(ov["exclude"]),
        }
    }


def test_exclude_drops_by_accession_version_insensitive():
    keep = _S("A", accession="NC_111.1")
    drop = _S("B", accession="NC_222.3")
    # List the unversioned form; must still match NC_222.3.
    kept, audit = apply_exclusions([keep, drop], _excl_cfg(["NC_222"]))
    assert kept == [keep]
    assert [(e["id"], e["action"]) for e in audit] == [("NC_222.3", "excluded")]
    assert audit[0]["detail"] == "matched on accession"


def test_exclude_matches_by_id_when_no_accession():
    s = _seq("WEIRD_ID", "ACGT")  # accession None
    kept, audit = apply_exclusions([s], _excl_cfg(["weird_id"]))
    assert kept == []
    assert audit[0]["action"] == "excluded"
    assert audit[0]["detail"] == "matched on id"


def test_exclude_does_not_match_isolate_id():
    # Deliberate: isolate_id isn't populated this early, so the blocklist is
    # header-id-only. An isolate_id-only match must NOT drop the sequence.
    s = _seq("X", "ACGT", accession="ACC1", isolate_id="isoZ")
    kept, audit = apply_exclusions([s], _excl_cfg(["isoZ"]))
    assert kept == [s]
    assert [e["action"] for e in audit] == ["unavailable"]


def test_exclude_reports_unavailable_id():
    s = _S("A", accession="ACC1")
    kept, audit = apply_exclusions([s], _excl_cfg(["GHOST"]))
    assert kept == [s]
    assert audit == [{
        "id": "GHOST", "action": "unavailable",
        "detail": "no input sequence matched (typo, or already absent)",
    }]


def test_exclude_noop_when_empty():
    s = _S("A", accession="ACC1")
    cfg = {"_overrides_runtime": {"exclude_ids": frozenset()}}
    kept, audit = apply_exclusions([s], cfg)
    assert kept == [s] and audit == []


def test_write_excluded_tsv_writes_and_skips(tmp_path):
    path = tmp_path / "x_excluded.tsv"
    assert write_excluded_tsv({}, path) is False
    assert not path.exists()
    cfg = {"_excluded_runtime": {"audit": [
        {"id": "NC_1.1", "action": "excluded", "detail": "matched on accession"},
        {"id": "GHOST", "action": "unavailable", "detail": "no input sequence matched"},
    ]}}
    assert write_excluded_tsv(cfg, path) is True
    lines = path.read_text().splitlines()
    assert lines[0] == "id\taction\tdetail"
    assert lines[1] == "NC_1.1\texcluded\tmatched on accession"
    assert lines[2] == "GHOST\tunavailable\tno input sequence matched"


# ---- validation: exclude block + contradiction guard ----

def test_validate_rejects_exclude_bad_types():
    cfg = copy.deepcopy(DEFAULTS)
    cfg["overrides"]["exclude"] = {"enabled": "yes", "ids": "NC_1", "ids_file": 3}
    errs = validate_config(cfg)
    assert any("overrides.exclude.enabled must be a boolean" in e for e in errs)
    assert any("overrides.exclude.ids must be a list" in e for e in errs)
    assert any("overrides.exclude.ids_file must be a path string" in e for e in errs)


def test_validate_rejects_exclude_protect_conflict():
    cfg = copy.deepcopy(DEFAULTS)
    cfg["overrides"]["ids"] = ["NC_045512.2"]
    cfg["overrides"]["protect_qc"] = True
    cfg["overrides"]["exclude"] = {"enabled": True, "ids": ["NC_045512"]}
    errs = validate_config(cfg)
    assert any("cannot be both" in e and "NC_045512" in e for e in errs)


def test_validate_no_conflict_when_keep_pin_off():
    # Same id in both lists is harmless when neither keep nor pin is active.
    cfg = copy.deepcopy(DEFAULTS)
    cfg["overrides"]["ids"] = ["NC_045512.2"]  # inert: flags off
    cfg["overrides"]["exclude"] = {"enabled": True, "ids": ["NC_045512.2"]}
    errs = validate_config(cfg)
    assert not any("cannot be both" in e for e in errs)


def test_validate_default_config_has_exclude_block():
    cfg = copy.deepcopy(DEFAULTS)
    assert cfg["overrides"]["exclude"] == {
        "enabled": False, "ids": [], "ids_file": None
    }
    assert validate_config(cfg) == []


# ---- --exclude-ids CLI flag (via _load_and_validate) ----

def test_exclude_ids_flag_merges_file_and_enables(tmp_path):
    from repseq.cli import _load_and_validate

    bad = tmp_path / "bad.txt"
    bad.write_text("# chimeric\nNC_9.1\n\nKJ642623\n")
    out = tmp_path / "out"
    cfg = _load_and_validate(
        config_path=None, output_dir=str(out), prefix="t",
        threads=None, seed=None, exclude_ids=str(bad),
    )
    assert cfg["overrides"]["exclude"]["enabled"] is True
    rt = cfg["_overrides_runtime"]
    assert "nc_9.1" in rt["exclude_ids"] and "kj642623" in rt["exclude_ids"]
    assert "NC_9.1" in rt["exclude_raw_ids"]


def test_resolve_overrides_skips_exclude_when_disabled():
    from repseq.cli import _resolve_overrides

    cfg = {"overrides": {"exclude": {"enabled": False, "ids": ["NC_1"]}}}
    _resolve_overrides(cfg)
    # Disabled => no ids resolved, so the blocklist is a no-op.
    assert cfg["_overrides_runtime"]["exclude_ids"] == frozenset()
