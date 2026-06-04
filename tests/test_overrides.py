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
    resolve_ids,
    resolve_stages,
)
from repseq.output.report import write_overrides_tsv
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
