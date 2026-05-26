"""Unit tests for repseq.polyprotein (specs + slicer).

Test fixtures must mirror the producer module's actual schema — every
synthetic HMM hit dict carries ``target`` (HMM profile name, NOT
``hmm_name``), ``ali_from``, ``ali_to``, ``dom_evalue``, and
``ali_span``. See ``feedback-test-fixtures-match-production`` in memory:
in v0.31.0 the conservation labeller silently misbehaved for two
releases because tests used ``hmm_name`` while the producer
(``hmm/hmmscan.py:_parse_domtblout``) writes ``target``.
"""

from __future__ import annotations

import pytest

from repseq.polyprotein import (
    PeptideSpec,
    PolyproteinSpec,
    collect_polyprotein_specs,
    slice_polyprotein,
)
from repseq.polyprotein.slicer import (
    compute_cuts,
    identify_parent_cds,
)


# ---------------------------------------------------------------------------
# Fixtures: producer-schema-shaped HMM hits + a synthetic polyprotein CDS.
# ---------------------------------------------------------------------------

def _hit(target: str, ali_from: int, ali_to: int, evalue: float = 1e-30) -> dict:
    """Build a synthetic HMM hit matching the real ``_parse_domtblout`` schema."""
    return {
        "target": target,
        "ali_from": ali_from,
        "ali_to": ali_to,
        "dom_evalue": evalue,
        "evalue": evalue,
        "ali_span": ali_to - ali_from + 1,
    }


def _polyprotein_cds(seq: str, hits: list[dict], pid: str = "YP_TEST.1") -> dict:
    return {
        "protein_id": pid,
        "sequence": seq,
        "length": len(seq),
        "hmm_hits": hits,
        "parent_accession": "NC_TEST",
    }


def _picornavirus_polyprotein() -> dict:
    """A 197-aa synthetic polyprotein with 4 peptides separated by LQ motifs."""
    seq = (
        "A" * 50 + "LQ" +
        "B" * 50 + "LQ" +
        "C" * 50 + "LQ" +
        "D" * 41
    )
    # seq length is 197; the final hit's ali_to must fit.
    assert len(seq) == 197
    hits = [
        _hit("P_VP4", 1,   50),
        _hit("P_VP2", 53,  103),
        _hit("P_VP3", 106, 156),
        _hit("P_VP1", 159, 197),
    ]
    return _polyprotein_cds(seq, hits)


def _picornavirus_spec(strategy: str) -> PolyproteinSpec:
    return PolyproteinSpec(
        name="P1",
        peptides=[
            PeptideSpec(name="VP4", hmm="P_VP4"),
            PeptideSpec(name="VP2", hmm="P_VP2", cleavage_motif="LQ"),
            PeptideSpec(name="VP3", hmm="P_VP3", cleavage_motif="LQ"),
            PeptideSpec(name="VP1", hmm="P_VP1", cleavage_motif="LQ"),
        ],
        cut_strategy=strategy,
        motif_window_aa=15,
        min_peptides_hit=2,
    )


# ---------------------------------------------------------------------------
# identify_parent_cds
# ---------------------------------------------------------------------------

class TestIdentifyParentCDS:
    def test_picks_cds_with_most_distinct_peptide_hits(self):
        parent = _picornavirus_polyprotein()
        decoy = _polyprotein_cds("X" * 100, [
            _hit("P_VP4", 1, 50),  # only one peptide HMM — doesn't clear threshold
        ], pid="YP_DECOY.1")
        spec = _picornavirus_spec("bisect")

        chosen = identify_parent_cds([decoy, parent], spec)
        assert chosen is not None
        assert chosen["protein_id"] == "YP_TEST.1"

    def test_returns_none_when_below_min_peptides_hit(self):
        parent = _polyprotein_cds("X" * 100, [_hit("P_VP4", 1, 50)])
        spec = _picornavirus_spec("bisect")  # min_peptides_hit=2
        assert identify_parent_cds([parent], spec) is None

    def test_returns_none_on_empty_list(self):
        spec = _picornavirus_spec("bisect")
        assert identify_parent_cds([], spec) is None


# ---------------------------------------------------------------------------
# Cut strategies
# ---------------------------------------------------------------------------

class TestCutStrategies:
    def test_boundary_uses_hit_coords_verbatim(self):
        parent = _picornavirus_polyprotein()
        spec = _picornavirus_spec("boundary")
        _, sliced = slice_polyprotein([parent], spec)

        # Each peptide should span its hit's ali_from..ali_to exactly.
        assert sliced[0].range_aa_from == 1 and sliced[0].range_aa_to == 50
        assert sliced[1].range_aa_from == 53 and sliced[1].range_aa_to == 103
        assert sliced[2].range_aa_from == 106 and sliced[2].range_aa_to == 156
        assert sliced[3].range_aa_from == 159 and sliced[3].range_aa_to == 197
        # Inter-hit residues (the LQ motifs) are dropped at boundaries.
        assert all(s.cut_method_actual == "boundary" for s in sliced)

    def test_bisect_no_residues_dropped(self):
        parent = _picornavirus_polyprotein()
        spec = _picornavirus_spec("bisect")
        _, sliced = slice_polyprotein([parent], spec)

        # First peptide starts at 1; last ends at protein C-term.
        assert sliced[0].range_aa_from == 1
        assert sliced[-1].range_aa_to == len(parent["sequence"])
        # Coverage: every residue accounted for exactly once.
        total = sum(s.length_aa for s in sliced)
        assert total == len(parent["sequence"])
        # Adjacent peptides are flush.
        for i in range(len(sliced) - 1):
            assert sliced[i].range_aa_to + 1 == sliced[i + 1].range_aa_from

    def test_motif_snaps_to_cleavage_site(self):
        parent = _picornavirus_polyprotein()
        spec = _picornavirus_spec("motif")
        _, sliced = slice_polyprotein([parent], spec)

        # VP2/VP3/VP1 each declare cleavage_motif="LQ". Their START
        # positions should snap to just after an LQ. The seq has LQ at
        # 51-52, 103-104, 155-156, so peptide starts are 53, 105, 157.
        assert sliced[1].range_aa_from == 53   # VP2
        assert sliced[2].range_aa_from == 105  # VP3
        assert sliced[3].range_aa_from == 157  # VP1
        # First peptide's start is the protein N-term, not motif-driven.
        assert sliced[0].cut_method_actual == "n-term"
        for i in (1, 2, 3):
            assert sliced[i].cut_method_actual == "motif:LQ"

    def test_motif_falls_back_to_bisect_when_motif_absent(self):
        # Same hits, but the parent has no LQ motifs anywhere.
        seq = "A" * 50 + "XX" + "B" * 50 + "XX" + "C" * 50 + "XX" + "D" * 41
        parent = _polyprotein_cds(seq, [
            _hit("P_VP4", 1, 50),
            _hit("P_VP2", 53, 103),
            _hit("P_VP3", 106, 156),
            _hit("P_VP1", 159, 199),
        ])
        spec = _picornavirus_spec("motif")
        _, sliced = slice_polyprotein([parent], spec)

        # No LQ in window → bisect fallback for each inter-peptide cut.
        for i in (1, 2, 3):
            assert sliced[i].cut_method_actual == "bisect"


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------

class TestSliceEdgeCases:
    def test_no_parent_cds_yields_audit_rows_only(self):
        seq = "X" * 100
        parent = _polyprotein_cds(seq, [_hit("P_VP4", 1, 50)])
        spec = _picornavirus_spec("bisect")  # needs ≥2 distinct hits
        chosen, sliced = slice_polyprotein([parent], spec)
        assert chosen is None
        assert len(sliced) == 4
        assert all(s.status == "no_parent_cds" for s in sliced)
        assert all(s.sequence == "" for s in sliced)

    def test_missing_peptide_skipped_neighbours_extend(self):
        # Drop the VP3 hit; VP2 and VP1 should still slice cleanly.
        seq = "A" * 50 + "LQ" + "B" * 50 + "LQ" + "C" * 50 + "LQ" + "D" * 41
        parent = _polyprotein_cds(seq, [
            _hit("P_VP4", 1, 50),
            _hit("P_VP2", 53, 103),
            # P_VP3 missing
            _hit("P_VP1", 159, 197),
        ])
        spec = _picornavirus_spec("bisect")
        _, sliced = slice_polyprotein([parent], spec)

        names = [s.peptide_name for s in sliced]
        assert names == ["VP4", "VP2", "VP3", "VP1"]
        assert sliced[2].status == "missing"
        assert sliced[2].sequence == ""
        # The other three should be ok.
        for i in (0, 1, 3):
            assert sliced[i].status == "ok"
            assert sliced[i].sequence
        # VP2 and VP1 are now adjacent (VP3 hole closed).
        assert sliced[1].range_aa_to + 1 == sliced[3].range_aa_from

    def test_out_of_order_hits_fail_spec(self):
        # P_VP2 hits AFTER P_VP3 — order violated.
        seq = "Z" * 200
        parent = _polyprotein_cds(seq, [
            _hit("P_VP4", 1, 50),
            _hit("P_VP3", 60, 100),  # earlier
            _hit("P_VP2", 110, 150),  # later but configured before VP3
            _hit("P_VP1", 160, 200),
        ])
        spec = _picornavirus_spec("bisect")
        _, sliced = slice_polyprotein([parent], spec)

        # Every peptide should land in the out_of_order failure state.
        assert all(s.status == "out_of_order" for s in sliced)
        assert all(s.sequence == "" for s in sliced)

    def test_overlap_flagged_but_sequences_still_emitted(self):
        # P_VP2 hit overlaps P_VP3 hit by 5 residues — should set overlap status,
        # still emit sequences (bisect falls back to midpoint of centres).
        seq = "Z" * 200
        parent = _polyprotein_cds(seq, [
            _hit("P_VP4", 1, 50),
            _hit("P_VP2", 60, 110),
            _hit("P_VP3", 100, 150),  # overlaps the previous by 11 aa
            _hit("P_VP1", 160, 199),
        ])
        spec = _picornavirus_spec("bisect")
        _, sliced = slice_polyprotein([parent], spec)
        # At least one peptide flagged overlap; sequences still produced.
        assert any(s.status == "overlap" for s in sliced)
        assert all(s.sequence for s in sliced if s.status in ("ok", "overlap"))


# ---------------------------------------------------------------------------
# Config → specs
# ---------------------------------------------------------------------------

class TestCollectPolyproteinSpecs:
    def test_non_segmented_simple(self):
        cfg = {
            "clustering": {
                "polyprotein": [
                    {
                        "name": "P1",
                        "peptides": [
                            {"name": "VP4", "hmm": "P_VP4"},
                            {"name": "VP2", "hmm": "P_VP2", "cleavage_motif": "LQ"},
                        ],
                        "min_peptides_hit": 2,
                    }
                ]
            }
        }
        specs = collect_polyprotein_specs(cfg)
        assert len(specs) == 1
        assert specs[0].name == "P1"
        assert specs[0].segment is None
        assert len(specs[0].peptides) == 2
        # cleavage_motif on one peptide → default strategy should be "motif"
        assert specs[0].cut_strategy == "motif"

    def test_segmented_per_segment_dict(self):
        cfg = {
            "segmented": {
                "enabled": True,
                "virus": "cov",
                "viruses": {
                    "cov": {
                        "segments": ["genome"],
                        "polyprotein": {
                            "genome": [
                                {
                                    "name": "ORF1ab",
                                    "peptides": [
                                        {"name": "NSP3", "hmm": "CoV_NSP3"},
                                        {"name": "NSP5", "hmm": "CoV_NSP5"},
                                    ],
                                }
                            ]
                        },
                    }
                },
            }
        }
        specs = collect_polyprotein_specs(cfg)
        assert len(specs) == 1
        assert specs[0].name == "ORF1ab"
        assert specs[0].segment == "genome"
        # No cleavage_motif declared → default strategy is "bisect"
        assert specs[0].cut_strategy == "bisect"

    def test_empty_returns_empty(self):
        assert collect_polyprotein_specs({}) == []
        assert collect_polyprotein_specs({"clustering": {"polyprotein": []}}) == []


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

class TestConfigValidation:
    def _validate(self, polyprotein_list):
        from repseq.config import _validate_polyprotein_list
        return _validate_polyprotein_list(polyprotein_list, "clustering.polyprotein")

    def test_valid_minimal_spec_passes(self):
        errs = self._validate([
            {
                "name": "P1",
                "peptides": [
                    {"name": "VP4", "hmm": "P_VP4"},
                    {"name": "VP2", "hmm": "P_VP2"},
                ],
            }
        ])
        assert errs == []

    def test_rejects_single_peptide(self):
        errs = self._validate([
            {"name": "P1", "peptides": [{"name": "VP4", "hmm": "P_VP4"}]}
        ])
        assert any("at least 2 peptide" in e for e in errs)

    def test_rejects_missing_hmm(self):
        errs = self._validate([
            {
                "name": "P1",
                "peptides": [
                    {"name": "VP4"},  # no hmm
                    {"name": "VP2", "hmm": "P_VP2"},
                ],
            }
        ])
        assert any("hmm must be a non-empty" in e for e in errs)

    def test_rejects_duplicate_spec_names(self):
        errs = self._validate([
            {"name": "P1", "peptides": [
                {"name": "VP4", "hmm": "P_VP4"},
                {"name": "VP2", "hmm": "P_VP2"},
            ]},
            {"name": "P1", "peptides": [
                {"name": "NSP3", "hmm": "CoV_NSP3"},
                {"name": "NSP5", "hmm": "CoV_NSP5"},
            ]},
        ])
        assert any("duplicate spec name" in e for e in errs)

    def test_rejects_duplicate_peptide_names(self):
        errs = self._validate([
            {"name": "P1", "peptides": [
                {"name": "VP4", "hmm": "P_VP4"},
                {"name": "VP4", "hmm": "P_VP4_alt"},
            ]}
        ])
        assert any("duplicate peptide name" in e for e in errs)

    def test_rejects_bad_cut_strategy(self):
        errs = self._validate([
            {"name": "P1", "cut_strategy": "magic", "peptides": [
                {"name": "VP4", "hmm": "P_VP4"},
                {"name": "VP2", "hmm": "P_VP2"},
            ]}
        ])
        assert any("cut_strategy 'magic' is not supported" in e for e in errs)
