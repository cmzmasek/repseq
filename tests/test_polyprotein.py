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
    _satisfying_span_for_token,
    _best_satisfying_alternative,
    slice_polyprotein,
)


# ---------------------------------------------------------------------------
# Fixtures: producer-schema-shaped HMM hits + a synthetic polyprotein CDS.
# ---------------------------------------------------------------------------

def _hit(
    target: str,
    ali_from: int,
    ali_to: int,
    evalue: float = 1e-30,
    passing: bool = True,
) -> dict:
    """Build a synthetic HMM hit matching the real ``_parse_domtblout`` schema.

    ``passing`` mirrors the per-hit flag set by ``cli.py:_run_hmm_qc`` in
    production. Defaults to True so existing fixtures keep their meaning;
    set False in regression tests that exercise the spurious-low-confidence-
    hit case the slicer now filters out.
    """
    return {
        "target": target,
        "ali_from": ali_from,
        "ali_to": ali_to,
        "dom_evalue": evalue,
        "evalue": evalue,
        "ali_span": ali_to - ali_from + 1,
        "passing": passing,
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
            PeptideSpec(name="VP4", hmms=["P_VP4"]),
            PeptideSpec(name="VP2", hmms=["P_VP2"], cleavage_motif="LQ"),
            PeptideSpec(name="VP3", hmms=["P_VP3"], cleavage_motif="LQ"),
            PeptideSpec(name="VP1", hmms=["P_VP1"], cleavage_motif="LQ"),
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

    def test_missing_peptide_leaves_unassigned_gap(self):
        # Drop the VP3 hit; VP2 and VP1 should still slice cleanly but
        # NOT extend across VP3's territory (v0.37.0+ behaviour). VP2's
        # C-side stays at its HMM hit's ali_to and VP1's N-side stays at
        # its HMM hit's ali_from — the gap between them is unassigned.
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
        # VP2's C-side stays at its HMM hit (103). VP1's N-side stays at
        # its HMM hit (159). The 104..158 gap is unassigned to either.
        assert sliced[1].range_aa_to == 103
        assert sliced[3].range_aa_from == 159
        # VP1's N-side method is hit-boundary because its predecessor in
        # spec (VP3) is missing.
        assert sliced[3].cut_method_actual == "hit-boundary"

    def test_missing_leading_peptide_does_not_let_successor_extend_to_n_term(self):
        # VP4 missing; VP2 (now the first surviving) must NOT extend to
        # AA 1 — it keeps its HMM-hit start (53). Pre-v0.37.0 the first
        # surviving peptide always started at AA 1, which on a real
        # polyprotein would absorb the missing peptide's territory.
        seq = "A" * 50 + "LQ" + "B" * 50 + "LQ" + "C" * 50 + "LQ" + "D" * 41
        parent = _polyprotein_cds(seq, [
            # P_VP4 missing
            _hit("P_VP2", 53, 103),
            _hit("P_VP3", 106, 156),
            _hit("P_VP1", 159, 197),
        ])
        spec = _picornavirus_spec("bisect")
        _, sliced = slice_polyprotein([parent], spec)

        assert sliced[0].status == "missing"  # VP4
        assert sliced[1].status == "ok"       # VP2
        assert sliced[1].range_aa_from == 53  # stays at HMM hit, NOT 1
        assert sliced[1].cut_method_actual == "hit-boundary"

    def test_missing_trailing_peptide_does_not_let_predecessor_extend_to_c_term(self):
        # VP1 missing; VP3 (now the last surviving) must NOT extend to
        # the C-terminus — it keeps its HMM-hit end (156).
        seq = "A" * 50 + "LQ" + "B" * 50 + "LQ" + "C" * 50 + "LQ" + "D" * 41
        parent = _polyprotein_cds(seq, [
            _hit("P_VP4", 1, 50),
            _hit("P_VP2", 53, 103),
            _hit("P_VP3", 106, 156),
            # P_VP1 missing
        ])
        spec = _picornavirus_spec("bisect")
        _, sliced = slice_polyprotein([parent], spec)

        assert sliced[3].status == "missing"  # VP1
        assert sliced[2].status == "ok"       # VP3
        assert sliced[2].range_aa_to == 156   # stays at HMM hit, NOT 197

    def test_missing_neighbour_cascade_sars_cov_2_shape(self):
        # Models the SARS-CoV-2 NSP1/NSP2-missing/NSP3 cascade that
        # motivated v0.37.0. With NSP2 (the middle peptide) missing,
        # the previous behaviour split the NSP1->NSP3 gap evenly,
        # inflating both flanking peptides by hundreds of residues.
        # New behaviour: NSP1 stops at its own HMM hit and NSP3 starts
        # at its own HMM hit — the missing peptide's territory is
        # unassigned to either flanker.
        seq = "Z" * 3000
        parent = _polyprotein_cds(seq, [
            _hit("P_NSP1", 1, 165),
            # P_NSP2 missing
            _hit("P_NSP3", 850, 2700),
        ])
        spec = PolyproteinSpec(
            name="ORF1a",
            peptides=[
                PeptideSpec(name="NSP1", hmms=["P_NSP1"]),
                PeptideSpec(name="NSP2", hmms=["P_NSP2"]),
                PeptideSpec(name="NSP3", hmms=["P_NSP3"]),
            ],
            cut_strategy="bisect",
            min_peptides_hit=2,
        )
        _, sliced = slice_polyprotein([parent], spec)

        assert sliced[0].status == "ok"
        assert sliced[1].status == "missing"
        assert sliced[2].status == "ok"
        # NSP1 stays at 1..165 — it does NOT inflate to 1..507 (the
        # midpoint between 165 and 850 that the old logic would have
        # picked).
        assert sliced[0].range_aa_from == 1
        assert sliced[0].range_aa_to == 165
        # NSP3 stays at 850..3000 — its N-side does NOT shift to 508.
        # C-side still extends to n_aa because NSP3 is the last declared.
        assert sliced[2].range_aa_from == 850
        assert sliced[2].range_aa_to == 3000
        assert sliced[2].cut_method_actual == "hit-boundary"

    def test_out_of_order_hits_rejected_at_parent_identification(self):
        # P_VP2 hits AFTER P_VP3 — declared order is violated on the CDS.
        # v0.36.1+: identify_parent_cds catches the order violation
        # before this CDS is elected as parent, so the spec fails with
        # no_parent_cds (more honest than the old out_of_order, which
        # implied "we elected a parent but then something went wrong").
        seq = "Z" * 200
        parent = _polyprotein_cds(seq, [
            _hit("P_VP4", 1, 50),
            _hit("P_VP3", 60, 100),  # earlier
            _hit("P_VP2", 110, 150),  # later but configured before VP3
            _hit("P_VP1", 160, 200),
        ])
        spec = _picornavirus_spec("bisect")
        _, sliced = slice_polyprotein([parent], spec)

        assert all(s.status == "no_parent_cds" for s in sliced)
        assert all(s.sequence == "" for s in sliced)
        assert all(s.parent_protein_id is None for s in sliced)

    def test_compute_cuts_out_of_order_defensive_guard(self):
        # Defense-in-depth: even if a parent CDS slipped through with
        # out-of-order spans (shouldn't happen post-v0.36.1, but the
        # check stays in place), compute_cuts itself flags it.
        from repseq.polyprotein.slicer import compute_cuts

        spec = _picornavirus_spec("bisect")
        peptide_spans = [
            (spec.peptides[0], (1, 50)),
            (spec.peptides[1], (110, 150)),  # VP2 at 110
            (spec.peptides[2], (60, 100)),   # VP3 at 60 — earlier than VP2
            (spec.peptides[3], (160, 200)),
        ]
        ranges, notes = compute_cuts("Z" * 200, spec, peptide_spans)
        assert all(r is None for r in ranges)
        assert any("out of N-to-C order" in n for n in notes)

    def test_compute_cuts_equal_starts_also_rejected(self):
        # Two peptides with identical span starts are biologically
        # nonsensical (two different peptides cannot start at the same
        # residue). v0.36.1+ uses `>=` instead of `>` so this fails too.
        from repseq.polyprotein.slicer import compute_cuts

        spec = _picornavirus_spec("bisect")
        peptide_spans = [
            (spec.peptides[0], (1, 50)),
            (spec.peptides[1], (60, 100)),
            (spec.peptides[2], (60, 120)),   # same start as VP2
            (spec.peptides[3], (160, 200)),
        ]
        ranges, notes = compute_cuts("Z" * 200, spec, peptide_spans)
        assert all(r is None for r in ranges)
        assert any("out of N-to-C order" in n for n in notes)

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

    def test_method_label_matches_cut_in_overlap_case(self):
        # Regression for the lockstep bug (pre-v0.36.1): when adjacent
        # peptide spans overlap, the first loop placed the cut at
        # midpoint-of-centres but the second loop recomputed the method
        # label using `(a_to + b_from)/2`, so a motif snap that fired in
        # one search window but not the other produced a method label
        # that disagreed with the actual cut location.
        from repseq.polyprotein.slicer import compute_cuts

        # Place an "LQ" motif at position 245 (inside the first-loop's
        # midpoint-of-centres window at 245) and NONE near the second
        # loop's old (200+180)/2=190 bisect point. Pre-fix: actual cut
        # snaps to 247, label says "bisect" (motif not found at 190).
        # Post-fix: label says "motif:LQ" matching the actual cut.
        seq = list("Z" * 500)
        seq[244:246] = list("LQ")  # 0-based 244-245 = 1-based 245-246
        seq_str = "".join(seq)
        spec = PolyproteinSpec(
            name="test",
            peptides=[
                PeptideSpec(name="P1", hmms=["A"]),
                PeptideSpec(name="P2", hmms=["B"], cleavage_motif="LQ"),
            ],
            cut_strategy="motif",
            motif_window_aa=15,
        )
        peptide_spans = [
            (spec.peptides[0], (100, 200)),
            (spec.peptides[1], (180, 500)),  # overlaps by 21 aa
        ]
        ranges, notes = compute_cuts(seq_str, spec, peptide_spans)
        assert ranges[0] is not None and ranges[1] is not None
        # The motif snap fired in the (midpoint-of-centres) window, so
        # the cut and the label must both reflect it.
        assert ranges[1][2] == "motif:LQ"
        # And the actual peptide 2 must start at the motif-snapped
        # position (just after "LQ" at 245-246, so position 247).
        assert ranges[1][0] == 247

    def test_every_overlap_gets_midpoint_of_centres(self):
        # Regression for the "only first overlap handled" bug: pre-v0.36.1
        # only the FIRST overlap in a spec used midpoint-of-centres;
        # subsequent overlaps fell through to (a_to + b_from)/2, which
        # for overlapping spans can produce a bisect point INSIDE
        # peptide A, chopping it short.
        from repseq.polyprotein.slicer import compute_cuts

        spec = PolyproteinSpec(
            name="test",
            peptides=[
                PeptideSpec(name="P1", hmms=["A"]),
                PeptideSpec(name="P2", hmms=["B"]),
                PeptideSpec(name="P3", hmms=["C"]),
            ],
            cut_strategy="bisect",
        )
        # Three peptides with TWO overlap pairs:
        # P1 = (100, 300), P2 = (250, 500), P3 = (450, 700)
        # Pair 1 overlaps by 51 aa; pair 2 overlaps by 51 aa.
        peptide_spans = [
            (spec.peptides[0], (100, 300)),
            (spec.peptides[1], (250, 500)),
            (spec.peptides[2], (450, 700)),
        ]
        ranges, _ = compute_cuts("Z" * 800, spec, peptide_spans)
        # First cut: midpoint(centres of P1, P2) = midpoint(200, 375) = 288
        # Second cut: midpoint(centres of P2, P3) = midpoint(375, 575) = 475
        # Pre-fix the second cut was (500+450)/2 = 475 — happens to match
        # here, so make the asymmetry more pronounced:
        peptide_spans = [
            (spec.peptides[0], (50, 400)),    # centre 225
            (spec.peptides[1], (350, 600)),   # centre 475 — overlaps P1
            (spec.peptides[2], (550, 700)),   # centre 625 — overlaps P2
        ]
        ranges, _ = compute_cuts("Z" * 800, spec, peptide_spans)
        # Cut 1: midpoint(225, 475) = 350
        # Cut 2: midpoint(475, 625) = 550  (vs pre-fix (600+550)/2 = 575)
        # Verify cut 2 used midpoint-of-centres, not (a_to+b_from)/2.
        # P2's end (cut2 - 1) and P3's start (cut2) should reflect 550,
        # not 575.
        assert ranges[1] is not None and ranges[2] is not None
        cut2 = ranges[2][0]
        assert cut2 == 550, (
            f"Cut 2 should be midpoint-of-centres 550 (every overlap "
            f"gets midpoint-of-centres post-v0.36.1), got {cut2}"
        )

    def test_walk_tiebreaker_prefers_best_evalue_when_ali_to_ties(self):
        # When two passing hits for the same HMM share the same ali_to
        # (the leftmost-by-position criterion), the walk should pick
        # the one with the better (lower) dom_evalue. Pre-v0.36.1 this
        # was untiebroken and depended on dict-iteration order.
        from repseq.polyprotein.slicer import _satisfying_span_for_token

        seq = "Z" * 500
        hits = [
            # Two A hits with the same ali_to=100 but different E-values.
            # Both are passing. The lower-E one should anchor the walk.
            _hit("A", 50, 100, evalue=1e-3),
            _hit("A", 60, 100, evalue=1e-30),  # much better E
            _hit("B", 200, 300),
        ]
        # Pick should give us span using the better-E A hit (ali_from=60).
        span = _satisfying_span_for_token(hits, "A--B")
        assert span is not None
        # span_from = min(ali_from of chosen hits) = min(60, 200) = 60.
        # Pre-fix: ali_to-only sort with no E-value tiebreak could pick
        # either A hit; post-fix, dom_evalue tiebreak prefers the 1e-30
        # one (ali_from=60).
        assert span[0] == 60, (
            f"Walk tiebreaker should prefer the better-E A hit "
            f"(ali_from=60), got span_from={span[0]}"
        )

    def test_identify_parent_cds_rejects_chimeric_out_of_order(self):
        # A CDS where every peptide HMM hits but they're scrambled
        # (e.g. a chimeric / contaminated assembly) used to be elected
        # as parent and then fail at compute_cuts with out_of_order.
        # Post-v0.36.1 identify_parent_cds catches the order violation
        # itself and returns None, so the slicer produces a clean
        # no_parent_cds audit row.
        seq = "Z" * 500
        # Spec declared order: VP4, VP2, VP3, VP1
        # CDS reality: VP4, VP3, VP2, VP1 (VP3 and VP2 swapped)
        parent = _polyprotein_cds(seq, [
            _hit("P_VP4", 1, 50),
            _hit("P_VP3", 60, 100),
            _hit("P_VP2", 110, 150),
            _hit("P_VP1", 160, 200),
        ])
        spec = _picornavirus_spec("bisect")
        result = identify_parent_cds([parent], spec)
        assert result is None

    def test_spurious_non_passing_hits_do_not_anchor_multidomain_walk(self):
        # Regression for the coronavirus-NSP15 case: a multidomain peptide
        # token like CoV_NSP15_N--CoV_NSP15_M--CoV_NSP15_C used to be
        # anchored by the LEFTMOST hit of the first HMM, even when that
        # hit was a spurious low-confidence match miles away from the real
        # peptide location. The walk would then pick real downstream
        # hits, and the union span (min ali_from .. max ali_to) would
        # engulf most of the polyprotein, fail the N->C order check
        # against earlier peptides, and mark the whole rep out_of_order.
        # Fix: the slicer filters by the per-hit `passing` flag (set by
        # cli.py:_run_hmm_qc) so non-passing spurious hits are ignored.
        seq = "Z" * 7000
        parent = _polyprotein_cds(seq, [
            # NSP14: real, passing
            _hit("CoV_ExoN", 5500, 6000),
            # NSP15: SPURIOUS early hit on the first HMM of the
            # multidomain token (non-passing — should be filtered out).
            _hit("CoV_NSP15_N", 130, 200, evalue=1.0, passing=False),
            # NSP15: REAL passing hits, all C-terminal.
            _hit("CoV_NSP15_N", 6010, 6100),
            _hit("CoV_NSP15_M", 6110, 6250),
            _hit("CoV_NSP15_C", 6260, 6400),
            # NSP16: real, passing
            _hit("CoV_Methyltr_2", 6450, 6700),
        ])
        spec = PolyproteinSpec(
            name="ORF1ab",
            peptides=[
                PeptideSpec(name="NSP14", hmms=["CoV_ExoN"]),
                PeptideSpec(
                    name="NSP15",
                    hmms=["CoV_NSP15_N--CoV_NSP15_M--CoV_NSP15_C"],
                ),
                PeptideSpec(name="NSP16", hmms=["CoV_Methyltr_2"]),
            ],
            cut_strategy="bisect",
            min_peptides_hit=2,
        )
        _, sliced = slice_polyprotein([parent], spec)
        # With the fix, NSP15's anchor is the real passing hit at 6010,
        # not the spurious one at 130 — order check passes, every
        # peptide slices ok.
        assert [s.status for s in sliced] == ["ok", "ok", "ok"]
        nsp15 = sliced[1]
        # Span starts inside the real NSP15 region, not at the spurious 130.
        assert nsp15.range_aa_from > 6000


# ---------------------------------------------------------------------------
# Multidomain peptide tokens (v0.34.0+)
# ---------------------------------------------------------------------------

def _multidomain_polyprotein() -> dict:
    """Synthetic NSP12-like polyprotein: a CoV_RPol_N domain followed by
    an RdRP_1 catalytic domain, plus a single CoV_NSP8 hit."""
    seq = (
        "M" * 20 +          # leader 1-20
        "N" * 80 +          # CoV_NSP8 region 21-100
        "LQ" +              # cleavage site 101-102
        "P" * 60 +          # CoV_RPol_N 103-162
        "R" * 100 +         # RdRP_1 163-262
        "T" * 10            # tail 263-272
    )
    hits = [
        _hit("CoV_NSP8",  21,  100),
        _hit("CoV_RPol_N", 103, 162),
        _hit("RdRP_1",    163, 262),
    ]
    return _polyprotein_cds(seq, hits)


def _multidomain_spec(cut_strategy: str = "motif") -> PolyproteinSpec:
    return PolyproteinSpec(
        name="ORF1ab",
        peptides=[
            PeptideSpec(name="NSP8", hmms=["CoV_NSP8"]),
            PeptideSpec(
                name="NSP12",
                hmms=["CoV_RPol_N--RdRP_1"],
                cleavage_motif="LQ",
            ),
        ],
        cut_strategy=cut_strategy,
        motif_window_aa=15,
        min_peptides_hit=2,
    )


class TestMultidomainPeptideTokens:
    def test_satisfying_span_single_hmm(self):
        hits = [_hit("CoV_NSP8", 21, 100)]
        span = _satisfying_span_for_token(hits, "CoV_NSP8")
        assert span == (21, 100)

    def test_satisfying_span_multidomain_in_order(self):
        hits = [
            _hit("CoV_RPol_N", 103, 162),
            _hit("RdRP_1",    163, 262),
        ]
        # The span should cover both domains: 103..262.
        span = _satisfying_span_for_token(
            hits, "CoV_RPol_N--RdRP_1", overlap_tolerance=10,
        )
        assert span == (103, 262)

    def test_satisfying_span_returns_none_when_one_domain_missing(self):
        hits = [_hit("CoV_RPol_N", 103, 162)]  # no RdRP_1 hit
        assert _satisfying_span_for_token(
            hits, "CoV_RPol_N--RdRP_1", overlap_tolerance=10,
        ) is None

    def test_satisfying_span_returns_none_when_domains_out_of_order(self):
        # RdRP_1 N-terminal of CoV_RPol_N → architecture violated.
        hits = [
            _hit("RdRP_1",    50,  120),
            _hit("CoV_RPol_N", 130, 200),
        ]
        assert _satisfying_span_for_token(
            hits, "CoV_RPol_N--RdRP_1", overlap_tolerance=10,
        ) is None

    def test_identify_parent_counts_satisfied_tokens(self):
        parent = _multidomain_polyprotein()
        # Decoy carries just CoV_NSP8 — only 1 token satisfied,
        # below min_peptides_hit=2.
        decoy = _polyprotein_cds(
            "X" * 200,
            [_hit("CoV_NSP8", 10, 90)],
            pid="YP_DECOY.1",
        )
        spec = _multidomain_spec()
        chosen = identify_parent_cds(
            [decoy, parent], spec, overlap_tolerance=30,
        )
        assert chosen is not None
        assert chosen["protein_id"] == "YP_TEST.1"

    def test_synthetic_cds_drops_boundary_spillover_hits(self):
        """Regression: neighbouring peptide's HMM that bleeds 1-2 aa
        across the cleavage boundary used to be clipped to a 1-aa
        synthetic domain box in the phyloXML output. The
        majority-residues rule drops those artifacts.
        """
        from repseq.polyprotein import PeptideSpec, PolyproteinSpec
        from repseq.polyprotein.slicer import SlicedPeptide
        from repseq.phylo.per_protein import _build_peptide_synthetic_cds

        parent = {
            "protein_id": "YP_TEST.1",
            "sequence": "X" * 300,
            "length": 300,
            "hmm_hits": [
                # Hit fully inside the peptide (100..200) — should be kept.
                _hit("RealDomain",     110, 190),
                # Hit ending at the peptide's start boundary (50..100) →
                # clips to ali_from=1, ali_to=1 — pure spillover, drop.
                _hit("SpilloverLeft",   50, 100),
                # Hit starting at the peptide's end boundary (200..280) →
                # clips to ali_from=101, ali_to=101 — pure spillover, drop.
                _hit("SpilloverRight", 200, 280),
                # Hit that straddles the start but >50% inside (95..150) →
                # 56 aa total, 51 aa inside — kept.
                _hit("MostlyIn",        95, 150),
            ],
        }
        sliced = SlicedPeptide(
            peptide_name="NSP_test",
            parent_protein_id="YP_TEST.1",
            parent_accession="NC_TEST",
            range_aa_from=100,
            range_aa_to=200,
            length_aa=101,
            sequence="X" * 101,
            cut_method_actual="boundary",
            status="ok",
        )
        synthetic = _build_peptide_synthetic_cds(parent, sliced, "NSP_test")
        kept_names = sorted(h.get("target") for h in synthetic["hmm_hits"])
        assert kept_names == ["MostlyIn", "RealDomain"]
        # And no domain of length ≤2 leaks through.
        assert all(
            int(h["ali_to"]) - int(h["ali_from"]) + 1 >= 3
            for h in synthetic["hmm_hits"]
        )

    def test_run_polyprotein_phylogeny_builds_one_tree_per_peptide(self, tmp_path):
        """End-to-end: 3 reps × 2 peptides → 2 trees, both produced
        because both peptides clear min_taxa=3."""
        import shutil as _shutil
        if not _shutil.which("mafft"):
            import pytest as _pytest
            _pytest.skip("mafft not on PATH")
        if not (_shutil.which("FastTree") or _shutil.which("fasttree")):
            import pytest as _pytest
            _pytest.skip("FastTree not on PATH")

        from repseq.models import RunResult, Sequence, SequenceType
        from repseq.phylo.per_protein import run_polyprotein_phylogeny

        # Build 3 reps each carrying a multidomain polyprotein with two
        # well-spaced mature peptides. The protein has CoV_NSP8 hits in
        # the N-terminal half and CoV_NSP12 ("RdRP") hits in the C-half.
        def _make_rep(rid: str, accent: str) -> Sequence:
            poly = (
                "M" + ("A" + accent) * 60   # NSP8 region (~120 aa) + filler
                + "LQ" +
                "G" + ("B" + accent) * 60   # NSP12 region (~120 aa)
            )
            n_len = len(poly)
            return Sequence(
                id=rid,
                header=rid,
                sequence="ACGT" * n_len,
                seq_type=SequenceType.NUCLEOTIDE,
                accession=rid,
                organism=f"Test virus {rid}",
                proteins=[{
                    "protein_id": f"YP_{rid}.1",
                    "product": "polyprotein",
                    "length": n_len,
                    "sequence": poly,
                    "hmm_hits": [
                        _hit("CoV_NSP8",  1,   121),
                        _hit("CoV_NSP12", 124, n_len),
                    ],
                }],
            )

        reps = [_make_rep("R001", "X"), _make_rep("R002", "Y"), _make_rep("R003", "Z")]
        cfg = {
            "_hmm_runtime": {"active": True},
            "clustering": {
                "polyprotein": [
                    {
                        "name": "P1",
                        "peptides": [
                            {"name": "NSP8",  "hmm": "CoV_NSP8"},
                            {"name": "NSP12", "hmm": "CoV_NSP12", "cleavage_motif": "LQ"},
                        ],
                        "min_peptides_hit": 2,
                    }
                ]
            },
            "segmented": {"enabled": False},
            "hmm": {"enabled": True, "multidomain_overlap_tolerance": 30},
            "phylo": {
                "tool": "fasttree",
                "per_protein": {"min_taxa": 3, "mafft": {"extra_args": []}},
                "coloring": {"enabled": False},  # speed up the test
                "lca": {"enabled": False},
                "rooting": {"method": "midpoint"},
            },
        }

        written = run_polyprotein_phylogeny(reps, cfg, tmp_path, "test")
        sub = tmp_path / "test_polyprotein"
        assert sub.exists()
        # Expect two trees built. _build_tree writes multiple files
        # (MSA, Newick, phyloXML, id_map); confirm one Newick per peptide.
        nwks = sorted(p.name for p in sub.glob("*.nwk"))
        assert nwks == ["P1_NSP12_tree.nwk", "P1_NSP8_tree.nwk"]
        assert any(p.name.endswith("_tree.xml") for p in written)

    def test_slice_multidomain_peptide_emits_full_footprint(self):
        from repseq.polyprotein.slicer import slice_polyprotein
        parent = _multidomain_polyprotein()
        spec = _multidomain_spec()
        chosen, sliced = slice_polyprotein(
            [parent], spec, overlap_tolerance=30,
        )
        assert chosen is not None
        # Two peptides, both ok. NSP12's range should cover the union
        # of its two domains (103..262 minus the inter-peptide cut).
        assert all(s.status == "ok" for s in sliced)
        # NSP12's start should snap to the residue after the LQ at 101-102,
        # so start = 103 (motif:LQ cut).
        nsp12 = sliced[1]
        assert nsp12.range_aa_from == 103
        assert nsp12.range_aa_to == len(parent["sequence"])
        assert nsp12.cut_method_actual == "motif:LQ"


# ---------------------------------------------------------------------------
# OR semantics across alternative peptide architectures (v0.34.0+)
# ---------------------------------------------------------------------------

class TestPeptideOrAlternatives:
    """A peptide whose architecture varies across virus genera carries
    multiple alternative tokens in ``hmms:``. Modelled here on the
    alpha- vs. beta-CoV NSP1 case: ``aCoV_NSP1`` and ``bCoV_NSP1`` are
    non-homologous HMMs; either one should locate NSP1 on the parent.
    """

    def test_best_alternative_returns_satisfied_token(self):
        # Only the beta-CoV NSP1 HMM hits.
        hits = [
            _hit("bCoV_NSP1", 1, 110),
            _hit("CoV_NSP2",  130, 200),
        ]
        out = _best_satisfying_alternative(
            hits, ["aCoV_NSP1", "bCoV_NSP1"],
        )
        assert out is not None
        assert out[2] == "bCoV_NSP1"
        assert out[0] == 1 and out[1] == 110

    def test_best_alternative_picks_lower_e_when_both_satisfied(self):
        # Both alternatives hit; the better-E one wins.
        hits = [
            _hit("aCoV_NSP1", 1, 100, evalue=1e-40),
            _hit("bCoV_NSP1", 1, 110, evalue=1e-10),
        ]
        out = _best_satisfying_alternative(
            hits, ["aCoV_NSP1", "bCoV_NSP1"],
        )
        assert out is not None
        assert out[2] == "aCoV_NSP1"

    def test_best_alternative_returns_none_when_no_token_satisfied(self):
        hits = [_hit("CoV_NSP2", 50, 100)]
        assert _best_satisfying_alternative(
            hits, ["aCoV_NSP1", "bCoV_NSP1"],
        ) is None

    def test_slice_records_matched_architecture(self):
        # Synthetic alpha-CoV ORF1ab fragment with two peptides; NSP1
        # has two alternative architectures, the alpha one fires.
        seq = "A" * 100 + "LQ" + "B" * 80
        parent = _polyprotein_cds(seq, [
            _hit("aCoV_NSP1", 1, 100),
            _hit("CoV_NSP2",  103, 180),
        ])
        spec = PolyproteinSpec(
            name="ORF1ab",
            peptides=[
                PeptideSpec(
                    name="NSP1",
                    hmms=["aCoV_NSP1", "bCoV_NSP1"],
                ),
                PeptideSpec(
                    name="NSP2",
                    hmms=["CoV_NSP2"],
                    cleavage_motif="LQ",
                ),
            ],
            cut_strategy="motif",
            motif_window_aa=10,
            min_peptides_hit=2,
        )
        _, sliced = slice_polyprotein([parent], spec)
        assert sliced[0].status == "ok"
        assert sliced[0].matched_token == "aCoV_NSP1"
        assert sliced[1].matched_token == "CoV_NSP2"


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
        # Legacy `hmm: <token>` wraps to a 1-element `hmms:` list.
        assert specs[0].peptides[0].hmms == ["P_VP4"]
        assert specs[0].peptides[1].hmms == ["P_VP2"]
        # cleavage_motif on one peptide → default strategy should be "motif"
        assert specs[0].cut_strategy == "motif"

    def test_hmms_list_round_trips(self):
        cfg = {
            "clustering": {
                "polyprotein": [
                    {
                        "name": "ORF1ab",
                        "peptides": [
                            {
                                "name": "NSP1",
                                "hmms": ["aCoV_NSP1", "bCoV_NSP1"],
                            },
                            {"name": "NSP2", "hmm": "CoV_NSP2"},
                        ],
                    }
                ]
            }
        }
        specs = collect_polyprotein_specs(cfg)
        assert specs[0].peptides[0].hmms == ["aCoV_NSP1", "bCoV_NSP1"]
        assert specs[0].peptides[1].hmms == ["CoV_NSP2"]

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
                    {"name": "VP4"},  # neither hmm nor hmms
                    {"name": "VP2", "hmm": "P_VP2"},
                ],
            }
        ])
        assert any(
            "must set either 'hmm'" in e or "must set either \"hmm\"" in e
            for e in errs
        )

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

    def test_accepts_multidomain_hmm_token(self):
        errs = self._validate([
            {
                "name": "ORF1ab",
                "peptides": [
                    {"name": "NSP8", "hmm": "CoV_NSP8"},
                    {"name": "NSP12", "hmm": "CoV_RPol_N--RdRP_1"},
                    {"name": "NSP13", "hmm": "CoV_ZBD--CoV_stalk--CoV_1B"},
                ],
            }
        ])
        assert errs == []

    def test_rejects_malformed_hmm_token(self):
        errs = self._validate([
            {
                "name": "P1",
                "peptides": [
                    {"name": "VP4", "hmm": "A----B"},  # empty middle
                    {"name": "VP2", "hmm": "P_VP2"},
                ],
            }
        ])
        assert any("empty component" in e for e in errs)

    def test_accepts_hmms_list_for_or_semantics(self):
        errs = self._validate([
            {
                "name": "ORF1ab",
                "peptides": [
                    {"name": "NSP1", "hmms": ["aCoV_NSP1", "bCoV_NSP1"]},
                    {"name": "NSP2", "hmm": "CoV_NSP2"},
                ],
            }
        ])
        assert errs == []

    def test_rejects_setting_both_hmm_and_hmms(self):
        errs = self._validate([
            {
                "name": "P1",
                "peptides": [
                    {"name": "X", "hmm": "A", "hmms": ["B", "C"]},
                    {"name": "Y", "hmm": "D"},
                ],
            }
        ])
        assert any(
            "set exactly one of 'hmm'" in e or "exactly one of \"hmm\"" in e
            for e in errs
        )

    def test_rejects_empty_hmms_list(self):
        errs = self._validate([
            {
                "name": "P1",
                "peptides": [
                    {"name": "X", "hmms": []},
                    {"name": "Y", "hmm": "D"},
                ],
            }
        ])
        assert any("hmms must be a non-empty list" in e for e in errs)

    def test_rejects_malformed_token_in_hmms_list(self):
        errs = self._validate([
            {
                "name": "P1",
                "peptides": [
                    {"name": "X", "hmms": ["A--B", "C----D"]},  # second is bad
                    {"name": "Y", "hmm": "E"},
                ],
            }
        ])
        assert any("hmms[1]" in e and "empty component" in e for e in errs)

    def test_rejects_bad_cut_strategy(self):
        errs = self._validate([
            {"name": "P1", "cut_strategy": "magic", "peptides": [
                {"name": "VP4", "hmm": "P_VP4"},
                {"name": "VP2", "hmm": "P_VP2"},
            ]}
        ])
        assert any("cut_strategy 'magic' is not supported" in e for e in errs)
