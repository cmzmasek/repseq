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
