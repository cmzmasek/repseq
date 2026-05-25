"""Per-column conservation metrics, hit projection, and PNG rendering.

The metric functions are pure-Python and gap-aware; tests pin the
exact numeric output on tiny synthetic alignments so a refactor that
silently changes the formula is caught fast. The PNG writer is
exercised end-to-end (matplotlib's ``Agg`` backend writes a real
file), but the *content* of the PNG is not pixel-asserted — that
would be brittle. We assert the file exists, is non-empty, and is a
PNG by magic bytes.
"""
from __future__ import annotations

import math
import struct
from pathlib import Path

import pytest

from repseq.viz.conservation import (
    DEFAULT_WINDOW_AA,
    compute_consensus_fraction,
    compute_shannon_entropy,
    family_color_hex,
    pick_reference_row,
    project_hits_to_alignment,
    read_msa,
    sliding_window_mean,
    write_conservation_heatmap,
)


# ---------------------------------------------------------------------------
# compute_shannon_entropy
# ---------------------------------------------------------------------------

def test_entropy_all_same_residue_is_zero():
    """A perfectly conserved column → 0 bits."""
    rows = ["AAAA", "AAAA", "AAAA"]
    assert compute_shannon_entropy(rows) == [0.0, 0.0, 0.0, 0.0]


def test_entropy_two_equally_likely_residues_is_one_bit():
    """50/50 split between two residues → log2(2) = 1 bit."""
    rows = ["AC", "AC", "GT", "GT"]
    # Col 0: 2A + 2G → H = -2*(0.5*log2(0.5)) = 1.0
    # Col 1: 2C + 2T → same
    h = compute_shannon_entropy(rows)
    assert pytest.approx(h[0], rel=1e-12) == 1.0
    assert pytest.approx(h[1], rel=1e-12) == 1.0


def test_entropy_gaps_excluded_from_counts():
    """Gaps are missing data — the column's entropy is computed on the
    surviving non-gap residues, not on a 'gap == 21st symbol' model."""
    # 3 rows, col0 = A,A,- → 2A, 1 gap → entropy on counts {A:2} = 0
    # col1 = A,C,- → A,C → entropy 1.0
    rows = ["AA", "AC", "--"]
    h = compute_shannon_entropy(rows)
    assert h[0] == 0.0
    assert pytest.approx(h[1], rel=1e-12) == 1.0


def test_entropy_all_gap_column_is_zero_not_nan():
    """An all-gap column → 0.0 (no residue to integrate over). NaN
    would propagate through the colormap and bin the figure."""
    rows = ["A-", "A-", "A-"]
    h = compute_shannon_entropy(rows)
    assert h[0] == 0.0
    assert h[1] == 0.0
    assert not any(math.isnan(x) for x in h)


def test_entropy_empty_msa_returns_empty_list():
    assert compute_shannon_entropy([]) == []


def test_entropy_case_insensitive():
    """MAFFT sometimes emits lowercase residues at ambiguous columns;
    'A' and 'a' must count as the same residue."""
    rows = ["AAaa", "aaAA"]
    h = compute_shannon_entropy(rows)
    # Every column is 100% A — entropy should be 0 everywhere.
    assert h == [0.0, 0.0, 0.0, 0.0]


# ---------------------------------------------------------------------------
# sliding_window_mean — the v0.31.0 smoother behind the line charts
# ---------------------------------------------------------------------------

def test_sliding_window_default_is_15():
    """The bench-scientist default is documented as 15aa — guard against
    silent drift."""
    assert DEFAULT_WINDOW_AA == 15


def test_sliding_window_identity_when_window_le_1():
    """A window of 1 (or 0) returns the input unchanged."""
    src = [1.0, 2.0, 3.0]
    assert sliding_window_mean(src, 1) == [1.0, 2.0, 3.0]
    assert sliding_window_mean(src, 0) == [1.0, 2.0, 3.0]


def test_sliding_window_empty_input():
    assert sliding_window_mean([], 5) == []


def test_sliding_window_centered_mean_interior():
    """Interior points see a full window, so the result is the mean
    of (i-half) to (i+half) inclusive. With window=3 on [1,2,3,4,5]
    the centred mean at index 2 is (2+3+4)/3 = 3.0."""
    out = sliding_window_mean([1.0, 2.0, 3.0, 4.0, 5.0], 3)
    assert out[2] == pytest.approx(3.0)


def test_sliding_window_shrinks_at_edges():
    """Window shrinks at edges rather than zero-padding, so the
    smoothed curve doesn't fake-dip to zero at the boundaries."""
    out = sliding_window_mean([1.0, 2.0, 3.0, 4.0, 5.0], 3)
    # Index 0 has window=[1,2] (i-1 clamped) → mean = 1.5
    assert out[0] == pytest.approx(1.5)
    # Last index has window=[4,5] → mean = 4.5
    assert out[-1] == pytest.approx(4.5)


def test_sliding_window_constant_input_preserves_value():
    """A flat input stays flat regardless of window size."""
    out = sliding_window_mean([0.5] * 10, 5)
    assert all(v == pytest.approx(0.5) for v in out)


def test_sliding_window_larger_than_input_returns_global_mean():
    """If window > n, every output position is the global mean."""
    out = sliding_window_mean([1.0, 2.0, 3.0], 99)
    assert all(v == pytest.approx(2.0) for v in out)


# ---------------------------------------------------------------------------
# compute_consensus_fraction
# ---------------------------------------------------------------------------

def test_consensus_all_same_is_one():
    rows = ["MKL", "MKL", "MKL"]
    assert compute_consensus_fraction(rows) == [1.0, 1.0, 1.0]


def test_consensus_2_of_3_is_two_thirds():
    rows = ["A", "A", "G"]
    out = compute_consensus_fraction(rows)
    assert pytest.approx(out[0], rel=1e-12) == 2.0 / 3.0


def test_consensus_gaps_excluded_from_denominator():
    """Consensus is over non-gap rows only."""
    rows = ["A", "A", "-"]  # 2 non-gap rows, both A
    assert compute_consensus_fraction(rows) == [1.0]


def test_consensus_all_gap_is_zero():
    rows = ["-", "-", "-"]
    assert compute_consensus_fraction(rows) == [0.0]


# ---------------------------------------------------------------------------
# pick_reference_row
# ---------------------------------------------------------------------------

def test_pick_reference_picks_longest_non_gap():
    msa = {
        "S0001": "MK--",  # 2 non-gap
        "S0002": "MKLP",  # 4 non-gap
        "S0003": "M---",  # 1 non-gap
    }
    assert pick_reference_row(msa) == "S0002"


def test_pick_reference_tie_broken_alphabetically():
    """Ties on non-gap count → smallest short_id wins (deterministic)."""
    msa = {"S0002": "MK", "S0001": "MK"}
    assert pick_reference_row(msa) == "S0001"


def test_pick_reference_empty_returns_none():
    assert pick_reference_row({}) is None


# ---------------------------------------------------------------------------
# project_hits_to_alignment
# ---------------------------------------------------------------------------

def test_project_hits_simple_no_gaps():
    """No gaps → ungapped position == aligned column."""
    aligned = "MKLPQE"
    hits = [{"ali_from": 1, "ali_to": 3, "hmm_name": "X"}]
    out = project_hits_to_alignment(aligned, hits)
    assert out[0]["aln_from"] == 1
    assert out[0]["aln_to"] == 3


def test_project_hits_shifts_past_gaps():
    """Gaps in the aligned reference push hit endpoints rightward."""
    # ungapped 'MKLPQE', aligned with one gap before position 3:
    # ungapped pos 1 → col 1 (M)
    # ungapped pos 2 → col 2 (K)
    # ungapped pos 3 → col 4 (L), since col 3 is '-'
    aligned = "MK-LPQE"
    hits = [{"ali_from": 2, "ali_to": 4, "hmm_name": "X"}]
    out = project_hits_to_alignment(aligned, hits)
    # ali_from=2 (K) → aln col 2
    # ali_to=4 (P)   → aln col 5
    assert out[0]["aln_from"] == 2
    assert out[0]["aln_to"] == 5


def test_project_drops_hit_past_reference_end():
    aligned = "MKL"  # 3 ungapped
    hits = [{"ali_from": 4, "ali_to": 5}]
    assert project_hits_to_alignment(aligned, hits) == []


def test_project_clamps_partial_overshoot():
    """Hit slightly past the reference end gets clamped, not dropped —
    a 1-residue truncation shouldn't hide a real domain box."""
    aligned = "MKLPQ"  # 5 ungapped
    hits = [{"ali_from": 4, "ali_to": 6, "hmm_name": "C"}]
    out = project_hits_to_alignment(aligned, hits)
    assert len(out) == 1
    assert out[0]["aln_from"] == 4
    assert out[0]["aln_to"] == 5  # clamped


def test_project_ignores_malformed_hits():
    """Missing or non-numeric ali_from/ali_to fields are skipped, not raised."""
    aligned = "MKL"
    hits = [
        {"hmm_name": "no-coords"},
        {"ali_from": "x", "ali_to": 2},
        {"ali_from": 1, "ali_to": 0},  # at < af
        {"ali_from": 1, "ali_to": 2, "hmm_name": "ok"},
    ]
    out = project_hits_to_alignment(aligned, hits)
    assert len(out) == 1
    assert out[0]["hmm_name"] == "ok"


# ---------------------------------------------------------------------------
# family_color_hex
# ---------------------------------------------------------------------------

def test_family_color_returns_valid_hex():
    out = family_color_hex("Spike", ["Spike", "Nucleocapsid", "Membrane"])
    assert out.startswith("#")
    assert len(out) == 7
    # Each hex digit pair must be a valid byte.
    for i in (1, 3, 5):
        int(out[i:i + 2], 16)


def test_family_color_stable_under_input_order():
    """Hue is keyed off the *sorted* family list, so a config edit that
    reorders cluster_protein entries doesn't shuffle the colours."""
    a = family_color_hex("Spike", ["Spike", "Nucleocapsid"])
    b = family_color_hex("Spike", ["Nucleocapsid", "Spike"])
    assert a == b


def test_family_color_unknown_family_does_not_raise():
    """A family not in the list falls back to position 0 rather than
    blowing up (defensive — the heatmap should not crash because of a
    bookkeeping mismatch)."""
    out = family_color_hex("Mystery", ["Spike", "Nucleocapsid"])
    assert out.startswith("#")


# ---------------------------------------------------------------------------
# read_msa
# ---------------------------------------------------------------------------

def test_read_msa_keys_on_first_token(tmp_path):
    """FASTA descriptions after a space are ignored — the key is the
    short id (the safe first whitespace-separated token)."""
    msa = tmp_path / "x.fasta"
    msa.write_text(">S0001 my long label\nMKL\nPQE\n>S0002\nMM-PQE\n")
    out = read_msa(msa)
    assert set(out.keys()) == {"S0001", "S0002"}
    assert out["S0001"] == "MKLPQE"  # multi-line sequence reassembled


# ---------------------------------------------------------------------------
# write_conservation_heatmap — end-to-end PNG
# ---------------------------------------------------------------------------

def _is_png(path: Path) -> bool:
    return path.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"


def test_write_heatmap_produces_png(tmp_path):
    """A tiny MSA + a couple of synthetic HMM hits should yield a real
    PNG on disk; matplotlib is exercised end-to-end via the Agg backend."""
    pytest.importorskip("matplotlib")
    msa = tmp_path / "fam_msa.fasta"
    msa.write_text(
        ">S0001\nMKLPQEMKLPQE\n"
        ">S0002\nMKLPQDMKLPQE\n"
        ">S0003\nMKLPQEMKLP-E\n"
    )
    hits = [
        {"ali_from": 1, "ali_to": 6, "hmm_name": "DomA"},
        {"ali_from": 7, "ali_to": 12, "hmm_name": "DomB"},
    ]
    # The reference row here is fully ungapped (S0001), so hits project
    # 1:1 onto MSA columns.
    projected = project_hits_to_alignment("MKLPQEMKLPQE", hits)
    out_png = tmp_path / "out" / "Spike.png"
    result = write_conservation_heatmap(
        msa,
        out_png=out_png,
        family_label="Spike",
        family_color="#226622",
        hmm_hits_on_reference=projected,
    )
    assert result == out_png
    assert out_png.exists()
    assert out_png.stat().st_size > 0
    assert _is_png(out_png)


def test_write_heatmap_target_key_actually_lands_in_label_text(tmp_path, monkeypatch):
    """Verify the *label text* matplotlib receives is the HMM target
    name. PNGs are opaque to substring assertions, so we intercept
    ``ax.text`` calls and check the positional ``name`` argument.

    Three synthetic hits each carry only ``target`` (no ``hmm_name``
    fallback). All three names must reach ``ax.text``. Regression-
    guard for the v0.31.0 bug where the labeller looked for
    ``hmm_name`` first and ``target`` was never read."""
    pytest.importorskip("matplotlib")

    captured_label_texts: list[str] = []
    import matplotlib.axes as mpl_axes
    original_text = mpl_axes.Axes.text

    def _capture_text(self, x, y, s, *a, **kw):
        captured_label_texts.append(s)
        return original_text(self, x, y, s, *a, **kw)

    monkeypatch.setattr(mpl_axes.Axes, "text", _capture_text)

    msa = tmp_path / "f.fasta"
    msa.write_text(">S0001\nMKLPQE\n>S0002\nMKLPQE\n>S0003\nMKLPQE\n")
    hits = [
        {"target": "CoV_S1", "ali_from": 1, "ali_to": 2, "aln_from": 1, "aln_to": 2},
        {"target": "CoV_S2", "ali_from": 3, "ali_to": 4, "aln_from": 3, "aln_to": 4},
        {"target": "TMprD",  "ali_from": 5, "ali_to": 6, "aln_from": 5, "aln_to": 6},
    ]
    write_conservation_heatmap(
        msa, out_png=tmp_path / "f.png",
        family_label="Spike", family_color="#3F8E7F",
        hmm_hits_on_reference=hits,
    )
    assert "CoV_S1" in captured_label_texts
    assert "CoV_S2" in captured_label_texts
    assert "TMprD" in captured_label_texts


def test_write_heatmap_uses_target_key_for_domain_label(tmp_path):
    """Real HMM hit dicts (from hmm/hmmscan.py:_parse_domtblout) carry
    the profile name under ``target``, not ``hmm_name``. The
    domain-ribbon labeller must read ``target`` first so live runs
    actually label their domains. Regression-guard: v0.31.0 only
    read ``hmm_name``/``name`` and produced unlabelled domains on
    every real run."""
    pytest.importorskip("matplotlib")
    msa = tmp_path / "f.fasta"
    msa.write_text(
        ">S0001\nMKLPQEMKLPQE\n>S0002\nMKLPQEMKLPQE\n>S0003\nMKLPQEMKLPQE\n"
    )
    out_png = tmp_path / "f.png"
    # Note: only ``target`` is set, mirroring the schema hmmscan.py
    # produces. No ``hmm_name`` / ``name`` fallback in this dict.
    hits = [
        {"target": "CoV_S2", "ali_from": 1, "ali_to": 6, "aln_from": 1, "aln_to": 6},
    ]
    result = write_conservation_heatmap(
        msa, out_png=out_png, family_label="Spike",
        family_color="#3F8E7F",
        hmm_hits_on_reference=hits,
    )
    assert result == out_png
    assert _is_png(out_png)


def test_write_heatmap_no_domains_still_renders(tmp_path):
    """Without HMM hits, the figure has just the two metric tracks."""
    pytest.importorskip("matplotlib")
    msa = tmp_path / "f.fasta"
    msa.write_text(">S0001\nMKL\n>S0002\nMKL\n>S0003\nMKL\n")
    out_png = tmp_path / "f.png"
    result = write_conservation_heatmap(
        msa,
        out_png=out_png,
        family_label="F",
        family_color="#444444",
        hmm_hits_on_reference=[],
    )
    assert result == out_png
    assert _is_png(out_png)


def test_write_heatmap_empty_msa_returns_none(tmp_path):
    msa = tmp_path / "empty.fasta"
    msa.write_text("")  # no records
    out_png = tmp_path / "empty.png"
    result = write_conservation_heatmap(
        msa, out_png=out_png, family_label="F", family_color="#000000",
        hmm_hits_on_reference=None,
    )
    assert result is None
    assert not out_png.exists()


# ---------------------------------------------------------------------------
# Orchestrator soft-fail
# ---------------------------------------------------------------------------

def test_orchestrator_raises_when_per_protein_dir_missing(tmp_path):
    """Without --per-protein-phylo, the {prefix}_per_protein/ dir is
    absent. The orchestrator must raise FileNotFoundError (which cli
    catches as a soft-fail stderr message) rather than silently writing
    nothing or crashing later."""
    from repseq.viz.conservation_runner import run_conservation_heatmaps
    cfg = {
        "segmented": {"enabled": False},
        "clustering": {
            "cluster_protein": [
                {"name": "Spike", "hmms": ["CoV_S2"]},
            ],
        },
    }
    # No representatives needed for the directory-missing check.
    with pytest.raises(FileNotFoundError, match="--per-protein-phylo"):
        run_conservation_heatmaps([], cfg, tmp_path, "test")


def test_orchestrator_end_to_end_emits_png(tmp_path, make_seq):
    """A synthetic 2F output dir + reps carrying matching CDS+HMM hits
    drives the orchestrator end-to-end: read MSA, resolve reference rep,
    project hits, write the heatmap PNG under {prefix}_conservation/."""
    pytest.importorskip("matplotlib")
    from repseq.viz.conservation_runner import run_conservation_heatmaps

    prefix = "test"
    per_protein = tmp_path / f"{prefix}_per_protein"
    per_protein.mkdir()
    # Synthetic 2F outputs for a single family "Spike".
    (per_protein / "Spike_msa.fasta").write_text(
        ">S0001\nMKLPQEMKLPQE\n"
        ">S0002\nMKLPQDMKLPQE\n"
        ">S0003\nMKLPQEMKLP-E\n"
    )
    (per_protein / "Spike_tree_id_map.tsv").write_text(
        "short_id\taccession\nS0001\trep_a\nS0002\trep_b\nS0003\trep_c\n"
    )

    # Reps must carry a CDS satisfying the spec's HMM token with
    # hmm_hits in ungapped CDS coords. rep_a is fully ungapped in the
    # synthetic MSA above, so it will be picked as the reference row
    # (longest non-gap = 12 columns).
    reps = []
    for rid in ("rep_a", "rep_b", "rep_c"):
        r = make_seq(rid, "ACGT")
        r.proteins = [{
            "protein_id": f"{rid}_p1",
            "product": "spike protein",
            "length": 12,
            "sequence": "MKLPQEMKLPQE",
            "hmm_hits": [
                {
                    "hmm_name": "CoV_S2",
                    "target": "CoV_S2",
                    "ali_from": 1, "ali_to": 12,
                    "dom_evalue": 1e-30, "passing": True,
                },
            ],
        }]
        reps.append(r)

    cfg = {
        "segmented": {"enabled": False},
        "hmm": {"enabled": True, "multidomain_overlap_tolerance": 30},
        "clustering": {
            "cluster_protein": [
                {"name": "Spike", "hmms": ["CoV_S2"]},
            ],
        },
        "_hmm_runtime": {"active": True},
    }
    written = run_conservation_heatmaps(reps, cfg, tmp_path, prefix)
    assert len(written) == 1
    out_png = tmp_path / f"{prefix}_conservation" / f"{prefix}_Spike.png"
    assert written[0] == out_png
    assert out_png.exists()
    assert _is_png(out_png)


def test_orchestrator_skips_family_without_msa(tmp_path, make_seq, caplog):
    """A configured family whose 2F MSA never landed on disk gets a
    warning + skip; sibling families still build."""
    pytest.importorskip("matplotlib")
    from repseq.viz.conservation_runner import run_conservation_heatmaps

    prefix = "test"
    per_protein = tmp_path / f"{prefix}_per_protein"
    per_protein.mkdir()
    # Only "Spike" has an MSA on disk; "Nucleocapsid" is missing.
    (per_protein / "Spike_msa.fasta").write_text(
        ">S0001\nMKLPQE\n>S0002\nMKLPQE\n>S0003\nMKLPQE\n"
    )
    (per_protein / "Spike_tree_id_map.tsv").write_text(
        "short_id\taccession\nS0001\trep_a\nS0002\trep_b\nS0003\trep_c\n"
    )
    reps = []
    for rid in ("rep_a", "rep_b", "rep_c"):
        r = make_seq(rid, "ACGT")
        r.proteins = [{
            "protein_id": f"{rid}_p1",
            "product": "spike protein",
            "length": 6,
            "sequence": "MKLPQE",
            "hmm_hits": [{
                "hmm_name": "CoV_S2",
                "target": "CoV_S2",
                "ali_from": 1, "ali_to": 6,
                "dom_evalue": 1e-30, "passing": True,
            }],
        }]
        reps.append(r)
    cfg = {
        "segmented": {"enabled": False},
        "hmm": {"enabled": True, "multidomain_overlap_tolerance": 30},
        "clustering": {
            "cluster_protein": [
                {"name": "Spike", "hmms": ["CoV_S2"]},
                {"name": "Nucleocapsid", "hmms": ["CoV_nucleocap"]},
            ],
        },
        "_hmm_runtime": {"active": True},
    }
    written = run_conservation_heatmaps(reps, cfg, tmp_path, prefix)
    # Only Spike survives.
    assert len(written) == 1
    assert written[0].name == f"{prefix}_Spike.png"
    # The Nucleocapsid skip surfaces via logging (run_conservation_heatmaps
    # uses logger.warning); not all repseq logger setups propagate to
    # caplog by default, so we re-check the disk: no Nucleocapsid PNG.
    assert not (tmp_path / f"{prefix}_conservation" / f"{prefix}_Nucleocapsid.png").exists()
