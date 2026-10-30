"""Threshold-search fill behaviour (``taxonomic1._binary_search_threshold``).

The search brackets the requested representative count between the largest
cluster count at-or-below it and the smallest above it. Cluster count is a
step function of identity, so the two brackets are often far apart and the
old code always kept the lower one — a group that stepped 1 → 7 returned 1
representative for a target of 5. On the H5N1 per-clade run (`custom
--field clade -n 5`) that left 32 of 45 clades short, clade 2.3.2.1c at 1 of
27 isolates.

Three changes (v0.80.0), pinned here:

* **resolution** — bisect to ``THRESHOLD_RESOLUTION`` (0.0005) instead of
  0.005, so a step hidden inside a narrow window near 99 % identity is
  found and the target can be hit exactly;
* **bracket choice** — ``overflow="trim"`` always takes the bracket above
  (only an above-target count can be pared to exactly the target), while
  ``overflow="keep"`` takes whichever bracket is closer, ties going up;
* **trim keeps the stratum whole** — the clusters the diversity sampler did
  not pick are folded into the nearest kept cluster instead of being
  dropped with their members.

The fakes below reproduce the step shapes seen in the real run logs rather
than a smooth count ~ threshold curve: that smoothness is what hid the bug.
Several tests pass ``resolution=`` explicitly to contrast the two stops —
without that contrast a revert of the default would keep the suite green.
"""
from __future__ import annotations

from unittest.mock import patch

from repseq.models import Cluster
from repseq.modes.taxonomic1 import THRESHOLD_RESOLUTION, _binary_search_threshold


def _cfg() -> dict:
    return {
        "seed": 42,
        "threads": 1,
        "clustering": {
            "backend": "mmseqs2",
            "mmseqs2_mode": "easy-linclust",
            "coverage": 0.8,
            "coverage_mode": 0,
            "extra_args": [],
            # Real runs cluster on the marker protein at k=5; the fixtures
            # below populate protein_sequence to match, so select_diverse
            # sees the same basis production does.
            "alphabet_for_clustering": "protein",
        },
        "representative": {"priority": ["refseq", "reviewed_uniprot", "longest"]},
    }


def _seqs(make_seq, n_seqs):
    """Sequences that are genuinely distinguishable under protein k=5.

    Each carries a distinct 60-residue protein block, so containment
    distances are well spread and MaxMin selection is driven by sequence
    content rather than by ``rng.choice`` over an all-ties pool.
    """
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    seqs = []
    for i in range(n_seqs):
        seq = make_seq(f"s{i}", "ACGT" * 30)
        residue = alphabet[i % len(alphabet)]
        seq.protein_sequence = (residue * 20) + alphabet * ((i % 3) + 1)
        seqs.append(seq)
    return seqs


def _step_clustering(steps):
    """A fake backend whose cluster count is a step function of identity.

    ``steps`` is ``[(threshold_at_or_above, n_clusters), ...]`` ascending;
    the highest matching entry wins. Members are distributed round-robin
    over the clusters, so cluster sizes cover the whole input the way a
    real backend's do.
    """
    def fake(sequences, threshold, cfg, tmp_dir=None):
        n = steps[0][1]
        for edge, count in steps:
            if threshold >= edge:
                n = count
        n = min(n, len(sequences))
        clusters = [
            Cluster(cluster_id=f"c{i}", representative=sequences[i])
            for i in range(n)
        ]
        for j, seq in enumerate(sequences[n:]):
            clusters[j % n].members.append(seq)
        return clusters
    return fake


def _run(make_seq, steps, n_target, overflow="keep", n_seqs=30, **kwargs):
    seqs = _seqs(make_seq, n_seqs)
    with patch("repseq.modes.taxonomic1.run_clustering",
               side_effect=_step_clustering(steps)):
        return _binary_search_threshold(
            seqs, n_target=n_target, cfg=_cfg(), overflow=overflow, **kwargs
        ), seqs


# ---------------------------------------------------------------------------
# Resolution: a step hidden inside the old 0.005 window
# ---------------------------------------------------------------------------

def test_a_step_inside_the_old_window_is_only_found_at_the_new_resolution(make_seq):
    """Clade 2.3.2.1c's shape, with the 5-cluster band narrow enough that
    the old 0.005 stop quit before entering it. This is the whole release
    in one assertion pair: same fake, two resolutions, different answers."""
    steps = [(0.0, 1), (0.9925, 5), (0.993, 7)]

    fine, _ = _run(make_seq, steps, n_target=5)
    assert len(fine.clusters) == 5
    assert fine.threshold >= 0.9925

    coarse, _ = _run(make_seq, steps, n_target=5, resolution=0.005)
    assert len(coarse.clusters) == 7, (
        "the 0.005 stop must miss the 5-cluster band — if it finds it, this "
        "fake no longer discriminates the two resolutions"
    )


def test_the_default_resolution_is_the_fine_one(make_seq):
    """Guards the default itself: a revert to 0.005 fails here even if
    every other test happens to break on an exact landing first."""
    assert THRESHOLD_RESOLUTION == 0.0005
    steps = [(0.0, 1), (0.9925, 5), (0.993, 7)]
    default, _ = _run(make_seq, steps, n_target=5)
    explicit, _ = _run(make_seq, steps, n_target=5,
                       resolution=THRESHOLD_RESOLUTION)
    assert len(default.clusters) == len(explicit.clusters) == 5


def test_the_search_still_terminates_when_no_exact_count_exists(make_seq):
    """No threshold yields 5; the search must stop at the resolution rather
    than run forever."""
    search, _ = _run(make_seq, [(0.0, 1), (0.99, 7)], n_target=5)
    assert len(search.clusters) == 7


def test_the_probe_budget_is_derived_from_window_and_resolution(make_seq):
    """The progress line's denominator must be the real number of probes:
    ceil(log2(0.7 / 0.0005)) = 11 for the mmseqs2 window."""
    seen: list[float] = []

    def counting(sequences, threshold, cfg, tmp_dir=None):
        seen.append(threshold)
        return [Cluster(cluster_id="c0", representative=sequences[0])]

    seqs = _seqs(make_seq, 30)
    with patch("repseq.modes.taxonomic1.run_clustering", side_effect=counting):
        _binary_search_threshold(seqs, n_target=5, cfg=_cfg(), overflow="keep")
    assert len(seen) == 11


# ---------------------------------------------------------------------------
# Bracket choice under overflow="keep"
# ---------------------------------------------------------------------------

def test_keep_fills_a_group_that_steps_past_the_target(make_seq):
    """1 → 7 for a target of 5: the upper bracket is closer (2 vs 4), so the
    group is filled instead of collapsing to a single representative."""
    search, _ = _run(make_seq, [(0.0, 1), (0.99, 7)], n_target=5)
    assert len(search.clusters) == 7


def test_keep_does_not_inflate_when_the_upper_bracket_is_far(make_seq):
    """4 → 200 for a target of 5: the lower bracket is closer, so the group
    keeps 4 rather than returning forty times the request."""
    search, _ = _run(make_seq, [(0.0, 4), (0.99, 200)], n_target=5,
                     n_seqs=250)
    assert len(search.clusters) == 4


def test_keep_breaks_a_tie_upward(make_seq):
    """4 → 6 for a target of 5: equal distance, so over-coverage wins —
    a panel short of a clade is worse than one carrying a spare."""
    search, _ = _run(make_seq, [(0.0, 4), (0.99, 6)], n_target=5)
    assert len(search.clusters) == 6


def test_keep_prefers_an_exact_landing_over_both_brackets(make_seq):
    search, _ = _run(make_seq, [(0.0, 2), (0.6, 5), (0.95, 40)], n_target=5)
    assert len(search.clusters) == 5


def test_an_only_above_bracket_is_taken_however_far_it_is(make_seq):
    """Nothing bounds this case, so it must be reported rather than passed
    off as the closeness rule: the message names the target."""
    search, _ = _run(make_seq, [(0.0, 40)], n_target=5, n_seqs=60)
    assert len(search.clusters) == 40
    assert search.note is not None and "40" in search.note


def test_an_empty_clustering_result_is_not_a_usable_bracket(make_seq):
    """A backend returning no clusters sets a count of 0, which must not
    beat a real bracket above the target — that path rebuilt one singleton
    per input sequence and called them representatives."""
    def sometimes_empty(sequences, threshold, cfg, tmp_dir=None):
        if threshold < 0.9:
            return []
        return [
            Cluster(cluster_id=f"c{i}", representative=sequences[i])
            for i in range(200)
        ]

    seqs = _seqs(make_seq, 250)
    with patch("repseq.modes.taxonomic1.run_clustering",
               side_effect=sometimes_empty):
        search = _binary_search_threshold(
            seqs, n_target=5, cfg=_cfg(), overflow="keep",
        )
    assert len(search.clusters) == 200


# ---------------------------------------------------------------------------
# Bracket choice under overflow="trim"
# ---------------------------------------------------------------------------

def test_trim_reaches_the_exact_target_through_the_upper_bracket(make_seq):
    """The old code only trimmed when NO threshold landed at or below the
    target, so this group returned 1. Rounding up first lets trim deliver
    exactly 5."""
    search, _ = _run(make_seq, [(0.0, 1), (0.99, 7)], n_target=5,
                     overflow="trim")
    assert len(search.clusters) == 5
    assert search.trimmed_from == 7


def test_trim_rounds_up_even_when_the_lower_bracket_is_closer(make_seq):
    """4 → 200 for a target of 5: 'keep' would take the 4, but a count below
    the target can never be grown, so trim must take the 200 and pare it."""
    search, _ = _run(make_seq, [(0.0, 4), (0.99, 200)], n_target=5,
                     overflow="trim", n_seqs=250)
    assert len(search.clusters) == 5


def test_trim_leaves_an_exact_landing_alone(make_seq):
    search, _ = _run(make_seq, [(0.0, 2), (0.6, 5), (0.95, 40)], n_target=5,
                     overflow="trim")
    assert len(search.clusters) == 5
    assert search.trimmed_from is None


def test_trim_keeps_every_input_sequence_inside_some_cluster(make_seq):
    """Trim used to delete the clusters it did not pick, members and all:
    _clusters.tsv then covered a fraction of the stratum while
    _group_counts.tsv still reported the full size."""
    search, seqs = _run(make_seq, [(0.0, 4), (0.99, 200)], n_target=5,
                        overflow="trim", n_seqs=250)
    covered = {id(c.representative) for c in search.clusters}
    for c in search.clusters:
        covered.update(id(m) for m in c.members)
    assert len(search.clusters) == 5
    assert covered == {id(s) for s in seqs}
    assert sum(c.size for c in search.clusters) == 250


def test_trim_honours_the_exact_count_when_representatives_share_an_id(make_seq):
    """Two records sharing an accession both survive QC, so filtering kept
    clusters by `.id` returned more than the target. Identity must decide."""
    seqs = _seqs(make_seq, 6)
    for seq in seqs:
        seq.id = "s0"
    fake = [
        Cluster(cluster_id=f"c{i}", representative=seqs[i]) for i in range(6)
    ]
    with patch("repseq.modes.taxonomic1.run_clustering", return_value=fake):
        search = _binary_search_threshold(
            seqs, n_target=1, cfg=_cfg(), overflow="trim",
        )
    assert len(search.clusters) == 1
    assert search.clusters[0].size == 6


# ---------------------------------------------------------------------------
# Genuinely unreachable targets still undershoot, honestly
# ---------------------------------------------------------------------------

def test_a_group_that_never_splits_keeps_its_single_representative(make_seq):
    """Two clonal H5 clades behaved this way at every identity. Neither fix
    can invent diversity, so the group still returns 1 — and says why."""
    search, _ = _run(make_seq, [(0.0, 1)], n_target=5)
    assert len(search.clusters) == 1


def test_the_undershoot_message_names_the_step(make_seq, capsys):
    _run(make_seq, [(0.0, 1)], n_target=5)
    out = capsys.readouterr().out
    assert "target 5 not reachable" in out
    assert "no threshold split this group past 1" in out


def test_the_overshoot_message_names_the_target(make_seq, capsys):
    """The silent case before v0.80.0: a group handed back more than asked
    for printed a bare count with no mention of the request."""
    _run(make_seq, [(0.0, 1), (0.99, 7)], n_target=5)
    out = capsys.readouterr().out
    assert "settled on 7 rep(s)" in out
    assert "target 5 overshot" in out
    assert "cluster count steps 1 → 7" in out


def test_the_trim_message_reports_the_final_count(make_seq, capsys):
    """The console must show what the user will actually get, not the
    pre-trim bracket."""
    _run(make_seq, [(0.0, 1), (0.99, 7)], n_target=5, overflow="trim")
    out = capsys.readouterr().out
    assert "settled on 5 rep(s)" in out
    assert "trimmed from 7" in out


def test_the_reported_threshold_is_the_one_that_produced_the_clusters(make_seq):
    """A clonal group probed up to 0.9995 used to report 0.6500 — the first
    probe to reach that count — and that value lands in the `cutoff`
    column of _group_counts.tsv."""
    search, _ = _run(make_seq, [(0.0, 1)], n_target=5)
    assert search.threshold > 0.99


def test_notes_travel_with_the_result_for_the_run_log(make_seq):
    under, _ = _run(make_seq, [(0.0, 1)], n_target=5)
    assert under.note is not None and "kept 1 of 5 requested" in under.note

    over, _ = _run(make_seq, [(0.0, 1), (0.99, 7)], n_target=5)
    assert over.note is not None and "kept 7 for 5 requested" in over.note

    exact, _ = _run(make_seq, [(0.0, 5)], n_target=5)
    assert exact.note is None


# ---------------------------------------------------------------------------
# Unchanged contracts
# ---------------------------------------------------------------------------

def test_cluster_members_survive_the_new_bracket_choice(make_seq):
    """Callers report cluster sizes off these objects, so the search must
    keep returning real Clusters with members intact — including on the
    above-bracket path, which is new and which `keep` now takes."""
    search, seqs = _run(make_seq, [(0.0, 1), (0.99, 7)], n_target=5)
    assert len(search.clusters) == 7
    assert all(isinstance(c, Cluster) for c in search.clusters)
    assert sum(c.size for c in search.clusters) == len(seqs)
    assert any(c.members for c in search.clusters)
