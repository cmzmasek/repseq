"""Mode smoke tests — verify dispatch and grouping logic without invoking MMseqs2."""
from __future__ import annotations

from unittest.mock import patch

from repseq.models import Cluster, RunResult, SequenceSource, SequenceType, TaxonomyInfo
from repseq.modes.global_mode import GlobalMode
from repseq.modes.taxonomic1 import TaxonomicMode1


def _base_cfg() -> dict:
    return {
        "seed": 42,
        "threads": 1,
        "clustering": {
            "backend": "mmseqs2",
            "mmseqs2_mode": "easy-linclust",
            "coverage": 0.8,
            "coverage_mode": 0,
            "extra_args": [],
        },
        "representative": {"priority": ["refseq", "reviewed_uniprot", "longest"]},
    }


# ---------------------------------------------------------------------------
# GlobalMode :: count sub-mode (uses select_diverse — no mmseqs needed)
# ---------------------------------------------------------------------------

def test_global_mode_count_returns_n_diverse(make_seq):
    seqs = [
        make_seq("a", "AAAAACCCCCGGGGGTTTTT"),
        make_seq("b", "AAAAACCCCCGGGGGTTTTC"),
        make_seq("c", "TTTTTGGGGGCCCCCAAAAA"),
    ]
    mode = GlobalMode(_base_cfg(), n_select=2)
    result = mode.run(seqs)
    assert result.mode == "global:count"
    assert len(result.representatives) == 2
    assert len(result.clusters) == 2


# ---------------------------------------------------------------------------
# GlobalMode :: threshold sub-mode (clustering mocked)
# ---------------------------------------------------------------------------

def test_global_mode_threshold_uses_mmseqs2(make_seq):
    cfg = _base_cfg()
    a = make_seq("rep_a", "A" * 50, is_refseq=True)
    b = make_seq("mem_b", "A" * 40)
    fake_clusters = [Cluster(cluster_id="c1", representative=a, members=[b])]

    with patch("repseq.modes.global_mode.run_clustering", return_value=fake_clusters) as m:
        result = GlobalMode(cfg, threshold=0.9).run([a, b])

    m.assert_called_once()
    assert result.mode == "global:threshold"
    assert [r.id for r in result.representatives] == ["rep_a"]


# ---------------------------------------------------------------------------
# TaxonomicMode1 :: groups by rank, skips clustering when group <= n_per_group
# ---------------------------------------------------------------------------

def _with_tax(seq, **ranks):
    seq.taxonomy = TaxonomyInfo(**ranks)
    return seq


def test_taxonomic1_groups_small_groups_skip_clustering(make_seq):
    """If a group already has <= n_per_group sequences, mode keeps them all
    without invoking the clustering backend."""
    s1 = _with_tax(make_seq("s1", "ACGT" * 20), genus="Foo")
    s2 = _with_tax(make_seq("s2", "ACGT" * 20), genus="Foo")
    s3 = _with_tax(make_seq("s3", "TGCA" * 20), genus="Bar")

    with patch("repseq.modes.taxonomic1.run_clustering") as m:
        result = TaxonomicMode1(_base_cfg(), rank="genus", n_per_group=5).run([s1, s2, s3])

    m.assert_not_called()
    assert {r.id for r in result.representatives} == {"s1", "s2", "s3"}


def test_taxonomic1_large_group_triggers_clustering(make_seq):
    """When a group exceeds n_per_group, clustering should be invoked
    (via the binary-search helper) until the count fits."""
    seqs = [_with_tax(make_seq(f"s{i}", "ACGT" * 20), genus="Foo") for i in range(10)]
    cfg = _base_cfg()

    fake_clusters = [
        Cluster(cluster_id=f"c{i}", representative=seqs[i]) for i in range(2)
    ]

    with patch("repseq.modes.taxonomic1.run_clustering", return_value=fake_clusters):
        result = TaxonomicMode1(cfg, rank="genus", n_per_group=2).run(seqs)

    assert len(result.representatives) == 2


def test_taxonomic1_records_group_stats_for_small_groups(make_seq):
    """Groups kept whole record clustered=False and no cutoff."""
    s1 = _with_tax(make_seq("s1", "ACGT" * 20), genus="Foo")
    s2 = _with_tax(make_seq("s2", "ACGT" * 20), genus="Foo")
    s3 = _with_tax(make_seq("s3", "TGCA" * 20), genus="Bar")

    with patch("repseq.modes.taxonomic1.run_clustering"):
        result = TaxonomicMode1(_base_cfg(), rank="genus", n_per_group=5).run([s1, s2, s3])

    stats = {gs.group: gs for gs in result.group_stats}
    assert set(stats) == {"Foo", "Bar"}
    assert all(gs.grouping == "genus" for gs in result.group_stats)
    assert (stats["Foo"].n_before, stats["Foo"].n_after) == (2, 2)
    assert stats["Foo"].clustered is False and stats["Foo"].cutoff is None


def test_taxonomic1_group_stats_record_cutoff_when_clustered(make_seq):
    """A clustered group records the binary-search threshold as its cutoff."""
    seqs = [_with_tax(make_seq(f"s{i}", "ACGT" * 20), genus="Foo") for i in range(10)]
    fake_clusters = [
        Cluster(cluster_id=f"c{i}", representative=seqs[i]) for i in range(2)
    ]

    with patch("repseq.modes.taxonomic1.run_clustering", return_value=fake_clusters):
        result = TaxonomicMode1(_base_cfg(), rank="genus", n_per_group=2).run(seqs)

    assert len(result.group_stats) == 1
    gs = result.group_stats[0]
    assert (gs.grouping, gs.group) == ("genus", "Foo")
    assert (gs.n_before, gs.n_after) == (10, 2)
    assert gs.clustered is True
    assert gs.cutoff is not None and 0.0 < gs.cutoff <= 1.0


def test_taxonomic1_preserves_cluster_members(make_seq):
    """Regression: binary-search modes must report real cluster sizes, not
    collapse every cluster to a singleton. _binary_search_threshold used to
    return only representatives, so the mode rebuilt member-less clusters and
    _clusters.tsv showed cluster_size=1 everywhere."""
    seqs = [_with_tax(make_seq(f"s{i}", "ACGT" * 20), genus="Foo") for i in range(10)]
    fake = [
        Cluster(cluster_id="c0", representative=seqs[0], members=list(seqs[1:5])),
        Cluster(cluster_id="c1", representative=seqs[5], members=list(seqs[6:10])),
    ]
    with patch("repseq.modes.taxonomic1.run_clustering", return_value=fake):
        result = TaxonomicMode1(_base_cfg(), rank="genus", n_per_group=2).run(seqs)

    assert len(result.clusters) == 2
    assert sorted(c.size for c in result.clusters) == [5, 5]
    assert all(c.members for c in result.clusters)
    # cluster_id relabeled with the group + representative id.
    assert all(c.cluster_id.startswith("Foo|") for c in result.clusters)


def test_binary_search_returns_clusters_with_members(make_seq):
    """The helper returns Cluster objects (members intact), not bare reps."""
    from repseq.modes.taxonomic1 import _binary_search_threshold

    seqs = [make_seq(f"s{i}", "ACGT" * 10) for i in range(10)]
    fake = [
        Cluster(cluster_id="c0", representative=seqs[0], members=list(seqs[1:6])),
        Cluster(cluster_id="c1", representative=seqs[6], members=list(seqs[7:10])),
    ]
    with patch("repseq.modes.taxonomic1.run_clustering", return_value=fake):
        clusters, _ = _binary_search_threshold(
            seqs, n_target=2, cfg=_base_cfg(), overflow="keep"
        )

    assert all(isinstance(c, Cluster) for c in clusters)
    assert sorted(c.size for c in clusters) == [4, 6]


def test_taxonomic2_accumulates_transitive_members(make_seq):
    """Hierarchical mode: a leaf representative stands for everything merged
    into it across levels, so leaf cluster_size must reflect the transitive
    member count, not 1."""
    from repseq.modes.taxonomic2 import TaxonomicMode2

    # 10 seqs, all genus=Foo, all species=Bar. Level 1 (genus) clusters
    # 10→2; level 2 (species) is a small group that passes through, so the
    # two level-1 reps reach the leaf carrying their members.
    seqs = [
        _with_tax(make_seq(f"s{i}", "ACGT" * 20), genus="Foo", species="Bar")
        for i in range(10)
    ]
    fake = [
        Cluster(cluster_id="c0", representative=seqs[0], members=list(seqs[1:5])),
        Cluster(cluster_id="c1", representative=seqs[5], members=list(seqs[6:10])),
    ]
    rank_levels = [
        {"rank": "genus", "n_per_group": 2},
        {"rank": "species", "n_per_group": 5},
    ]
    # _binary_search_threshold lives in taxonomic1 and resolves run_clustering
    # in that module's namespace, so patch it there (not in taxonomic2).
    with patch("repseq.modes.taxonomic1.run_clustering", return_value=fake):
        result = TaxonomicMode2(_base_cfg(), rank_levels=rank_levels).run(seqs)

    assert len(result.clusters) == 2
    assert sorted(c.size for c in result.clusters) == [5, 5]


def test_global_mode_records_group_stats(make_seq):
    seqs = [
        make_seq("a", "AAAAACCCCCGGGGGTTTTT"),
        make_seq("b", "AAAAACCCCCGGGGGTTTTC"),
        make_seq("c", "TTTTTGGGGGCCCCCAAAAA"),
    ]
    result = GlobalMode(_base_cfg(), n_select=2).run(seqs)
    assert len(result.group_stats) == 1
    gs = result.group_stats[0]
    assert (gs.grouping, gs.group) == ("global", "(all)")
    assert (gs.n_before, gs.n_after) == (3, 2)
    # Diversity selection, not threshold clustering — no cutoff.
    assert gs.clustered is False and gs.cutoff is None


def test_taxonomic1_missing_rank_grouped_as_unknown(make_seq):
    s1 = _with_tax(make_seq("s1", "ACGT" * 20), genus="Foo")
    s2 = make_seq("s2", "ACGT" * 20)  # no taxonomy at all
    with patch("repseq.modes.taxonomic1.run_clustering"):
        result = TaxonomicMode1(_base_cfg(), rank="genus", n_per_group=5).run([s1, s2])
    # Both kept because both groups are small
    assert {r.id for r in result.representatives} == {"s1", "s2"}


# ---------------------------------------------------------------------------
# Binary-search threshold direction (regression: the search used to walk
# AWAY from the target because the MMseqs2 threshold→count relationship
# was inverted).
# ---------------------------------------------------------------------------

def test_binary_search_climbs_toward_target(make_seq):
    from repseq.modes.taxonomic1 import _binary_search_threshold

    seqs = [make_seq(f"s{i}", "ACGT" * 10) for i in range(20)]

    def fake_run_clustering(sequences, threshold, cfg, tmp_dir=None):
        # Correct MMseqs2 semantics: a higher identity threshold yields
        # MORE clusters. count ~ threshold * 20.
        n = max(1, min(len(sequences), round(threshold * 20)))
        return [
            Cluster(cluster_id=f"c{i}", representative=sequences[i])
            for i in range(n)
        ]

    with patch("repseq.modes.taxonomic1.run_clustering", side_effect=fake_run_clustering):
        reps, threshold = _binary_search_threshold(
            seqs, n_target=10, cfg=_base_cfg(), overflow="keep"
        )

    # count ≈ threshold*20, so the target of 10 sits near threshold 0.5.
    # The search must land close to 10 — not undershoot to 1-3.
    assert len(reps) <= 10
    assert len(reps) >= 8


def test_binary_search_respects_cdhit_floor(make_seq):
    # When the backend is cd-hit-est, the wrapper refuses identity < 0.80.
    # The binary search must clamp ``lo`` to that floor so it never
    # presents the backend with a value it would reject. We give the
    # cluster mock a function that depends on threshold so we can see
    # the lo-clamp in the (lo, hi) endpoints.
    from repseq.modes.taxonomic1 import _binary_search_threshold
    from repseq.models import SequenceType

    seqs = [make_seq(f"s{i}", "ACGT" * 10, seq_type=SequenceType.NUCLEOTIDE)
            for i in range(20)]
    cfg = _base_cfg()
    cfg["clustering"]["backend"] = "cdhit"

    seen_thresholds: list[float] = []

    def fake_run_clustering(sequences, threshold, cfg_, tmp_dir=None):
        seen_thresholds.append(threshold)
        # Always produce enough clusters to force the search to lower
        # the threshold; the floor should stop it at 0.80.
        return [
            Cluster(cluster_id=f"c{i}", representative=sequences[i])
            for i in range(len(sequences))
        ]

    with patch("repseq.modes.taxonomic1.run_clustering", side_effect=fake_run_clustering):
        _binary_search_threshold(seqs, n_target=2, cfg=cfg, overflow="keep")

    # Mid is always >= lo, and lo is clamped to 0.80 for cd-hit-est.
    assert seen_thresholds, "binary search did not iterate"
    assert min(seen_thresholds) >= 0.80


def test_binary_search_trim_enforces_exact_count(make_seq):
    from repseq.modes.taxonomic1 import _binary_search_threshold

    seqs = [make_seq(f"s{i}", "ACGT" * 10 + "A" * i) for i in range(20)]

    def fake_run_clustering(sequences, threshold, cfg, tmp_dir=None):
        # Always returns more than the target so 'trim' has to act.
        return [
            Cluster(cluster_id=f"c{i}", representative=sequences[i])
            for i in range(len(sequences))
        ]

    with patch("repseq.modes.taxonomic1.run_clustering", side_effect=fake_run_clustering):
        reps, _ = _binary_search_threshold(
            seqs, n_target=5, cfg=_base_cfg(), overflow="trim"
        )

    assert len(reps) == 5


def test_time_mode_rejects_nonpositive_or_garbage_window():
    import pytest

    from repseq.modes.time_mode import TimeMode

    for bad in ("0", "-3", "banana"):
        with pytest.raises(ValueError):
            TimeMode(_base_cfg(), n_per_window=5, window=bad)

    # valid windows construct fine
    TimeMode(_base_cfg(), n_per_window=5, window="year")
    TimeMode(_base_cfg(), n_per_window=5, window="decade")
    TimeMode(_base_cfg(), n_per_window=5, window="5")
