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


def test_taxonomic1_missing_rank_grouped_as_unknown(make_seq):
    s1 = _with_tax(make_seq("s1", "ACGT" * 20), genus="Foo")
    s2 = make_seq("s2", "ACGT" * 20)  # no taxonomy at all
    with patch("repseq.modes.taxonomic1.run_clustering"):
        result = TaxonomicMode1(_base_cfg(), rank="genus", n_per_group=5).run([s1, s2])
    # Both kept because both groups are small
    assert {r.id for r in result.representatives} == {"s1", "s2"}
