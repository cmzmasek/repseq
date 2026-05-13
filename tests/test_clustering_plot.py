"""Smoke test for the optional clustering UMAP plot.

Skipped automatically when the ``[viz]`` extras aren't installed.
"""

from __future__ import annotations

import pytest

pytest.importorskip("umap")
pytest.importorskip("matplotlib")

import matplotlib

matplotlib.use("Agg")  # noqa: E402

from repseq.models import Cluster, RunResult, Sequence, SequenceType, TaxonomyInfo
from repseq.viz.clustering_plot import write_clustering_plot


def _mkseq(sid: str, sequence: str, genus: str) -> Sequence:
    return Sequence(
        id=sid,
        header=sid,
        sequence=sequence,
        seq_type=SequenceType.NUCLEOTIDE,
        accession=sid,
        organism=f"virus_{genus}",
        taxonomy=TaxonomyInfo(genus=genus),
    )


def test_writes_png_for_small_clustering(tmp_path):
    # Three clusters with distinct k-mer profiles so UMAP separates them.
    a = [_mkseq(f"A{i}", "ACGT" * 25 + "A" * i, "Alphavirus") for i in range(5)]
    b = [_mkseq(f"B{i}", "GGCC" * 25 + "G" * i, "Betavirus") for i in range(4)]
    c = [_mkseq(f"C{i}", "TTAA" * 25 + "T" * i, "Alphavirus") for i in range(3)]

    clusters = [
        Cluster(cluster_id="c1", representative=a[0], members=a[1:]),
        Cluster(cluster_id="c2", representative=b[0], members=b[1:]),
        Cluster(cluster_id="c3", representative=c[0], members=c[1:]),
    ]
    result = RunResult(
        mode="test",
        representatives=[a[0], b[0], c[0]],
        clusters=clusters,
    )

    out = tmp_path / "clustering.png"
    written = write_clustering_plot(result, out, seed=0)
    assert written == out
    assert out.exists()
    assert out.stat().st_size > 2000  # not an empty/error PNG


def test_returns_none_when_no_clusters(tmp_path):
    result = RunResult(mode="test", representatives=[], clusters=[])
    out = tmp_path / "nothing.png"
    assert write_clustering_plot(result, out) is None
    assert not out.exists()


def test_returns_none_when_too_few_points(tmp_path):
    s = _mkseq("X", "ACGT" * 20, "Alphavirus")
    result = RunResult(
        mode="test",
        representatives=[s],
        clusters=[Cluster(cluster_id="c1", representative=s, members=[])],
    )
    out = tmp_path / "tiny.png"
    assert write_clustering_plot(result, out) is None
    assert not out.exists()


def test_subsampling_keeps_all_representatives(tmp_path):
    # Build a result with more points than max_points so subsampling fires.
    # Two reps + many members; both reps must survive.
    rep1 = _mkseq("R1", "ACGT" * 60, "Alphavirus")
    rep2 = _mkseq("R2", "TGCA" * 60, "Betavirus")
    members1 = [
        _mkseq(f"m1_{i}", "ACGT" * 60 + "A" * (i % 9), "Alphavirus")
        for i in range(30)
    ]
    members2 = [
        _mkseq(f"m2_{i}", "TGCA" * 60 + "T" * (i % 9), "Betavirus")
        for i in range(30)
    ]

    result = RunResult(
        mode="test",
        representatives=[rep1, rep2],
        clusters=[
            Cluster(cluster_id="c1", representative=rep1, members=members1),
            Cluster(cluster_id="c2", representative=rep2, members=members2),
        ],
    )

    out = tmp_path / "subsampled.png"
    written = write_clustering_plot(result, out, seed=0, max_points=10)
    assert written == out
    assert out.exists()
