"""Smoke test for the optional clustering scatter plot.

matplotlib is the only hard requirement (the ``[viz]`` extra); the
embedding prefers UMAP when ``umap-learn`` imports cleanly and otherwise
falls back to a numpy-only classical MDS. These tests therefore skip only
when matplotlib is missing — they exercise whichever embedding the
environment supports, plus the MDS path explicitly via monkeypatch.
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("matplotlib")

import matplotlib

matplotlib.use("Agg")  # noqa: E402

from repseq.models import Cluster, RunResult, Sequence, SequenceType, TaxonomyInfo
from repseq.viz import clustering_plot
from repseq.viz.clustering_plot import (
    _classical_mds,
    _embed,
    write_clustering_plot,
)


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


def _force_mds(monkeypatch):
    """Make _umap_status report unusable so _embed must use classical MDS."""
    monkeypatch.setattr(
        clustering_plot, "_umap_status",
        lambda: (False, "umap-learn not installed (test)"),
    )


def test_writes_png_for_small_clustering(tmp_path):
    # Three clusters with distinct k-mer profiles so the embedding separates them.
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


# ---------------------------------------------------------------------------
# Embedding: classical MDS fallback (no umap dependency)
# ---------------------------------------------------------------------------

def test_classical_mds_returns_2d_coords():
    # Build a small symmetric distance matrix and confirm Nx2 output.
    d = np.array(
        [
            [0.0, 0.2, 0.9, 0.95],
            [0.2, 0.0, 0.85, 0.9],
            [0.9, 0.85, 0.0, 0.15],
            [0.95, 0.9, 0.15, 0.0],
        ]
    )
    coords = _classical_mds(d)
    assert coords.shape == (4, 2)
    # Near-neighbours (0,1) should embed closer than the cross-block pair (0,2).
    near = np.linalg.norm(coords[0] - coords[1])
    far = np.linalg.norm(coords[0] - coords[2])
    assert near < far


def test_classical_mds_is_deterministic():
    rng = np.random.default_rng(0)
    base = rng.random((6, 6))
    d = np.abs(base + base.T)
    np.fill_diagonal(d, 0.0)
    a = _classical_mds(d)
    b = _classical_mds(d)
    assert np.allclose(a, b)


def test_embed_falls_back_to_mds_when_umap_unusable(monkeypatch):
    _force_mds(monkeypatch)
    d = np.array([[0.0, 0.3, 0.8], [0.3, 0.0, 0.7], [0.8, 0.7, 0.0]])
    coords, method = _embed(d, n=3, seed=0)
    assert method == "MDS"
    assert coords.shape == (3, 2)


def test_embed_uses_umap_when_available():
    pytest.importorskip("umap")
    d = np.array([[0.0, 0.3, 0.8], [0.3, 0.0, 0.7], [0.8, 0.7, 0.0]])
    coords, method = _embed(d, n=3, seed=0)
    assert method == "UMAP"
    assert coords.shape == (3, 2)


def test_plot_renders_via_mds_fallback(tmp_path, monkeypatch):
    # The whole figure must render with matplotlib alone (no umap).
    _force_mds(monkeypatch)
    a = [_mkseq(f"A{i}", "ACGT" * 25 + "A" * i, "Alphavirus") for i in range(5)]
    b = [_mkseq(f"B{i}", "GGCC" * 25 + "G" * i, "Betavirus") for i in range(4)]
    result = RunResult(
        mode="test",
        representatives=[a[0], b[0]],
        clusters=[
            Cluster(cluster_id="c1", representative=a[0], members=a[1:]),
            Cluster(cluster_id="c2", representative=b[0], members=b[1:]),
        ],
    )
    out = tmp_path / "mds.png"
    written = write_clustering_plot(result, out, seed=0)
    assert written == out
    assert out.stat().st_size > 2000


def test_require_matplotlib_passes_when_installed():
    # matplotlib is importorskip-guaranteed at module load — no raise.
    clustering_plot._require_matplotlib()
