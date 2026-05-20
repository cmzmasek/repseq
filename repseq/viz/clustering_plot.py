"""Two-panel scatter of repseq clustering results.

Rendered from a finished ``RunResult``:

* **Left**  — every clustered sequence embedded on a k-mer Jaccard
  distance, colored by genus (top 10 + "Other").
* **Right** — same coordinates, colored by cluster, with point size scaling
  with cluster size, faint member→representative lines, and an inset
  cluster-size histogram.

The only hard requirement is ``matplotlib`` (``pip install 'repseq[viz]'``),
which ships wheels for every platform and installs without compiling. The
embedding prefers **UMAP** when ``umap-learn`` imports cleanly (the
``[viz-umap]`` extra — nicer separation on large diverse sets) and otherwise
falls back to a **numpy-only classical MDS** (PCoA) on the same distance
matrix. ``umap-learn`` pins narrow numpy/scipy/numba ranges and frequently
fails to import in a shared conda env, so the fallback is what keeps
``--plot`` working everywhere; the figure honestly labels which method ran.

For large runs the embedding is subsampled (representatives always kept)
to keep distance-matrix cost bounded; lines are auto-suppressed past a
density threshold to avoid spaghetti.
"""

from __future__ import annotations

import importlib
import logging
import random
from collections import Counter
from pathlib import Path
from typing import Optional

from ..clustering.diversity import _jaccard_distance, _kmer_set
from ..models import RunResult, Sequence

logger = logging.getLogger(__name__)

DEFAULT_MAX_POINTS = 2000
DEFAULT_MAX_LINES = 500


def _require_matplotlib() -> None:
    """Verify matplotlib — the one hard requirement for ``--plot`` — imports.

    Distinguishes "not installed" from "installed but failing to import".
    The latter is almost always a NumPy/SciPy version clash in the
    environment, not a repseq problem; surface the real error rather than
    telling the user to reinstall a package they already have. Raises
    ``ImportError`` with a plain-English next step on either.
    """
    try:
        importlib.import_module("matplotlib")
    except ModuleNotFoundError as exc:
        if exc.name == "matplotlib":
            raise ImportError(
                "Plotting requires matplotlib — install it with: "
                "pip install 'repseq[viz]'"
            ) from exc
        raise ImportError(
            f"matplotlib is installed but a dependency is missing ({exc}); "
            "the cleanest fix is a dedicated environment for repseq"
        ) from exc
    except ImportError as exc:
        raise ImportError(
            f"matplotlib is installed but failed to import [{exc}] — this is "
            "almost always a NumPy/SciPy version mismatch in the environment, "
            "not a repseq problem; the cleanest fix is a dedicated "
            "environment for repseq"
        ) from exc


def _umap_status() -> tuple[bool, str]:
    """Return ``(usable, detail)`` for the optional UMAP upgrade.

    "Usable" means ``import umap`` actually succeeds. umap-learn pins narrow
    numpy/scipy/numba ranges and a mismatch surfaces as an ``ImportError``
    (or, from numba/llvmlite, occasionally a non-Import exception) at import
    time. Any of those means "fall back to MDS" — we never want a broken
    optional dependency to abort a plot the fallback can produce.
    """
    try:
        importlib.import_module("umap")
        return True, "umap-learn available"
    except ModuleNotFoundError as exc:
        if exc.name == "umap":
            return False, "umap-learn not installed (install 'repseq[viz-umap]')"
        return False, f"umap-learn dependency missing: {exc}"
    except ImportError as exc:
        return False, f"umap-learn failed to import ({exc})"
    except Exception as exc:  # numba/llvmlite can raise non-ImportError here
        return False, f"umap-learn failed to import ({type(exc).__name__}: {exc})"


def _classical_mds(dist) -> "object":
    """Classical MDS (PCoA) into 2-D from a precomputed distance matrix.

    numpy-only — no scipy, no umap. Double-centre the squared-distance
    matrix ``B = -½ J D² J`` (``J`` the centring matrix) and take the two
    largest-eigenvalue eigenvectors, scaling each axis by ``sqrt(λ)``.
    Jaccard distance is not Euclidean, so ``B`` can carry small negative
    eigenvalues; we simply keep the top-two *positive* ones (an axis with a
    non-positive eigenvalue collapses to zero), which is the standard PCoA
    handling for non-Euclidean inputs. Deterministic up to the arbitrary
    sign of each eigenvector — no ``random_state`` lottery.
    """
    import numpy as np

    d = np.asarray(dist, dtype=np.float64)
    n = d.shape[0]
    d2 = d ** 2
    j = np.eye(n) - np.full((n, n), 1.0 / n)
    b = -0.5 * (j @ d2 @ j)
    # B is symmetric → eigh gives real eigenpairs in ascending order.
    vals, vecs = np.linalg.eigh(b)
    order = np.argsort(vals)[::-1][:2]
    coords = np.zeros((n, 2), dtype=np.float64)
    for axis, idx in enumerate(order):
        lam = vals[idx]
        if lam > 0:
            coords[:, axis] = vecs[:, idx] * np.sqrt(lam)
    return coords


def _embed(dist, *, n: int, seed: int) -> tuple["object", str]:
    """Embed the precomputed distance matrix into 2-D.

    Prefers UMAP when ``umap-learn`` imports cleanly (better separation on
    large diverse sets); otherwise — or if the UMAP fit itself raises —
    falls back to numpy-only classical MDS. Returns ``(coords, method)``
    where ``method`` is ``"UMAP"`` or ``"MDS"`` so the figure can label
    honestly which embedding produced it.
    """
    usable, detail = _umap_status()
    if usable:
        try:
            import umap

            n_neighbors = max(2, min(15, n - 1))
            reducer = umap.UMAP(
                metric="precomputed",
                n_neighbors=n_neighbors,
                min_dist=0.1,
                random_state=seed,
                n_jobs=1,
            )
            return reducer.fit_transform(dist), "UMAP"
        except Exception as exc:
            logger.warning(
                "UMAP embedding failed (%s) — falling back to classical MDS.",
                exc,
            )
    else:
        logger.info(
            "%s — using classical MDS (PCoA) for the clustering plot.", detail
        )
    return _classical_mds(dist), "MDS"


def write_clustering_plot(
    result: RunResult,
    path: Path,
    *,
    k: int = 5,
    seed: int = 42,
    max_points: int = DEFAULT_MAX_POINTS,
    max_lines: int = DEFAULT_MAX_LINES,
) -> Optional[Path]:
    """Render the before/after clustering scatter and write a PNG.

    Returns the written path, or ``None`` if the run had nothing to plot
    (no clusters, or fewer than 3 sequences).
    """
    _require_matplotlib()

    if not result.clusters:
        logger.warning("No clusters in result — skipping clustering plot.")
        return None

    seqs: list[Sequence] = []
    cluster_idx: list[int] = []
    is_rep: list[bool] = []
    for cid, cluster in enumerate(result.clusters):
        seqs.append(cluster.representative)
        cluster_idx.append(cid)
        is_rep.append(True)
        for m in cluster.members:
            seqs.append(m)
            cluster_idx.append(cid)
            is_rep.append(False)

    n = len(seqs)
    if n < 3:
        logger.warning("Too few points (%d) for clustering plot — skipping.", n)
        return None

    n_input = n
    if n > max_points:
        rng = random.Random(seed)
        rep_positions = [i for i, r in enumerate(is_rep) if r]
        non_rep_positions = [i for i, r in enumerate(is_rep) if not r]
        keep_extra = max(0, max_points - len(rep_positions))
        sampled = rng.sample(
            non_rep_positions, min(keep_extra, len(non_rep_positions))
        )
        keep = sorted(set(rep_positions) | set(sampled))
        seqs = [seqs[i] for i in keep]
        cluster_idx = [cluster_idx[i] for i in keep]
        is_rep = [is_rep[i] for i in keep]
        n = len(seqs)
        logger.info("Subsampled %d → %d points for plot.", n_input, n)

    import numpy as np

    kmer_sets = [_kmer_set(s.sequence, k) for s in seqs]
    dist = np.zeros((n, n), dtype=np.float32)
    for i in range(n):
        for j in range(i + 1, n):
            d = _jaccard_distance(kmer_sets[i], kmer_sets[j])
            dist[i, j] = dist[j, i] = d

    coords, method = _embed(dist, n=n, seed=seed)
    xs = coords[:, 0]
    ys = coords[:, 1]

    genera = [
        (s.taxonomy.genus if s.taxonomy and s.taxonomy.genus else "Unknown")
        for s in seqs
    ]
    size_per_cid = Counter(cluster_idx)
    cid_to_rep_pos = {
        cid: i for i, (cid, r) in enumerate(zip(cluster_idx, is_rep)) if r
    }

    _render(
        xs=xs, ys=ys,
        cluster_idx=cluster_idx, is_rep=is_rep,
        genera=genera, size_per_cid=size_per_cid,
        cid_to_rep_pos=cid_to_rep_pos,
        n_total_clusters=len(result.clusters),
        n_total_reps=len(result.representatives),
        n_plotted=n,
        n_input=n_input,
        max_lines=max_lines,
        method=method,
        path=path,
    )
    return path


def _render(
    *,
    xs,
    ys,
    cluster_idx,
    is_rep,
    genera,
    size_per_cid,
    cid_to_rep_pos,
    n_total_clusters,
    n_total_reps,
    n_plotted,
    n_input,
    max_lines,
    method: str,
    path: Path,
) -> None:
    import numpy as np
    import matplotlib.pyplot as plt

    xs = np.asarray(xs)
    ys = np.asarray(ys)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # ---- Left: by genus (top 10 + Other) ----
    cmap_g = plt.get_cmap("tab10")
    genus_counts = Counter(genera)
    top_genera = [g for g, _ in genus_counts.most_common(10)]
    genus_color = {g: cmap_g(i % 10) for i, g in enumerate(top_genera)}
    other_color = "#bbbbbb"
    has_other = any(g not in genus_color for g in genera)

    point_colors = [genus_color.get(g, other_color) for g in genera]
    ax1.scatter(
        xs, ys, c=point_colors, s=46, alpha=0.75,
        edgecolor="white", linewidth=0.5,
    )
    for g in top_genera:
        ax1.scatter([], [], c=[genus_color[g]], label=g, s=46)
    if has_other:
        ax1.scatter([], [], c=other_color, label="Other", s=46)
    ax1.legend(
        loc="best", fontsize=8, title="Genus",
        title_fontsize=9, frameon=True,
    )
    sub_left = (
        f"({n_plotted} of {n_input} shown — subsampled)"
        if n_plotted < n_input
        else f"({n_plotted} sequences)"
    )
    ax1.set_title(
        f"Before clustering — colored by genus  {sub_left}",
        fontsize=12, pad=10,
    )
    ax1.set_xlabel(f"{method} 1"); ax1.set_ylabel(f"{method} 2")
    ax1.set_xticks([]); ax1.set_yticks([])
    for s in ax1.spines.values():
        s.set_edgecolor("#cccccc")

    # ---- Right: clusters ----
    cmap_c = plt.get_cmap("tab10")

    def cluster_color(cid: int):
        return cmap_c(cid % 10)

    n_non_rep = sum(1 for r in is_rep if not r)
    draw_lines = n_non_rep <= max_lines

    if draw_lines:
        for i, cid in enumerate(cluster_idx):
            if is_rep[i]:
                continue
            r = cid_to_rep_pos.get(cid)
            if r is None:
                continue
            ax2.plot(
                [xs[i], xs[r]], [ys[i], ys[r]],
                color=cluster_color(cid), alpha=0.25,
                linewidth=0.6, zorder=1,
            )

    def point_size(cid: int) -> float:
        return 28.0 + 12.0 * float(size_per_cid[cid]) ** 0.5

    sizes = [point_size(cid) for cid in cluster_idx]
    member_colors = [cluster_color(cid) for cid in cluster_idx]
    ax2.scatter(
        xs, ys, c=member_colors, s=sizes, alpha=0.70,
        edgecolor="white", linewidth=0.5, zorder=2,
    )

    rep_pos = [cid_to_rep_pos[cid] for cid in sorted(cid_to_rep_pos)]
    ax2.scatter(
        xs[rep_pos], ys[rep_pos],
        facecolor="none", edgecolor="black", s=240,
        linewidth=1.8, marker="o", zorder=3, label="representative",
    )

    subtitle_parts = ["point size ∝ √(cluster size)"]
    if draw_lines:
        subtitle_parts.append("lines = cluster membership")
    ax2.set_title(
        f"After clustering — {n_total_reps} reps / {n_total_clusters} clusters\n"
        + " · ".join(subtitle_parts),
        fontsize=12, pad=10,
    )
    ax2.set_xlabel(f"{method} 1"); ax2.set_ylabel(f"{method} 2")
    ax2.set_xticks([]); ax2.set_yticks([])
    for s in ax2.spines.values():
        s.set_edgecolor("#cccccc")
    ax2.legend(loc="upper right", fontsize=9, frameon=True)

    # Inset cluster-size histogram
    inset = ax2.inset_axes([0.04, 0.04, 0.30, 0.26])
    sizes_arr = np.array(list(size_per_cid.values()))
    max_s = int(sizes_arr.max())
    if max_s <= 1:
        bins = [1, 2]
    else:
        bins = np.unique(
            np.geomspace(1, max_s + 1, num=8).astype(int)
        )
        if len(bins) < 2:
            bins = np.array([1, max_s + 1])
    inset.hist(sizes_arr, bins=bins, color="#444444",
               edgecolor="white", linewidth=0.6)
    if max_s > 1:
        inset.set_xscale("log")
    inset.tick_params(axis="both", labelsize=7)
    inset.set_title("cluster size distribution", fontsize=8, pad=2)
    inset.set_facecolor("#fafafa")
    for s in inset.spines.values():
        s.set_edgecolor("#cccccc")

    embed_label = (
        "UMAP" if method == "UMAP"
        else "classical MDS (umap-learn unavailable)"
    )
    fig.suptitle(
        f"repseq clustering — sequence embedding ({embed_label}, k-mer Jaccard)",
        fontsize=13, y=1.02,
    )
    plt.tight_layout()

    path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(path, dpi=130, bbox_inches="tight", facecolor="white")
    plt.close(fig)
