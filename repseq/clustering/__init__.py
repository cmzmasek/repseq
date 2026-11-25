"""Clustering backends.

``run_clustering`` is a thin dispatcher that selects the backend from
``cfg['clustering']['backend']`` and forwards to the matching wrapper.

Modes should import ``run_clustering`` from this package — not from
``mmseqs2`` directly — so swapping backends is a config change, not a
code change. Tests that need to mock clustering should patch the symbol
in the *mode* module's namespace (e.g. ``repseq.modes.global_mode.
run_clustering``), since Python re-exports the dispatcher by name there.
"""

from __future__ import annotations

import time
from typing import Any, Optional

import click

from ..models import Cluster, Sequence
from . import cdhit, mmseqs2


def run_clustering(
    sequences: list[Sequence],
    threshold: float,
    cfg: dict[str, Any],
    tmp_dir: Optional[str] = None,
) -> list[Cluster]:
    """Dispatch to the clustering backend named in ``cfg['clustering']['backend']``."""
    backend = cfg.get("clustering", {}).get("backend", "mmseqs2")
    if backend == "mmseqs2":
        clusters = mmseqs2.run_clustering(sequences, threshold, cfg, tmp_dir=tmp_dir)
    elif backend == "cdhit":
        clusters = cdhit.run_clustering(sequences, threshold, cfg, tmp_dir=tmp_dir)
    else:
        raise ValueError(f"Unknown clustering backend: {backend!r}")
    return canonical_order(clusters)


def _order_key(seq: Sequence) -> str:
    return str(seq.id or "")


def canonical_order(clusters: list[Cluster]) -> list[Cluster]:
    """Put clusters, and the members inside each, in a fixed order.

    A multi-threaded backend does not return its clusters in a stable
    order: clustering one 215-sequence file six times with
    ``mmseqs easy-cluster --threads 10`` gave the same clusters every time
    in FOUR different orders (2026-11-25). Everything downstream inherits
    that order — the representative tables, every FASTA, the tie-break
    when two members score equally in representative election, the input
    order of the diversity sampler under ``--overflow trim`` — so two
    identical runs wrote files that differed line by line.

    Sorted by sequence id: clusters by their representative's, members by
    their own; cluster ids are then renumbered ``cluster_000001`` ... in
    that order. In place and stable; the ``Sequence`` objects themselves are
    untouched (``force_select`` relies on their identity), and each
    cluster keeps the backend's representative, so who wins a tie between a
    representative and an equally good member does not change.
    """
    for cluster in clusters:
        cluster.members.sort(key=_order_key)
    clusters.sort(key=lambda c: _order_key(c.representative))
    # The backends number their clusters in the order they happened to
    # return them, so the ids moved between runs too (`global` mode writes
    # them to _clusters.tsv as they are). Renumber in the fixed order, same
    # format.
    for i, cluster in enumerate(clusters):
        cluster.cluster_id = f"cluster_{i + 1:06d}"
    return clusters


def min_threshold(cfg: dict[str, Any], sequences: list[Sequence]) -> float:
    """Lowest identity threshold the active backend will accept.

    The binary-search wrappers in modes use this as a stopping floor to
    avoid asking cd-hit for thresholds it would refuse (e.g. below 0.80
    for cd-hit-est).
    """
    backend = cfg.get("clustering", {}).get("backend", "mmseqs2")
    if backend == "cdhit":
        return cdhit.min_threshold(sequences, cfg)
    return 0.0


def compute_diversity_curve(
    sequences: list[Sequence],
    cfg: dict[str, Any],
    label: Optional[str] = None,
) -> Optional[dict[float, Optional[int]]]:
    """Cluster ``sequences`` at each configured "standard" identity threshold
    and return ``{cutoff: n_clusters}``. Reporting-only — does not influence
    representative selection.

    Cutoffs below the active backend's identity floor (``min_threshold``)
    map to ``None`` and the backend is not invoked for them. The TSV
    writer renders ``None`` as ``NA``.

    Returns ``None`` (not an empty dict) when the feature is disabled —
    either because ``clustering.diversity_curve_cutoffs`` is missing /
    empty, or because the caller passed no sequences. ``None`` lets the
    TSV writer distinguish "feature off" from "feature on, all cells
    below floor".

    Configured cutoffs are de-duplicated and sorted high → low before use,
    so a cutoff listed twice costs one clustering pass (not two) and the
    console ordering matches the ``_group_counts.tsv`` column ordering,
    which ``output/report.py`` sorts the same way.

    Progress is echoed to the console in the same shape as the modes'
    binary search (``[label] `` tag, per-pass timing). This step runs one
    full clustering pass per cutoff *after* selection has already settled,
    so on a large group it can take minutes — without the echo it reads as
    a hang. ``label`` is the group name the caller is working on (None for
    a single ungrouped pass). Thresholds print with ``:g`` — the form the
    user wrote in YAML and the form the TSV column names use
    (``n_clusters_0.99``) — so the two surfaces cross-reference directly.
    (The binary search prints ``:.4f`` instead because its thresholds are
    continuous search results, not config constants.)
    """
    raw_cutoffs = cfg.get("clustering", {}).get("diversity_curve_cutoffs", []) or []
    if not raw_cutoffs or not sequences:
        return None
    # De-dup + sort descending. Without the de-dup a cutoff repeated in the
    # config would run the backend twice yet still collapse to ONE key in
    # `out` (and one TSV column), so the announced cutoff count would
    # overstate what the run actually delivers.
    cutoffs = sorted({float(c) for c in raw_cutoffs}, reverse=True)
    if not cutoffs:
        return None
    floor = min_threshold(cfg, sequences)
    tag = f"[{label}] " if label else ""
    click.echo(
        f"    {tag}diversity curve (report only): {len(sequences):,} sequence(s) "
        f"at {len(cutoffs)} standard cutoff(s) ..."
    )
    out: dict[float, Optional[int]] = {}
    for i, c in enumerate(cutoffs):
        if c < floor:
            out[c] = None
            click.echo(
                f"      {tag}cutoff {i + 1}/{len(cutoffs)}: threshold={c:g} "
                f"→ skipped (below the backend's {floor:g} floor)"
            )
            continue
        t0 = time.perf_counter()
        clusters = run_clustering(sequences, c, cfg)
        out[c] = len(clusters)
        click.echo(
            f"      {tag}cutoff {i + 1}/{len(cutoffs)}: threshold={c:g} "
            f"→ {len(clusters):,} cluster(s) [{time.perf_counter() - t0:.1f}s]"
        )
    return out


__all__ = [
    "run_clustering", "min_threshold", "compute_diversity_curve",
    "canonical_order",
    "mmseqs2", "cdhit",
]
