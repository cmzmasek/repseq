"""Taxonomic Lineage Aware Mode 1: N representatives per taxonomic group.

For each group at the specified rank, binary-search for the clustering
threshold whose cluster count lands on ``n_per_group`` — or, where the step
function offers no exact landing, on the closest attainable count, which may
sit above the target (see :func:`_binary_search_threshold` for which side
wins). A group therefore yields *about* ``n_per_group`` representatives, not
at most ``n_per_group``.
"""

from __future__ import annotations

import logging
import math
import time
from typing import Any, Literal, NamedTuple, Optional

import click

from ..clustering import compute_diversity_curve, min_threshold, run_clustering
from ..clustering.diversity import nearest_by_kmer, select_diverse
from ..models import Cluster, GroupStat, RunResult, Sequence
from ..representative.selector import apply_representative_selection
from .base import BaseMode

logger = logging.getLogger(__name__)

OverflowStrategy = Literal["keep", "trim"]

#: Identity-unit width the per-group search bisects down to before giving up
#: on an exact landing. Fine (0.0005) because a whole step — 1 → 7 clusters —
#: can sit inside a 0.005-wide window near 99 % identity, which is how a
#: 27-isolate H5N1 clade came back with one representative. THE single
#: source for this number: ``output/summary.py`` imports it for the Methods
#: prose rather than restating the literal.
THRESHOLD_RESOLUTION = 0.0005


class ThresholdSearch(NamedTuple):
    """What one per-group threshold search produced.

    ``trimmed_from`` and ``note`` exist so a caller's :class:`GroupStat` can
    carry the search's own account of itself into
    ``{prefix}_group_counts.tsv`` and the run log, rather than leaving it in
    console output that scrolls away.
    """

    clusters: list[Cluster]
    threshold: float
    trimmed_from: Optional[int] = None
    note: Optional[str] = None


def _group_by_rank(sequences: list[Sequence], rank: str) -> dict[str, list[Sequence]]:
    groups: dict[str, list[Sequence]] = {}
    for seq in sequences:
        if seq.taxonomy:
            label = seq.taxonomy.get_rank(rank)
        else:
            label = None
        key = label or "Unknown"
        groups.setdefault(key, []).append(seq)
    return groups


def _binary_search_threshold(
    sequences: list[Sequence],
    n_target: int,
    cfg: dict[str, Any],
    overflow: OverflowStrategy,
    lo: float = 0.3,
    hi: float = 1.0,
    label: str = "",
    resolution: float = THRESHOLD_RESOLUTION,
) -> ThresholdSearch:
    """Binary-search for the clustering threshold that yields n_target clusters.

    Identity-threshold semantics (same direction for both supported
    backends): a *higher* threshold makes clustering stricter and produces
    *more*, smaller clusters; a *lower* threshold merges more aggressively
    into *fewer* clusters. The search moves accordingly.

    Cluster count is a step function of the threshold, so an exact landing
    on ``n_target`` is not always possible. The search bisects until the
    window is narrower than ``resolution`` (see
    :data:`THRESHOLD_RESOLUTION`), keeping the largest count it saw
    at-or-below the target and the smallest it saw above it. The result is
    then:

    * an exact landing, when some threshold gave exactly ``n_target``;
    * otherwise, when **both** brackets exist, under ``overflow="trim"``
      the one above — which :func:`select_diverse` pares to exactly
      ``n_target``, a count below the target being impossible to grow — and
      under ``overflow="keep"`` whichever is closer, ties going up. So
      1 → 7 for a target of 5 yields 7, while 4 → 200 yields 4;
    * otherwise the only bracket there is, **whatever its distance from the
      target**: a group whose count never reaches ``n_target`` keeps fewer
      representatives (a clonal clade that never splits returns 1), and one
      whose count never falls to ``n_target`` keeps more — possibly far
      more, since nothing bounds how large the smallest attainable count
      can be. Both cases are echoed, logged as warnings and recorded in the
      returned ``note``, because neither is something to discover later
      from a file.

    The bracket rules changed in v0.80.0: before it a count at-or-below the
    target always won, so a group stepping 1 → 7 kept the 1 and asking for
    5 per influenza HA clade returned 1 of clade 2.3.2.1c's 27 isolates.
    Representative counts therefore shift upward for every grouped mode on
    datasets with tight groups.

    The lower bound is clamped to ``min_threshold(cfg, sequences)`` so the
    search never asks the backend for a threshold it would refuse — cd-hit
    rejects ``-c`` below 0.40 (protein) or 0.80 (nucleotide); mmseqs2's
    floor is 0.0, so this is a no-op there.

    Returns a :class:`ThresholdSearch`. Its ``clusters`` are full ``Cluster``
    objects with members intact (representative via ``c.representative``),
    so callers report accurate cluster sizes instead of treating every
    representative as a singleton, and every input sequence stays inside
    some returned cluster even when ``trim`` dropped the cluster it first
    landed in.
    """
    floor = min_threshold(cfg, sequences)
    if floor > lo:
        lo = floor
    if lo >= hi:
        raise ValueError(
            f"Binary search cannot run: clustering backend floor ({floor}) "
            f"is at or above the search upper bound ({hi}). Lower the "
            f"backend's threshold floor or switch backends."
        )

    # The two candidates the bisection brackets the target with: the
    # largest count at-or-below it and the smallest count above it. Both
    # are kept, because which one wins depends on `overflow` (see the
    # docstring) and a count below the target can never be grown. `None`
    # means "never set" — not a count, so no real n_target can collide
    # with it and an empty clustering result cannot pose as a bracket.
    below_clusters: list[Cluster] = []
    below_count: Optional[int] = None
    below_threshold = hi
    above_clusters: list[Cluster] = []
    above_count: Optional[int] = None
    above_threshold = hi

    # The window halves exactly once per probe, so the budget needed to
    # narrow it below `resolution` follows from the window itself. Derived
    # rather than capped at a guess, so the progress line's denominator is
    # the real number of probes a user can plan their wait against.
    span = hi - lo
    n_probes = (
        max(1, math.ceil(math.log2(span / resolution)))
        if span > resolution
        else 1
    )

    tag = f"[{label}] " if label else ""
    click.echo(
        f"  {tag}clustering {len(sequences)} sequences "
        f"(target = {n_target} reps) ..."
    )

    for i in range(n_probes):
        mid = (lo + hi) / 2
        t0 = time.perf_counter()
        clusters = run_clustering(sequences, mid, cfg)
        clusters = apply_representative_selection(clusters, cfg)
        n_reps = len(clusters)
        dt = time.perf_counter() - t0
        click.echo(
            f"    {tag}iter {i + 1}/{n_probes}: threshold={mid:.4f} "
            f"→ {n_reps} cluster(s) [{dt:.1f}s]"
        )

        # `>=` / `<=`, not `>` / `<`: a later probe repeating a count has
        # walked closer to the step's edge, and its threshold is the one
        # worth reporting as the identity this group was clustered at. On
        # `>` the reported cutoff stayed at the FIRST probe to reach that
        # count — 0.65 for a clade that is clonal up to 0.9995 — and no
        # longer matched the clusters handed back.
        if n_reps <= n_target:
            # Few enough — keep the count closest to the target from below,
            # then raise the threshold to split further.
            if below_count is None or n_reps >= below_count:
                below_clusters = clusters
                below_count = n_reps
                below_threshold = mid
            lo = mid
        else:
            # Too many clusters — keep the closest from above, then lower
            # the threshold to merge more.
            if above_count is None or n_reps <= above_count:
                above_clusters = clusters
                above_count = n_reps
                above_threshold = mid
            hi = mid

        if n_reps == n_target or (hi - lo) < resolution:
            break

    # --- pick between the two brackets -----------------------------------
    # A bracket counts as available only if it actually holds clusters: a
    # backend that returns an empty result sets a count of 0, and an empty
    # "below" bracket must not beat a usable one above (mmseqs2 warns about
    # exactly this round trip in its own wrapper).
    have_below = below_count is not None and bool(below_clusters)
    have_above = above_count is not None and bool(above_clusters)

    if have_below and below_count == n_target:
        # Exact landing.
        best_clusters, best_threshold = below_clusters, below_threshold
    elif not have_above:
        # The target was never exceeded at any threshold.
        best_clusters, best_threshold = below_clusters, below_threshold
    elif not have_below:
        # Every threshold overshot: only the above bracket exists.
        best_clusters, best_threshold = above_clusters, above_threshold
    elif overflow == "trim" or (
        (above_count - n_target) <= (n_target - below_count)
    ):
        best_clusters, best_threshold = above_clusters, above_threshold
    else:
        best_clusters, best_threshold = below_clusters, below_threshold

    if not best_clusters:
        # Nothing clustered at all — every sequence is its own singleton.
        best_clusters = [
            Cluster(cluster_id=f"cluster_{i + 1:06d}", representative=s)
            for i, s in enumerate(sequences)
        ]
        best_threshold = hi

    # If overflow == "trim", diversity-select exactly n_target clusters by
    # picking the most diverse representatives and keeping their clusters.
    trimmed_from: Optional[int] = None
    if overflow == "trim" and len(best_clusters) > n_target:
        seed = cfg.get("seed", 42)
        reps = [c.representative for c in best_clusters]
        chosen = select_diverse(reps, n_target, seed=seed, cfg=cfg)
        # Object identity, not `.id`: two records sharing an accession but
        # differing by a base both survive QC (dedup is exact-sequence
        # based, and segmented runs skip segment-level dedup entirely), so
        # an id set would keep a cluster the sampler did not choose and
        # quietly break trim's "exactly n_target" promise.
        chosen_ids = {id(s) for s in chosen}
        kept = [c for c in best_clusters if id(c.representative) in chosen_ids]
        dropped = [
            c for c in best_clusters if id(c.representative) not in chosen_ids
        ]
        # Fold each dropped cluster into the kept cluster whose
        # representative is closest, so the stratum stays fully accounted
        # for: without this, trim silently deletes the members of every
        # cluster it did not pick — on a 250-isolate group pared to 5 that
        # is 240-odd sequences gone from _clusters.tsv, its cluster_size
        # and purity columns, and (in taxonomic2) from the next rank level.
        if kept and dropped:
            nearest = nearest_by_kmer(
                [c.representative for c in dropped],
                [c.representative for c in kept],
                cfg,
            )
            for cluster, target_i in zip(dropped, nearest):
                kept[target_i].members.append(cluster.representative)
                kept[target_i].members.extend(cluster.members)
        trimmed_from = len(best_clusters)
        best_clusters = kept

    # Report the FINAL count, after any trim — the number the user will see
    # in the output — rather than the pre-trim bracket.
    final_count = len(best_clusters)
    where = f" for '{label}'" if label else ""
    trim_note = f" (trimmed from {trimmed_from})" if trimmed_from else ""
    if have_below and have_above:
        detail = (
            f"cluster count steps {below_count} → {above_count} between "
            f"adjacent identities"
        )
    elif have_below:
        detail = f"no threshold split this group past {below_count}"
    elif have_above:
        detail = f"no threshold merged this group below {above_count}"
    else:
        detail = "clustering produced no usable result at any threshold"

    note: Optional[str] = None
    if final_count < n_target:
        note = (
            f"kept {final_count} of {n_target} requested — {detail}"
            f"{trim_note}"
        )
        click.echo(
            f"    {tag}settled on {final_count} rep(s) at threshold "
            f"{best_threshold:.4f}{trim_note} (target {n_target} not "
            f"reachable — {detail})"
        )
        logger.warning(
            "Binary search%s: requested %d representatives but cluster count "
            "is a step function of identity; closest achievable is %d at "
            "threshold %.3f (%s).",
            where, n_target, final_count, best_threshold, detail,
        )
    elif final_count > n_target:
        # The other side of the same coin, and the one that used to pass in
        # silence: nothing bounds how far above the target the nearest
        # attainable count sits, so say so instead of printing a bare count.
        note = (
            f"kept {final_count} for {n_target} requested — {detail}"
            f"{trim_note}"
        )
        click.echo(
            f"    {tag}settled on {final_count} rep(s) at threshold "
            f"{best_threshold:.4f}{trim_note} (target {n_target} overshot "
            f"— {detail})"
        )
        logger.warning(
            "Binary search%s: requested %d representatives but the smallest "
            "attainable cluster count above the target is %d at threshold "
            "%.3f (%s); returning it. Use --overflow trim to pare a group "
            "like this to exactly %d.",
            where, n_target, final_count, best_threshold, detail, n_target,
        )
    else:
        if trimmed_from:
            note = f"trimmed from {trimmed_from} to {final_count}"
        click.echo(
            f"    {tag}settled on {final_count} rep(s) at threshold "
            f"{best_threshold:.4f}{trim_note}"
        )
        logger.info(
            "Binary search%s: %d representatives at threshold %.3f.",
            where, final_count, best_threshold,
        )

    return ThresholdSearch(best_clusters, best_threshold, trimmed_from, note)


class TaxonomicMode1(BaseMode):
    def __init__(
        self,
        cfg: dict[str, Any],
        rank: str,
        n_per_group: int,
        overflow: OverflowStrategy = "keep",
    ) -> None:
        super().__init__(cfg)
        self.rank = rank
        self.n_per_group = n_per_group
        self.overflow = overflow

    def run(self, sequences: list[Sequence]) -> RunResult:
        groups = _group_by_rank(sequences, self.rank)
        all_reps: list[Sequence] = []
        all_clusters: list[Cluster] = []
        group_stats: list[GroupStat] = []

        for group_label, group_seqs in groups.items():
            if len(group_seqs) <= self.n_per_group:
                # Keep all — no clustering needed
                for seq in group_seqs:
                    all_clusters.append(
                        Cluster(
                            cluster_id=f"{group_label}|{seq.id}",
                            representative=seq,
                        )
                    )
                all_reps.extend(group_seqs)
                group_stats.append(GroupStat(
                    grouping=self.rank, group=group_label,
                    n_before=len(group_seqs), n_after=len(group_seqs),
                    clustered=False,
                ))
            else:
                search = _binary_search_threshold(
                    group_seqs,
                    self.n_per_group,
                    self.cfg,
                    self.overflow,
                    label=group_label,
                )
                clusters = search.clusters
                for c in clusters:
                    c.cluster_id = f"{group_label}|{c.representative.id}"
                all_clusters.extend(clusters)
                all_reps.extend(c.representative for c in clusters)
                group_stats.append(GroupStat(
                    grouping=self.rank, group=group_label,
                    n_before=len(group_seqs), n_after=len(clusters),
                    clustered=True, cutoff=search.threshold,
                    trimmed_from=search.trimmed_from, note=search.note,
                    cutoff_counts=compute_diversity_curve(
                        group_seqs, self.cfg, label=group_label),
                ))

        return RunResult(
            mode="taxonomic1",
            representatives=all_reps,
            clusters=all_clusters,
            group_stats=group_stats,
            config_snapshot={
                "rank": self.rank,
                "n_per_group": self.n_per_group,
                "overflow": self.overflow,
            },
        )
