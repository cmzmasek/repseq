"""Pairwise topological incongruence between trees (unrooted Robinson-Foulds).

The per-protein step (2F) builds one tree per HMM domain-architecture
marker. Different markers can carry different evolutionary histories —
reassortment in segmented viruses, recombination in non-segmented ones —
and that shows up as *topological incongruence* between the marker trees.
This module turns "squint at the Spike tree next to the N tree" into a
number: the **Robinson-Foulds (RF) distance** between every pair of trees.

Design choices:

* **Unrooted RF** (bipartition symmetric difference). Marker trees often
  root differently; counting that as incongruence would be misleading, so
  the root is ignored. ``rf = 0`` iff the two trees have the same
  unrooted topology.
* **Common taxa only.** Marker trees have different leaf sets (a
  representative may lack an architecture). Each pair is scored on the
  intersection of its taxa — equivalent to pruning both trees to the
  shared leaves before computing RF.
* **Dependency-free.** RF is computed directly from Bio.Phylo clade
  membership; no dendropy / ete3. Leaves are relabelled from each tree's
  local short ids (``S0001…``) back to the representative id via that
  tree's ``_tree_id_map.tsv`` so taxa are comparable across trees.

``norm_rf`` is ``rf / (2·(n−3))`` for ``n`` common taxa ≥ 4 (the maximum
RF between two binary unrooted trees on ``n`` leaves); with fewer than 4
common taxa every unrooted topology is identical, so it is reported as
``NA``.
"""

from __future__ import annotations

import itertools
import logging
from pathlib import Path
from typing import Optional

from Bio import Phylo

logger = logging.getLogger(__name__)


def _read_id_map(path: Path) -> dict[str, str]:
    """Parse a ``short_id<TAB>accession`` map (header row skipped)."""
    mapping: dict[str, str] = {}
    with open(path) as fh:
        for i, line in enumerate(fh):
            parts = line.rstrip("\n").split("\t")
            if i == 0 and parts[:1] == ["short_id"]:
                continue
            if len(parts) >= 2 and parts[0]:
                mapping[parts[0]] = parts[1]
    return mapping


def load_tree_clusters(
    nwk_path: Path, id_map_path: Path
) -> Optional[tuple[frozenset[str], list[frozenset[str]]]]:
    """Load a Newick tree and return ``(taxa, clusters)`` over rep ids.

    ``taxa`` is the set of leaf ids (relabelled from short ids via the
    id-map). ``clusters`` is the descendant-leaf set of every internal
    node with at least two leaves — the raw material for unrooted
    bipartitions (the root's all-taxa cluster and singletons are
    dropped). Returns None when the tree can't be parsed or has < 2
    leaves (nothing to compare).
    """
    try:
        tree = Phylo.read(str(nwk_path), "newick")
    except Exception as exc:  # malformed Newick → caller skips this tree
        logger.warning("[incongruence] could not parse %s: %s", nwk_path, exc)
        return None

    id_map = _read_id_map(id_map_path) if id_map_path.exists() else {}

    def rep_id(short: Optional[str]) -> Optional[str]:
        if short is None:
            return None
        return id_map.get(short, short)

    taxa = {rep_id(t.name) for t in tree.get_terminals() if t.name}
    taxa.discard(None)
    if len(taxa) < 2:
        return None

    clusters: list[frozenset[str]] = []
    for clade in tree.get_nonterminals():
        members = {rep_id(t.name) for t in clade.get_terminals() if t.name}
        members.discard(None)
        # Drop the whole-tree cluster and singletons; everything in
        # between is a candidate bipartition, filtered per-pair on the
        # common-taxa restriction below.
        if 2 <= len(members) <= len(taxa) - 1:
            clusters.append(frozenset(members))
    return frozenset(taxa), clusters


def rf_distance(
    taxa_a: frozenset[str],
    clusters_a: list[frozenset[str]],
    taxa_b: frozenset[str],
    clusters_b: list[frozenset[str]],
) -> tuple[int, Optional[float], int]:
    """Unrooted RF on the common taxa → ``(rf, norm_rf, n_common)``.

    ``norm_rf`` is None when fewer than 4 taxa are shared (RF is then
    trivially 0 and normalisation is undefined).
    """
    common = taxa_a & taxa_b
    n = len(common)
    if n == 0:
        return 0, None, 0
    pivot = min(common)

    def induced(clusters: list[frozenset[str]]) -> set[frozenset[str]]:
        out: set[frozenset[str]] = set()
        for x in clusters:
            side = x & common
            # Non-trivial on the restricted taxon set: both parts ≥ 2.
            if not (2 <= len(side) <= n - 2):
                continue
            # Canonicalise: the side NOT containing the pivot, so a
            # bipartition and its mirror collapse to one key.
            key = side if pivot not in side else (common - side)
            out.add(frozenset(key))
        return out

    rf = len(induced(clusters_a) ^ induced(clusters_b))
    if n >= 4:
        max_rf = 2 * (n - 3)
        norm: Optional[float] = (rf / max_rf) if max_rf > 0 else 0.0
    else:
        norm = None
    return rf, norm, n


def compute_incongruence(
    trees: list[tuple[str, Path, Path]],
) -> list[dict]:
    """Pairwise RF over a list of ``(label, nwk_path, id_map_path)``.

    Trees that fail to load are dropped. Returns one row dict per pair of
    successfully-loaded trees, in input order (so genome comparisons,
    appended last by the caller, sort after the marker-vs-marker pairs).
    """
    loaded: list[tuple[str, frozenset[str], list[frozenset[str]]]] = []
    for label, nwk_path, id_map_path in trees:
        result = load_tree_clusters(nwk_path, id_map_path)
        if result is None:
            continue
        taxa, clusters = result
        loaded.append((label, taxa, clusters))

    rows: list[dict] = []
    for (la, ta, ca), (lb, tb, cb) in itertools.combinations(loaded, 2):
        rf, norm, n = rf_distance(ta, ca, tb, cb)
        rows.append(
            {
                "tree_a": la,
                "tree_b": lb,
                "rf": rf,
                "norm_rf": norm,
                "n_common_taxa": n,
            }
        )
    return rows


def write_incongruence_tsv(rows: list[dict], path: Path) -> None:
    """Write the pairwise RF table. ``norm_rf`` renders as ``NA`` when None."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write("tree_a\ttree_b\trf\tnorm_rf\tn_common_taxa\n")
        for r in rows:
            norm = r["norm_rf"]
            norm_str = "NA" if norm is None else f"{norm:.4f}"
            fh.write(
                f"{r['tree_a']}\t{r['tree_b']}\t{r['rf']}\t"
                f"{norm_str}\t{r['n_common_taxa']}\n"
            )
