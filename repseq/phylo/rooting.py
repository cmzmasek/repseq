"""Tree rooting for the phyloXML output.

Both FastTree and IQ-TREE produce *unrooted* trees by default — the
root the parser picks is arbitrary and rarely biologically meaningful.
A bench scientist staring at a tree where Hantaviridae is split across
the root has lost the main thing the tree is supposed to tell them.

This module roots the tree using a priority chain, first success wins:

1. **Taxonomy-guided** (``_taxonomy_guided_root``): for every branch,
   try rooting there and score the result by mean LCA specificity of
   the internal nodes weighted by clade size. Higher = the rooting
   agrees better with the NCBI taxonomy. Skipped if too few leaves
   carry lineage data (the score becomes meaningless).

2. **MAD** (``_mad_root`` — Minimal Ancestor Deviation, Tria et al.
   2017): for every branch, find the split point that minimises
   Σ ρ(i,j)² across all leaf pairs, where
   ρ(i,j) = (d(i, root) − d(j, root)) / d(i,j). Pure Python, no
   external dependency. Robust when taxonomy is sparse.

3. **Midpoint** (Bio.Phylo's ``root_at_midpoint``): final fallback;
   always succeeds on a tree with at least one branch.

The user can pin a single method via ``phylo.rooting.method`` (``auto``
runs the chain; ``none`` disables rooting entirely — useful when the
input tree is already rooted by an outgroup).
"""

from __future__ import annotations

import copy
import logging
from typing import Optional

from Bio.Phylo.BaseTree import Clade, Tree

from ..models import Sequence

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Lineage extraction
# ---------------------------------------------------------------------------

def _lineage_for(seq: Sequence) -> list[str]:
    """Return the seq's taxonomy lineage as a root→tip list of names.

    Built from ``TaxonomyInfo.lineage`` (the full rank→name map from
    NCBI's ``LineageEx``) when available; falls back to the standard
    fields in a coarse-to-fine order. Returns ``[]`` if there's no
    taxonomy at all — those leaves are excluded from the LCA vote.
    """
    if seq.taxonomy is None:
        return []
    # Prefer the lineage map (NCBI's authoritative root→tip list);
    # values come out in insertion order, which matches NCBI's
    # ordering when the cache was built. Append species/genus/...
    # standard fields only if they aren't already present.
    out: list[str] = []
    seen: set[str] = set()
    for name in seq.taxonomy.lineage.values():
        if name and name not in seen:
            out.append(name)
            seen.add(name)
    for field in (
        seq.taxonomy.superkingdom, seq.taxonomy.kingdom, seq.taxonomy.phylum,
        seq.taxonomy.class_, seq.taxonomy.order, seq.taxonomy.family,
        seq.taxonomy.genus, seq.taxonomy.species,
    ):
        if field and field not in seen:
            out.append(field)
            seen.add(field)
    return out


def _build_short_id_lineage_map(
    reps_by_short_id: dict[str, Sequence],
) -> dict[str, list[str]]:
    """Map short_id → lineage list (root→tip)."""
    return {sid: _lineage_for(seq) for sid, seq in reps_by_short_id.items()}


# ---------------------------------------------------------------------------
# Taxonomy-guided rooting
# ---------------------------------------------------------------------------

def _lca_prefix(lineages: list[list[str]]) -> list[str]:
    """Longest common prefix of a list of root→tip lineages."""
    if not lineages:
        return []
    prefix: list[str] = []
    for tokens in zip(*lineages):
        if all(t == tokens[0] for t in tokens):
            prefix.append(tokens[0])
        else:
            break
    return prefix


def _mean_lca_specificity(
    tree: Tree,
    short_id_lineage: dict[str, list[str]],
    coverage_threshold: float = 0.5,
) -> float:
    """Score a rooting by how well its internals agree with taxonomy.

    For every internal node, take the LCA prefix of its terminals'
    lineages, weighted by terminal count. Internal nodes whose
    terminals lack lineage data above ``coverage_threshold`` are
    skipped. Higher score = the rooting groups taxa the way NCBI
    classifies them.
    """
    total_weight = 0
    total_score = 0
    for node in tree.get_nonterminals():
        terminals = node.get_terminals()
        if not terminals:
            continue
        lineages = [
            short_id_lineage.get(t.name, []) for t in terminals
        ]
        covered = [lg for lg in lineages if lg]
        if not covered or len(covered) / len(lineages) < coverage_threshold:
            continue
        prefix = _lca_prefix(covered)
        total_score += len(prefix) * len(terminals)
        total_weight += len(terminals)
    if total_weight == 0:
        return 0.0
    return total_score / total_weight


def _candidate_clades(tree: Tree) -> list[Clade]:
    """Every clade that can serve as an outgroup target (i.e. every
    non-root clade reachable from the root). Used by both the
    taxonomy-guided and MAD searches."""
    return [c for c in tree.find_clades() if c is not tree.root]


def _taxonomy_guided_root(
    tree: Tree, short_id_lineage: dict[str, list[str]],
) -> Optional[Tree]:
    """Root by maximum mean LCA specificity. Returns the best rooted
    copy of the tree, or None if no candidate has any lineage signal.
    """
    candidates = _candidate_clades(tree)
    if not candidates:
        return None
    if not any(short_id_lineage.values()):
        return None

    best_tree: Optional[Tree] = None
    best_score = -1.0
    best_branch_len = float("inf")
    for cand in candidates:
        attempt = copy.deepcopy(tree)
        # Re-find the matching clade in the deep-copied tree by name +
        # branch length (clades have no stable identity across copies).
        matched = _find_matching_clade(attempt, cand)
        if matched is None:
            continue
        try:
            attempt.root_with_outgroup(matched)
        except Exception:
            continue
        # Re-extract lineages on the rooted copy.
        score = _mean_lca_specificity(attempt, short_id_lineage)
        bl = matched.branch_length or 0.0
        if (score > best_score) or (score == best_score and bl < best_branch_len):
            best_score = score
            best_branch_len = bl
            best_tree = attempt
    if best_score <= 0:
        return None
    return best_tree


def _find_matching_clade(tree: Tree, target: Clade) -> Optional[Clade]:
    """Find ``target``'s equivalent in a deep-copied tree.

    Identity isn't preserved by ``copy.deepcopy``, so we match on
    (terminal-name set, branch_length) which is unique enough for the
    trees the orchestrator produces (every leaf has a distinct short
    id and branch lengths are real numbers from MAFFT/MSA).
    """
    target_terms = frozenset(
        t.name for t in target.get_terminals() if t.name
    )
    target_bl = target.branch_length
    for clade in tree.find_clades():
        terms = frozenset(t.name for t in clade.get_terminals() if t.name)
        if terms == target_terms and clade.branch_length == target_bl:
            return clade
    return None


# ---------------------------------------------------------------------------
# MAD rooting (Minimal Ancestor Deviation)
# ---------------------------------------------------------------------------

def _pairwise_leaf_distances(tree: Tree) -> dict[tuple[str, str], float]:
    """All-pairs leaf distances by BFS through the tree.

    Returns a dict keyed by sorted (leaf_a_name, leaf_b_name) tuples
    for stable lookup. Cost is O(L²) which is fine for the rep counts
    repseq produces (typically dozens to a few hundred).
    """
    terminals = list(tree.get_terminals())
    out: dict[tuple[str, str], float] = {}
    for i, a in enumerate(terminals):
        for b in terminals[i + 1:]:
            d = tree.distance(a, b)
            key = tuple(sorted([a.name, b.name]))
            out[key] = d
    return out


def _leaves_under(clade: Clade) -> set[str]:
    return {t.name for t in clade.get_terminals() if t.name}


def _mad_score_for_branch(
    tree: Tree, branch: Clade,
    pair_d: dict[tuple[str, str], float],
    leaf_dist_to_branch: dict[str, float],
) -> tuple[float, float]:
    """Return (score, split_position) for placing the root on ``branch``.

    For each branch length L, the root sits at position x ∈ [0, L]
    measured from the branch's distal (child) endpoint. For each leaf
    pair (i, j) we compute the deviation
    ρ(i,j) = (d(i,r) − d(j,r)) / d(i,j) and minimise Σ ρ² over x.

    Solved analytically: it's a sum of squares, quadratic in x.
    """
    L = branch.branch_length or 0.0
    if L <= 0:
        return float("inf"), 0.0
    leaves_below = _leaves_under(branch)
    if not leaves_below:
        return float("inf"), 0.0

    # For a given x, d(i,r) = d(i,branch) - x for i below the branch,
    # else d(i,branch) + x (the rest of the tree).
    # ρ(i,j) for i below, j above:
    #   ρ = ((d_i - x) - (d_j + x)) / d_ij = (d_i - d_j - 2x) / d_ij
    # For both below or both above, x drops out — those pairs add a
    # constant Σ((d_i - d_j) / d_ij)².
    leaves = list(leaf_dist_to_branch.keys())
    below = leaves_below
    sq_sum_const = 0.0
    a = 0.0
    b = 0.0
    c = 0.0
    for i_idx, li in enumerate(leaves):
        di = leaf_dist_to_branch[li]
        for lj in leaves[i_idx + 1:]:
            dj = leaf_dist_to_branch[lj]
            key = tuple(sorted([li, lj]))
            d_ij = pair_d.get(key)
            if not d_ij or d_ij == 0:
                continue
            li_below = li in below
            lj_below = lj in below
            if li_below == lj_below:
                # x drops out.
                num = di - dj if li_below else dj - di
                sq_sum_const += (num / d_ij) ** 2
            else:
                # One above, one below; sign depends on which is which.
                if li_below:
                    num0 = di - dj
                else:
                    num0 = dj - di
                # ρ = (num0 - 2x) / d_ij  (after fixing the sign).
                # ρ² = (num0² - 4·num0·x + 4x²) / d_ij²
                c += (num0 / d_ij) ** 2
                b += -4 * num0 / (d_ij ** 2)
                a += 4 / (d_ij ** 2)
    if a == 0:
        return sq_sum_const + c, 0.0
    # Minimum of ax² + bx + c is at x* = -b / (2a); clamp to [0, L].
    x_star = max(0.0, min(L, -b / (2 * a)))
    score = a * x_star ** 2 + b * x_star + c + sq_sum_const
    return score, x_star


def _mad_root(tree: Tree) -> Optional[Tree]:
    """MAD root the tree. Returns a rooted deep copy, or None on
    failure (e.g. tree has no branch lengths)."""
    terminals = list(tree.get_terminals())
    if len(terminals) < 3:
        return None
    pair_d = _pairwise_leaf_distances(tree)
    candidates = _candidate_clades(tree)
    if not candidates:
        return None

    best_tree: Optional[Tree] = None
    best_score = float("inf")
    for branch in candidates:
        if branch.branch_length in (None, 0):
            continue
        # Distance from each leaf to the *child* (distal) endpoint of
        # the branch.
        leaf_dist: dict[str, float] = {}
        for leaf in terminals:
            try:
                leaf_dist[leaf.name] = tree.distance(branch, leaf)
            except Exception:
                continue
        if not leaf_dist:
            continue
        score, _x = _mad_score_for_branch(tree, branch, pair_d, leaf_dist)
        if score < best_score:
            best_score = score
            best_tree = (branch, _x)

    if best_tree is None:
        return None

    branch, _x = best_tree
    rooted = copy.deepcopy(tree)
    matched = _find_matching_clade(rooted, branch)
    if matched is None:
        return None
    try:
        # Bio.Phylo's root_with_outgroup splits the branch at its
        # midpoint by default; passing outgroup_branch_length=x would
        # let us honour the optimal split, but Bio.Phylo's API only
        # accepts that as a length; we use it here.
        rooted.root_with_outgroup(matched, outgroup_branch_length=_x)
    except TypeError:
        # Older Bio.Phylo without the keyword argument — accept the
        # default midpoint split on the chosen branch.
        rooted.root_with_outgroup(matched)
    except Exception:
        return None
    return rooted


# ---------------------------------------------------------------------------
# Top-level entry point
# ---------------------------------------------------------------------------

VALID_METHODS = ("auto", "taxonomy", "mad", "midpoint", "none")


def root_tree(
    tree: Tree,
    reps_by_short_id: dict[str, Sequence],
    method: str = "auto",
    *,
    coverage_threshold: float = 0.5,
) -> tuple[Tree, str]:
    """Root ``tree`` and return ``(rooted_tree, method_used)``.

    ``method`` is one of ``"auto"`` (default — taxonomy → MAD →
    midpoint chain), ``"taxonomy"``, ``"mad"``, ``"midpoint"``, or
    ``"none"`` (return the tree unchanged). When a specific method is
    requested but fails, the function falls through to midpoint rather
    than raising — a usable tree is more valuable than a perfect one.
    """
    if method not in VALID_METHODS:
        method = "auto"
    if method == "none":
        return tree, "none"

    short_id_lineage = _build_short_id_lineage_map(reps_by_short_id)

    if method in ("auto", "taxonomy"):
        rooted = _taxonomy_guided_root(tree, short_id_lineage)
        if rooted is not None:
            return rooted, "taxonomy"
        if method == "taxonomy":
            logger.info(
                "[phylo] taxonomy rooting failed (no lineage signal); "
                "falling back to midpoint"
            )
    if method in ("auto", "mad"):
        rooted = _mad_root(tree)
        if rooted is not None:
            return rooted, "mad"
        if method == "mad":
            logger.info("[phylo] MAD rooting failed; falling back to midpoint")
    # Midpoint is the final fallback — it cannot fail on a tree with
    # any non-zero branches.
    try:
        rooted = copy.deepcopy(tree)
        rooted.root_at_midpoint()
        return rooted, "midpoint"
    except Exception:
        return tree, "none"
