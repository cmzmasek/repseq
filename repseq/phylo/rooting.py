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
import sys
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

    # For a given x measured from the branch's child (distal) endpoint:
    #   i below the branch:  d(i,r) = d(i,branch) + x   (ascend i→child→r)
    #   j above the branch:  d(j,r) = d(j,branch) - x   (d(j,branch) already
    #     includes the full branch length L because the path j→child crosses
    #     it, so descending the remaining L-x from the parent nets to -x)
    # ρ(i,j) for i below, j above:
    #   ρ = ((d_i + x) - (d_j - x)) / d_ij = (d_i - d_j + 2x) / d_ij
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
                # num0 = d_below - d_above (sign fixed above), so
                # ρ = (num0 + 2x) / d_ij
                # ρ² = (num0² + 4·num0·x + 4x²) / d_ij²
                c += (num0 / d_ij) ** 2
                b += 4 * num0 / (d_ij ** 2)
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
# Outgroup rooting (user-specified)
# ---------------------------------------------------------------------------

def _resolve_outgroup_short_ids(
    outgroup_spec,
    outgroup_rank,
    reps_by_short_id: dict[str, Sequence],
) -> set[str]:
    """Resolve the user's outgroup spec to a set of short_ids on the tree.

    ``outgroup_spec`` is a single accession string, a list of accession
    strings, or ``None``; matched (case-insensitively) against each rep's
    ``accession``, ``id``, or ``isolate_id``. ``outgroup_rank`` is an
    optional ``{rank: taxon}`` dict; every rep whose taxonomy carries that
    taxon at that rank is added. The two sources are unioned — a user can
    name one explicit accession AND a clade by rank.
    """
    result: set[str] = set()

    accs: list[str] = []
    if isinstance(outgroup_spec, str):
        accs = [outgroup_spec]
    elif isinstance(outgroup_spec, list):
        accs = [a for a in outgroup_spec if isinstance(a, str) and a.strip()]
    if accs:
        wanted = {a.strip().lower() for a in accs}
        for sid, rep in reps_by_short_id.items():
            candidates = (rep.accession, rep.id, rep.isolate_id)
            if any(c and c.lower() in wanted for c in candidates):
                result.add(sid)

    if isinstance(outgroup_rank, dict):
        for rank, taxon in outgroup_rank.items():
            if not isinstance(rank, str) or not isinstance(taxon, str):
                continue
            tlower = taxon.strip().lower()
            if not tlower:
                continue
            for sid, rep in reps_by_short_id.items():
                if rep.taxonomy is None:
                    continue
                v = rep.taxonomy.get_rank(rank)
                if v and v.lower() == tlower:
                    result.add(sid)

    return result


def _best_outgroup_branch(tree: Tree, wanted: set) -> tuple[Clade, int]:
    """Pick the branch that best separates the outgroup from the rest.

    The tree a builder hands us is rooted arbitrarily (IQ-TREE roots on
    its first sequence), so "the MRCA of the outgroup" is meaningless
    as parsed: if the outgroup holds that first sequence its MRCA is the
    root and re-rooting on it does nothing. A branch, by contrast, is a
    property of the UNROOTED tree. Every branch splits the leaves in two;
    we take the branch whose outgroup-bearing side holds every outgroup
    leaf and the fewest others. Returns ``(clade_below_branch, extra)``
    where ``extra`` is the number of non-outgroup leaves on the outgroup
    side (0 when the outgroup is a clean group).
    """
    n_total = 0
    n_out_total = 0
    counts: dict[int, tuple[int, int]] = {}
    for clade in tree.find_clades(order="postorder"):
        if clade.is_terminal():
            n, k = 1, (1 if clade.name in wanted else 0)
            n_total += 1
            n_out_total += k
        else:
            n = sum(counts[id(c)][0] for c in clade.clades)
            k = sum(counts[id(c)][1] for c in clade.clades)
        counts[id(clade)] = (n, k)
    best: Optional[Clade] = None
    best_extra = n_total + 1
    for clade in tree.find_clades():
        if clade is tree.root:
            continue
        n, k = counts[id(clade)]
        if k == n_out_total:          # outgroup entirely below this branch
            extra = n - k
        elif k == 0:                  # outgroup entirely above this branch
            extra = (n_total - n) - n_out_total
        else:
            continue
        if extra < best_extra:
            best, best_extra = clade, extra
    return best, best_extra


def _root_on_branch(tree: Tree, clade: Clade) -> None:
    """Place the root at the middle of the branch above ``clade``.

    Passing ``outgroup_branch_length`` matters: without it Bio.Phylo makes
    an internal ``clade`` itself the new root, which yields a three-way
    root with the outgroup drawn split across two root branches.
    """
    root = tree.root
    if len(root.clades) == 2 and clade in root.clades:
        return  # already rooted on this very branch
    half = (clade.branch_length or 0.0) / 2.0
    tree.root_with_outgroup(clade, outgroup_branch_length=half)
    # Bio.Phylo leaves the FULL original length on the far side when the
    # branch hung directly off the old root, lengthening the tree by half
    # a branch. Both halves are known, so set them outright.
    if len(tree.root.clades) == 2:
        for child in tree.root.clades:
            child.branch_length = half


def _outgroup_root(
    tree: Tree,
    reps_by_short_id: dict[str, Sequence],
    outgroup_spec,
    outgroup_rank,
) -> Optional[Tree]:
    """Root by a user-specified outgroup (single accession or a clade).

    Returns the rooted copy on success, ``None`` if the spec resolves to
    no leaves on the tree, or to ALL of them (the caller then falls back
    to midpoint). When several leaves resolve, the root goes on the
    branch separating them from every other leaf (see
    ``_best_outgroup_branch``), giving an ordinary two-way root. If the
    outgroup leaves are not a single group on this tree, the root goes on
    the branch enclosing them with the fewest other leaves and a stderr
    note says so.
    """
    wanted = _resolve_outgroup_short_ids(
        outgroup_spec, outgroup_rank, reps_by_short_id,
    )
    if not wanted:
        return None
    rooted = copy.deepcopy(tree)
    terminals = rooted.get_terminals()
    targets = [t for t in terminals if t.name in wanted]
    if not targets:
        return None
    if len(targets) == len(terminals):
        print(
            "[phylo] outgroup rooting: the configured outgroup matches "
            "EVERY leaf of this tree, so there is nothing to root "
            "against; falling back to midpoint.",
            file=sys.stderr,
        )
        return None
    try:
        if len(targets) == 1:
            rooted.root_with_outgroup(targets[0])
        else:
            branch, extra = _best_outgroup_branch(
                rooted, {t.name for t in targets},
            )
            if branch is None:
                return None
            _root_on_branch(rooted, branch)
            if extra:
                print(
                    f"[phylo] outgroup rooting: the {len(targets)} outgroup "
                    "sequences do not form a single group on this tree. "
                    "Rooted on the smallest group that holds all of them, "
                    f"which also holds {extra} other sequence(s) — check "
                    "the tree before reading its root.",
                    file=sys.stderr,
                )
    except Exception as exc:
        logger.info("[phylo] outgroup rooting raised: %s", exc)
        return None
    return rooted


def _split_key(leaves: frozenset, all_leaves: frozenset) -> frozenset:
    """One canonical name for a branch: the smaller side of its split
    (ties broken by sorted names), so a branch has the same key however
    the tree is rooted."""
    other = all_leaves - leaves
    return min(leaves, other, key=lambda s: (len(s), sorted(s)))


def _restore_supports(original: Tree, rooted: Tree) -> None:
    """Put every support value back on the branch it was estimated for.

    Bio.Phylo stores a support value on the NODE below a branch. When a
    tree is re-rooted, the branches between the old and the new root
    reverse direction, but the values stay on their nodes — so each one
    ends up describing the neighbouring branch. Supports belong to
    branches (splits of the leaf set), which re-rooting cannot change,
    so we read split → support off the tree as built and write it back.
    Both root children of a two-way root sit on the same branch and get
    the same value.
    """
    all_leaves = frozenset(t.name for t in original.get_terminals())
    below: dict[int, frozenset] = {}
    support: dict[frozenset, object] = {}
    for clade in original.find_clades(order="postorder"):
        if clade.is_terminal():
            below[id(clade)] = frozenset([clade.name])
            continue
        leaves = frozenset().union(*(below[id(c)] for c in clade.clades))
        below[id(clade)] = leaves
        if clade is original.root or clade.confidence is None:
            continue
        support.setdefault(_split_key(leaves, all_leaves), clade.confidence)
    below = {}
    for clade in rooted.find_clades(order="postorder"):
        if clade.is_terminal():
            below[id(clade)] = frozenset([clade.name])
            continue
        leaves = frozenset().union(*(below[id(c)] for c in clade.clades))
        below[id(clade)] = leaves
        if clade is rooted.root:
            clade.confidence = None
            continue
        clade.confidence = support.get(_split_key(leaves, all_leaves))


# ---------------------------------------------------------------------------
# Top-level entry point
# ---------------------------------------------------------------------------

VALID_METHODS = ("auto", "taxonomy", "mad", "midpoint", "outgroup", "none")


def root_tree(
    tree: Tree,
    reps_by_short_id: dict[str, Sequence],
    method: str = "auto",
    *,
    coverage_threshold: float = 0.5,
    outgroup=None,
    outgroup_rank=None,
) -> tuple[Tree, str]:
    """Root ``tree`` and return ``(rooted_tree, method_used)``.

    ``method`` is one of ``"auto"`` (default — taxonomy → MAD →
    midpoint chain), ``"taxonomy"``, ``"mad"``, ``"midpoint"``,
    ``"outgroup"`` (user-specified, needs ``outgroup`` and/or
    ``outgroup_rank``), or ``"none"`` (return the tree unchanged). When
    a specific method is requested but fails, the function falls
    through to midpoint rather than raising — a usable tree is more
    valuable than a perfect one.
    """
    rooted, used = _root_tree(
        tree, reps_by_short_id, method,
        coverage_threshold=coverage_threshold,
        outgroup=outgroup, outgroup_rank=outgroup_rank,
    )
    if rooted is not tree:
        # Re-rooting leaves support values one branch off along the path
        # between the old and new root; put them back (see the helper).
        try:
            _restore_supports(tree, rooted)
        except Exception as exc:  # never lose a tree over its supports
            print(
                f"[phylo] WARNING: could not re-check branch supports "
                f"after rooting ({exc}); support values near the root may "
                "be misplaced.",
                file=sys.stderr,
            )
    return rooted, used


def _root_tree(
    tree: Tree,
    reps_by_short_id: dict[str, Sequence],
    method: str = "auto",
    *,
    coverage_threshold: float = 0.5,
    outgroup=None,
    outgroup_rank=None,
) -> tuple[Tree, str]:
    """The rooting chain itself; ``root_tree`` wraps it."""
    if method not in VALID_METHODS:
        method = "auto"
    if method == "none":
        return tree, "none"

    short_id_lineage = _build_short_id_lineage_map(reps_by_short_id)

    if method == "outgroup":
        rooted = _outgroup_root(tree, reps_by_short_id, outgroup, outgroup_rank)
        if rooted is not None:
            return rooted, "outgroup"
        logger.info(
            "[phylo] outgroup rooting failed: no representative matched "
            "the configured outgroup (`phylo.rooting.outgroup` / "
            "`outgroup_rank`); falling back to midpoint"
        )

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
