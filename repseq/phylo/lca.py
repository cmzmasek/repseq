"""Internal-node LCA annotation for phylogenetic trees.

After rooting (:mod:`repseq.phylo.rooting`), every internal node is
labelled with the taxonomic LCA of its descendants. The LCA's rank
(family, genus, species, …) is read from the resolver's lineage map
when available; falls back to NCBI's virus-ICTV suffix convention
(``…viridae`` = family, ``…virales`` = order, etc.) when the rank
isn't on the taxonomy object.

Three behaviours worth pinning:

1. **Coverage gate.** An internal node is annotated only if
   ≥ ``coverage_threshold`` of its terminals (default 0.5) carry
   lineage data. A handful of well-annotated leaves inside a mostly
   bare clade would otherwise dictate a too-specific label.

2. **min_rank filter.** Leaves whose lineage doesn't reach
   ``min_rank`` (default ``"genus"``) are *excluded from the LCA
   vote* — they remain on the tree, but their thin lineage doesn't
   blur internal labels. ``"none"`` keeps every leaf in the vote.

3. **keep_deepest_labels cleanup.** The same taxon often appears at
   multiple nested internal nodes (e.g. every clade inside
   Hantaviridae has Hantaviridae as its LCA). The cleanup walks
   internals largest-first and clears every nested duplicate label.

A separate :func:`suppress_same_species_pairs` removes the LCA name
from internals whose only children are two leaves of the same
species — the species name is already on the leaves, so duplicating
it on the parent is just noise.
"""

from __future__ import annotations

import logging
from typing import Optional

from Bio.Phylo.BaseTree import Clade, Tree

from ..models import Sequence

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Rank ordering (used for the min_rank gate)
# ---------------------------------------------------------------------------

# Coarse → fine. Anything not in this list compares as "more specific
# than the most-specific known rank", i.e. never gates a leaf out.
_RANK_ORDER = [
    "superkingdom", "realm", "kingdom", "subkingdom",
    "phylum", "subphylum", "class", "subclass",
    "order", "suborder", "family", "subfamily",
    "genus", "subgenus", "species",
]


def _rank_index(rank: str) -> int:
    rank = rank.lower()
    if rank in _RANK_ORDER:
        return _RANK_ORDER.index(rank)
    # Unknown ranks (e.g. "no rank", "clade") never gate anything out
    # — return a value past species so any min_rank comparison passes.
    return len(_RANK_ORDER)


# ICTV suffix → rank, used when the lineage map doesn't carry the
# rank for a given taxon. Order matters: longer suffixes first so
# "viricetes" wins over "virae".
_ICTV_SUFFIXES: list[tuple[str, str]] = [
    ("viricotina", "subphylum"),
    ("viricetes", "class"),
    ("viricota", "phylum"),
    ("virinae", "subfamily"),
    ("viridae", "family"),
    ("virales", "order"),
    ("virae", "kingdom"),
    ("viria", "realm"),
    ("vira", "subrealm"),
]


def _infer_rank_from_name(name: str) -> str:
    """Return the most likely rank for ``name`` based on ICTV suffix.

    Falls back to ``"genus"`` for single-word ``…virus`` names (the
    ICTV convention) and ``"no rank"`` otherwise.
    """
    if not name:
        return "no rank"
    lower = name.lower()
    for suffix, rank in _ICTV_SUFFIXES:
        if lower.endswith(suffix):
            return rank
    if " " not in name and lower.endswith("virus"):
        return "genus"
    return "no rank"


# ---------------------------------------------------------------------------
# Lineage extraction (with ranks)
# ---------------------------------------------------------------------------

def _ranked_lineage_for(seq: Sequence) -> list[tuple[str, str]]:
    """Return the seq's lineage as a root→tip list of (rank, name).

    Built from ``TaxonomyInfo.lineage`` (rank → name from NCBI's
    LineageEx) plus the standard fields. Used both to gate leaves by
    ``min_rank`` and to compute the LCA prefix.
    """
    if seq.taxonomy is None:
        return []
    out: list[tuple[str, str]] = []
    seen_names: set[str] = set()
    for rank, name in seq.taxonomy.lineage.items():
        if name and name not in seen_names:
            out.append((rank, name))
            seen_names.add(name)
    # Append standard fields if not already present.
    for rank, name in (
        ("superkingdom", seq.taxonomy.superkingdom),
        ("kingdom", seq.taxonomy.kingdom),
        ("phylum", seq.taxonomy.phylum),
        ("class", seq.taxonomy.class_),
        ("order", seq.taxonomy.order),
        ("family", seq.taxonomy.family),
        ("genus", seq.taxonomy.genus),
        ("species", seq.taxonomy.species),
    ):
        if name and name not in seen_names:
            out.append((rank, name))
            seen_names.add(name)
    return out


def _reaches_min_rank(lineage: list[tuple[str, str]], min_rank: str) -> bool:
    """True if the lineage extends to ``min_rank`` or finer.

    ``min_rank="none"`` always returns True (the gate is disabled).
    """
    if min_rank == "none":
        return True
    min_idx = _rank_index(min_rank)
    return any(_rank_index(r) >= min_idx for r, _ in lineage)


# ---------------------------------------------------------------------------
# LCA computation
# ---------------------------------------------------------------------------

def _lca_prefix(
    lineages: list[list[tuple[str, str]]],
) -> list[tuple[str, str]]:
    """Longest common (rank, name) prefix across a list of lineages."""
    if not lineages:
        return []
    prefix: list[tuple[str, str]] = []
    for entries in zip(*lineages):
        first = entries[0]
        if all(e == first for e in entries):
            prefix.append(first)
        else:
            break
    return prefix


def _resolve_rank(name: str, lineage: list[tuple[str, str]]) -> str:
    """Return the rank for ``name`` — from the lineage if present,
    else inferred from the ICTV suffix."""
    for rank, taxon in lineage:
        if taxon == name and rank:
            return rank
    return _infer_rank_from_name(name)


# ---------------------------------------------------------------------------
# Tree annotation
# ---------------------------------------------------------------------------

def annotate_internal_nodes(
    tree: Tree,
    reps_by_short_id: dict[str, Sequence],
    *,
    min_rank: str = "genus",
    coverage_threshold: float = 0.5,
) -> None:
    """Set ``_lca_name`` / ``_lca_rank`` on each internal clade in place.

    Internals whose terminals fail the coverage gate are left
    unlabelled (no attributes set). The writer reads these attributes
    to emit ``<name>`` and ``<taxonomy>`` on internal nodes; absence
    of the attributes means "no label".

    Always run after rooting — an unrooted tree has no biologically
    meaningful internal-node clade membership.
    """
    short_id_lineage = {
        sid: _ranked_lineage_for(seq) for sid, seq in reps_by_short_id.items()
    }
    # Filter by min_rank: leaves with short lineages stay in the
    # tree but don't vote on internals' LCA.
    voting_lineages = {
        sid: lg for sid, lg in short_id_lineage.items()
        if lg and _reaches_min_rank(lg, min_rank)
    }

    for node in tree.get_nonterminals():
        terminals = node.get_terminals()
        if not terminals:
            continue
        leaf_names = [t.name for t in terminals if t.name]
        votes = [voting_lineages[n] for n in leaf_names if n in voting_lineages]
        # Coverage measured against the *total* terminal count, not
        # just voting leaves — so a clade of 100 leaves with only 10
        # well-annotated entries falls below the gate.
        if not terminals or len(votes) / len(terminals) < coverage_threshold:
            continue
        prefix = _lca_prefix(votes)
        if not prefix:
            continue
        lca_rank, lca_name = prefix[-1]
        node._lca_name = lca_name
        node._lca_rank = lca_rank or _infer_rank_from_name(lca_name)


# ---------------------------------------------------------------------------
# Label cleanup: keep deepest only
# ---------------------------------------------------------------------------

def keep_deepest_labels(tree: Tree) -> None:
    """Keep each LCA label on the crown of every maximal same-name run.

    A labelled internal node's name is cleared only when an **ancestor**
    already carries the same name — i.e. it is a nested repeat inside one
    clade. Two *disjoint* clades that resolve to the same taxon (a
    non-monophyletic taxon) therefore EACH keep their label, so the split is
    visible on the tree rather than collapsed to a single occurrence. A
    monophyletic taxon is still labelled exactly once, on its crown
    (the topmost node of the run, which is also the largest — unchanged from
    the previous size-based behaviour for that case).

    Iterative depth-first walk carrying the set of label names seen on the
    path from the root, so it is safe on deep / ladderized trees (no Python
    recursion limit) and deterministic (the keep/clear decision depends only
    on ancestry, never on visit order).
    """
    stack: list[tuple[Clade, frozenset[str]]] = [(tree.root, frozenset())]
    while stack:
        clade, ancestor_names = stack.pop()
        name = getattr(clade, "_lca_name", None)
        child_ancestors = ancestor_names
        if name:
            if name in ancestor_names:
                clade._lca_name = None
                clade._lca_rank = None
            else:
                child_ancestors = ancestor_names | {name}
        for child in clade.clades:
            stack.append((child, child_ancestors))


# ---------------------------------------------------------------------------
# Same-species pair suppression
# ---------------------------------------------------------------------------

def _terminal_species(terminal: Clade, reps_by_short_id: dict[str, Sequence]) -> Optional[str]:
    seq = reps_by_short_id.get(terminal.name)
    if seq is None or seq.taxonomy is None:
        return None
    return seq.taxonomy.species


def suppress_same_species_pairs(
    tree: Tree,
    reps_by_short_id: dict[str, Sequence],
) -> None:
    """Clear the LCA label from any internal whose only children are
    two leaves of the same species — the species is already on the
    leaves, so duplicating it on the parent is just noise."""
    for node in tree.get_nonterminals():
        if not getattr(node, "_lca_name", None):
            continue
        children = node.clades
        if len(children) != 2:
            continue
        if not all(c.is_terminal() for c in children):
            continue
        species = [
            _terminal_species(c, reps_by_short_id) for c in children
        ]
        if species[0] and species[0] == species[1]:
            node._lca_name = None
            node._lca_rank = None


# ---------------------------------------------------------------------------
# PhyloXML rank validation
# ---------------------------------------------------------------------------

# PhyloXML 1.10 schema's permitted rank enumeration. Anything outside
# this list must be written as "other" (or the file fails to validate).
_PHYLOXML_VALID_RANKS = frozenset({
    "domain", "kingdom", "subkingdom", "branch",
    "infrakingdom", "superphylum", "phylum", "subphylum", "infraphylum",
    "microphylum", "superdivision", "division", "subdivision", "infradivision",
    "superclass", "class", "subclass", "infraclass",
    "superlegion", "legion", "sublegion", "infralegion",
    "supercohort", "cohort", "subcohort", "infracohort",
    "superorder", "order", "suborder", "superfamily", "family", "subfamily",
    "supertribe", "tribe", "subtribe", "infratribe",
    "genus", "subgenus", "superspecies", "species", "subspecies",
    "variety", "varietas", "subvariety",
    "form", "subform", "cultivar",
    "strain", "section", "subsection",
    "unknown", "other",
})


def phyloxml_rank(rank: Optional[str]) -> str:
    """Return ``rank`` if PhyloXML accepts it, else ``"other"``."""
    if rank and rank.lower() in _PHYLOXML_VALID_RANKS:
        return rank.lower()
    return "other"
