"""Per-taxon monophyly assessment over the annotated phyloXML trees (2J).

For every taxon at every rank present among a tree's leaves, classify how its
members sit on the rooted tree:

* **monophyletic** — the taxon's leaves form exactly one clade (the MRCA of
  its members subtends only members).
* **paraphyletic** — non-monophyletic, and the foreign leaves inside the
  members' MRCA span form a *single* excluded clade (the taxon is an
  ancestral grade with one derived group removed).
* **polyphyletic** — non-monophyletic, and the foreign leaves form *two or
  more* separate excluded clades (the taxon is interrupted in several places —
  the recombination / reassortment / misannotation signal).

Each rank is assessed **only among leaves classified at that rank** — a leaf
with no genus (say) is neither a member nor an intruder when genus monophyly
is judged, so an annotation gap is never mistaken for non-monophyly. A taxon
with fewer than two classified leaves has no meaningful status and is omitted.

The monophyletic call is an exact MRCA test; the **para-vs-poly split is a
documented HEURISTIC** — it keys on how the *intruders* are distributed
(``intruder_clusters``: 1 → para, ≥2 → poly), the standard topology-only
convention. The strict cladistic definition lets a taxon with several nested
carve-outs still be paraphyletic, so a taxon flagged polyphyletic here whose
``intruder_clusters`` is small may be paraphyletic under that stricter reading.
The unambiguous supporting numbers — ``n_clusters`` (separate in-group
blocks), ``n_intruders``, ``intruder_clusters`` — are reported alongside the
tag so the reader can always judge for themselves.

Computed by sweeping the retained ``*_tree.xml`` files at end of run (every
tree repseq builds — whole-genome 2E, per-protein / extra / segment /
polyprotein 2F, pre-cluster — ends in ``_tree.xml``), reading each leaf's
taxonomy from the ``repseq:<rank>`` ``<property>`` ladder the phyloXML writer
emits. Fully decoupled (no threading), same posture as the PDF and
conservation sweeps: soft-fails to no file when there are no trees, and
per-tree on an unparseable XML. Running it on every tree is the point —
a taxon monophyletic on the whole-genome tree but polyphyletic on a marker
tree is the per-marker reassortment/recombination signal.

The most informative power-pairing is with ``{prefix}_incongruence.tsv``
(tree-vs-tree RF distance): this report is the taxon-resolved companion.
"""

from __future__ import annotations

import logging
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

_PXNS = "{http://www.phyloxml.org}"

# Ranks assessed, fine→coarse. Species is deliberately skipped — viral species
# monophyly is too often violated to be a useful signal (same rationale as the
# taxonomic reports).
_RANKS: tuple[str, ...] = (
    "subgenus", "genus", "subfamily", "family",
    "suborder", "order", "subclass", "class",
)
_RANK_SET = frozenset(_RANKS)
_RANK_INDEX = {r: i for i, r in enumerate(_RANKS)}


class _Node:
    __slots__ = ("children", "leaves")

    def __init__(self) -> None:
        self.children: list["_Node"] = []
        self.leaves: frozenset[int] = frozenset()


def _parse_tree(path: Path) -> Optional[tuple[_Node, list[dict[str, str]]]]:
    """Parse a phyloXML into ``(root, leaf_taxa)``.

    ``leaf_taxa[i]`` is leaf *i*'s ``{rank: taxon}`` map (read from its
    ``repseq:<rank>`` properties); ``root.leaves`` etc. carry the set of leaf
    indices below each node. Returns None on a malformed / non-phyloXML file.
    """
    try:
        root_el = ET.parse(path).getroot()
    except Exception as exc:  # malformed XML → caller skips this tree
        logger.warning("[monophyly] could not parse %s: %s", path, exc)
        return None
    phylogeny = root_el.find(f"{_PXNS}phylogeny")
    if phylogeny is None:
        return None
    root_clade = phylogeny.find(f"{_PXNS}clade")
    if root_clade is None:
        return None

    leaf_taxa: list[dict[str, str]] = []

    def build(elem: ET.Element) -> _Node:
        node = _Node()
        child_els = elem.findall(f"{_PXNS}clade")
        if child_els:
            for ce in child_els:
                node.children.append(build(ce))
        else:
            idx = len(leaf_taxa)
            taxa: dict[str, str] = {}
            for prop in elem.findall(f"{_PXNS}property"):
                ref = prop.get("ref", "")
                if ref.startswith("repseq:"):
                    rank = ref.split(":", 1)[1]
                    text = (prop.text or "").strip()
                    if rank in _RANK_SET and text:
                        taxa[rank] = text
            leaf_taxa.append(taxa)
            node.leaves = frozenset({idx})
        return node

    root = build(root_clade)

    def fill(node: _Node) -> None:
        if node.children:
            acc: set[int] = set()
            for child in node.children:
                fill(child)
                acc |= child.leaves
            node.leaves = frozenset(acc)

    fill(root)
    return root, leaf_taxa


def _all_nodes(root: _Node) -> list[_Node]:
    out: list[_Node] = []
    stack = [root]
    while stack:
        node = stack.pop()
        out.append(node)
        stack.extend(node.children)
    return out


def _count_clusters(
    root: _Node, members: frozenset[int], labeled: frozenset[int]
) -> int:
    """Number of maximal clades whose *labeled* leaves are all members.

    Restricted to ``labeled`` (leaves classified at the rank in question), so
    leaves with no label at that rank don't split an in-group block — they're
    treated as absent for this rank, not as foreign.
    """
    count = 0
    stack = [root]
    while stack:
        node = stack.pop()
        induced = node.leaves & labeled
        if induced and induced <= members:
            count += 1            # maximal block — don't descend further
        else:
            stack.extend(node.children)
    return count


def assess_tree(root: _Node, leaf_taxa: list[dict[str, str]]) -> list[dict]:
    """Return one status row per taxon (n_leaves ≥ 2) per assessed rank."""
    nodes = _all_nodes(root)
    rows: list[dict] = []
    for rank in _RANKS:
        members_by_taxon: dict[str, set[int]] = {}
        leaf_taxon: list[Optional[str]] = []
        for taxa in leaf_taxa:
            t = taxa.get(rank)
            leaf_taxon.append(t)
            if t:
                members_by_taxon.setdefault(t, set()).add(len(leaf_taxon) - 1)
        # Only leaves classified at this rank participate; unlabeled leaves are
        # neither members nor intruders (an annotation gap must not read as
        # non-monophyly).
        labeled = frozenset(i for i, t in enumerate(leaf_taxon) if t)

        for taxon, member_set in members_by_taxon.items():
            if len(member_set) < 2:
                continue  # a single leaf has no meaningful monophyly status
            members = frozenset(member_set)
            # MRCA = the smallest node whose leaf set is a superset of members.
            mrca_leaves = root.leaves
            for node in nodes:
                if members <= node.leaves and len(node.leaves) < len(mrca_leaves):
                    mrca_leaves = node.leaves
            # Intruders = leaves classified at this rank that fall inside the
            # members' MRCA span but belong to another taxon (unlabeled leaves
            # excluded).
            intruders = (mrca_leaves & labeled) - members
            n_clusters = _count_clusters(root, members, labeled)
            if not intruders:
                status = "monophyletic"
                intruder_clusters = 0
                intruder_taxa: list[str] = []
            else:
                intruder_taxa = sorted(
                    {leaf_taxon[j] for j in intruders if leaf_taxon[j]}
                )
                # Heuristic para/poly call by how the INTRUDERS are
                # distributed inside the members' MRCA span: one excluded
                # clade → the taxon is an ancestral grade with a single
                # derived group removed (paraphyletic); two or more separate
                # excluded clades → the taxon is interrupted in multiple
                # places (polyphyletic). This (unlike an MRCA-span test)
                # correctly flags a taxon shattered across the whole tree.
                intruder_clusters = _count_clusters(root, intruders, labeled)
                status = (
                    "paraphyletic" if intruder_clusters == 1 else "polyphyletic"
                )
            rows.append({
                "rank": rank,
                "taxon": taxon,
                "n_leaves": len(member_set),
                "status": status,
                "n_clusters": n_clusters,
                "n_intruders": len(intruders),
                "intruder_clusters": intruder_clusters,
                "intruder_taxa": ";".join(intruder_taxa),
            })
    return rows


def write_monophyly_report(out_dir: Path, prefix: str) -> Optional[Path]:
    """Sweep every ``*_tree.xml`` under *out_dir* into ``{prefix}_monophyly.tsv``.

    Returns the written path, or None when there are no trees / nothing to
    report. Soft by construction — the caller wraps it so a parse failure
    never voids the trees that were already built.
    """
    xmls = sorted(out_dir.rglob("*_tree.xml"))
    if not xmls:
        return None

    rows_out: list[dict] = []
    for path in xmls:
        rel = path.relative_to(out_dir)
        parsed = _parse_tree(path)
        if parsed is None:
            continue
        root, leaf_taxa = parsed
        for row in assess_tree(root, leaf_taxa):
            row["tree"] = str(rel)
            rows_out.append(row)

    if not rows_out:
        return None

    rows_out.sort(
        key=lambda r: (r["tree"], _RANK_INDEX.get(r["rank"], 99), r["taxon"])
    )
    path = out_dir / f"{prefix}_monophyly.tsv"
    with open(path, "w") as fh:
        fh.write(
            "tree\trank\ttaxon\tn_leaves\tstatus\t"
            "n_clusters\tn_intruders\tintruder_clusters\tintruder_taxa\n"
        )
        for r in rows_out:
            fh.write(
                f"{r['tree']}\t{r['rank']}\t{r['taxon']}\t{r['n_leaves']}\t"
                f"{r['status']}\t{r['n_clusters']}\t{r['n_intruders']}\t"
                f"{r['intruder_clusters']}\t{r['intruder_taxa']}\n"
            )
    return path
