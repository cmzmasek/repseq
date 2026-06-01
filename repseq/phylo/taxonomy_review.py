"""Phylogeny-based taxonomy review (v0.39.0).

Some input records carry **missing** or **wrong** taxonomic ranks
(blank genus, "unclassified", an out-of-date genus that pre-dates an
ICTV reclassification, …). Once repseq builds a well-supported,
LCA-labelled tree of the representatives, the tree itself is evidence
about those ranks: a leaf nested deep inside a strongly-supported,
taxonomically *pure* clade is very likely the same taxon as its
neighbours.

This module turns that intuition into an auditable report. For each
representative leaf and each configured rank (default family → genus →
subgenus; species is deliberately excluded — viral species monophyly
is too often violated), it walks from the leaf up to the **smallest
enclosing clade** that is trustworthy:

  * branch support ≥ ``min_support`` (normalised to 0-100; FastTree
    SH-like values are rescaled from [0,1] like the phyloXML writer),
  * a single value ``M`` holds a strong majority among the clade's
    *labelled* leaves at that rank — purity ≥ ``min_purity``,
  * at least ``min_agreeing`` labelled neighbours back ``M``, ideally
    including a RefSeq / reviewed **anchor** (required when
    ``require_refseq_anchor``).

Two operations fall out of the same walk:

  * **impute_missing** — the leaf is blank at the rank → fill from
    ``M`` (high-confidence only is written into the corrected output
    copies; everything ≥ medium lands in the review TSV).
  * **conflict_flag** — the leaf is populated but disagrees with a
    confident ``M`` → reported as a suggestion only, **never**
    auto-applied.

Ranks are processed coarse → fine and imputations are kept
**hierarchy-consistent**: a finer rank is only imputed from neighbours
that agree with the leaf's effective coarser ranks (original or
already-imputed), and a blank coarser rank is back-filled when the
agreeing neighbours share one value there. So we never impute a
subgenus whose implied genus contradicts the leaf's genus.

The tree, its colouring, and its LCA labels are **not** modified — this
is a side report (v1 decision: no feedback loop). The whole step is
opt-in via ``phylo.taxonomy_review.enabled`` and soft-fails (a stderr
note, no file) on any error, like the other phylo extras.
"""

from __future__ import annotations

import copy
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

from Bio.Phylo.BaseTree import Tree

from ..models import Sequence, TaxonomyInfo
from .phyloxml_writer import _confidence_type_for, _normalize_confidence

logger = logging.getLogger(__name__)

DEFAULT_RANKS = ["family", "genus", "subgenus"]

# Coarse → fine, so a configured rank list is always evaluated in the
# right order regardless of how the user wrote it.
_RANK_ORDER = ["class", "subclass", "order", "suborder", "family",
               "subfamily", "genus", "subgenus", "species"]

# Values that mean "no real label" — treated as blank (imputation
# targets), mirroring the colouring module's missing set.
_MISSING = {"", "unknown", "na", "n/a", "none", "null", "?", "unclassified"}


def _clean(value: Optional[str]) -> Optional[str]:
    """Return a real taxon label or None (blank/placeholder → None)."""
    if value is None:
        return None
    v = str(value).strip()
    return None if v.lower() in _MISSING else v


def _rank_value(seq: Sequence, rank: str) -> Optional[str]:
    if seq.taxonomy is None:
        return None
    return _clean(seq.taxonomy.get_rank(rank))


@dataclass
class Verdict:
    """One (leaf, rank) finding written to the review TSV."""

    seq_id: str
    accession: str
    organism: str
    rank: str
    current_value: str          # "" when blank
    suggested_value: str
    action: str                 # "impute_missing" | "conflict_flag"
    clade_support: int
    clade_purity: float
    n_agreeing: int
    anchor_refseq: bool
    confidence: str             # "high" | "medium"


def _ancestors_small_to_large(tree: Tree, leaf) -> list:
    """Enclosing internal clades of ``leaf``, smallest (parent) first.

    ``Tree.get_path`` returns ``[child_of_root, …, parent, leaf]`` (root
    excluded). The enclosing internals are that list minus the leaf,
    reversed, with the root appended last.
    """
    path = tree.get_path(leaf)
    internals = list(reversed(path[:-1])) if len(path) > 1 else []
    internals.append(tree.root)
    return internals


def run_taxonomy_review(
    tree: Tree,
    reps_by_short_id: dict[str, Sequence],
    *,
    tree_tool: str,
    cfg: dict[str, Any],
    out_dir: Path,
    file_prefix: str,
) -> dict[str, Any]:
    """Walk the rooted tree and write ``{prefix}_taxonomy_review.tsv``.

    Returns a dict ``{"path": Path|None, "imputations": {seq_id:
    {rank: value}}, "verdicts": [Verdict, …]}``. ``imputations`` carries
    only the **high-confidence impute_missing** calls (what the caller
    applies to the corrected output copies); conflicts and medium calls
    are reported in the TSV but never auto-applied. Returns ``{}`` when
    the step is disabled.
    """
    rcfg = (cfg.get("phylo", {}) or {}).get("taxonomy_review", {}) or {}
    if not rcfg.get("enabled", False):
        return {}

    ranks_cfg = [str(r).lower() for r in (rcfg.get("ranks") or DEFAULT_RANKS)]
    ranks = sorted(
        {r for r in ranks_cfg if r in _RANK_ORDER},
        key=_RANK_ORDER.index,
    )
    if not ranks:
        return {}

    min_support = float(rcfg.get("min_support", 90))
    min_purity = float(rcfg.get("min_purity", 0.9))
    min_agreeing = int(rcfg.get("min_agreeing", 3))
    require_anchor = bool(rcfg.get("require_refseq_anchor", True))
    # Relaxed bar for the "medium" tier (reported but not auto-applied).
    med_support = min(min_support, 70.0)
    med_purity = min(min_purity, 0.8)

    conf_type = _confidence_type_for(tree_tool, None)

    # Precompute per-internal: enclosed leaf short-ids and 0-100 support.
    node_leaf_ids: dict[int, list[str]] = {}
    node_support: dict[int, Optional[int]] = {}
    for node in tree.get_nonterminals():
        node_leaf_ids[id(node)] = [t.name for t in node.get_terminals() if t.name]
        node_support[id(node)] = _normalize_confidence(
            getattr(node, "confidence", None), conf_type
        )

    verdicts: list[Verdict] = []
    imputations: dict[str, dict[str, str]] = {}

    for leaf in tree.get_terminals():
        sid = leaf.name
        seq = reps_by_short_id.get(sid)
        if seq is None:
            continue
        ancestors = _ancestors_small_to_large(tree, leaf)
        # Effective coarser context: original populated ranks, updated
        # with imputations as we move coarse → fine (hierarchy gate).
        effective: dict[str, Optional[str]] = {
            r: _rank_value(seq, r) for r in ranks
        }
        for rank in ranks:
            current = _rank_value(seq, rank)
            coarser = [r for r in ranks if _RANK_ORDER.index(r) < _RANK_ORDER.index(rank)]
            constraint = {r: effective[r] for r in coarser if effective[r]}

            ev = _find_evidence(
                sid, rank, ancestors, node_leaf_ids, node_support,
                reps_by_short_id, constraint,
                med_support=med_support, med_purity=med_purity,
                min_agreeing=min_agreeing,
            )
            if ev is None:
                continue
            value, support, purity, n_agree, anchor, agreeing_seqs = ev

            is_high = (
                support >= min_support
                and purity >= min_purity
                and (anchor or not require_anchor)
            )
            confidence = "high" if is_high else "medium"

            if current is None:
                action = "impute_missing"
                effective[rank] = value
                if confidence == "high":
                    imputations.setdefault(seq.id, {})[rank] = value
                    # Back-fill any still-blank coarser rank the agreeing
                    # neighbours unanimously share (keeps the imputed
                    # lineage internally consistent).
                    _backfill_coarser(
                        seq, ranks, rank, effective, agreeing_seqs,
                        imputations, verdicts, support, purity, n_agree,
                        anchor,
                    )
            elif value != current:
                action = "conflict_flag"
            else:
                continue  # consistent — no row

            verdicts.append(Verdict(
                seq_id=seq.id,
                accession=seq.accession or seq.isolate_id or seq.id,
                organism=seq.organism or "",
                rank=rank,
                current_value=current or "",
                suggested_value=value,
                action=action,
                clade_support=support,
                clade_purity=round(purity, 3),
                n_agreeing=n_agree,
                anchor_refseq=anchor,
                confidence=confidence,
            ))

    path = _write_review_tsv(verdicts, out_dir, file_prefix) if verdicts else None
    return {"path": path, "imputations": imputations, "verdicts": verdicts}


def _find_evidence(
    sid: str,
    rank: str,
    ancestors: list,
    node_leaf_ids: dict[int, list[str]],
    node_support: dict[int, Optional[int]],
    reps_by_short_id: dict[str, Sequence],
    constraint: dict[str, str],
    *,
    med_support: float,
    med_purity: float,
    min_agreeing: int,
) -> Optional[tuple[str, int, float, int, bool, list[Sequence]]]:
    """Smallest enclosing clade that confidently names ``rank``.

    Returns ``(value, support, purity, n_agreeing, anchor, agreeing_seqs)``
    or None. Neighbours are restricted to those matching ``constraint``
    (the leaf's effective coarser ranks) so the call is hierarchy-safe.
    """
    for node in ancestors:
        support = node_support.get(id(node))
        if support is None or support < med_support:
            continue
        # Labelled neighbours (exclude the leaf itself), honouring the
        # coarser-rank constraint.
        labelled: list[tuple[str, Sequence]] = []
        for nid in node_leaf_ids.get(id(node), []):
            if nid == sid:
                continue
            nseq = reps_by_short_id.get(nid)
            if nseq is None:
                continue
            if any(_rank_value(nseq, cr) != cv for cr, cv in constraint.items()):
                continue
            val = _rank_value(nseq, rank)
            if val is not None:
                labelled.append((val, nseq))
        if not labelled:
            continue
        # Modal value.
        counts: dict[str, int] = {}
        for val, _ in labelled:
            counts[val] = counts.get(val, 0) + 1
        modal = max(counts, key=lambda k: counts[k])
        n_agree = counts[modal]
        purity = n_agree / len(labelled)
        if n_agree < min_agreeing or purity < med_purity:
            continue
        agreeing_seqs = [s for v, s in labelled if v == modal]
        anchor = any(s.is_refseq or s.is_reviewed for s in agreeing_seqs)
        return modal, support, purity, n_agree, anchor, agreeing_seqs
    return None


def _backfill_coarser(
    seq: Sequence,
    ranks: list[str],
    rank: str,
    effective: dict[str, Optional[str]],
    agreeing_seqs: list[Sequence],
    imputations: dict[str, dict[str, str]],
    verdicts: list[Verdict],
    support: int,
    purity: float,
    n_agree: int,
    anchor: bool,
) -> None:
    """Fill still-blank coarser ranks the agreeing neighbours all share.

    Guarantees the imputed lineage is internally consistent: imputing
    ``subgenus = Sarbecovirus`` from a set of leaves that are all
    ``genus = Betacoronavirus`` also fills a blank genus.
    """
    for cr in ranks:
        if _RANK_ORDER.index(cr) >= _RANK_ORDER.index(rank):
            continue
        if effective.get(cr):
            continue  # already populated/imputed
        vals = {_rank_value(s, cr) for s in agreeing_seqs}
        vals.discard(None)
        if len(vals) != 1:
            continue
        value = next(iter(vals))
        effective[cr] = value
        imputations.setdefault(seq.id, {})[cr] = value
        verdicts.append(Verdict(
            seq_id=seq.id,
            accession=seq.accession or seq.isolate_id or seq.id,
            organism=seq.organism or "",
            rank=cr,
            current_value="",
            suggested_value=value,
            action="impute_missing",
            clade_support=support,
            clade_purity=round(purity, 3),
            n_agreeing=n_agree,
            anchor_refseq=anchor,
            confidence="high",
        ))


def _impute_tax(tax: Optional[TaxonomyInfo], imp: dict[str, str]) -> TaxonomyInfo:
    """Return a copy of ``tax`` with the imputed ranks filled in.

    Standard ranks land on their dataclass field; everything else
    (subgenus, subfamily, …) goes into the ``lineage`` map, matching
    where ``TaxonomyInfo.get_rank`` reads each rank from.
    """
    out = copy.deepcopy(tax) if tax is not None else TaxonomyInfo()
    for rank, value in imp.items():
        if rank == "genus":
            out.genus = value
        elif rank == "family":
            out.family = value
        elif rank == "order":
            out.order = value
        elif rank == "class":
            out.class_ = value
        elif rank == "species":
            out.species = value
        else:
            out.lineage[rank] = value
    return out


def apply_imputations(
    representatives: list[Sequence],
    complete_isolates: Optional[dict[str, list[Sequence]]],
    imputations: dict[str, dict[str, str]],
) -> tuple[list[Sequence], Optional[dict[str, list[Sequence]]]]:
    """Build imputation-corrected copies of the reps (+ segmented segments).

    Returns ``(corrected_reps, corrected_complete_isolates)``. Originals
    are never mutated — each touched Sequence is shallow-copied with a
    fresh, imputed ``TaxonomyInfo``. ``imputations`` is keyed by rep
    ``Sequence.id``; in segmented mode the same rank fills propagate to
    every segment of the rep's isolate so the corrected protein FASTA
    carries them too.
    """
    iso_to_imp: dict[str, dict[str, str]] = {}
    corrected_reps: list[Sequence] = []
    for rep in representatives:
        imp = imputations.get(rep.id)
        if not imp:
            corrected_reps.append(rep)
            continue
        rep2 = copy.copy(rep)
        rep2.taxonomy = _impute_tax(rep.taxonomy, imp)
        corrected_reps.append(rep2)
        iso = rep.isolate_id
        if not iso and rep.id.startswith("CONCAT|") and "|" in rep.id:
            iso = rep.id.split("|", 1)[1]
        if iso:
            iso_to_imp[iso] = imp

    corrected_ci: Optional[dict[str, list[Sequence]]] = None
    if complete_isolates:
        corrected_ci = {}
        for iso, segs in complete_isolates.items():
            imp = iso_to_imp.get(iso)
            if not imp:
                corrected_ci[iso] = segs
                continue
            new_segs = []
            for s in segs:
                s2 = copy.copy(s)
                s2.taxonomy = _impute_tax(s.taxonomy, imp)
                new_segs.append(s2)
            corrected_ci[iso] = new_segs
    return corrected_reps, corrected_ci


_REVIEW_COLUMNS = [
    "accession", "organism", "rank", "current_value", "suggested_value",
    "action", "clade_support", "clade_purity", "n_agreeing",
    "anchor_refseq", "confidence",
]


def _write_review_tsv(
    verdicts: list[Verdict], out_dir: Path, file_prefix: str
) -> Path:
    """Write the verdicts to ``{prefix}_taxonomy_review.tsv``.

    Sorted impute_missing first, then conflict_flag, then by accession /
    rank so a reviewer reads the safe fills before the flags.
    """
    path = out_dir / f"{file_prefix}_taxonomy_review.tsv"
    out_dir.mkdir(parents=True, exist_ok=True)
    order = {"impute_missing": 0, "conflict_flag": 1}
    rows = sorted(
        verdicts,
        key=lambda v: (order.get(v.action, 9), v.accession,
                       _RANK_ORDER.index(v.rank) if v.rank in _RANK_ORDER else 99),
    )
    with open(path, "w") as fh:
        fh.write("\t".join(_REVIEW_COLUMNS) + "\n")
        for v in rows:
            fh.write("\t".join([
                v.accession, v.organism, v.rank, v.current_value,
                v.suggested_value, v.action, str(v.clade_support),
                f"{v.clade_purity:.3f}", str(v.n_agreeing),
                "TRUE" if v.anchor_refseq else "FALSE", v.confidence,
            ]) + "\n")
    return path
