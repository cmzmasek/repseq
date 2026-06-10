"""Cross-tree segment-status matrix (Ext-2) — the taxon-resolved
localisation of reassortment / recombination.

``{prefix}_monophyly.tsv`` reports one row per (taxon, rank, tree): a taxon's
mono/para/poly status on each tree a run built. That is the raw signal, but
the *biological* question is cross-tree: **is a taxon clean on the whole-genome
tree yet broken on exactly one marker tree?** — the signature of a single
segment having a different history (the per-marker reassortment call). The
incongruence table (``{prefix}_incongruence.tsv``) answers "do two trees
disagree?" with one global number; this matrix answers "*which taxon*, on
*which tree*?".

It is a pure pivot of ``{prefix}_monophyly.tsv`` — no trees are re-read, so it
inherits that report's support-aware collapse and (when
``phylo.monophyly.include_species`` is on) its species rows, which is the rank
where intra-genus reassortment lives. Written right after the monophyly sweep;
soft-fails to no file when there is no monophyly report.

Columns (one row per assessed (rank, taxon)):

* ``rank``, ``taxon``
* ``n_leaves`` — classified members (max across the trees the taxon appears on)
* ``n_trees`` — how many trees assessed this taxon
* ``n_nonmono`` — on how many of them it is para/poly
* ``genome_status`` — its status on the whole-genome tree (blank if that tree
  was not built, e.g. a ``--per-protein-phylo``-only run)
* ``nonmono_trees`` — ``;``-joined short labels of every tree it breaks on
* ``single_marker_break`` — **the localised reassortment call**: the short
  label of the *one* marker tree it breaks on, populated only when the taxon is
  monophyletic on the whole-genome tree AND breaks on exactly one (marker)
  tree. Blank otherwise. Filtering this column to non-empty yields the clean
  single-segment-discordant candidate list.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

from .flags import _read_tsv, _short_tree
from ..phylo.monophyly import _RANK_INDEX

logger = logging.getLogger(__name__)


def write_segment_status_matrix(out_dir: Path, prefix: str) -> Optional[Path]:
    """Pivot ``{prefix}_monophyly.tsv`` into ``{prefix}_segment_status_matrix.tsv``.

    Returns the written path, or None when there is no monophyly report / no
    rows. Soft by construction — a synthesis bug never voids the trees or the
    monophyly report it reads.
    """
    rows = _read_tsv(out_dir / f"{prefix}_monophyly.tsv")
    if not rows:
        return None

    genome_file = f"{prefix}_tree.xml"
    # (rank, taxon) -> {tree_relpath: status_row}
    by_taxon: dict[tuple[str, str], dict[str, dict]] = {}
    for r in rows:
        by_taxon.setdefault((r["rank"], r["taxon"]), {})[r["tree"]] = r

    out_rows: list[dict] = []
    for (rank, taxon), trees in by_taxon.items():
        genome = trees.get(genome_file)
        genome_status = genome["status"] if genome is not None else ""
        n_leaves = max(int(rr["n_leaves"]) for rr in trees.values())

        nonmono = {
            t: rr for t, rr in trees.items() if rr["status"] != "monophyletic"
        }
        nonmono_labels = sorted(_short_tree(t, prefix) for t in nonmono)

        # The localised call: clean on the whole-genome tree, broken on exactly
        # one (necessarily marker) tree. ``nonmono`` here can only hold marker
        # trees because the genome tree is monophyletic in this branch.
        single_marker_break = ""
        if genome is not None and genome_status == "monophyletic" and len(nonmono) == 1:
            single_marker_break = next(iter(nonmono_labels))

        out_rows.append({
            "rank": rank,
            "taxon": taxon,
            "n_leaves": n_leaves,
            "n_trees": len(trees),
            "n_nonmono": len(nonmono),
            "genome_status": genome_status,
            "nonmono_trees": ";".join(nonmono_labels),
            "single_marker_break": single_marker_break,
        })

    if not out_rows:
        return None

    # Sort fine→coarse rank, then localised-discordance rows first within a rank
    # (so the candidate calls surface at the top), then by taxon.
    out_rows.sort(
        key=lambda r: (
            _RANK_INDEX.get(r["rank"], 99),
            r["single_marker_break"] == "",
            r["taxon"],
        )
    )

    path = out_dir / f"{prefix}_segment_status_matrix.tsv"
    with open(path, "w") as fh:
        fh.write(
            "rank\ttaxon\tn_leaves\tn_trees\tn_nonmono\tgenome_status\t"
            "nonmono_trees\tsingle_marker_break\n"
        )
        for r in out_rows:
            fh.write(
                f"{r['rank']}\t{r['taxon']}\t{r['n_leaves']}\t{r['n_trees']}\t"
                f"{r['n_nonmono']}\t{r['genome_status']}\t{r['nonmono_trees']}\t"
                f"{r['single_marker_break']}\n"
            )
    return path
