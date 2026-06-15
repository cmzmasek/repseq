"""Pre-cluster overview tree (2H).

A rough, single-pass phylogeny of **every post-QC sequence** (one leaf
per CONCAT isolate in segmented mode), built BEFORE clustering would
have collapsed redundancy. The point is diagnostic, not analytical:
the bench scientist opens this tree alongside the post-cluster tree
(2E) and can see at a glance where the elected representatives land
in the broader diversity of the input pool. Representative leaves get
a ``[repr] `` prefix on their phyloXML ``<name>`` for visual
identification; non-rep leaves carry the same formatted label without
the prefix. Every leaf also carries a machine-readable
``repseq:is_representative`` boolean ``<property>`` (``true``/``false``)
so a viewer (Archaeopteryx) can filter / select / colour by it — this
is the only tree that emits it, since every other repseq tree contains
only representatives.

The pipeline here is intentionally **hard-coded for speed**, regardless
of the rest of ``phylo:``:

* **Leaf cap** ``phylo.pre_cluster_tree.max_leaves`` (default 5000).
  FastTree memory scales ~linearly with leaf count, so an uncapped pool
  of tens of thousands of leaves needs hundreds of GB and gets
  OOM-killed — and a tree that large isn't legible anyway. Above the cap
  the tree is built on **all representatives + a random sample of the
  non-representative background** up to ``max_leaves`` (reps are never
  dropped — they're the point of the overview); ``0`` = no cap. This
  mirrors the ``--plot`` subsample (``viz/clustering_plot.py``: cap 2000,
  seed 42, reps always kept).
* MAFFT ``--retree 1`` (single-pass FFT-NS-1) — no ``--auto``,
  no L-INS-i. **Above** ``phylo.pre_cluster_tree.parttree_threshold``
  (default 10000) leaves it switches to ``--retree 2 --parttree``
  (MAFFT's PartTree guide), which skips the O(N²) distance matrix that
  OOMs ``--retree 1`` on huge pools and scales to 10⁵⁺ sequences.
* **FastTree** (no IQ-TREE / ModelFinder / UFBoot) regardless of
  alphabet — the user explicitly asked for "rough", and any model
  selection would dominate the runtime budget.
* **Midpoint rooting** only — no taxonomy-guided rooting, no MAD.
* **No LCA annotation, no trimAl, no bootstrap.** The internal nodes
  carry only the branch lengths FastTree wrote.

Outputs land alongside the rest of the run's phylo files:

* ``{prefix}_pre_cluster_tree.nwk`` — Newick, short-id leaves.
* ``{prefix}_pre_cluster_tree.xml`` — phyloXML with taxonomy
  ``<property>`` enrichment + per-leaf taxonomy colouring (the same
  colour palette 2E and 2F use), the ``[repr] `` prefix on rep leaves,
  and a ``repseq:is_representative`` boolean ``<property>`` on every
  leaf.
* ``{prefix}_pre_cluster_tree_id_map.tsv`` — three columns
  (``short_id``, ``accession``, ``is_rep``) so a user grepping for
  reps without opening the XML can find them.

Soft-fails ``PhyloError`` like the rest of the phylo step: missing
MAFFT/FastTree, fewer than 3 sequences, or any subprocess error
surfaces as a stderr line and the rest of the run continues.
"""

from __future__ import annotations

import logging
import random
import tempfile
from pathlib import Path
from typing import Any, Optional

from Bio import Phylo

from ..models import Sequence, SequenceType
from .coloring import build_color_scheme
from .fasttree import FastTreeError, run_fasttree
from .mafft import MafftError, run_mafft
from .phyloxml_writer import write_phyloxml
from .pipeline import (
    PhyloError,
    _SHORT_ID_FMT,
    _write_short_id_fasta,
)
from .rooting import root_tree
from . import fasttree as fasttree_mod
from . import mafft as mafft_mod

logger = logging.getLogger(__name__)


def _use_protein_sequence(
    sequences: list[Sequence], cfg: dict[str, Any],
) -> bool:
    """True when the pre-cluster tree should be built on ``protein_sequence``.

    Mirrors ``pipeline._use_protein_sequence``: the alphabet that fed
    clustering is what the overview tree should depict, so the user
    can compare like with like. A pool where any sequence lacks
    ``protein_sequence`` falls back to NT — we never want a
    half-protein, half-NT alignment.
    """
    if cfg.get("clustering", {}).get("alphabet_for_clustering") != "protein":
        return False
    return all(s.protein_sequence for s in sequences)


def _build_id_map_pre(sequences: list[Sequence]) -> dict[str, str]:
    """``S0001`` → ``seq.id`` map for the pre-cluster pool.

    Insertion order matches the input order, so re-running on the
    same post-QC pool produces identical ids — convenient for diffing.
    """
    return {
        _SHORT_ID_FMT.format(i + 1): seq.id
        for i, seq in enumerate(sequences)
    }


def _subsample_pool(
    sequences: list[Sequence],
    rep_ids: set[str],
    max_leaves: int,
    *,
    seed: int = 42,
) -> tuple[list[Sequence], int]:
    """Cap the pre-cluster pool at ``max_leaves`` leaves for tractability.

    Mirrors the ``--plot`` subsample (``viz/clustering_plot.py``): **all
    representatives are always kept**; the remaining budget is filled with
    a deterministic random (``seed``) sample of the non-representative
    background. FastTree memory scales ~linearly with leaf count, so an
    uncapped pool of tens of thousands of leaves needs hundreds of GB and
    gets OOM-killed — and a tree that large isn't legible anyway.

    A pool at or below the cap (or ``max_leaves <= 0``) is returned
    unchanged. When the representatives alone meet or exceed the cap, all
    of them are kept and no background is added (the cap bounds only the
    background — reps are never dropped). Input order is preserved among
    the kept sequences so re-runs are diff-stable.

    Returns ``(kept_sequences, n_background_kept)``.
    """
    n = len(sequences)
    n_background_total = sum(1 for s in sequences if s.id not in rep_ids)
    if max_leaves <= 0 or n <= max_leaves:
        return sequences, n_background_total

    reps = [s for s in sequences if s.id in rep_ids]
    background = [s for s in sequences if s.id not in rep_ids]
    budget = max(0, max_leaves - len(reps))
    if budget < len(background):
        background = random.Random(seed).sample(background, budget)
    keep = {id(s) for s in reps} | {id(s) for s in background}
    kept = [s for s in sequences if id(s) in keep]
    return kept, len(background)


def _write_id_map_with_rep_flag(
    id_map: dict[str, str],
    rep_ids: set[str],
    path: Path,
) -> None:
    """Three-column TSV: short_id, accession, is_rep.

    Extends the two-column 2E format with an ``is_rep`` boolean
    (``TRUE``/``FALSE``) so a user grepping the id_map without
    opening the phyloXML can identify representative leaves.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write("short_id\taccession\tis_rep\n")
        for short, original in id_map.items():
            flag = "TRUE" if original in rep_ids else "FALSE"
            fh.write(f"{short}\t{original}\t{flag}\n")


def run_pre_cluster_phylogeny(
    sequences: list[Sequence],
    representatives: list[Sequence],
    cfg: dict[str, Any],
    out_dir: Path,
    prefix: str,
) -> list[Path]:
    """Build the pre-cluster overview tree of every post-QC sequence.

    ``sequences`` is the post-QC pool fed to the mode (CONCAT isolates
    in segmented mode). ``representatives`` is the elected subset;
    leaves with an id in this set get a ``[repr] `` prefix in the
    phyloXML ``<name>`` AND a ``repseq:is_representative`` boolean
    ``<property>`` (every leaf carries the property — ``true`` for the
    elected subset, ``false`` otherwise).

    Returns the list of files written (Newick, phyloXML, id_map).
    Raises ``PhyloError`` when the step cannot proceed — soft-failed
    upstream by the caller exactly like 2E.
    """
    if len(sequences) < 3:
        raise PhyloError(
            f"pre-cluster tree needs >= 3 sequences, got {len(sequences)}"
        )

    rep_ids = {r.id for r in representatives}

    # Leaf cap (phylo.pre_cluster_tree.max_leaves, default 5000): FastTree
    # memory scales ~linearly with leaf count, so subsample a huge pool to
    # all reps + a random background sample before building anything.
    # Mirrors the --plot subsample; 0 = no cap. n below is the CAPPED size,
    # so it also drives the parttree decision (a capped pool is back under
    # parttree_threshold and uses the fast --retree 1 path).
    pc_cfg = (cfg.get("phylo", {}) or {}).get("pre_cluster_tree", {}) or {}
    max_leaves = pc_cfg.get("max_leaves", 5000)
    n_pool = len(sequences)
    sequences, n_background = _subsample_pool(sequences, rep_ids, max_leaves)
    n = len(sequences)
    if n < n_pool:
        logger.warning(
            "[pre-cluster] pool of %d isolates capped to %d leaves (all %d "
            "representatives + %d background sample; set "
            "phylo.pre_cluster_tree.max_leaves: 0 to disable)",
            n_pool, n, n - n_background, n_background,
        )

    use_protein = _use_protein_sequence(sequences, cfg)

    id_map = _build_id_map_pre(sequences)
    # The Newick keeps the short ids (safe across phylo tools); the
    # human-readable label lands in the phyloXML <name> via the
    # writer's existing label-formatting path. We pass [repr] prefixes
    # so reps stand out at a glance in the XML.
    label_prefix_by_id = {
        seq.id: ("[repr] " if seq.id in rep_ids else "")
        for seq in sequences
    }

    bodies = {
        s.id: ((s.protein_sequence if use_protein else s.sequence) or "")
        for s in sequences
    }

    # MAFFT strategy, hard-coded regardless of phylo.mafft. The single-pass
    # FFT-NS-1 path (--retree 1) builds an O(N^2) distance matrix and gets
    # OOM-killed on huge pools, so at/above parttree_threshold we switch to
    # the PartTree guide (--retree 2 --parttree), which skips the full matrix
    # and scales to 10^5+ sequences. (pc_cfg was read above for max_leaves.)
    threshold = pc_cfg.get("parttree_threshold", 10000)
    if (
        isinstance(threshold, int)
        and not isinstance(threshold, bool)
        and n >= threshold
    ):
        mafft_args = ["--retree", "2", "--parttree"]
        logger.warning(
            "[pre-cluster] %d sequences >= parttree_threshold (%d) — using "
            "MAFFT --parttree (large-input mode)", n, threshold,
        )
    else:
        mafft_args = ["--retree", "1"]

    out_dir.mkdir(parents=True, exist_ok=True)
    newick_path = out_dir / f"{prefix}_pre_cluster_tree.nwk"
    phyloxml_path = out_dir / f"{prefix}_pre_cluster_tree.xml"
    id_map_path = out_dir / f"{prefix}_pre_cluster_tree_id_map.tsv"

    # MSA is intermediate-only: MAFFT + FastTree both need a real file,
    # but the user's locked decision is "Newick + phyloXML only", so
    # we write the MSA into the run's temp dir and clean it up at
    # exit.
    work_root = Path(cfg.get("temp_dir") or "/tmp/repseq")
    work_root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        dir=work_root, prefix="pre_cluster_",
    ) as td:
        td = Path(td)
        input_fasta = td / "input.fasta"
        msa_fasta = td / "msa.fasta"
        _write_short_id_fasta(
            sequences, id_map, input_fasta, bodies=bodies,
        )

        # use_auto=False drops --auto so the chosen extra args alone govern
        # (mafft_args picked above: --retree 1, or --retree 2 --parttree for
        # a large pool).
        try:
            run_mafft(
                input_fasta, msa_fasta, cfg,
                extra_args=mafft_args,
                use_auto=False,
            )
        except MafftError as exc:
            raise PhyloError(
                f"pre-cluster MAFFT failed: {exc}"
            ) from exc

        # FastTree regardless of clustering alphabet — IQ-TREE's
        # ModelFinder would dominate runtime on a many-thousand-leaf
        # rough tree. `is_protein` selects the FastTree model
        # (defaults to JTT for protein; -nt -gtr for nucleotide).
        try:
            run_fasttree(
                msa_fasta, newick_path, cfg,
                is_protein=use_protein,
            )
        except FastTreeError as exc:
            raise PhyloError(
                f"pre-cluster FastTree failed: {exc}"
            ) from exc

    # Parse the Newick, midpoint root, write phyloXML.
    try:
        parsed = Phylo.read(str(newick_path), "newick")
    except Exception as exc:
        raise PhyloError(
            f"pre-cluster Newick parse failed for {newick_path}: {exc}"
        ) from exc

    reps_by_short_id = dict(zip(id_map.keys(), sequences))
    try:
        parsed, rooting_used = root_tree(
            parsed, reps_by_short_id, method="midpoint",
        )
    except Exception as exc:
        logger.warning(
            "[pre-cluster] midpoint rooting failed: %s; leaving as-is", exc,
        )
        rooting_used = "none"

    # Shared colour palette is built over the FULL pool so the same
    # taxon gets the same colour in 2E / 2F / 2H of the same run.
    color_scheme = build_color_scheme(sequences, cfg)

    alphabet_label = "protein" if use_protein else "nucleotide"
    try:
        write_phyloxml(
            None,
            phyloxml_path,
            sequences,
            id_map,
            cfg=cfg or {},
            prefix=f"{prefix}_pre_cluster",
            alphabet=alphabet_label,
            msa_tool="MAFFT",
            msa_version=mafft_mod.tool_version(),
            tree_tool="FastTree",
            tree_version=fasttree_mod.tool_version(),
            model=("JTT" if use_protein else "GTR"),
            ufboot=None,
            extra_msa_args=mafft_args,
            extra_tree_args=[],
            tree=parsed,
            rooting_method=rooting_used,
            color_scheme=color_scheme,
            leaf_protein_ids=None,
            domain_architecture=False,
            label_prefix_by_id=label_prefix_by_id,
            representative_ids=rep_ids,
            basis_role="pre_cluster",
        )
    except Exception as exc:
        raise PhyloError(
            f"pre-cluster phyloXML write failed: {exc}"
        ) from exc

    _write_id_map_with_rep_flag(id_map, rep_ids, id_map_path)

    return [newick_path, phyloxml_path, id_map_path]
