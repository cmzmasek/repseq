"""Phylogeny orchestrator: short-id rename → MAFFT → FastTree → phyloXML.

The pipeline:

1. Skip with a warning if there are fewer than 3 representatives —
   FastTree cannot build a tree on a pair or singleton.
2. Decide protein vs nucleotide from the first representative's
   ``seq_type`` (the upstream pipeline guarantees a single alphabet).
3. Assign deterministic short ids (``S0001``, ``S0002``, …) so that
   downstream tools — many of which choke on long names, pipes, or
   whitespace — see only clean tokens. The mapping is written to
   ``{prefix}_tree_id_map.tsv`` so the MSA and Newick (which keep the
   short ids) remain readable.
4. Run MAFFT → MSA FASTA (short ids retained).
5. Run FastTree → Newick (short ids retained).
6. Convert Newick → phyloXML via ``Bio.Phylo``, restoring every
   terminal clade's ``name`` to the original ``seq.id``.

Outputs (all under ``{prefix}_*``):
    {prefix}_msa.fasta           aligned MSA, short-id headers
    {prefix}_tree.nwk            FastTree Newick, short-id leaves
    {prefix}_tree.xml            phyloXML, original seq.id leaves
    {prefix}_tree_id_map.tsv     short_id<TAB>original_id

The orchestrator never raises out of the click command — it catches its
own subprocess and conversion errors and reports them to stderr, matching
the existing ``--plot`` behaviour.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Optional

from Bio import Phylo

from ..models import Sequence, SequenceType
from .fasttree import FastTreeError, run_fasttree
from .mafft import MafftError, run_mafft

logger = logging.getLogger(__name__)


class PhyloError(RuntimeError):
    """Raised when the phylogeny step cannot proceed — surfaced to stderr,
    not propagated, so a tree failure never destroys an otherwise good run."""


_SHORT_ID_FMT = "S{:04d}"


def _use_protein_sequence(reps: list[Sequence], cfg: Optional[dict[str, Any]]) -> bool:
    """True when the phylo input should be each rep's protein_sequence.

    Active when ``clustering.alphabet="protein"`` AND every rep actually
    carries a protein_sequence (a no-resolve fallback or missing-marker
    drop could leave some empty; we never want a half-protein, half-NT
    alignment).
    """
    if cfg is None:
        return False
    if cfg.get("clustering", {}).get("alphabet") != "protein":
        return False
    return all(r.protein_sequence for r in reps)


def _is_protein(reps: list[Sequence]) -> bool:
    """Pick a single alphabet for the whole rep set when reading seq.sequence.

    The upstream pipeline never mixes protein and nucleotide reps in a
    single run, so checking the first one is enough; unknown defaults
    to protein since FastTree's protein model is the default and works
    on anything alphabetic.
    """
    return reps[0].seq_type != SequenceType.NUCLEOTIDE


def _build_id_map(reps: list[Sequence]) -> dict[str, str]:
    """Assign ``S0001``, ``S0002``, … to each representative.

    Returns a short_id → original seq.id mapping. Order is stable
    (the input order), so re-running on the same reps produces
    identical ids — convenient for diff-based comparison.
    """
    return {_SHORT_ID_FMT.format(i + 1): rep.id for i, rep in enumerate(reps)}


def _write_short_id_fasta(
    reps: list[Sequence],
    id_map: dict[str, str],
    path: Path,
    use_protein: bool = False,
) -> None:
    """Write a FASTA whose header is the short id and nothing else.

    Mirrors ``clustering.mmseqs2._write_id_fasta`` — sole-token headers
    so the alignment/tree tools cannot lose track of identity at any
    whitespace or pipe boundary. With ``use_protein=True`` the body comes
    from ``rep.protein_sequence`` (the marker / per-isolate AA concat).
    """
    reverse = {v: k for k, v in id_map.items()}
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        for rep in reps:
            short = reverse[rep.id]
            fh.write(f">{short}\n")
            seq = (rep.protein_sequence if use_protein else rep.sequence) or ""
            for i in range(0, len(seq), 70):
                fh.write(seq[i:i + 70] + "\n")


def _write_id_map(id_map: dict[str, str], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write("short_id\toriginal_id\n")
        for short, original in id_map.items():
            fh.write(f"{short}\t{original}\n")


def _newick_to_phyloxml(
    newick_path: Path,
    phyloxml_path: Path,
    id_map: dict[str, str],
) -> None:
    """Convert Newick → phyloXML, restoring each terminal's original id.

    FastTree's Newick uses internal-node labels for support values; Bio.Phylo
    parses those as clade ``confidence`` rather than ``name``, so renaming
    only the terminals leaves internal-node support intact.
    """
    tree = Phylo.read(str(newick_path), "newick")
    for terminal in tree.get_terminals():
        if terminal.name in id_map:
            terminal.name = id_map[terminal.name]
    phyloxml_path.parent.mkdir(parents=True, exist_ok=True)
    Phylo.write([tree], str(phyloxml_path), "phyloxml")


def run_phylogeny(
    representatives: list[Sequence],
    cfg: dict[str, Any],
    out_dir: Path,
    prefix: str,
) -> list[Path]:
    """Run the full phylogeny pipeline on a list of representative sequences.

    Returns the list of files written so the caller can append them to its
    output manifest. Raises :class:`PhyloError` when the step cannot
    proceed (under 3 reps, missing binary, or subprocess failure) — the
    caller is expected to catch and report rather than crash the run.
    """
    n = len(representatives)
    if n < 3:
        raise PhyloError(
            f"need >= 3 representatives to build a tree, got {n}"
        )

    use_protein = _use_protein_sequence(representatives, cfg)
    is_protein = use_protein or _is_protein(representatives)
    id_map = _build_id_map(representatives)

    input_fasta = out_dir / f"{prefix}_phylo_input.fasta"
    msa_fasta = out_dir / f"{prefix}_msa.fasta"
    newick_path = out_dir / f"{prefix}_tree.nwk"
    phyloxml_path = out_dir / f"{prefix}_tree.xml"
    id_map_path = out_dir / f"{prefix}_tree_id_map.tsv"

    _write_short_id_fasta(representatives, id_map, input_fasta, use_protein=use_protein)
    _write_id_map(id_map, id_map_path)

    try:
        run_mafft(input_fasta, msa_fasta, cfg)
    except MafftError as exc:
        raise PhyloError(str(exc)) from exc

    try:
        run_fasttree(msa_fasta, newick_path, cfg, is_protein=is_protein)
    except FastTreeError as exc:
        raise PhyloError(str(exc)) from exc

    try:
        _newick_to_phyloxml(newick_path, phyloxml_path, id_map)
    except Exception as exc:
        raise PhyloError(f"Newick → phyloXML conversion failed: {exc}") from exc

    # The temp input is redundant once the MSA is written.
    try:
        input_fasta.unlink()
    except OSError:
        pass

    return [msa_fasta, newick_path, phyloxml_path, id_map_path]
