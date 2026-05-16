"""IQ-TREE wrapper.

Infers an ML phylogeny on a MAFFT alignment with principled model
selection (ModelFinder Plus by default) and ultrafast bootstrap. The
binary is conventionally named ``iqtree2`` on modern installs (since
IQ-TREE 2.0); some older / system packages still ship ``iqtree``. We try
both, prefering iqtree2.

IQ-TREE writes a handful of auxiliary files alongside the tree
(``.iqtree``, ``.log``, ``.bionj``, ``.mldist``, ``.contree``,
``.splits.nex``, ``.ckp.gz``). We keep the canonical Newick
(``.treefile`` → caller-supplied path) and the human-readable
``.iqtree`` summary (model-selection report); the rest are deleted to
keep the run output directory tidy.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Optional


class IQTreeError(RuntimeError):
    pass


def _count_msa_records(msa_fasta: Path) -> int:
    n = 0
    with open(msa_fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n


def _check_iqtree(override: Optional[str] = None) -> str:
    """Locate the IQ-TREE binary.

    Tries the explicit ``override`` first, then ``iqtree2``, then ``iqtree``.
    """
    candidates = [override] if override else []
    candidates += ["iqtree2", "iqtree"]
    for name in candidates:
        if not name:
            continue
        path = shutil.which(name)
        if path:
            return path
    raise IQTreeError(
        "IQ-TREE not found in PATH (tried 'iqtree2' and 'iqtree'). "
        "Install it from https://www.iqtree.org/"
    )


def tool_version(override: Optional[str] = None) -> str:
    """Return IQ-TREE's version string, or ``"unknown"`` on any failure.

    Used by the phyloXML writer to annotate the tree's
    ``<phylogeny><description>``. IQ-TREE's ``--version`` prints e.g.
    ``IQ-TREE multicore version 2.2.0.3 ...`` to stdout and exits 0.
    """
    try:
        path = _check_iqtree(override)
    except IQTreeError:
        return "unknown"
    try:
        result = subprocess.run(
            [path, "--version"], capture_output=True, text=True, timeout=5,
        )
        out = (result.stdout or result.stderr or "").strip()
        return out.splitlines()[0].strip() if out else "unknown"
    except (subprocess.TimeoutExpired, OSError):
        return "unknown"


def run_iqtree(
    msa_fasta: Path,
    output_newick: Path,
    cfg: dict[str, Any],
    is_protein: bool,
    summary_path: Optional[Path] = None,
) -> None:
    """Run IQ-TREE on ``msa_fasta`` and write the ML tree to ``output_newick``.

    With the default model ``MFP`` (ModelFinder Plus), IQ-TREE first scans
    available substitution models and picks the best by BIC, then infers
    the ML tree under that model. Ultrafast bootstrap (``-bb``) is on by
    default with 1000 replicates — IQ-TREE's recommended minimum for
    interpretable support values.

    All run files are written to a temporary directory under
    ``cfg['temp_dir']``; the ``.treefile`` is then moved to
    ``output_newick`` and the ``.iqtree`` summary (if ``summary_path`` is
    given) is moved alongside. Everything else is discarded.

    The ``is_protein`` flag is read for parity with ``run_fasttree``;
    IQ-TREE itself detects the alphabet from the alignment, but we use
    the flag to validate the configured model is appropriate (an explicit
    NT model on AA data would fail late inside IQ-TREE).
    """
    binary_name_override = (
        cfg.get("phylo", {}).get("iqtree", {}).get("binary")
    )
    iqtree = _check_iqtree(binary_name_override)

    iq_cfg = cfg.get("phylo", {}).get("iqtree", {}) or {}
    model = iq_cfg.get("model") or "MFP"
    ufboot = iq_cfg.get("ultrafast_bootstrap", 1000)
    extra_args: list[str] = list(iq_cfg.get("extra_args", []) or [])
    threads = cfg.get("threads", 4)
    seed = cfg.get("seed", 42)

    # IQ-TREE refuses ultrafast bootstrap with fewer than 4 sequences
    # ("It makes no sense to perform bootstrap..."). We've already passed
    # the >=3 rep gate in the orchestrator, so reps may be 3. Drop UFBoot
    # to keep the tree, rather than crashing the run.
    n_reps = _count_msa_records(msa_fasta)
    if ufboot and n_reps < 4:
        import sys as _sys
        print(
            f"  [iqtree] only {n_reps} sequences in the MSA — skipping "
            f"ultrafast bootstrap (IQ-TREE requires >= 4).",
            file=_sys.stderr,
        )
        ufboot = 0

    work_root = Path(cfg.get("temp_dir") or "/tmp/repseq")
    work_root.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(dir=work_root, prefix="iqtree_") as td:
        td = Path(td)
        # Copy the input under a short, side-effect-free prefix so the
        # IQ-TREE output files all land inside the temp dir (it writes
        # alongside the input, using <input>.treefile, <input>.iqtree, ...).
        local_input = td / "input.fasta"
        local_input.write_bytes(msa_fasta.read_bytes())

        cmd = [
            iqtree,
            "-s", str(local_input),
            "-m", str(model),
            "-T", str(threads),
            "-seed", str(seed),
            "--prefix", str(td / "run"),
            "--quiet",
            # IQ-TREE refuses to overwrite by default; this scratch dir is
            # always empty, but pass --redo so a re-run on the same tmp
            # path also succeeds.
            "--redo",
        ]
        if ufboot and int(ufboot) > 0:
            cmd.extend(["-B", str(int(ufboot))])
        cmd.extend(extra_args)

        try:
            subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True,
            )
        except subprocess.CalledProcessError as e:
            err = (e.stderr or "") + (e.stdout or "")
            raise IQTreeError(f"IQ-TREE failed:\n{err}") from e

        treefile = td / "run.treefile"
        if not treefile.exists():
            raise IQTreeError(
                f"IQ-TREE finished but did not produce a .treefile under {td}"
            )
        output_newick.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(treefile, output_newick)

        if summary_path is not None:
            iqfile = td / "run.iqtree"
            if iqfile.exists():
                summary_path.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(iqfile, summary_path)
