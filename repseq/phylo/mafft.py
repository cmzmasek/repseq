"""MAFFT MSA wrapper.

Shells out to ``mafft`` and writes the aligned FASTA. Default mode is
``--auto`` which picks the algorithm by input size (FFT-NS-2 for large
inputs, L-INS-i for small, etc.) — fine for most cases. Power users can
override via ``phylo.mafft.extra_args``. The per-protein-tree path
(2F) passes its own ``extra_args`` (default L-INS-i:
``--maxiterate 1000 --localpair``) with ``use_auto=False``, since those
single-gene alignments are small enough to afford high accuracy.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any


class MafftError(RuntimeError):
    pass


def _check_mafft() -> str:
    path = shutil.which("mafft")
    if not path:
        raise MafftError(
            "mafft not found in PATH. Install it from https://mafft.cbrc.jp/alignment/software/"
        )
    return path


def tool_version() -> str:
    """Return MAFFT's version string, or ``"unknown"`` if it cannot be
    determined.

    Used by the phyloXML writer to annotate the tree's
    ``<phylogeny><description>`` with the alignment-tool provenance.
    MAFFT prints its version to stderr (e.g. ``v7.520 (2023/Mar/16)``)
    and exits non-zero on ``--version``, so we capture stderr and
    ignore the exit code.
    """
    try:
        path = _check_mafft()
    except MafftError:
        return "unknown"
    try:
        result = subprocess.run(
            [path, "--version"], capture_output=True, text=True, timeout=5,
        )
        out = (result.stderr or result.stdout or "").strip()
        return out.splitlines()[0].strip() if out else "unknown"
    except (subprocess.TimeoutExpired, OSError):
        return "unknown"


def run_mafft(
    input_fasta: Path,
    output_fasta: Path,
    cfg: dict[str, Any],
    *,
    extra_args: list[str] | None = None,
    use_auto: bool = True,
) -> None:
    """Align ``input_fasta`` with MAFFT, writing the result to ``output_fasta``.

    MAFFT prints the alignment to stdout, so we capture it directly into
    the output file rather than parsing structured output.

    ``extra_args`` overrides ``phylo.mafft.extra_args`` when given (the
    per-protein path passes its own L-INS-i args). ``use_auto`` controls
    whether ``--auto`` is prepended: it must be **False** when the caller
    supplies an explicit pairwise strategy (``--localpair`` /
    ``--globalpair`` / ``--genafpair``), because ``--auto`` overrides
    those and the chosen strategy would silently not take effect.
    """
    mafft = _check_mafft()
    phylo_cfg = cfg.get("phylo", {}) or {}
    mafft_cfg = phylo_cfg.get("mafft", {}) or {}
    threads = cfg.get("threads", 4)
    if extra_args is None:
        extra_args = list(mafft_cfg.get("extra_args", []) or [])

    cmd = [mafft]
    if use_auto:
        cmd.append("--auto")
    cmd.extend(["--thread", str(threads)])
    cmd.extend(extra_args)
    cmd.append(str(input_fasta))

    # Bench-scientist progress: the MSA step can run for minutes on a
    # large input, and a silent terminal makes the user wonder if the
    # pipeline froze. Echo the args (without the binary path or the
    # input file) before the run, plus elapsed time on success.
    display_args = " ".join(cmd[1:-1])
    print(f"[phylo] starting MAFFT ({display_args})", file=sys.stderr)
    t0 = time.time()

    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    with open(output_fasta, "w") as fh:
        try:
            subprocess.run(
                cmd,
                check=True,
                stdout=fh,
                stderr=subprocess.PIPE,
                text=True,
            )
        except subprocess.CalledProcessError as e:
            raise MafftError(f"mafft failed:\n{e.stderr}") from e

    print(
        f"[phylo] MAFFT finished ({time.time() - t0:.1f}s)",
        file=sys.stderr,
    )
