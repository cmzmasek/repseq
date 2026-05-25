"""FastTree wrapper.

Infers an approximate-ML phylogeny from a MAFFT alignment. Auto-picks the
model from the sequence alphabet — FastTree's default is protein (JTT);
for nucleotide we pass ``-nt -gtr``. The binary is conventionally named
``FastTree`` (capital F) on macOS/conda but ``fasttree`` (lower-case) on
some Linux packages, so we look for both.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

from ..utils.subprocess_stream import StreamedProcessError, run_streaming


class FastTreeError(RuntimeError):
    pass


def _check_fasttree() -> str:
    for name in ("FastTree", "fasttree"):
        path = shutil.which(name)
        if path:
            return path
    raise FastTreeError(
        "FastTree not found in PATH (tried 'FastTree' and 'fasttree'). "
        "Install it from http://www.microbesonline.org/fasttree/"
    )


def tool_version() -> str:
    """Return FastTree's version string, or ``"unknown"`` on any failure.

    FastTree prints a banner including the version when invoked with no
    input (the first line looks like
    ``FastTree Version 2.1.11 Double precision (No SSE3)``); it then
    exits non-zero waiting for stdin. We capture stderr and ignore the
    exit code.
    """
    try:
        path = _check_fasttree()
    except FastTreeError:
        return "unknown"
    try:
        result = subprocess.run(
            [path], input="", capture_output=True, text=True, timeout=5,
        )
        out = (result.stderr or result.stdout or "").strip()
        if not out:
            return "unknown"
        for line in out.splitlines():
            if "version" in line.lower():
                return line.strip()
        return out.splitlines()[0].strip()
    except (subprocess.TimeoutExpired, OSError):
        return "unknown"


def run_fasttree(
    msa_fasta: Path,
    output_newick: Path,
    cfg: dict[str, Any],
    is_protein: bool,
) -> None:
    """Run FastTree on ``msa_fasta`` and write a Newick tree to ``output_newick``.

    FastTree prints the tree to stdout; we capture it straight into the
    output file. ``is_protein`` selects the substitution model: default
    JTT for protein, ``-gtr`` for nucleotide.
    """
    fasttree = _check_fasttree()
    phylo_cfg = cfg.get("phylo", {}) or {}
    ft_cfg = phylo_cfg.get("fasttree", {}) or {}
    extra_args: list[str] = list(ft_cfg.get("extra_args", []) or [])

    cmd = [fasttree]
    if not is_protein:
        cmd.extend(["-nt", "-gtr"])
    cmd.extend(extra_args)
    cmd.append(str(msa_fasta))

    # Bench-scientist progress message: echo args (without binary path
    # or input file) before the run, plus elapsed time on success.
    display_args = " ".join(cmd[1:-1])
    print(f"[phylo] starting FastTree ({display_args})", file=sys.stderr)
    t0 = time.time()

    output_newick.parent.mkdir(parents=True, exist_ok=True)
    # Stream FastTree's progress lines (it prints "Iteration N" / "ML…"
    # to stderr) so a long run shows life.
    try:
        run_streaming(
            cmd,
            stdout_file=output_newick,
            stream_prefix="[fasttree] ",
        )
    except StreamedProcessError as e:
        raise FastTreeError(f"FastTree failed:\n{e.stderr}") from e

    print(
        f"[phylo] FastTree finished ({time.time() - t0:.1f}s)",
        file=sys.stderr,
    )
