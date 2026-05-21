"""trimAl wrapper — optional alignment trimming between MAFFT and the tree.

trimAl (Capella-Gutiérrez et al. 2009) removes poorly-aligned / gap-rich
columns from an MSA before tree inference. It is **off by default** and
configured independently for the whole-genome tree (``phylo.trimal``) and
the per-protein trees (``phylo.per_protein.trimal``); the default method is
``-automated1`` (trimAl's heuristic that picks the best column-trimming for
maximum-likelihood reconstruction).

The public entry point is :func:`maybe_trim`, which embodies the soft-fail
contract the user chose: if trimAl is disabled it does nothing; if enabled
but the binary is missing, the run fails, or the trimmed alignment is
degenerate, it emits a **loud** stderr warning and returns ``False`` so the
caller builds the tree on the UNTRIMMED alignment rather than losing the
tree. It only trims columns (``-automated1`` never drops sequences), so the
leaf set / id-map stay valid in every path.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Optional

# trimAl's column-trimming heuristic methods (no sequence removal). Threshold
# trimming (-gt/-st/-cons) is reachable via the `extra_args` escape hatch.
KNOWN_MODES = (
    "automated1", "gappyout", "strict", "strictplus", "nogaps", "noallgaps",
)

# A trimmed alignment narrower than this is treated as degenerate (trimAl
# stripped essentially everything) → fall back to the untrimmed alignment.
_MIN_COLS = 1


class TrimalError(RuntimeError):
    pass


def _check_trimal(override: Optional[str] = None) -> str:
    """Locate the trimal binary, or raise :class:`TrimalError`."""
    for name in ([override] if override else []) + ["trimal"]:
        if not name:
            continue
        path = shutil.which(name)
        if path:
            return path
    raise TrimalError(
        "trimal not found in PATH. Install it from "
        "http://trimal.cgenomics.org/ (conda: `conda install -c bioconda trimal`)."
    )


def tool_version(override: Optional[str] = None) -> str:
    """trimAl version string, or ``"unknown"``. Used by the phyloXML
    description and ``repseq doctor``. trimAl prints e.g.
    ``trimAl v1.4.rev15 build[2013-12-17] ...``."""
    try:
        path = _check_trimal(override)
    except TrimalError:
        return "unknown"
    try:
        result = subprocess.run(
            [path, "--version"], capture_output=True, text=True, timeout=5,
        )
        out = (result.stdout or result.stderr or "").strip()
        for line in out.splitlines():
            if line.strip():
                return line.strip()
        return "unknown"
    except (subprocess.TimeoutExpired, OSError):
        return "unknown"


def _alignment_width(path: Path) -> int:
    """Number of columns in a FASTA alignment (length of the first record)."""
    width = 0
    with open(path) as fh:
        seen_header = False
        for line in fh:
            if line.startswith(">"):
                if seen_header:
                    break
                seen_header = True
            elif seen_header:
                width += len(line.strip())
    return width


def trim_note(settings: Optional[dict[str, Any]], binary: Optional[str] = None) -> str:
    """A short provenance string for the phyloXML description, e.g.
    ``"trimAl v1.4.rev15 -automated1"``."""
    mode = (settings or {}).get("mode", "automated1")
    extra = list((settings or {}).get("extra_args", []) or [])
    flags = (f"-{mode}" if mode else "") + (" " + " ".join(extra) if extra else "")
    return f"{tool_version(binary)} {flags}".strip()


def maybe_trim(
    input_fasta: Path,
    output_fasta: Path,
    cfg: dict[str, Any],
    settings: Optional[dict[str, Any]],
    *,
    label: str = "",
) -> bool:
    """Trim ``input_fasta`` → ``output_fasta`` if ``settings.enabled``.

    Returns ``True`` when ``output_fasta`` is a usable trimmed alignment;
    ``False`` (with a loud stderr warning) when trimming was disabled,
    skipped (binary missing), failed, or produced a degenerate alignment —
    in which case the caller builds the tree on the untrimmed input. Never
    raises: a trimming problem must never destroy an otherwise-good tree.
    """
    if not settings or not settings.get("enabled"):
        return False

    mode = settings.get("mode", "automated1")
    extra = list(settings.get("extra_args", []) or [])
    binary = settings.get("binary")
    where = f" [{label}]" if label else ""

    try:
        trimal = _check_trimal(binary)
    except TrimalError as exc:
        print(
            f"[phylo] trimal requested{where} but unavailable: {exc} "
            f"Building the tree on the UNTRIMMED alignment.",
            file=sys.stderr,
        )
        return False

    cmd = [trimal, "-in", str(input_fasta), "-out", str(output_fasta)]
    if mode:
        cmd.append(f"-{mode}")
    cmd.extend(extra)
    print(
        f"[phylo] trimming alignment{where} with trimal "
        f"({('-' + mode) if mode else ''}{(' ' + ' '.join(extra)) if extra else ''})",
        file=sys.stderr,
    )
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as exc:
        err = ((exc.stderr or "") + (exc.stdout or "")).strip()
        print(
            f"[phylo] trimal failed{where} — using the UNTRIMMED alignment:\n{err}",
            file=sys.stderr,
        )
        return False

    # Surface stderr on the sanity-check failure rather than swallowing it
    # (trimAl can exit 0 yet warn / strip everything).
    if not output_fasta.exists() or _alignment_width(output_fasta) < _MIN_COLS:
        err = (result.stderr or "").strip()
        print(
            f"[phylo] trimal{where} produced a degenerate alignment "
            f"(< {_MIN_COLS} columns) — using the UNTRIMMED alignment."
            + (f"\n{err}" if err else ""),
            file=sys.stderr,
        )
        return False
    return True
