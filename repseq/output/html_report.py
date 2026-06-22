"""Single-file HTML run report — a shareable overview of a run.

Bundles, into one ``{prefix}_report.html`` a bench scientist can open in any
browser or e-mail to a collaborator:

* the run's **analysis flags** (from :mod:`repseq.output.flags` — the
  reassortment / non-monophyly / taxonomy-conflict synthesis);
* a **gallery of every tree figure** (the ``*_tree.png`` rendered by the PDF
  sweep), embedded as base64 so the visual payload travels with the file;
* an **index of all output files** (with sizes and relative links) so the full
  machine-readable detail is one click away when the report sits in its
  output directory.

Pure-Python and dependency-free (string templating + ``base64`` +
``html.escape``); a decoupled post-hoc sweep like the PDF / conservation /
monophyly steps, so it never affects the analysis and soft-fails to nothing.
"""

from __future__ import annotations

import base64
import html
import logging
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from .flags import _SECTIONS, Flag, collect_flags

logger = logging.getLogger(__name__)

_CSS = """
body { font-family: -apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif;
       max-width: 1000px; margin: 2rem auto; padding: 0 1rem; color: #222;
       line-height: 1.45; }
h1 { border-bottom: 2px solid #444; padding-bottom: .3rem; }
h2 { margin-top: 2rem; border-bottom: 1px solid #ccc; padding-bottom: .2rem; }
.meta { color: #666; font-size: .9rem; }
.flag { padding: .4rem .6rem; margin: .3rem 0; border-radius: 4px; }
.flag.warn { background: #fff3cd; border-left: 4px solid #e0a800; }
.flag.info { background: #e7f1ff; border-left: 4px solid #4a90d9; }
.clean { background: #e6f4ea; border-left: 4px solid #34a853; padding: .5rem .7rem;
         border-radius: 4px; }
details { margin: .5rem 0; }
summary { cursor: pointer; font-weight: 600; }
img.tree { max-width: 100%; height: auto; border: 1px solid #ddd; margin: .4rem 0; }
table { border-collapse: collapse; width: 100%; font-size: .9rem; }
th, td { border: 1px solid #ddd; padding: .3rem .5rem; text-align: left; }
th { background: #f4f4f4; }
td.size { text-align: right; color: #666; white-space: nowrap; }
"""


def _e(text: object) -> str:
    return html.escape(str(text))


def _human_size(n: int) -> str:
    f = float(n)
    for unit in ("B", "KB", "MB", "GB"):
        if f < 1024 or unit == "GB":
            return f"{f:.0f} {unit}" if unit == "B" else f"{f:.1f} {unit}"
        f /= 1024
    return f"{n} B"


def _flags_html(flags: list[Flag]) -> str:
    if not flags:
        return '<div class="clean">No taxonomy / tree conflicts flagged. ✓</div>'
    parts: list[str] = []
    # Iterate the canonical section order from flags._SECTIONS so the HTML
    # gallery and _flags.txt never drift (e.g. the "Taxa eliminated by QC"
    # section appears in both).
    for category, title in _SECTIONS:
        group = [f for f in flags if f.category == category]
        if not group:
            continue
        parts.append(f"<h3>{_e(title)}</h3>")
        for f in group:
            parts.append(f'<div class="flag {_e(f.severity)}">{_e(f.message)}</div>')
    return "\n".join(parts)


def _tree_gallery_html(out_dir: Path, prefix: str) -> str:
    pngs = sorted(p for p in out_dir.rglob("*_tree.png"))
    if not pngs:
        return "<p class='meta'>No tree figures were rendered this run.</p>"
    parts: list[str] = []
    for png in pngs:
        rel = png.relative_to(out_dir)
        try:
            b64 = base64.b64encode(png.read_bytes()).decode("ascii")
        except OSError:
            continue
        # The whole-genome tree is the headline — show it open, the rest folded.
        is_genome = png.name == f"{prefix}_tree.png"
        open_attr = " open" if is_genome else ""
        parts.append(
            f"<details{open_attr}><summary>{_e(rel)}</summary>"
            f'<img class="tree" alt="{_e(rel)}" '
            f'src="data:image/png;base64,{b64}"></details>'
        )
    return "\n".join(parts)


def _file_index_html(out_dir: Path, prefix: str, self_name: str) -> str:
    rows: list[str] = []
    for p in sorted(out_dir.rglob("*")):
        if not p.is_file():
            continue
        rel = p.relative_to(out_dir)
        if str(rel) == self_name:
            continue  # don't list the report in its own index
        rows.append(
            f"<tr><td><a href='{_e(rel)}'>{_e(rel)}</a></td>"
            f"<td class='size'>{_e(_human_size(p.stat().st_size))}</td></tr>"
        )
    if not rows:
        return ""
    return (
        "<table><tr><th>File</th><th>Size</th></tr>\n"
        + "\n".join(rows)
        + "\n</table>"
    )


def write_html_report(
    out_dir: Path, prefix: str, cfg: Optional[dict] = None
) -> Optional[Path]:
    """Write ``{prefix}_report.html``. Returns the path, or None if the output
    directory has nothing worth bundling (no flags sources and no figures)."""
    try:
        from .. import __version__ as version
    except Exception:  # pragma: no cover - defensive
        version = "?"

    have_flag_sources = (
        (out_dir / f"{prefix}_monophyly.tsv").exists()
        or (out_dir / f"{prefix}_per_protein" / f"{prefix}_incongruence.tsv").exists()
        or (out_dir / f"{prefix}_taxonomy_review.tsv").exists()
    )
    have_figures = any(out_dir.rglob("*_tree.png"))
    # Collect the flags once and reuse for both the qc-drop guard below and the
    # flags section in the body (re-parsing the source TSVs twice for one render
    # is wasted I/O).
    flags = collect_flags(out_dir, prefix)
    # A QC-elimination flag (genus+ wiped out by QC) is worth a report even on
    # a plain clustering run with no conflict tables or tree figures.
    have_qc_drop = any(f.category == "qc_drop" for f in flags)
    if not have_flag_sources and not have_figures and not have_qc_drop:
        return None

    self_name = f"{prefix}_report.html"
    when = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    summary_link = (
        f"<p>Full Methods-style write-up: "
        f"<a href='{prefix}_summary.md'>{prefix}_summary.md</a></p>"
        if (out_dir / f"{prefix}_summary.md").exists() else ""
    )

    body = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<title>{_e(prefix)} — repseq report</title>
<style>{_CSS}</style></head><body>
<h1>{_e(prefix)} — repseq run report</h1>
<p class="meta">repseq {_e(version)} · generated {_e(when)}</p>
{summary_link}

<h2>Analysis flags</h2>
{_flags_html(flags)}

<h2>Tree figures</h2>
{_tree_gallery_html(out_dir, prefix)}

<h2>Output files</h2>
{_file_index_html(out_dir, prefix, self_name)}
</body></html>
"""
    path = out_dir / self_name
    path.write_text(body)
    return path
