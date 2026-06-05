"""Render a graphical PDF + PNG of every phyloXML tree built this run (2H).

This is the sibling of the plain-text Newick / annotated phyloXML outputs:
for every ``{prefix}..._tree.xml`` a phylo step writes, draw a ladderized
rectangular phylogram and save it as ``{prefix}..._tree.pdf`` (vector) plus
``{prefix}..._tree.png`` (150 dpi raster). The look is ported from the sister
project ``vfam_trees`` (matplotlib + ``Bio.Phylo.draw``): taxonomy-coloured
leaf labels, a genus/subfamily colour legend, internal-node labels for the
genus→family ranks, and branch-support labels for nodes with confidence ≥ 50.

The renderer is fully decoupled from the tree builders — it reconstructs
**everything** (leaf labels, leaf colours, the genus legend, internal ranks,
support values) from the phyloXML itself, which already carries it all:

* leaf ``<name>``                                   → leaf display label
* ``<property ref="style:font_color">``             → leaf label colour
* ``<property ref="repseq:genus" / ":subfamily">``  → legend reconstruction
* internal ``<taxonomy><rank>`` + ``<name>``        → internal-node labels
* ``<confidence>``                                  → branch-support labels

So a single end-of-run sweep over the tracked ``*_tree.xml`` files renders
figures for every tree type — whole-genome (2E), per-protein / extra / segment
(2F), pre-cluster, and partition — with no per-builder threading.

Hard-requires only matplotlib (the ``[viz]`` extra) and biopython (already a
core dependency). Missing matplotlib soft-fails with a single message and
emits no figures — the same posture as ``--plot``. A single tree failing to
render warns and is skipped; the rest still draw.
"""
from __future__ import annotations

import importlib
from pathlib import Path
from typing import Optional


# Internal-node labels are shown only for these ranks (species-level and
# unranked internal annotations are suppressed to cut clutter; leaf labels
# are always shown in full). Mirrors vfam_trees' _SHOW_INTERNAL_RANKS.
_SHOW_INTERNAL_RANKS = frozenset({"genus", "subgenus", "subfamily", "family"})

_BRANCH_LINEWIDTH = 0.5
_PNG_DPI = 150
# Confidence (branch-support) values below this are not drawn — low support
# would only clutter the figure.
_MIN_SHOWN_CONFIDENCE = 50


# ---------------------------------------------------------------------------
# Dependency check (mirrors viz/clustering_plot._require_matplotlib)
# ---------------------------------------------------------------------------

def matplotlib_unavailable_reason() -> Optional[str]:
    """Return a plain-English reason matplotlib can't be used, or None.

    Distinguishes "not installed" from "installed but failing to import"
    (almost always a NumPy/SciPy clash in the environment, not a repseq
    problem). Returns None when matplotlib imports cleanly.
    """
    try:
        importlib.import_module("matplotlib")
    except ModuleNotFoundError as exc:
        if exc.name == "matplotlib":
            return (
                "tree PDF rendering requires matplotlib — install it with: "
                "pip install 'repseq[viz]'"
            )
        return (
            f"matplotlib is installed but a dependency is missing ({exc}); "
            "the cleanest fix is a dedicated environment for repseq"
        )
    except ImportError as exc:
        return (
            f"matplotlib is installed but failed to import [{exc}] — this is "
            "almost always a NumPy/SciPy version mismatch in the environment, "
            "not a repseq problem; the cleanest fix is a dedicated "
            "environment for repseq"
        )
    return None


# ---------------------------------------------------------------------------
# phyloXML clade readers (Bio.Phylo PhyloXML objects)
# ---------------------------------------------------------------------------

def _clade_property(clade, ref: str) -> Optional[str]:
    """Return the value of the first ``<property ref=...>`` on a clade, or None."""
    for prop in getattr(clade, "properties", None) or []:
        if getattr(prop, "ref", None) == ref:
            value = getattr(prop, "value", None)
            return str(value) if value is not None else None
    return None


def _clade_rank(clade) -> str:
    """Return the taxonomy rank stored on a clade's first ``<taxonomy>``, or ''."""
    taxes = getattr(clade, "taxonomies", None) or []
    if taxes:
        return (getattr(taxes[0], "rank", None) or "").strip()
    return ""


def _confidence_value(clade) -> Optional[float]:
    """Return a clade's numeric branch-support value, or None.

    ``Clade.confidence`` is a plain float on Newick-parsed trees but a property
    (or a ``Confidence`` object) on phyloXML-parsed trees, so read defensively.
    """
    try:
        c = clade.confidence
    except (ValueError, AttributeError):
        c = None
    if c is None:
        confs = getattr(clade, "confidences", None) or []
        if confs:
            c = getattr(confs[0], "value", None)
    if hasattr(c, "value"):
        c = c.value
    return c if isinstance(c, (int, float)) else None


def _internal_label(clade) -> str:
    """label_func for Phylo.draw: full name on leaves, rank-gated on internals."""
    if clade.is_terminal():
        return clade.name or ""
    if _clade_rank(clade) in _SHOW_INTERNAL_RANKS:
        return clade.name or ""
    return ""


def _build_color_maps(tree):
    """Reconstruct leaf-label and legend colours from the phyloXML leaves.

    Returns ``(display_to_color, genus_to_color, subfamily_to_genera)`` mirroring
    the dicts vfam_trees computes separately — here derived from each leaf's
    ``style:font_color`` plus its ``repseq:genus`` / ``repseq:subfamily``
    properties.
    """
    display_to_color: dict[str, str] = {}
    genus_to_color: dict[str, str] = {}
    subfamily_to_genera: dict[str, set] = {}

    for leaf in tree.get_terminals():
        color = _clade_property(leaf, "style:font_color")
        name = (leaf.name or "").strip()
        if name and color:
            display_to_color[name] = color
        genus = (_clade_property(leaf, "repseq:genus") or "").strip()
        if genus and color:
            genus_to_color.setdefault(genus, color)
            subfamily = (_clade_property(leaf, "repseq:subfamily") or "").strip()
            if subfamily:
                subfamily_to_genera.setdefault(subfamily, set()).add(genus)

    return display_to_color, genus_to_color, subfamily_to_genera


# ---------------------------------------------------------------------------
# Drawing (ported from vfam_trees/report.py)
# ---------------------------------------------------------------------------

def _wrap(text: str, width: int) -> str:
    """Wrap a long string at word boundaries."""
    if len(text) <= width:
        return text
    words = text.split()
    lines: list[str] = []
    current = ""
    for word in words:
        if len(current) + len(word) + 1 > width:
            lines.append(current)
            current = word
        else:
            current = f"{current} {word}".strip()
    if current:
        lines.append(current)
    return "\n".join(lines)


def _place_figure_header(fig, title: str, info: str) -> float:
    """Place a bold title + optional grey caption at the top of *fig*.

    Both are written with ``fig.text`` in figure coordinates (so they sit
    close together, unlike ``suptitle`` + ``set_title``). Returns the
    figure-fraction y below which the axes should start. Ported from
    vfam_trees.
    """
    pts_per_frac = fig.get_figheight() * 72.0

    TITLE_Y = 0.99
    TITLE_SIZE = 11
    INFO_SIZE = 8
    INFO_LS = 1.4
    GAP_PTS = 3

    title_frac = TITLE_SIZE / pts_per_frac
    fig.text(0.5, TITLE_Y, title, ha="center", va="top",
             fontsize=TITLE_SIZE, fontweight="bold")

    if info:
        n_lines = info.count("\n") + 1
        info_y = TITLE_Y - title_frac - GAP_PTS / pts_per_frac
        info_frac = INFO_SIZE * INFO_LS * n_lines / pts_per_frac
        fig.text(0.5, info_y, info, ha="center", va="top",
                 fontsize=INFO_SIZE, color="#555555", linespacing=INFO_LS)
        axes_top = info_y - info_frac - GAP_PTS / pts_per_frac
    else:
        axes_top = TITLE_Y - title_frac - GAP_PTS / pts_per_frac

    return float(max(0.50, min(0.97, axes_top)))


def _thin_tree_lines(ax, linewidth: float) -> None:
    """Set an absolute linewidth on every branch line (Phylo.draw uses ~1.5)."""
    for line in ax.get_lines():
        line.set_linewidth(linewidth)


def _draw_taxonomy_legend(ax, genus_to_color: dict, subfamily_to_genera: dict) -> None:
    """Add a genus/subfamily colour legend to the axes. Ported from vfam_trees."""
    import matplotlib.patches as mpatches

    handles = []
    placed: set = set()
    if subfamily_to_genera:
        for subfamily, genera in sorted(subfamily_to_genera.items()):
            handles.append(mpatches.Patch(color="none", label=f"  {subfamily}"))
            for genus in sorted(genera):
                color = genus_to_color.get(genus, "#888888")
                handles.append(mpatches.Patch(color=color, label=f"    {genus}"))
                placed.add(genus)
    for genus, color in sorted(genus_to_color.items()):
        if genus not in placed:
            handles.append(mpatches.Patch(color=color, label=genus))

    if not handles:
        return

    ncol = max(1, len(handles) // 30)
    ax.legend(
        handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1.0),
        fontsize=6, frameon=True, framealpha=0.8, ncol=ncol,
        title="Genus (by subfamily)", title_fontsize=7,
        borderpad=0.5, handlelength=1.0, handletextpad=0.4,
    )


def _draw_tree_figure(tree, title_name: str, info: str):
    """Return a matplotlib Figure of *tree*, or None on error. Ported from vfam_trees."""
    import copy
    import matplotlib.pyplot as plt
    import matplotlib.text as _mtext
    from Bio import Phylo

    # Work on a copy — ladderize() would otherwise mutate the parsed tree.
    tree = copy.deepcopy(tree)
    n_leaves = sum(1 for _ in tree.get_terminals())
    if n_leaves < 2:
        return None

    display_to_color, genus_to_color, subfamily_to_genera = _build_color_maps(tree)

    fig_h = max(6, n_leaves * 0.18)
    has_legend = bool(genus_to_color)
    fig_w = 14 if has_legend else 11
    fig, ax = plt.subplots(figsize=(fig_w, min(fig_h, 24)))

    def _confidence_label(clade):
        conf = _confidence_value(clade)
        if conf is None or conf < _MIN_SHOWN_CONFIDENCE:
            return None
        return str(int(conf))

    tree.ladderize(reverse=True)
    with plt.rc_context({"lines.linewidth": _BRANCH_LINEWIDTH}):
        Phylo.draw(
            tree, axes=ax, do_show=False,
            label_func=_internal_label,
            branch_labels=_confidence_label,
        )
    ax.set_title("")  # Phylo.draw sets ax.title = tree.name; clear for our header
    _thin_tree_lines(ax, _BRANCH_LINEWIDTH)
    ax.axis("off")

    axes_top = _place_figure_header(
        fig, f"{title_name} ({n_leaves} external nodes)", info
    )
    fig.subplots_adjust(top=axes_top)

    font_size = max(4, min(8, int(200 / max(n_leaves, 1))))
    for artist in ax.get_children():
        if not isinstance(artist, _mtext.Text) or artist is ax.title:
            continue
        artist.set_fontsize(font_size)
        color = display_to_color.get(artist.get_text().strip())
        if color:
            artist.set_color(color)

    if has_legend:
        _draw_taxonomy_legend(ax, genus_to_color, subfamily_to_genera)

    return fig


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def render_one(xml_path: Path, want_png: bool = True) -> list[Path]:
    """Render one phyloXML tree to a sibling PDF (+ PNG). Returns created paths.

    Raises on a hard failure (unreadable XML, draw error); the caller decides
    whether one tree failing aborts the sweep (it does not — see
    ``render_tree_pdfs``).
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from Bio import Phylo

    xml_path = Path(xml_path)
    tree = Phylo.read(str(xml_path), "phyloxml")

    title_name = (getattr(tree, "name", None) or xml_path.stem).strip()
    description = (getattr(tree, "description", None) or "").strip()
    info = _wrap(description, 110) if description else ""

    fig = _draw_tree_figure(tree, title_name, info)
    if fig is None:
        return []

    created: list[Path] = []
    try:
        targets = [(xml_path.with_suffix(".pdf"), {})]
        if want_png:
            targets.append((xml_path.with_suffix(".png"), {"dpi": _PNG_DPI}))
        for path, kwargs in targets:
            fig.savefig(str(path), bbox_inches="tight", **kwargs)
            created.append(path)
    finally:
        plt.close(fig)
    return created


def render_tree_pdfs(xml_paths, want_png: bool = True):
    """Render PDFs/PNGs for a batch of phyloXML trees.

    Checks matplotlib once. Returns ``(created, skipped_reason, failures)``:
    ``created`` is the flat list of written PDF/PNG paths, ``skipped_reason``
    is a single message when matplotlib is unavailable (nothing rendered), and
    ``failures`` is a list of ``(xml_path, error_message)`` for individual
    trees that errored (the rest still render). The caller handles all
    user-facing reporting.
    """
    xml_paths = [Path(p) for p in xml_paths]
    if not xml_paths:
        return [], None, []

    reason = matplotlib_unavailable_reason()
    if reason is not None:
        return [], reason, []

    created: list[Path] = []
    failures: list[tuple] = []
    for xml_path in xml_paths:
        try:
            created.extend(render_one(xml_path, want_png=want_png))
        except Exception as exc:  # noqa: BLE001 — one bad tree must not abort the rest
            failures.append((xml_path, str(exc)))
    return created, None, failures
