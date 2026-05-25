"""Per-marker conservation heatmaps.

One PNG per declared marker spec, written to ``{prefix}_conservation/`` —
the visual companion to the per-protein trees (2F) and per-protein
FASTAs. Two stacked tracks on top, one underneath:

1. **Shannon entropy** per column, in bits — high = variable, low =
   conserved. Gaps are treated as missing data and excluded from the
   per-column counts (otherwise the metric collapses to "how
   gappy is this column").
2. **Fraction matching consensus** per column — fraction of non-gap
   rows whose residue equals the column's mode. Complementary to
   entropy: entropy weighs *all* residues by frequency, %-consensus
   answers "how often does the most common residue dominate?".
3. **Domain architecture** ribbon — coloured boxes drawn from the HMM
   hits on the *longest satisfying CDS across all reps*, projected
   from ungapped CDS coordinates into MSA-column coordinates so the
   box edges land on the columns the user sees in the heatmap above.

Heatmap cells use a sequential ``viridis`` colormap (continuous
magnitude data); domain boxes are coloured by family with a stable
golden-angle hue (the same HSV scheme that drives the taxonomy palette
in ``phylo/coloring.py``), so a marker that is "the green tree" in 2F
is also "the green domain ribbon" here. Tasks that aren't matplotlib's
strong suit (computing per-column metrics, mapping coordinates) are
done in pure Python so the math is unit-testable without spinning up
the renderer.

Triggered by ``--conservation-heatmap`` (independent of ``--phylo`` and
``--per-protein-phylo``), but **requires ``--per-protein-phylo`` to
have run first** so the per-family MSA is already on disk under
``{prefix}_per_protein/<family>_msa.fasta``. Without that MSA the step
soft-fails with a stderr note rather than triggering a fresh MAFFT run
(MAFFT is the slow part — making conservation a free piggyback on 2F
is the whole point).
"""

from __future__ import annotations

import colorsys
import importlib
import logging
import math
from pathlib import Path
from typing import Any, Optional

logger = logging.getLogger(__name__)

# Match the HSV defaults used by phylo/coloring.py so family-coloured
# boxes here look visually continuous with the taxonomy-coloured leaves
# in the matching 2F tree.
_HSV_SATURATION = 0.65
_HSV_VALUE = 0.90
_GOLDEN_ANGLE = 137.5077640500378  # 360 / phi^2

# Default smoothing window for the entropy / consensus line charts.
# 15 residues is a single-residue-resolution-smoothing default — wide
# enough to suppress single-column spikes (especially from columns with
# few non-gap rows) but narrow enough to keep real conservation peaks
# in their actual structural location. Override via the
# ``window`` kwarg to :func:`write_conservation_heatmap`.
DEFAULT_WINDOW_AA = 15


def _require_matplotlib() -> None:
    """Import matplotlib eagerly so a missing/broken install fails fast
    with a plain-English next step (same pattern as ``--plot``)."""
    try:
        importlib.import_module("matplotlib")
    except ModuleNotFoundError as exc:
        raise ImportError(
            "Conservation heatmaps require matplotlib — install with: "
            "pip install 'repseq[viz]'"
        ) from exc


# ---------------------------------------------------------------------------
# MSA I/O
# ---------------------------------------------------------------------------

def read_msa(path: Path) -> dict[str, str]:
    """Parse a FASTA MSA into ``{first_header_token: aligned_seq}``.

    Mirrors :func:`repseq.phylo.partition.read_msa` so the conservation
    step reads the *same* records out of a 2F per-family MSA as the
    partitioned-supermatrix builder would. Kept as a local copy to
    avoid pulling the whole partition module into the viz path.
    """
    records: dict[str, str] = {}
    cur: Optional[str] = None
    buf: list[str] = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None:
                    records[cur] = "".join(buf)
                cur = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if cur is not None:
        records[cur] = "".join(buf)
    return records


# ---------------------------------------------------------------------------
# Per-column metrics — pure Python, gap-aware, testable in isolation
# ---------------------------------------------------------------------------

_GAP_CHARS = frozenset("-.")


def _columns(rows: list[str]) -> int:
    """Number of MSA columns (the common length; rows are assumed
    aligned, so any row's length is the truth)."""
    return len(rows[0]) if rows else 0


def compute_shannon_entropy(rows: list[str]) -> list[float]:
    """Per-column Shannon entropy in bits, gaps excluded.

    ``H = -Σ p_i log₂ p_i`` over the residue frequencies in the column
    (each residue uppercased; gap characters ``-`` and ``.`` dropped
    before the frequency calculation). An all-gap column or an empty
    MSA contributes ``0.0`` rather than NaN — the renderer expects a
    finite value per column.
    """
    n_cols = _columns(rows)
    out: list[float] = [0.0] * n_cols
    for col in range(n_cols):
        counts: dict[str, int] = {}
        total = 0
        for row in rows:
            ch = row[col].upper()
            if ch in _GAP_CHARS:
                continue
            counts[ch] = counts.get(ch, 0) + 1
            total += 1
        if total == 0:
            continue
        h = 0.0
        for c in counts.values():
            p = c / total
            if p > 0:
                h -= p * math.log2(p)
        out[col] = h
    return out


def sliding_window_mean(values: list[float], window: int) -> list[float]:
    """Centered sliding-window mean over ``values``.

    For each position ``i``, returns the mean of
    ``values[max(0, i-w//2) : min(n, i+w//2+1)]``. The window shrinks
    at the edges rather than zero-padding, so the smoothed series
    starts and ends at meaningful values rather than dipping to zero.
    ``window <= 1`` returns a copy of the input (identity smoothing).
    An empty input returns an empty list.
    """
    n = len(values)
    if n == 0 or window <= 1:
        return list(values)
    half = window // 2
    out: list[float] = []
    for i in range(n):
        lo = max(0, i - half)
        hi = min(n, i + half + 1)
        slc = values[lo:hi]
        out.append(sum(slc) / len(slc) if slc else 0.0)
    return out


def compute_consensus_fraction(rows: list[str]) -> list[float]:
    """Per-column fraction of non-gap rows that carry the column's mode.

    Complement to Shannon entropy: a column where 99/100 rows show ``M``
    and one shows ``L`` has *low* entropy (residues are nearly
    monomorphic) but the entropy number alone doesn't say *what*
    fraction matched. This metric does.

    All-gap columns return ``0.0`` (no consensus residue to match).
    """
    n_cols = _columns(rows)
    out: list[float] = [0.0] * n_cols
    for col in range(n_cols):
        counts: dict[str, int] = {}
        total = 0
        for row in rows:
            ch = row[col].upper()
            if ch in _GAP_CHARS:
                continue
            counts[ch] = counts.get(ch, 0) + 1
            total += 1
        if total == 0:
            continue
        out[col] = max(counts.values()) / total
    return out


# ---------------------------------------------------------------------------
# Reference-row pick + ungapped→aligned coordinate projection
# ---------------------------------------------------------------------------

def pick_reference_row(msa: dict[str, str]) -> Optional[str]:
    """Return the short_id of the row with the **most non-gap columns**.

    Per the design spec: "Longest satisfying CDS across all reps".
    Tie-broken alphabetically by short_id for determinism. Returns
    ``None`` on an empty MSA.
    """
    if not msa:
        return None
    best_id: Optional[str] = None
    best_score = -1
    for short_id in sorted(msa.keys()):
        non_gap = sum(1 for ch in msa[short_id] if ch not in _GAP_CHARS)
        if non_gap > best_score:
            best_score = non_gap
            best_id = short_id
    return best_id


def project_hits_to_alignment(
    aligned_seq: str, hits: list[dict],
) -> list[dict]:
    """Map HMM hits from ungapped CDS coords into MSA-column coords.

    ``hits`` is the list as stored in ``protein["hmm_hits"]`` — each dict
    carries 1-based ``ali_from`` / ``ali_to`` in the **ungapped** CDS.
    The aligned sequence carries the same residues with gaps inserted;
    we walk it once to build the ungapped→aligned column map and then
    project every hit. Hits that fall outside the ungapped sequence
    (which shouldn't happen but defensively guard against) are dropped.

    Returns hit dicts shallow-copied with ``aln_from`` / ``aln_to``
    populated (1-based, inclusive). The original ``ali_from`` / ``ali_to``
    are preserved on the copy.
    """
    # 1-based mapping from ungapped position → aligned column.
    # ungapped[0] is unused; ungapped[k] for 1 <= k <= len(non-gap) maps
    # to the 1-based column index in the aligned sequence.
    ungapped_to_aligned: list[int] = [0]
    for i, ch in enumerate(aligned_seq, start=1):
        if ch not in _GAP_CHARS:
            ungapped_to_aligned.append(i)
    max_ungapped = len(ungapped_to_aligned) - 1

    projected: list[dict] = []
    for hit in hits or []:
        af = hit.get("ali_from")
        at = hit.get("ali_to")
        if af is None or at is None:
            continue
        try:
            af = int(af)
            at = int(at)
        except (TypeError, ValueError):
            continue
        if af < 1 or at < af:
            continue
        if af > max_ungapped:
            # Hit lies past the reference CDS — can't project (would
            # require choosing a different reference row).
            continue
        # Clamp at_to_unkapped to the reference's last non-gap column
        # rather than dropping the hit when it slightly overshoots a
        # truncated reference (cosmetic, not load-bearing).
        at_eff = min(at, max_ungapped)
        copy = dict(hit)
        copy["aln_from"] = ungapped_to_aligned[af]
        copy["aln_to"] = ungapped_to_aligned[at_eff]
        projected.append(copy)
    return projected


# ---------------------------------------------------------------------------
# Family colour — golden-angle hue for stable, distinct boxes
# ---------------------------------------------------------------------------

def family_color_hex(family_label: str, families: list[str]) -> str:
    """Return a stable ``#RRGGBB`` for ``family_label``.

    Hue is the golden-angle multiple of the family's position in the
    *sorted* families list, so a given family's colour is the same
    across runs as long as the marker config is unchanged. Saturation
    and value match the taxonomy palette defaults in
    :mod:`repseq.phylo.coloring`, so a family-coloured ribbon here
    visually rhymes with the taxonomy-coloured leaves in its 2F tree.
    Unknown families (not in ``families``) get position 0.
    """
    try:
        idx = sorted(families).index(family_label)
    except ValueError:
        idx = 0
    hue_deg = (idx * _GOLDEN_ANGLE) % 360.0
    r, g, b = colorsys.hsv_to_rgb(hue_deg / 360.0, _HSV_SATURATION, _HSV_VALUE)
    return "#{:02X}{:02X}{:02X}".format(
        int(round(r * 255)), int(round(g * 255)), int(round(b * 255)),
    )


# ---------------------------------------------------------------------------
# PNG writer
# ---------------------------------------------------------------------------

def write_conservation_heatmap(
    msa_path: Path,
    *,
    out_png: Path,
    family_label: str,
    family_color: str,
    hmm_hits_on_reference: Optional[list[dict]] = None,
    title_suffix: Optional[str] = None,
    window: int = DEFAULT_WINDOW_AA,
) -> Optional[Path]:
    """Render the three-track conservation figure to ``out_png``.

    Two metric tracks (Shannon entropy + fraction matching consensus)
    drawn as line charts smoothed by a ``window``-residue centered
    sliding mean, over a domain-architecture ribbon with one labelled
    rectangle per HMM hit on the reference CDS. The line-chart
    rendering supersedes the v0.29.0 viridis heatmap: at MSA scales
    of hundreds-to-thousands of columns the heatmap reduces to a
    near-uniform smear, while a line chart shows the peaks and
    valleys that matter biologically.

    Returns the path on success, ``None`` when the MSA was empty or
    degenerate (no columns, no rows) — the caller treats ``None`` as
    "nothing to draw, skip this family".

    ``hmm_hits_on_reference`` are the projected hits returned by
    :func:`project_hits_to_alignment`; pass ``None`` or ``[]`` to omit
    the domain ribbon entirely (the figure then has just the two
    metric tracks). ``family_color`` is the hex string from
    :func:`family_color_hex` and is used for the domain rectangles.
    """
    _require_matplotlib()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    msa = read_msa(msa_path)
    if not msa:
        return None
    rows = list(msa.values())
    n_cols = _columns(rows)
    n_rows = len(rows)
    if n_cols == 0 or n_rows == 0:
        return None

    entropy_raw = compute_shannon_entropy(rows)
    consensus_raw = compute_consensus_fraction(rows)
    entropy = sliding_window_mean(entropy_raw, window)
    consensus = sliding_window_mean(consensus_raw, window)

    # Layout: two stacked line-chart panels + a thin domain ribbon. Width
    # scales gently with MSA length so a 2000-column protein doesn't
    # produce a postage-stamp figure; floor + ceiling so the artefact
    # stays printable.
    width = max(7.0, min(20.0, n_cols / 90.0))
    has_domains = bool(hmm_hits_on_reference)
    height = 4.6 if has_domains else 3.6

    fig = plt.figure(figsize=(width, height))
    if has_domains:
        gs = fig.add_gridspec(
            nrows=3, ncols=1, height_ratios=[5, 5, 1], hspace=0.30,
        )
    else:
        gs = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[1, 1], hspace=0.40)

    x_vals = list(range(1, n_cols + 1))

    # Shannon entropy track — line chart in bits.
    ax_h = fig.add_subplot(gs[0])
    ax_h.plot(x_vals, entropy, color="#2c3e50", linewidth=1.1)
    ax_h.fill_between(x_vals, entropy, color="#2c3e50", alpha=0.12)
    ax_h.set_ylabel(
        f"Shannon entropy (bits)\n[{window}-aa window]",
        fontsize=8,
    )
    ax_h.set_xlim(0.5, n_cols + 0.5)
    ax_h.set_ylim(bottom=0.0)
    ax_h.grid(True, axis="y", linewidth=0.4, alpha=0.4)
    ax_h.tick_params(axis="both", labelsize=7)

    # Fraction-matching-consensus track — line chart in [0, 1].
    ax_c = fig.add_subplot(gs[1], sharex=ax_h)
    ax_c.plot(x_vals, consensus, color="#1f5e8c", linewidth=1.1)
    ax_c.fill_between(x_vals, consensus, color="#1f5e8c", alpha=0.12)
    ax_c.set_ylabel(
        f"Fraction matching consensus\n[{window}-aa window]",
        fontsize=8,
    )
    ax_c.set_xlim(0.5, n_cols + 0.5)
    ax_c.set_ylim(0.0, 1.05)
    ax_c.grid(True, axis="y", linewidth=0.4, alpha=0.4)
    ax_c.tick_params(axis="both", labelsize=7)

    # Domain-architecture ribbon. Every hit gets a labelled rectangle —
    # very short hits have their label overflow the box rather than
    # vanish, because the user explicitly asked for every domain to be
    # labelled (long labels on narrow domains read fine in practice
    # since the ribbon row has no other content competing for space).
    if has_domains:
        ax_d = fig.add_subplot(gs[2], sharex=ax_h)
        ax_d.set_ylim(0, 1)
        ax_d.set_yticks([])
        ax_d.set_ylabel("Domains", fontsize=8, rotation=0,
                        ha="right", va="center", labelpad=18)
        ax_d.set_xlim(0.5, n_cols + 0.5)
        ax_d.set_xlabel("MSA column", fontsize=8)
        ax_d.tick_params(axis="x", labelsize=7)
        ax_d.set_frame_on(False)
        for hit in hmm_hits_on_reference:
            af = hit.get("aln_from")
            at = hit.get("aln_to")
            if af is None or at is None:
                continue
            w = max(1, at - af + 1)
            rect = Rectangle(
                (af - 0.5, 0.15), w, 0.7,
                facecolor=family_color, edgecolor="#202020",
                linewidth=0.6, alpha=0.85,
            )
            ax_d.add_patch(rect)
            name = hit.get("hmm_name") or hit.get("name") or ""
            if name:
                ax_d.text(
                    af + w / 2.0 - 0.5, 0.5, name,
                    ha="center", va="center", fontsize=7, color="#101010",
                    clip_on=False,
                )
    else:
        ax_c.set_xlabel("MSA column", fontsize=8)

    title = f"Conservation — {family_label}  ({n_rows} sequences, {n_cols} columns)"
    if title_suffix:
        title = f"{title} — {title_suffix}"
    fig.suptitle(title, fontsize=10, y=0.995)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return out_png
