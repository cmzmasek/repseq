"""One conservation number per MSA, collected for a whole run.

Every tree path in repseq writes its alignment to disk as a
``*_msa*.fasta`` under the run's output directory — the whole-genome
tree (2E), the partitioned per-family alignments + supermatrix, the
per-protein marker trees (2F), the extra-protein trees, the polyprotein
peptide trees, and the per-segment NT trees. After the phylo steps
finish, :func:`write_msa_conservation_report` sweeps the output tree,
scores each alignment, and collects the results into a single
``{prefix}_msa_conservation.tsv`` (one row per MSA). Being a post-hoc
sweep, it is fully decoupled — a scoring bug can never affect a tree,
and it captures every MSA regardless of which code path produced it.

The score is the **mean per-column Jensen-Shannon divergence to a
residue background** (Capra & Singh, *Bioinformatics* 2007), with two
standard corrections:

* **Henikoff & Henikoff (1994) position-based sequence weighting**, so
  a clade of near-identical sequences plus one outlier doesn't read as
  artificially conserved.
* a **gap penalty**: each column's JSD is multiplied by
  ``(1 - weighted_gap_fraction)``, so a mostly-gap column contributes
  little.

JSD is computed with log base 2 and the symmetric mixture
(``m = (p + q) / 2``), so a column score is bounded ``[0, 1]``. The
run-level number is the mean over all columns. Interpretation:

* ``~0``     — columns look like the background (unrelated sequences);
* ``~0.6-0.8`` — typical for a real protein family of close homologs;
* ``~0.85-0.95`` — a perfectly conserved **protein** column. JSD-to-
  background cannot reach exactly 1 even for an invariant column; the
  ceiling depends on the conserved residue's background frequency (a
  rarer residue lands closer to 1). This is inherent to the metric and
  is the reason published JSD conservation values never saturate at 1.

For **nucleotide** alignments the ceiling is lower still — a perfectly
conserved column scores ~0.55 (exactly 0.549 against the uniform
4-letter background), because there are only four symbols and the
background is flat. So conservation numbers are only comparable
*within* an alphabet; don't compare a protein MSA's score against a
nucleotide MSA's. The ``alphabet`` column in the output flags which is
which.

The module is pure-Python (no numpy dependency) and reads aligned FASTA
directly, so it works on any ``_msa.fasta`` repseq writes.
"""

from __future__ import annotations

import logging
import math
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

# Robinson & Robinson amino-acid background frequencies — the
# ``blosum62.distribution`` Capra & Singh (2007) used as the JSD
# background. Re-normalised in code (the published numbers sum to
# ~1.002), so small rounding here is harmless.
_AA_BACKGROUND: dict[str, float] = {
    "A": 0.078, "R": 0.051, "N": 0.041, "D": 0.052, "C": 0.024,
    "Q": 0.034, "E": 0.059, "G": 0.083, "H": 0.025, "I": 0.062,
    "L": 0.092, "K": 0.056, "M": 0.024, "F": 0.044, "P": 0.043,
    "S": 0.059, "T": 0.055, "W": 0.014, "Y": 0.034, "V": 0.072,
}

# Uniform nucleotide background (U folded onto T upstream).
_NT_BACKGROUND: dict[str, float] = {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}

_GAP_CHARS = set("-.~")


def _normalise(dist: dict[str, float]) -> dict[str, float]:
    total = sum(dist.values())
    if total <= 0:
        return dict(dist)
    return {k: v / total for k, v in dist.items()}


_AA_BG = _normalise(_AA_BACKGROUND)
_NT_BG = _normalise(_NT_BACKGROUND)


def detect_alphabet(rows: list[str]) -> str:
    """Return ``"nucleotide"`` or ``"protein"`` from the residue content.

    Counts non-gap characters; if ≥90% are in ``{A,C,G,T,U,N}`` the
    alignment is nucleotide, else protein. Robust to the occasional
    ``N`` / ambiguity code.
    """
    nt_like = 0
    total = 0
    nt_set = set("ACGTUN")
    for row in rows:
        for ch in row.upper():
            if ch in _GAP_CHARS:
                continue
            total += 1
            if ch in nt_set:
                nt_like += 1
    if total == 0:
        return "protein"
    return "nucleotide" if nt_like / total >= 0.90 else "protein"


def henikoff_weights(rows: list[str]) -> list[float]:
    """Global Henikoff & Henikoff (1994) position-based sequence weights.

    For each column, a sequence carrying residue *r* contributes
    ``1 / (k * n_r)`` where *k* is the number of distinct symbols in the
    column (gaps counted as a symbol) and *n_r* the count of *r*. Summed
    across columns and **normalised to sum to 1**, so the weights are
    used directly as a probability mass over sequences. Falls back to
    uniform weights for a degenerate alignment (no columns, or every
    column monomorphic).
    """
    n = len(rows)
    if n == 0:
        return []
    if n == 1:
        return [1.0]
    ncols = max((len(r) for r in rows), default=0)
    weights = [0.0] * n
    for c in range(ncols):
        col = [row[c].upper() if c < len(row) else "-" for row in rows]
        counts: dict[str, int] = {}
        for ch in col:
            counts[ch] = counts.get(ch, 0) + 1
        k = len(counts)
        if k <= 0:
            continue
        for i, ch in enumerate(col):
            weights[i] += 1.0 / (k * counts[ch])
    total = sum(weights)
    if total <= 0:
        return [1.0 / n] * n
    return [w / total for w in weights]


def _column_distribution(
    col: list[str], weights: list[float], background: dict[str, float]
) -> tuple[dict[str, float], float]:
    """Weighted residue distribution for one column + its gap weight.

    Returns ``(p, gap_weight)`` where *p* is the weight-normalised
    distribution over the background alphabet (residues outside it —
    gaps, ``X``, ``N``, ambiguity codes — are excluded and folded into
    *gap_weight*), and *gap_weight* is the summed weight of the excluded
    rows (already on the [0,1] scale since the weights sum to 1). When
    every row is excluded, *p* is empty and *gap_weight* is 1.0.
    """
    alphabet = set(background)
    p: dict[str, float] = {}
    in_weight = 0.0
    gap_weight = 0.0
    for ch, w in zip(col, weights):
        c = ch.upper()
        if c == "U" and "T" in alphabet:  # RNA → DNA fold
            c = "T"
        if c in alphabet:
            p[c] = p.get(c, 0.0) + w
            in_weight += w
        else:
            gap_weight += w
    if in_weight > 0:
        p = {k: v / in_weight for k, v in p.items()}
    return p, gap_weight


def _jsd(p: dict[str, float], q: dict[str, float]) -> float:
    """Jensen-Shannon divergence (λ=0.5, log base 2) of *p* against *q*.

    ``m = (p + q) / 2``; ``JSD = ½·KL(p‖m) + ½·KL(q‖m)``. Bounded
    ``[0, 1]`` in base 2. ``0·log0`` terms are taken as 0; *m* is
    strictly positive wherever *q* is, so the logs are always defined.
    """
    keys = set(p) | set(q)
    kl_pm = 0.0
    kl_qm = 0.0
    for k in keys:
        pk = p.get(k, 0.0)
        qk = q.get(k, 0.0)
        mk = 0.5 * (pk + qk)
        if mk <= 0:
            continue
        if pk > 0:
            kl_pm += pk * math.log2(pk / mk)
        if qk > 0:
            kl_qm += qk * math.log2(qk / mk)
    jsd = 0.5 * kl_pm + 0.5 * kl_qm
    # Clamp tiny negative / >1 drift from floating point.
    if jsd < 0.0:
        return 0.0
    if jsd > 1.0:
        return 1.0
    return jsd


def score_rows(rows: list[str]) -> Optional[dict]:
    """Score an aligned set of sequences.

    Returns ``{"alphabet", "n_seqs", "n_sites", "mean_conservation",
    "mean_conservation_core"}`` or ``None`` when there is nothing to
    score (no rows or zero-width alignment).

    * ``mean_conservation`` — mean over all columns of the
      gap-penalised JSD (the headline number).
    * ``mean_conservation_core`` — mean *raw* JSD (no gap penalty) over
      columns with occupancy ≥ 0.5, or ``None`` when no column clears
      that occupancy. A view of conservation in the well-aligned core,
      insensitive to ragged ends.
    """
    rows = [r for r in rows if r != ""]
    if not rows:
        return None
    ncols = max((len(r) for r in rows), default=0)
    if ncols == 0:
        return None
    alphabet = detect_alphabet(rows)
    background = _NT_BG if alphabet == "nucleotide" else _AA_BG
    weights = henikoff_weights(rows)

    penalised: list[float] = []
    core_raw: list[float] = []
    for c in range(ncols):
        col = [row[c] if c < len(row) else "-" for row in rows]
        p, gap_weight = _column_distribution(col, weights, background)
        occupancy = 1.0 - gap_weight
        if occupancy <= 0 or not p:
            penalised.append(0.0)
            continue
        jsd = _jsd(p, background)
        penalised.append(jsd * occupancy)
        if occupancy >= 0.5:
            core_raw.append(jsd)

    mean_cons = sum(penalised) / len(penalised) if penalised else 0.0
    mean_core = sum(core_raw) / len(core_raw) if core_raw else None
    return {
        "alphabet": alphabet,
        "n_seqs": len(rows),
        "n_sites": ncols,
        "mean_conservation": mean_cons,
        "mean_conservation_core": mean_core,
    }


def _read_msa_rows(path: Path) -> list[str]:
    """Read an aligned FASTA into a list of (upper-case-preserving) rows."""
    rows: list[str] = []
    buf: list[str] = []
    started = False
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if started:
                    rows.append("".join(buf))
                buf = []
                started = True
            elif started:
                buf.append(line.strip())
    if started:
        rows.append("".join(buf))
    return rows


def _classify(rel: Path, prefix: str) -> tuple[str, str, bool]:
    """Infer ``(role, label, trimmed)`` from an MSA path relative to out_dir.

    * Top-level ``{prefix}_msa.fasta`` → ``("genome", "", True)``;
      ``{prefix}_msa_{family}.fasta`` (the partitioned per-family
      alignments) → ``("partition_family", family, True)``.
    * Files inside a ``{prefix}_<kind>/`` subdir are named
      ``<label>_msa.fasta``; the subdir maps to a role
      (``per_protein`` → ``marker``, ``extra_protein`` →
      ``extra_protein``, ``polyprotein`` → ``peptide``, ``per_segment``
      → ``segment_nt``; anything else uses the subdir name).
    * An ``_untrimmed`` infix flags the retained pre-trimAl companion
      (``trimmed=False``).
    """
    name = rel.name
    stem = name[:-len(".fasta")] if name.endswith(".fasta") else name
    trimmed = True
    if stem.endswith("_untrimmed"):
        trimmed = False
        stem = stem[: -len("_untrimmed")]

    parts = rel.parts
    subdir = parts[-2] if len(parts) >= 2 else None

    def _strip_prefix(s: str) -> str:
        pfx = f"{prefix}_"
        return s[len(pfx):] if s.startswith(pfx) else s

    if subdir is not None:
        # Files written under a per-step subdirectory: "<label>_msa".
        label = stem[: -len("_msa")] if stem.endswith("_msa") else stem
        label = _strip_prefix(label)
        sub = subdir
        for suffix, role in (
            ("_per_protein", "marker"),
            ("_extra_protein", "extra_protein"),
            ("_polyprotein", "peptide"),
            ("_per_segment", "segment_nt"),
        ):
            if sub.endswith(suffix):
                return role, label, trimmed
        return _strip_prefix(sub), label, trimmed

    # Top-level: "{prefix}_msa" or "{prefix}_msa_{family}".
    base = f"{prefix}_msa"
    if stem == base:
        return "genome", "", trimmed
    if stem.startswith(base + "_"):
        return "partition_family", stem[len(base) + 1:], trimmed
    return "genome", _strip_prefix(stem), trimmed


def _fmt(v: Optional[float]) -> str:
    return "" if v is None else f"{v:.4f}"


def write_msa_conservation_report(
    out_dir: Path, prefix: str, cfg: Optional[dict] = None
) -> Optional[Path]:
    """Score every ``*_msa*.fasta`` under *out_dir* into one TSV.

    Writes ``{prefix}_msa_conservation.tsv`` (one row per MSA, sorted
    genome-first then by role/label) and returns its path, or ``None``
    when conservation scoring is disabled, no MSA files exist, or none
    could be scored. Soft by construction — callers wrap in try/except
    so a scoring failure never voids the trees that were already built.
    """
    if cfg is not None:
        cons_cfg = (cfg.get("phylo", {}) or {}).get("conservation", {}) or {}
        if not cons_cfg.get("enabled", True):
            return None

    msa_paths = sorted(out_dir.rglob("*_msa*.fasta"))
    if not msa_paths:
        return None

    # Stable role ordering: whole-genome tree first, then its partition
    # families, then the per-protein / extra / peptide / segment trees.
    role_order = {
        "genome": 0, "partition_family": 1, "marker": 2,
        "extra_protein": 3, "peptide": 4, "segment_nt": 5,
    }
    rows_out: list[tuple] = []
    for path in msa_paths:
        rel = path.relative_to(out_dir)
        role, label, trimmed = _classify(rel, prefix)
        try:
            seqs = _read_msa_rows(path)
            metrics = score_rows(seqs)
        except Exception as exc:  # pragma: no cover - defensive
            logger.warning("conservation: could not score %s: %s", rel, exc)
            continue
        if metrics is None:
            continue
        sort_key = (role_order.get(role, 99), not trimmed, label, str(rel))
        rows_out.append((
            sort_key,
            str(rel),
            role,
            label,
            metrics["alphabet"],
            "TRUE" if trimmed else "FALSE",
            metrics["n_seqs"],
            metrics["n_sites"],
            _fmt(metrics["mean_conservation"]),
            _fmt(metrics["mean_conservation_core"]),
        ))

    if not rows_out:
        return None
    rows_out.sort(key=lambda r: r[0])

    path = out_dir / f"{prefix}_msa_conservation.tsv"
    header = (
        "msa\trole\tlabel\talphabet\ttrimmed\tn_seqs\tn_sites\t"
        "mean_conservation\tmean_conservation_core"
    )
    with open(path, "w") as fh:
        fh.write(header + "\n")
        for r in rows_out:
            fh.write("\t".join(str(x) for x in r[1:]) + "\n")
    return path
