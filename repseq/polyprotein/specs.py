"""Parse polyprotein-cutting declarations into typed specs.

Mirrors :func:`repseq.phylo.per_protein.collect_extra_specs`: reads from
``clustering.polyprotein`` (non-segmented) or
``virus.polyprotein`` (per-segment dict, segmented) and returns one
:class:`PolyproteinSpec` per declaration, paired with the segment scope
(``None`` for non-segmented). The config has already passed validation
by the time these helpers run, so the parser is allowed to be terse —
unrecognised entries are dropped silently rather than raising.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any, Optional

from ..config import get_virus_config


_SANITIZE_RE = re.compile(r"[^A-Za-z0-9._-]+")


def _sanitize(name: str) -> str:
    """Filesystem-safe component (mirrors :mod:`phylo.per_protein`)."""
    return _SANITIZE_RE.sub("_", (name or "").strip()).strip("_") or "polyprotein"


@dataclass
class PeptideSpec:
    """One mature peptide within a polyprotein.

    ``hmms`` is a list of **alternative HMM tokens** (v0.34.0+) —
    each token is either a single HMM profile name
    (``"CoV_NSP8"``) or a multidomain architecture in N-to-C order
    (``"CoV_RPol_N--RdRP_1"``). The peptide is located when *any one* of
    the listed tokens is satisfied on the parent CDS (OR across the
    list); for a multidomain token, the peptide is satisfied only when
    *all* of the named HMMs hit in the declared order (AND within the
    token). This matches the existing ``cluster_protein.hmms[]``
    semantics: same ``parse_hmm_token`` validator, same N-to-C
    enforcement, AND-within-token, OR-across-list.

    The biological reason for OR support: a mature peptide's
    architecture can differ across virus genera. Coronavirus NSP1 has
    distinct alpha- vs. beta-CoV folds (the alpha form has no homologue
    in beta and vice-versa), so writing one `hmms:` list with both
    architectures lets each rep's CDS match whichever variant it
    carries.

    Config users may supply either the singular ``hmm: <token>``
    (legacy v0.33.0/early-v0.34.0 form, kept for backwards
    compatibility — wrapped into a 1-element list internally) or the
    plural ``hmms: [<token>, ...]`` (the new OR form). Exactly one of
    the two must be set.

    When multiple alternatives are satisfied on one CDS, the slicer
    picks the alternative with the best (lowest) worst-domain E-value,
    mirroring how :func:`repseq.phylo.per_protein._best_satisfying_cds_any`
    chooses across alternative architectures for the
    ``cluster_protein`` HMM gate. The selected token is recorded on
    each :class:`SlicedPeptide` for audit.

    ``cleavage_motif`` (optional) is the residue motif found just
    N-terminal of the cut that liberates this peptide (3CL: ``"LQ"``,
    picornavirus 3C: ``"Q"``, etc.). The motif refines the cut between
    the *previous* peptide and this one when the spec's cut strategy is
    ``"motif"``.
    """

    name: str
    hmms: list[str]  # OR-list of tokens, each "Name" or "A--B--C"
    cleavage_motif: Optional[str] = None


@dataclass
class PolyproteinSpec:
    """One polyprotein declaration.

    ``peptides`` is ordered N-to-C; ``cut_strategy`` chooses how
    inter-peptide cuts are placed:

    * ``boundary`` — each peptide spans its HMM hit's ``ali_from..ali_to``
      verbatim (lossy at the seams; deterministic).
    * ``bisect`` — each peptide spans midpoint(prev_hit, this_hit) ..
      midpoint(this_hit, next_hit); endpoints extend to the N-/C-term.
      No residues dropped; cuts are artificial.
    * ``motif`` — like bisect, but for each inter-peptide cut, if the
      *downstream* peptide carries ``cleavage_motif``, the cut snaps to
      the last occurrence of the motif within ``±motif_window_aa`` of
      the bisect point. Falls back to bisect for that cut if no motif is
      found in the window.

    ``min_peptides_hit`` is the parent-CDS identification threshold: a
    CDS qualifies as the polyprotein when its proteins carry hits from
    at least this many of the declared peptide HMMs.
    """

    name: str
    peptides: list[PeptideSpec]
    cut_strategy: str = "bisect"  # boundary | bisect | motif
    motif_window_aa: int = 10
    min_peptides_hit: int = 2
    segment: Optional[str] = None  # segmented mode: the scope segment

    @property
    def file_basename(self) -> str:
        """Filesystem-safe basename component (sanitised, no segment prefix).

        Spec names are required to be unique across all segments (validated
        upstream), so no segment prefix is needed — the name alone keys the
        output filenames.
        """
        return _sanitize(self.name)


def _parse_peptide(entry: dict) -> Optional[PeptideSpec]:
    """Build one :class:`PeptideSpec` from a validated peptide config dict.

    Accepts either ``hmm: <token>`` (legacy singular form, wrapped to a
    1-element list internally) or ``hmms: [<token>, ...]`` (the new OR
    form). Validation upstream guarantees exactly one of the two is set.
    """
    name = (entry.get("name") or "").strip()
    if not name:
        return None

    tokens: list[str] = []
    if isinstance(entry.get("hmms"), list):
        tokens = [
            t.strip() for t in entry["hmms"]
            if isinstance(t, str) and t.strip()
        ]
    elif isinstance(entry.get("hmm"), str):
        single = entry["hmm"].strip()
        if single:
            tokens = [single]
    if not tokens:
        return None

    motif = entry.get("cleavage_motif")
    if isinstance(motif, str):
        motif = motif.strip() or None
    else:
        motif = None
    return PeptideSpec(name=name, hmms=tokens, cleavage_motif=motif)


def _parse_polyprotein_entry(
    entry: dict, segment: Optional[str]
) -> Optional[PolyproteinSpec]:
    """Build one :class:`PolyproteinSpec` from a validated config dict."""
    name = (entry.get("name") or "").strip()
    if not name:
        return None
    peptides_raw = entry.get("peptides") or []
    peptides = [_parse_peptide(p) for p in peptides_raw if isinstance(p, dict)]
    peptides = [p for p in peptides if p is not None]
    if len(peptides) < 2:
        return None

    has_motif = any(p.cleavage_motif for p in peptides)
    # Per the design: if the user didn't pin a strategy, default to motif
    # when any peptide declared one, else bisect. Boundary is opt-in only.
    cut_strategy = entry.get("cut_strategy") or ("motif" if has_motif else "bisect")
    if cut_strategy not in ("boundary", "bisect", "motif"):
        cut_strategy = "bisect"

    return PolyproteinSpec(
        name=name,
        peptides=peptides,
        cut_strategy=cut_strategy,
        motif_window_aa=int(entry.get("motif_window_aa", 10) or 10),
        min_peptides_hit=max(1, int(entry.get("min_peptides_hit", 2) or 2)),
        segment=segment,
    )


def collect_polyprotein_specs(cfg: dict[str, Any]) -> list[PolyproteinSpec]:
    """Return one :class:`PolyproteinSpec` per declaration.

    Reads ``clustering.polyprotein`` for non-segmented runs and
    ``virus.polyprotein`` (per-segment dict) for segmented runs. Order
    of declaration is preserved; in segmented mode the order is
    (segment-declaration-order, within-segment-list-order).
    """
    specs: list[PolyproteinSpec] = []

    segmented = bool((cfg.get("segmented", {}) or {}).get("enabled"))
    if segmented:
        virus = get_virus_config(cfg) or {}
        per_seg = virus.get("polyprotein") or {}
        # Honour the configured segment order for deterministic output
        # filenames; fall back to insertion order for unknown segments.
        segments = list(virus.get("segments") or [])
        for extra in per_seg.keys():
            if extra not in segments:
                segments.append(extra)
        for seg in segments:
            for entry in per_seg.get(seg, []) or []:
                if not isinstance(entry, dict):
                    continue
                parsed = _parse_polyprotein_entry(entry, seg)
                if parsed is not None:
                    specs.append(parsed)
    else:
        for entry in (cfg.get("clustering", {}) or {}).get("polyprotein", []) or []:
            if not isinstance(entry, dict):
                continue
            parsed = _parse_polyprotein_entry(entry, None)
            if parsed is not None:
                specs.append(parsed)

    return specs
