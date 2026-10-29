"""2F — one phylogenetic tree per declared HMM marker family.

Where :func:`repseq.phylo.pipeline.run_phylogeny` (2E) builds a *single*
tree from each representative's whole sequence (or marker-protein concat),
this module builds *one tree per protein family*, where a "family" is a
marker **spec** the user already declared for QC
(``virus.segment_markers[seg]`` / ``virus.cluster_protein[seg]`` for
segmented runs, an entry of ``clustering.cluster_protein`` otherwise).

A spec's ``hmms:`` list holds **alternative domain architectures (OR)** —
e.g. coronavirus Spike as ``["CoV_S1--CoV_S2",
"bCoV_S1_N--bCoV_S1_RBD--CoV_S2"]`` so alpha- and beta-CoV Spikes land in
*one* tree. For each representative we take the CDS that satisfies *any*
of the spec's tokens (via :func:`repseq.hmm.runner.cds_satisfies_token` —
intervening/extra domains tolerated exactly as the QC gate tolerates
them), and build an MSA + tree + phyloXML on those protein translations.
The heavy lifting (MAFFT → IQ-TREE/FastTree → root → LCA → phyloXML) is
the same :func:`repseq.phylo.pipeline._build_tree` engine 2E uses, so
rooting, LCA labelling, ``phylo.tool``, and per-leaf annotation all carry
over.

Requirements (soft-fail, mirroring ``--plot`` / ``--phylo``):
  * The HMM tier must have run this session (``hmm.enabled`` + configured
    ``hmms:`` + populated ``hmm_hits`` on the proteins). Without it no CDS
    can satisfy any token; we raise :class:`PhyloError` with that reason.
  * Each family needs at least ``phylo.per_protein.min_taxa`` (default 3,
    the tree-builder floor) representatives carrying the architecture;
    sparser families are skipped with a log note.

Outputs land in a ``{prefix}_per_protein/`` subdirectory, one set per
built family: ``<family>_msa.fasta``, ``<family>_tree.nwk``,
``<family>_tree.xml``, ``<family>_tree_id_map.tsv`` (+
``<family>_iqtree_summary.txt`` for IQ-TREE). The family name is the
spec's ``name:`` when given (``Spike``), else the single token, else the
first token + ``_altN`` for an unnamed multi-architecture spec; in
segmented mode it is prefixed with the segment (``S_Bunya_nucleocap``,
``M_Spike``) so two segments declaring the same marker never collide.
"""

from __future__ import annotations

import logging
import re
from functools import partial
from pathlib import Path
from typing import Any, Callable, Optional, Union

from ..clustering.marker import (
    LengthBounds,
    format_length_bounds,
    parse_length_bounds,
    marker_candidate_ok,
    require_met_start_active,
    starts_with_met,
    within_length_bounds,
)
from ..config import get_virus_config
from ..hmm.runner import cds_satisfies_token, parse_hmm_token
from ..models import Sequence
from ..overrides import ProtectionPolicy
from ..segmented.completeness import resolve_segment_marker
from .coloring import build_color_scheme
from .pipeline import PhyloError, _build_tree

logger = logging.getLogger(__name__)


def _sanitize(name: str) -> str:
    """Filesystem-safe family component.

    Keeps alphanumerics and ``-`` / ``_`` / ``.`` (so the multidomain
    ``--`` separator survives intact), collapses any other run — notably
    whitespace and path separators — to a single ``_``.
    """
    return re.sub(r"[^A-Za-z0-9._-]+", "_", (name or "").strip()).strip("_") or "family"


def _dedup(items: list[str]) -> list[str]:
    """Order-preserving de-duplication of a token list."""
    seen: set[str] = set()
    out: list[str] = []
    for it in items:
        if it not in seen:
            seen.add(it)
            out.append(it)
    return out


def _resolve_segment_marker(
    seg_name: str,
    segment_markers: dict,
    cluster_protein_per_seg: dict,
) -> tuple[Optional[str], list[str]]:
    """``(spec_name, tokens)`` gating a segment — mirrors ``cli._resolve_segment_hmms``.

    Replicated here (rather than imported) because ``cli`` imports the
    phylo package, so importing back would be circular. Lookup order:
    ``segment_markers[seg]`` first (v0.13+ HMM-aware form). If it carries
    ``hmms:``, those tokens win. If it exists but is alias-only (empty
    ``hmms:``), fall through to any dict-form ``cluster_protein[seg]``
    entry's ``hmms`` so a user who split aliases into ``segment_markers``
    and HMMs into ``cluster_protein`` still gets per-protein trees built
    for the segment. The spec name — when present on the marker dict —
    is used for the family label; ``segment_markers`` entries usually
    have none (the segment is the identity), so the segment prefix
    carries the meaning.
    """
    name_from_sm: Optional[str] = None
    if seg_name in segment_markers:
        spec = segment_markers[seg_name] or {}
        tokens = list(spec.get("hmms") or [])
        if tokens:
            return spec.get("name"), tokens
        name_from_sm = spec.get("name")
        # Fall through: segment_markers is alias-only here — try cluster_protein.
    if seg_name in cluster_protein_per_seg:
        tokens = []
        name: Optional[str] = name_from_sm
        for entry in cluster_protein_per_seg[seg_name] or []:
            if isinstance(entry, dict):
                name = name or entry.get("name")
                tokens.extend(entry.get("hmms") or [])
        return name, tokens
    return name_from_sm, []


def _family_label(
    name: Optional[str], tokens: list[str], segment: Optional[str]
) -> str:
    """Filesystem-safe family label, name-first.

    Priority: the spec ``name:`` when given (``Spike``), else the single
    token verbatim (``CoV_S1--CoV_S2``), else — for an unnamed spec with
    several alternative architectures — the first token plus an ``_altN``
    suffix (``CoV_S1--CoV_S2_alt2``), so we never glue whole architectures
    into one unwieldy filename. Segmented families are prefixed with the
    segment so two segments declaring the same marker never collide.

    Alias-only specs (no tokens) require a ``name:`` — if neither name
    nor tokens are present we fall back to ``"marker"`` (sanitised), but
    the segment prefix in segmented mode still disambiguates.
    """
    if name and name.strip():
        base = _sanitize(name)
    elif len(tokens) == 1:
        base = _sanitize(tokens[0])
    elif tokens:
        base = f"{_sanitize(tokens[0])}_alt{len(tokens)}"
    else:
        base = "marker"
    if segment is not None:
        return f"{_sanitize(segment)}_{base}"
    return base


def _resolve_segment_marker_full(
    seg_name: str,
    segment_markers: dict,
    cluster_protein_per_seg: dict,
) -> tuple[Optional[str], list[str], list[str]]:
    """``(name, hmm_tokens, aliases)`` for one segment's marker.

    Delegates to :func:`repseq.segmented.completeness.resolve_segment_marker`,
    the one resolver clustering, the QC gates and these output pickers share
    (``segment_markers`` wins; an alias-only entry pulls ``hmms:`` from
    ``cluster_protein[seg]``).
    """
    return resolve_segment_marker(seg_name, segment_markers, cluster_protein_per_seg)


def _iter_marker_entries(
    cfg: dict[str, Any],
) -> list[tuple[str, list[str], list[str], Optional[str], Optional[LengthBounds]]]:
    """Every usable cluster-marker spec as ``(label, tokens, aliases,
    segment_or_None, length_bounds)``, in declaration order.

    Single source for :func:`collect_marker_specs` and
    :func:`collect_marker_length_bounds`, so the two can never disagree on
    which specs exist or how they are labelled. ``length_bounds`` comes from
    ``segment_markers[seg].length_aa`` (segmented; the legacy per-segment
    ``cluster_protein`` form takes no bound) or the
    ``clustering.cluster_protein`` entry's ``length_aa`` (non-segmented).
    """
    out: list[tuple[str, list[str], list[str], Optional[str], Optional[LengthBounds]]] = []

    segmented = bool((cfg.get("segmented", {}) or {}).get("enabled"))
    if segmented:
        virus = get_virus_config(cfg) or {}
        segment_markers = virus.get("segment_markers") or {}
        cluster_protein = virus.get("cluster_protein") or {}
        segments = list(virus.get("segments") or [])
        for extra in (*segment_markers.keys(), *cluster_protein.keys()):
            if extra not in segments:
                segments.append(extra)
        for seg in segments:
            name, tokens, aliases = _resolve_segment_marker_full(
                seg, segment_markers, cluster_protein
            )
            tokens = _dedup(tokens)
            aliases = _dedup(aliases)
            if not tokens and not aliases:
                continue
            raw = (segment_markers.get(seg) or {}).get("length_aa")
            out.append((
                _family_label(name, tokens, seg), tokens, aliases, seg,
                parse_length_bounds(raw),
            ))
    else:
        for entry in (cfg.get("clustering", {}) or {}).get("cluster_protein", []) or []:
            if not isinstance(entry, dict):
                continue
            tokens = _dedup(list(entry.get("hmms") or []))
            aliases = _dedup(list(entry.get("aliases") or []))
            if not tokens and not aliases:
                continue
            out.append((
                _family_label(entry.get("name"), tokens, None), tokens, aliases,
                None, parse_length_bounds(entry.get("length_aa")),
            ))
    return out


def collect_marker_specs(
    cfg: dict[str, Any],
) -> list[tuple[str, list[str], list[str], Optional[str]]]:
    """Return ``[(family_label, hmm_tokens, aliases, segment_or_None), …]``.

    Like :func:`collect_family_specs` but also includes alias-only specs
    (those declaring ``aliases:`` without ``hmms:``) — used by the
    always-on per-protein FASTA writer (``{prefix}_per_protein_fasta/``).
    Specs with neither aliases nor tokens are skipped.

    For each spec, downstream callers can pick the satisfying CDS via the
    HMM gate (when ``hmm_tokens`` is non-empty and the HMM tier ran) and
    fall back to alias matching against ``/product`` otherwise — same
    chain :func:`repseq.clustering.marker.select_marker_protein` uses.
    The spec's ``length_aa`` is not in this tuple; pick through
    :class:`MarkerBoundsContext`, which applies it.
    """
    return [(l, t_, a, s) for l, t_, a, s, _b in _iter_marker_entries(cfg)]


def collect_marker_length_bounds(
    cfg: dict[str, Any],
) -> dict[tuple[str, Optional[str]], LengthBounds]:
    """``{(family_label, segment_or_None): (min, max)}`` for every cluster
    marker spec that declares ``length_aa``."""
    return {
        (l, s): b for l, _t, _a, s, b in _iter_marker_entries(cfg) if b is not None
    }


#: ``exempt`` arguments take a bool or a zero-argument callable. The callable
#: form defers the protection lookup (``ProtectionPolicy.protects_any`` over a
#: rep and all its segments) until a bounded pick actually misses, which is
#: rare; the coverage report calls the pickers millions of times.
Exempt = Union[bool, Callable[[], bool]]


def _is_exempt(exempt: Exempt) -> bool:
    return bool(exempt()) if callable(exempt) else bool(exempt)


class MarkerBoundsContext:
    """What an output-side picker needs to take EXACTLY the CDS clustering
    took for a cluster marker: each spec's ``length_aa``, the global
    5'-completeness requirement, and the ``marker_quality`` protection
    policy.

    Clustering (``clustering.marker.select_marker_protein``) picks the
    longest satisfying CDS inside the bound that also starts with Met when
    ``qc.marker_completeness.require_met_start`` is on; an isolate protected
    against the ``marker_quality`` QC stage is picked as if unrestricted.
    Every per-protein FASTA / tree / supermatrix / coverage pick goes
    through :meth:`pick` (or :func:`pick_marker_cds` with :meth:`bounds_for`
    and :meth:`exempt`) so the same rule holds there.
    """

    def __init__(
        self,
        bounds: dict[tuple[str, Optional[str]], LengthBounds],
        policy: Optional[ProtectionPolicy] = None,
        require_met_start: bool = False,
    ) -> None:
        self.bounds = bounds
        self.policy = policy
        self.require_met_start = require_met_start

    def bounds_for(self, label: str, segment: Optional[str]) -> Optional[LengthBounds]:
        return self.bounds.get((label, segment))

    def exempt(self, rep: Sequence) -> bool:
        if not self.bounds and not self.require_met_start:
            return False
        if self.policy is None:
            return False
        return self.policy.protects_any(
            [rep, *(rep.concat_segments or [])], "marker_quality"
        )

    def pick(
        self,
        rep: Sequence,
        proteins: list[dict],
        label: str,
        tokens: list[str],
        aliases: list[str],
        segment: Optional[str],
        hmm_active: bool,
        overlap_tolerance: int,
    ) -> Optional[dict]:
        return pick_marker_cds(
            proteins, tokens, aliases, hmm_active,
            overlap_tolerance=overlap_tolerance,
            bounds=self.bounds_for(label, segment),
            exempt=lambda: self.exempt(rep),
            require_met_start=self.require_met_start,
        )


def build_marker_bounds_context(cfg: dict[str, Any]) -> MarkerBoundsContext:
    """Cluster-marker restrictions + the protection policy, once per writer."""
    return MarkerBoundsContext(
        collect_marker_length_bounds(cfg),
        ProtectionPolicy.from_cfg(cfg or {}),
        require_met_start_active(cfg),
    )


def _iter_extra_entries(
    cfg: dict[str, Any],
) -> list[tuple[str, list[str], list[str], Optional[str], dict]]:
    """Every usable ``extra_protein`` entry as ``(label, tokens, aliases,
    segment_or_None, raw_entry)``, in declaration order.

    Single source for :func:`collect_extra_specs` and
    :func:`collect_extra_length_bounds`, so the two can never disagree on
    which entries exist or how they are labelled.
    """
    def _one(entry: Any, seg: Optional[str]):
        if not isinstance(entry, dict):
            return None
        name = (entry.get("name") or "").strip()
        if not name:
            return None
        tokens = _dedup(list(entry.get("hmms") or []))
        aliases = _dedup(list(entry.get("aliases") or []))
        if not tokens and not aliases:
            return None
        return (_sanitize(name), tokens, aliases, seg, entry)

    out = []
    segmented = bool((cfg.get("segmented", {}) or {}).get("enabled"))
    if segmented:
        virus = get_virus_config(cfg) or {}
        extra = virus.get("extra_protein") or {}
        for seg_name, entries in extra.items():
            for entry in entries or []:
                item = _one(entry, seg_name)
                if item is not None:
                    out.append(item)
    else:
        for entry in (cfg.get("clustering", {}) or {}).get("extra_protein", []) or []:
            item = _one(entry, None)
            if item is not None:
                out.append(item)
    return out


def collect_extra_specs(
    cfg: dict[str, Any],
) -> list[tuple[str, list[str], list[str], Optional[str]]]:
    """Return ``[(family_label, hmm_tokens, aliases, segment_or_None), …]``
    for ``extra_protein`` declarations.

    Like :func:`collect_marker_specs` but reads exclusively from the
    ``extra_protein`` schema:

    * non-segmented: ``clustering.extra_protein`` (list of dicts);
    * segmented: ``virus.extra_protein`` (per-segment dict of lists).

    Validation already required every entry to be a dict with a non-empty
    ``name:`` and at least one of ``aliases:`` / ``hmms:``, so the label
    is taken from the spec name verbatim (sanitised for the filesystem) —
    no segment prefix even in segmented mode, since the user's filenames
    use the name alone.

    The optional ``length_aa`` bounds are NOT part of this tuple (its
    4-element shape is shared with the marker specs); read them with
    :func:`collect_extra_length_bounds`, or pick CDSes through
    :func:`pick_extra_cds`, which applies them.
    """
    return [
        (label, tokens, aliases, seg)
        for label, tokens, aliases, seg, _entry in _iter_extra_entries(cfg)
    ]


# ---------------------------------------------------------------------------
# extra_protein CDS picking (v0.75.0): "HMM specs choose first" + length_aa
# ---------------------------------------------------------------------------
#
# An alias is a case-insensitive SUBSTRING of /product, and GenBank /product
# labels are sometimes wrong — e.g. a tospovirus M segment whose 1,129 aa
# glycoprotein precursor is labelled "NSm" while its real 303 aa NSm is
# labelled "Gn/Gc" (submitter label swap). The alias picker keeps the LONGEST
# match, so the mislabelled glycoprotein won and landed in the NSm file.
#
# Two guards for ``extra_protein`` specs. (Cluster-marker specs take a
# ``length_aa`` too since v0.76.0, but theirs is a QC requirement enforced in
# every picker — clustering included — so a marker's per-protein FASTA keeps
# showing exactly the CDS clustering used; see :class:`MarkerBoundsContext`.)
#
# 1. Claiming. Within one scope (a segment; the whole sequence when
#    non-segmented), a CDS that a declared HMM-bearing spec PICKS — a cluster
#    marker or an HMM extra_protein spec — is withheld from alias matching.
#    HMM evidence of identity outranks a free-text label. Only the CDS a spec
#    actually picks is claimed (not every CDS that merely satisfies its
#    tokens), so a shorter protein sharing a domain with the marker (the
#    influenza PA-X / PA shape) stays available. Claiming needs the HMM tier;
#    without it nothing is claimed and picking is exactly as before.
# 2. ``length_aa: [min, max]`` (inclusive; either end may be null). Candidate
#    CDSes outside the bounds are dropped BEFORE the pick, so an out-of-range
#    longer match cannot hide an in-range one. Applies to alias and HMM
#    extra_protein specs alike.

# The bound helpers are shared with the clustering picker (one definition,
# so clustering and every output-side picker read a bound identically). The
# private aliases keep this module's historical names working.
_length_bounds = parse_length_bounds
_within_bounds = within_length_bounds


def collect_extra_length_bounds(
    cfg: dict[str, Any],
) -> dict[tuple[str, Optional[str]], LengthBounds]:
    """``{(family_label, segment_or_None): (min, max)}`` for every
    ``extra_protein`` entry that declares ``length_aa``."""
    out: dict[tuple[str, Optional[str]], LengthBounds] = {}
    for label, _t, _a, seg, entry in _iter_extra_entries(cfg):
        bounds = _length_bounds(entry.get("length_aa"))
        if bounds is not None:
            out[(label, seg)] = bounds
    return out


def _matches_any_alias(prot: dict, aliases: list[str]) -> bool:
    product = (prot.get("product") or "").lower()
    return any(
        (a or "").lower().strip() and (a or "").lower().strip() in product
        for a in aliases
    )


class ExtraPickContext:
    """Per-run state for :func:`pick_extra_cds`: the HMM claimants per scope
    and the ``length_aa`` bounds. Build once per writer with
    :func:`build_extra_pick_context`; it memoises claimed sets, so it must
    not outlive the Sequence objects it was used on."""

    def __init__(
        self,
        *,
        hmm_active: bool,
        overlap_tolerance: int,
        bounds: dict[tuple[str, Optional[str]], LengthBounds],
        claimants: dict[
            Optional[str],
            list[tuple[str, list[str], Optional[LengthBounds], bool]],
        ],
        marker_ctx: Optional["MarkerBoundsContext"] = None,
        require_met_start: bool = False,
    ) -> None:
        self.hmm_active = hmm_active
        self.overlap_tolerance = overlap_tolerance
        self.bounds = bounds
        self.require_met_start = require_met_start
        # (label, tokens, bounds, is_cluster_marker) per scope.
        self.claimants = claimants
        self.marker_ctx = marker_ctx
        self._claim_cache: dict[tuple, dict[int, str]] = {}

    def exempt(self, rep: Sequence) -> bool:
        """True when ``rep`` is protected against the ``marker_quality`` QC
        stage, so its cluster markers are picked as if unrestricted."""
        return self.marker_ctx is not None and self.marker_ctx.exempt(rep)

    def claimed(
        self, proteins: list[dict], segment: Optional[str], exempt: Exempt = False,
    ) -> dict[int, str]:
        """``{id(cds): claimant_label}`` for CDSes an HMM spec picks in scope.

        Keyed on object identity: the claimant's pick comes from the same
        protein dicts the alias spec sees. Empty when the HMM tier is off.
        A cluster-marker claimant claims exactly the CDS the marker picker
        takes — inside its ``length_aa`` and starting with Met when the run
        requires it, or (``exempt``: the rep is protected against
        ``marker_quality``) its unrestricted pick when nothing qualifies.
        ``exempt`` is only evaluated when a restricted marker claimant misses.
        """
        if not self.hmm_active:
            return {}
        with_seq = [p for p in (proteins or []) if p.get("sequence")]
        ids = tuple(id(p) for p in with_seq)
        restricted = self._claims(segment, ids, with_seq, unbounded_markers=False)
        if not restricted[1] or not _is_exempt(exempt):
            return restricted[0]
        return self._claims(segment, ids, with_seq, unbounded_markers=True)[0]

    def _claims(
        self, segment, ids, with_seq, *, unbounded_markers: bool,
    ) -> tuple[dict[int, str], bool]:
        """``(claims, a_restricted_marker_missed)``, memoised per scope."""
        key = (segment, unbounded_markers, ids)
        cached = self._claim_cache.get(key)
        if cached is not None:
            return cached
        out: dict[int, str] = {}
        missed = False
        met = self.require_met_start and not unbounded_markers
        for label, tokens, bounds, is_marker in self.claimants.get(segment, []):
            pick = pick_marker_cds(
                with_seq, tokens, [], True,
                overlap_tolerance=self.overlap_tolerance,
                bounds=bounds, exempt=(unbounded_markers and is_marker),
                # The exemption is marker-only, exactly as the bounds waiver
                # above is: an HMM extra_protein spec claims the CDS its own
                # `pick_extra_cds` call will take, or the claim lands on a
                # CDS that spec withholds and leaves its real pick free for
                # another alias spec to take as well.
                require_met_start=(met if is_marker else self.require_met_start),
            )
            if pick is None and is_marker and (bounds or self.require_met_start):
                missed = True
            if pick is not None:
                out.setdefault(id(pick), label)
        self._claim_cache[key] = (out, missed)
        return out, missed


def build_extra_pick_context(cfg: dict[str, Any], hmm_active: bool) -> ExtraPickContext:
    """Collect every HMM-bearing spec (cluster markers + HMM extra_protein
    specs) per scope, plus the extra_protein ``length_aa`` bounds. Cluster
    markers claim with their own ``length_aa`` so the claimed CDS is exactly
    the one the marker picked."""
    bounds = collect_extra_length_bounds(cfg)
    marker_ctx = build_marker_bounds_context(cfg)
    claimants: dict[
        Optional[str], list[tuple[str, list[str], Optional[LengthBounds], bool]]
    ] = {}
    for label, tokens, _aliases, seg in collect_marker_specs(cfg):
        if tokens:
            claimants.setdefault(seg, []).append(
                (label, tokens, marker_ctx.bounds_for(label, seg), True)
            )
    for label, tokens, _aliases, seg in collect_extra_specs(cfg):
        if tokens:
            claimants.setdefault(seg, []).append(
                (label, tokens, bounds.get((label, seg)), False)
            )
    return ExtraPickContext(
        hmm_active=hmm_active,
        overlap_tolerance=overlap_tolerance_from_cfg(cfg),
        bounds=bounds,
        claimants=claimants,
        marker_ctx=marker_ctx,
        require_met_start=require_met_start_active(cfg),
    )


def pick_extra_cds(
    proteins: list[dict],
    family_label: str,
    hmm_tokens: list[str],
    aliases: list[str],
    segment: Optional[str],
    ctx: ExtraPickContext,
    exempt: Exempt = False,
) -> tuple[Optional[dict], list[tuple[dict, str]]]:
    """Pick one CDS for an ``extra_protein`` spec; return ``(cds, set_aside)``.

    ``exempt`` (from :meth:`ExtraPickContext.exempt`; a bool, or a
    zero-argument callable evaluated only when needed) only changes which
    CDS the cluster markers claim — see :meth:`ExtraPickContext.claimed`.

    Same chain as :func:`pick_marker_cds` (HMM-first when the spec has
    ``hmms:`` and the tier ran, else alias) plus the two guards described
    above. ``set_aside`` lists ``(cds, reason)`` for every CDS that WOULD
    have been a candidate for this spec (satisfies its tokens / matches an
    alias) but was withheld — so a caller can surface likely GenBank
    mislabels to the user. Callers that only need the pick ignore it.
    """
    with_seq = [p for p in (proteins or []) if p.get("sequence")]
    if not with_seq:
        return None, []
    bounds = ctx.bounds.get((family_label, segment))
    met_required = ctx.require_met_start

    def _withheld(p: dict) -> Optional[str]:
        """Why this CDS may not be picked for the spec, or None.

        Both restrictions apply to accessory proteins exactly as to cluster
        markers (parity; the designer's call 2026-10-29) — only the
        CONSEQUENCE differs, since a missing accessory protein is simply
        absent rather than a drop.
        """
        if not _within_bounds(p, bounds):
            return f"outside length_aa {format_length_bounds(bounds)}"
        if met_required and not starts_with_met(p):
            return "5'-partial CDS (translation does not start with Met)"
        return None

    eligible = [p for p in with_seq if _withheld(p) is None]
    elig_ids = {id(p) for p in eligible}
    set_aside: list[tuple[dict, str]] = []

    if hmm_tokens and ctx.hmm_active:
        cds = pick_marker_cds(
            eligible, hmm_tokens, [], True, overlap_tolerance=ctx.overlap_tolerance,
        )
        if bounds or met_required:
            for p in with_seq:
                if id(p) in elig_ids:
                    continue
                if pick_marker_cds(
                    [p], hmm_tokens, [], True,
                    overlap_tolerance=ctx.overlap_tolerance,
                ) is not None:
                    set_aside.append((p, _withheld(p) or ""))
        return cds, set_aside

    # Alias path: an alias-only spec, or an HMM spec on a run where the HMM
    # tier did not run (then nothing is claimed — behaviour as before).
    claimed = ctx.claimed(with_seq, segment, exempt)
    pool = [p for p in eligible if id(p) not in claimed]
    cds = pick_marker_cds(pool, [], aliases, False)
    for p in with_seq:
        if not _matches_any_alias(p, aliases):
            continue
        if id(p) in claimed:
            set_aside.append((p, f"already selected by HMM spec {claimed[id(p)]}"))
        elif id(p) not in elig_ids:
            set_aside.append((p, _withheld(p) or ""))
    return cds, set_aside


def collect_family_specs(
    cfg: dict[str, Any],
) -> list[tuple[str, list[str], Optional[str]]]:
    """Return ``[(family_label, tokens, segment_or_None), …]``, one per spec.

    A "family" is a marker **spec**, not a single token: the tokens in one
    spec's ``hmms:`` list are alternative architectures (OR) that all feed
    one tree. ``segment`` scopes the CDS search: in segmented mode a spec
    declared on segment ``S`` only ever pulls CDS from the ``S`` segment of
    each isolate; ``None`` (non-segmented) searches the representative's own
    proteins. Declaration order is preserved; tokens are de-duplicated
    within a spec.
    """
    specs: list[tuple[str, list[str], Optional[str]]] = []

    segmented = bool((cfg.get("segmented", {}) or {}).get("enabled"))
    if segmented:
        virus = get_virus_config(cfg) or {}
        segment_markers = virus.get("segment_markers") or {}
        cluster_protein = virus.get("cluster_protein") or {}
        segments = list(virus.get("segments") or [])
        # Include any segment that only appears under the marker maps.
        for extra in (*segment_markers.keys(), *cluster_protein.keys()):
            if extra not in segments:
                segments.append(extra)
        for seg in segments:
            name, tokens = _resolve_segment_marker(
                seg, segment_markers, cluster_protein
            )
            tokens = _dedup(tokens)
            if not tokens:
                continue
            specs.append((_family_label(name, tokens, seg), tokens, seg))
    else:
        for entry in (cfg.get("clustering", {}) or {}).get("cluster_protein", []) or []:
            if not isinstance(entry, dict):
                continue
            tokens = _dedup(list(entry.get("hmms") or []))
            if not tokens:
                continue
            specs.append((_family_label(entry.get("name"), tokens, None), tokens, None))
    return specs


def _segment_proteins(rep: Sequence, segment: Optional[str]) -> list[dict]:
    """Proteins to search for ``rep`` under family ``segment`` scope.

    Non-segmented (``segment is None``): the rep's own ``proteins``.
    Segmented: the matching entry in ``concat_segments`` (the per-segment
    Sequences carried on a CONCAT rep); falls back to the rep itself when
    a non-CONCAT rep happens to carry that segment label.
    """
    if segment is None:
        return rep.proteins or []
    for seg in rep.concat_segments or []:
        if seg.segment == segment:
            return seg.proteins or []
    if rep.segment == segment:
        return rep.proteins or []
    return []


def _best_satisfying_cds_any(
    proteins: list[dict],
    parsed_tokens: list[list[str]],
    *,
    overlap_tolerance: int = 0,
) -> Optional[dict]:
    """The CDS that best satisfies **any** of the alternative tokens (OR).

    ``parsed_tokens`` is a list of already-parsed token HMM-lists (the
    spec's alternative architectures). A CDS qualifies when it satisfies at
    least one of them; its score is the best (lowest) worst-domain E-value
    across the tokens it satisfies. "Best" CDS = longest translation
    (mirrors ``select_marker_protein``'s default), tie-broken by that
    score. Proteins without a translation can't seed a tree leaf and are
    skipped.

    ``overlap_tolerance`` is forwarded to :func:`cds_satisfies_token` so
    Pfam-boundary fuzz at multidomain seams (e.g. coronavirus S1/S2 across
    the furin site) doesn't reject biologically valid CDSes.
    """
    candidates: list[tuple[dict, float]] = []
    for prot in proteins:
        seq = prot.get("sequence")
        if not seq:
            continue
        hits = prot.get("hmm_hits") or []
        best_e: Optional[float] = None
        for token_hmms in parsed_tokens:
            worst_e = cds_satisfies_token(
                hits, token_hmms, overlap_tolerance=overlap_tolerance,
            )
            if worst_e is not None:
                best_e = worst_e if best_e is None else min(best_e, worst_e)
        if best_e is not None:
            candidates.append((prot, best_e))
    if not candidates:
        return None
    candidates.sort(
        key=lambda pe: (-(pe[0].get("length") or len(pe[0]["sequence"])), pe[1])
    )
    return candidates[0][0]


def _best_satisfying_cds(
    proteins: list[dict],
    token_hmms: list[str],
    *,
    overlap_tolerance: int = 0,
) -> Optional[dict]:
    """Single-token convenience wrapper over :func:`_best_satisfying_cds_any`.

    Kept for direct callers / unit tests that pass one parsed token.
    """
    return _best_satisfying_cds_any(
        proteins, [token_hmms], overlap_tolerance=overlap_tolerance,
    )


def overlap_tolerance_from_cfg(cfg: dict[str, Any]) -> int:
    """Read ``hmm.multidomain_overlap_tolerance`` from cfg, defaulting to 30.

    Default matches ``DEFAULTS`` so a hand-crafted cfg (e.g. in a test or
    legacy entry point) that omits the key still gets a sensible v0.22+
    tolerance rather than the pre-v0.22 strict-zero behaviour. Callers
    that genuinely want strict-zero should set the key to 0 explicitly.
    """
    return int(((cfg or {}).get("hmm", {}) or {}).get(
        "multidomain_overlap_tolerance", 30,
    ))


def pick_marker_cds(
    proteins: list[dict],
    hmm_tokens: list[str],
    aliases: list[str],
    hmm_active: bool,
    *,
    overlap_tolerance: int = 0,
    bounds: Optional[LengthBounds] = None,
    exempt: Exempt = False,
    require_met_start: bool = False,
) -> Optional[dict]:
    """Pick one CDS to represent a marker spec on one sequence/segment.

    ``bounds`` (the spec's ``length_aa``) and ``require_met_start``
    (``qc.marker_completeness``) restrict the candidates before the pick,
    exactly as the clustering picker does. ``exempt`` (the rep is protected
    against the ``marker_quality`` QC stage) falls back to the unrestricted
    pick when nothing survives both — again matching clustering, where one
    protect token waives both restrictions.
    """
    with_seq = [p for p in (proteins or []) if p.get("sequence")]
    if not with_seq:
        return None
    if bounds or require_met_start:
        pick = _pick_marker_cds_core(
            [
                p for p in with_seq
                if marker_candidate_ok(
                    p, bounds, require_met_start=require_met_start,
                )
            ],
            hmm_tokens, aliases, hmm_active, overlap_tolerance,
        )
        if pick is None and _is_exempt(exempt):
            pick = _pick_marker_cds_core(
                with_seq, hmm_tokens, aliases, hmm_active, overlap_tolerance,
            )
        return pick
    return _pick_marker_cds_core(
        with_seq, hmm_tokens, aliases, hmm_active, overlap_tolerance,
    )


def _pick_marker_cds_core(
    with_seq: list[dict],
    hmm_tokens: list[str],
    aliases: list[str],
    hmm_active: bool,
    overlap_tolerance: int,
) -> Optional[dict]:
    """The unbounded marker pick (body of :func:`pick_marker_cds`).

    Mirrors :func:`repseq.clustering.marker.select_marker_protein`'s chain:
    if the spec declares ``hmms:`` AND the HMM tier ran, prefer the longest
    CDS satisfying any token; otherwise fall back to the longest CDS whose
    ``/product`` substring-matches any alias (case-insensitive). Proteins
    without a translation are skipped — they cannot seed a FASTA record
    or a tree leaf.

    Used by both the always-on per-protein/extra-protein FASTA writers and
    the ``--per-protein-phylo`` tree builder so the artifact set is
    coherent (the FASTA's CDS is the tree's CDS).
    """
    if not with_seq:
        return None

    if hmm_tokens and hmm_active:
        parsed: list[list[str]] = []
        for token in hmm_tokens:
            try:
                parsed.append(parse_hmm_token(token))
            except ValueError:
                continue
        if parsed:
            hit = _best_satisfying_cds_any(
                with_seq, parsed, overlap_tolerance=overlap_tolerance,
            )
            if hit is not None:
                return hit
        # HMM tier active but nothing satisfied — DO NOT fall through to
        # aliases (matches select_marker_protein: HMM specs are authoritative).
        return None

    # Alias-only spec, or HMM tier inactive for this run.
    for alias in aliases:
        needle = (alias or "").lower().strip()
        if not needle:
            continue
        matches = [
            p for p in with_seq
            if needle in (p.get("product") or "").lower()
        ]
        if matches:
            return max(matches, key=lambda p: len(p["sequence"]))
    return None


def _min_taxa(cfg: dict[str, Any]) -> int:
    """Per-family minimum, never below the tree-builder floor of 3."""
    raw = ((cfg or {}).get("phylo", {}) or {}).get("per_protein", {}) or {}
    try:
        return max(3, int(raw.get("min_taxa", 3)))
    except (TypeError, ValueError):
        return 3


def _per_protein_mafft(cfg: dict[str, Any]) -> tuple[Optional[list[str]], bool]:
    """MAFFT args for the per-protein alignments → ``(extra_args, use_auto)``.

    These single-gene alignments are small enough to afford high-accuracy
    L-INS-i, so the default is ``--maxiterate 1000 --localpair``. When the
    per-protein ``mafft.extra_args`` list is non-empty it is used verbatim
    with ``--auto`` OFF (so the explicit strategy takes effect); when
    empty, returns ``(None, True)`` to fall back to the 2E behaviour
    (``--auto`` + ``phylo.mafft.extra_args``).
    """
    m = (
        ((cfg or {}).get("phylo", {}) or {}).get("per_protein", {}) or {}
    ).get("mafft", {}) or {}
    args = list(m.get("extra_args", []) or [])
    if args:
        return args, False
    return None, True


def _hmm_tier_ran(cfg: dict[str, Any], representatives: list[Sequence]) -> bool:
    """Did the HMM scan actually annotate any protein this session?

    Prefer the runtime flag the pipeline stashes (``_hmm_runtime.active``);
    fall back to scanning for any populated ``hmm_hits`` so a warm cache
    or an unusual entry point is still recognised.
    """
    if ((cfg or {}).get("_hmm_runtime") or {}).get("active"):
        return True
    for rep in representatives:
        pools = [rep.proteins or []]
        pools.extend((seg.proteins or []) for seg in rep.concat_segments or [])
        for pool in pools:
            if any(p.get("hmm_hits") for p in pool):
                return True
    return False


def _build_specs_trees(
    specs: list[tuple[str, list[str], list[str], Optional[str]]],
    representatives: list[Sequence],
    cfg: dict[str, Any],
    sub_dir: Path,
    prefix: str,
    label: str,
    *,
    color_scheme,
    pp_mafft_args,
    pp_mafft_auto: bool,
    emit_domains: bool,
    pp_cfg: dict[str, Any],
    min_taxa: int,
    hmm_active: bool,
    basis_role: str = "marker",
    extra_ctx: Optional[ExtraPickContext] = None,
) -> tuple[list[Path], list[str], list[str]]:
    """Build one tree per spec into ``sub_dir``; shared by cluster + extra.

    ``basis_role`` is the :func:`repseq.phylo.basis.describe_tree_basis`
    role for these trees — ``"marker"`` for the cluster-driving markers
    (2F) and ``"extra_protein"`` for the accessory-protein trees.

    Spec form: ``(family_label, hmm_tokens, aliases, segment_or_None)`` —
    same shape ``collect_marker_specs`` / ``collect_extra_specs`` return.
    CDS picking honours :func:`pick_marker_cds` (HMM-first, alias-fallback)
    so a single tree leaf matches its FASTA-writer record. For the
    extra_protein group pass ``extra_ctx``: picks then go through
    :func:`pick_extra_cds` (claiming + ``length_aa``), exactly as the
    extra FASTA writer and the protein report do.

    Returns ``(files_written, built_labels, sparse_summaries)``. ``label``
    appears in log messages so cluster vs. extra runs are distinguishable.
    """
    written: list[Path] = []
    built_labels: list[str] = []
    sparse: list[str] = []
    tol = overlap_tolerance_from_cfg(cfg)
    marker_ctx = build_marker_bounds_context(cfg) if extra_ctx is None else None

    for family_label, tokens, aliases, segment in specs:
        leaf_reps: list[Sequence] = []
        bodies: dict[str, str] = {}
        leaf_protein_ids: dict[str, set[str]] = {}
        for rep in representatives:
            if extra_ctx is not None:
                cds, _set_aside = pick_extra_cds(
                    _segment_proteins(rep, segment),
                    family_label, tokens, aliases, segment, extra_ctx,
                    partial(extra_ctx.exempt, rep),
                )
            else:
                cds = marker_ctx.pick(
                    rep, _segment_proteins(rep, segment),
                    family_label, tokens, aliases, segment, hmm_active, tol,
                )
            if cds is None:
                continue
            leaf_reps.append(rep)
            bodies[rep.id] = cds["sequence"]
            pid = cds.get("protein_id")
            if pid:
                leaf_protein_ids[rep.id] = {pid}

        if len(leaf_reps) < min_taxa:
            sparse.append(f"{family_label} ({len(leaf_reps)})")
            logger.info(
                "[%s] family %s: %d representative(s) carry it "
                "(< %d) — skipped", label, family_label, len(leaf_reps), min_taxa,
            )
            continue

        logger.info(
            "[%s] building %s from %d representative(s)…",
            label, family_label, len(leaf_reps),
        )
        try:
            files = _build_tree(
                leaf_reps,
                bodies,
                is_protein=True,
                cfg=cfg,
                out_dir=sub_dir,
                file_prefix=family_label,
                xml_name_prefix=f"{prefix}_{family_label}",
                color_scheme=color_scheme,
                leaf_protein_ids=leaf_protein_ids,
                mafft_extra_args=pp_mafft_args,
                mafft_use_auto=pp_mafft_auto,
                domain_architecture=emit_domains,
                trimal_settings=pp_cfg.get("trimal"),
                basis_role=basis_role,
                basis_family=family_label,
                basis_segment=segment,
                basis_architecture=(" OR ".join(tokens) if tokens else None),
            )
        except PhyloError as exc:
            logger.warning("[%s] family %s failed: %s", label, family_label, exc)
            continue
        written.extend(files)
        built_labels.append(family_label)

    return written, built_labels, sparse


def run_per_protein_phylogeny(
    representatives: list[Sequence],
    cfg: dict[str, Any],
    out_dir: Path,
    prefix: str,
) -> list[Path]:
    """Build one tree per declared marker spec AND extra_protein spec AND
    polyprotein peptide.

    Three output groups, distinct subdirectories:

    * marker (cluster_protein / segment_markers) trees → ``{prefix}_per_protein/``
    * extra_protein trees → ``{prefix}_extra_protein/``
    * polyprotein peptide trees → ``{prefix}_polyprotein/`` (alongside the
      unaligned peptide FASTAs from :func:`write_polyprotein_outputs`)

    The whole step soft-fails (``PhyloError``) only when there is genuinely
    nothing to do across all three groups (no specs configured anywhere,
    or no group has any family/peptide clearing ``min_taxa``). The groups
    are independent — a sparse extra_protein family does not stop the
    marker trees from being emitted, and a sparse peptide does not stop
    the marker or extra trees.
    """
    from ..polyprotein import collect_polyprotein_specs as _collect_poly_specs

    marker_specs = collect_family_specs(cfg)
    # Convert marker specs to the 4-tuple form pick_marker_cds expects.
    # Marker specs are HMM-only in collect_family_specs by design (the gate
    # is authoritative for clustering markers); pass empty aliases.
    marker_specs_full: list[tuple[str, list[str], list[str], Optional[str]]] = [
        (lab, tokens, [], seg) for lab, tokens, seg in marker_specs
    ]
    extra_specs = collect_extra_specs(cfg)
    poly_specs = _collect_poly_specs(cfg)

    if not marker_specs_full and not extra_specs and not poly_specs:
        raise PhyloError(
            "nothing to build: no cluster_protein / segment_markers / "
            "extra_protein / polyprotein specs declare aliases, HMM tokens, "
            "or peptides"
        )

    hmm_active = _hmm_tier_ran(cfg, representatives)
    # Marker specs are HMM-only — without the HMM tier there's no way to
    # assign CDS for them. Extra specs may be alias-only (built without HMM).
    needs_hmm_for_anything = (
        bool(marker_specs_full)
        or any(t and not a for _, t, a, _ in extra_specs)
    )
    has_any_alias_extra = any(a for _, _, a, _ in extra_specs)
    if not hmm_active and needs_hmm_for_anything and not has_any_alias_extra:
        raise PhyloError(
            "the HMM tier did not run this session (need hmm.enabled plus "
            "configured hmms:), so no CDS can be assigned to a protein family"
        )

    min_taxa = _min_taxa(cfg)
    written: list[Path] = []

    # Build the leaf-colour palette once over the FULL representative set
    # so a taxon keeps the same colour across every per-protein tree, every
    # extra-protein tree, AND the whole-genome tree (2E) — what makes
    # cross-tree incongruence legible at a glance.
    color_scheme = build_color_scheme(representatives, cfg)

    # High-accuracy MAFFT (L-INS-i by default) for these single-gene
    # alignments; resolved once and applied to every family in either group.
    pp_mafft_args, pp_mafft_auto = _per_protein_mafft(cfg)
    # Render each leaf protein's HMM hits as a phyloXML <domain_architecture>
    # (Archaeopteryx draws the domain boxes; its E-value slider filters them).
    pp_cfg = ((cfg or {}).get("phylo", {}) or {}).get("per_protein", {}) or {}
    emit_domains = bool(pp_cfg.get("domain_architecture", True))

    marker_files: list[Path] = []
    marker_built: list[str] = []
    marker_sparse: list[str] = []
    if marker_specs_full:
        marker_sub_dir = out_dir / f"{prefix}_per_protein"
        marker_files, marker_built, marker_sparse = _build_specs_trees(
            marker_specs_full, representatives, cfg, marker_sub_dir, prefix,
            "per-protein",
            color_scheme=color_scheme,
            pp_mafft_args=pp_mafft_args,
            pp_mafft_auto=pp_mafft_auto,
            emit_domains=emit_domains,
            pp_cfg=pp_cfg,
            min_taxa=min_taxa,
            hmm_active=hmm_active,
        )
        written.extend(marker_files)
        if marker_built:
            inc_path = _write_incongruence_table(
                marker_built, marker_sub_dir, out_dir, prefix, cfg,
            )
            if inc_path is not None:
                written.append(inc_path)

    extra_files: list[Path] = []
    extra_built: list[str] = []
    extra_sparse: list[str] = []
    if extra_specs:
        extra_sub_dir = out_dir / f"{prefix}_extra_protein"
        extra_files, extra_built, extra_sparse = _build_specs_trees(
            extra_specs, representatives, cfg, extra_sub_dir, prefix,
            "extra-protein",
            color_scheme=color_scheme,
            pp_mafft_args=pp_mafft_args,
            pp_mafft_auto=pp_mafft_auto,
            emit_domains=emit_domains,
            pp_cfg=pp_cfg,
            min_taxa=min_taxa,
            hmm_active=hmm_active,
            basis_role="extra_protein",
            extra_ctx=build_extra_pick_context(cfg, hmm_active),
        )
        written.extend(extra_files)
        # extra_protein trees deliberately do NOT participate in the marker
        # incongruence table — they live in a separate sub_dir and may be
        # sparse / accessory by design. If the user wants RF distances
        # against them, we can add a separate table in a later iteration.

    # Polyprotein peptide trees — third independent group.
    poly_files: list[Path] = []
    poly_built = False
    poly_skip_reason: Optional[str] = None
    if poly_specs:
        try:
            poly_files = run_polyprotein_phylogeny(
                representatives, cfg, out_dir, prefix,
            )
            poly_built = bool(poly_files)
            written.extend(poly_files)
        except PhyloError as exc:
            poly_skip_reason = str(exc)
            logger.warning("[polyprotein] %s", exc)

    if not marker_built and not extra_built and not poly_built:
        detail_parts: list[str] = []
        if marker_sparse:
            detail_parts.append("marker families: " + ", ".join(marker_sparse))
        if extra_sparse:
            detail_parts.append("extra_protein families: " + ", ".join(extra_sparse))
        if poly_skip_reason:
            detail_parts.append(f"polyprotein peptides: {poly_skip_reason}")
        detail = f" ({'; '.join(detail_parts)})" if detail_parts else ""
        raise PhyloError(
            f"no protein family or peptide had >= {min_taxa} representatives "
            f"carrying the declared architecture/aliases — nothing built{detail}"
        )

    return written


def _build_peptide_synthetic_cds(
    parent_cds: dict, sliced, peptide_name: str,
) -> dict:
    """Build a synthetic CDS dict representing one mature peptide.

    Hits that lie **mostly within** the peptide range are kept and
    **shifted to peptide-local coordinates** (1-based, range
    1..length_aa) so the phyloXML ``<domain_architecture>`` draws
    domain boxes on the peptide's own span rather than on the
    polyprotein's. Multidomain peptide tokens (e.g. NSP12 =
    ``CoV_RPol_N--RdRP_1``) get both subdomain boxes visible;
    single-HMM peptides get one.

    "Mostly within" = more than 50% of the original ``ali_span`` falls
    inside the peptide range. Hits where ≤50% lies in the peptide are
    dropped — those are Pfam-boundary fuzz from the *neighbouring*
    peptide's domain that the slicer would otherwise emit as
    1-or-2-aa artifact boxes at the peptide's start or end. A hit truly
    belonging to this peptide always passes (100% inside ≫ 50%), and a
    genuinely 50/50-straddling hit gets dropped from both adjacent
    peptides — honest, since such a hit can't be biologically assigned
    to either mature peptide.
    """
    pep_from = sliced.range_aa_from
    pep_to = sliced.range_aa_to
    shifted_hits: list[dict] = []
    for h in (parent_cds.get("hmm_hits") or []):
        try:
            h_from = int(h.get("ali_from", 0))
            h_to = int(h.get("ali_to", 0))
        except (TypeError, ValueError):
            continue
        if h_to < pep_from or h_from > pep_to:
            continue
        original_span = h_to - h_from + 1
        if original_span <= 0:
            continue
        clipped_from = max(h_from, pep_from)
        clipped_to = min(h_to, pep_to)
        clipped_span = clipped_to - clipped_from + 1
        # Drop boundary-spillover artifacts: a hit must have a majority
        # of its residues inside this peptide to be considered a real
        # domain of the peptide. Strict ">"; equal/below the half-mark
        # is ambiguous and dropped from both sides.
        if clipped_span * 2 <= original_span:
            continue
        new_from = clipped_from - (pep_from - 1)
        new_to = clipped_to - (pep_from - 1)
        shifted = dict(h)
        shifted["ali_from"] = new_from
        shifted["ali_to"] = new_to
        shifted["ali_span"] = clipped_span
        shifted_hits.append(shifted)
    parent_pid = parent_cds.get("protein_id") or "polyprotein"
    return {
        "protein_id": f"{parent_pid}:{peptide_name}",
        "product": peptide_name,
        "length": sliced.length_aa,
        "sequence": sliced.sequence,
        "hmm_hits": shifted_hits,
    }


def run_polyprotein_phylogeny(
    representatives: list[Sequence],
    cfg: dict[str, Any],
    out_dir: Path,
    prefix: str,
) -> list[Path]:
    """2F-equivalent for mature peptides: one tree per polyprotein peptide.

    For each declared ``clustering.polyprotein:`` / ``virus.polyprotein:``
    spec, and each peptide within it, slice every representative's
    polyprotein CDS and build a phylogenetic tree from the spliced
    peptide AA sequences (``ok``-status slices only — peptides whose
    token wasn't satisfied are simply absent from their tree's input).
    All knobs come from ``phylo.per_protein`` so the peptide trees
    inherit the same MAFFT / trimAl / IQ-TREE-or-FastTree / rooting /
    LCA / colouring / domain-architecture settings as 2F.

    Outputs land alongside the unaligned peptide FASTAs in
    ``{prefix}_polyprotein/`` with a matching basename schema:

    * ``{prefix}_<spec>_<peptide>_msa.fasta`` (MAFFT alignment)
    * ``{prefix}_<spec>_<peptide>_tree.nwk`` (Newick)
    * ``{prefix}_<spec>_<peptide>_tree.xml`` (phyloXML)
    * ``{prefix}_<spec>_<peptide>_tree_id_map.tsv`` (id map)

    When ``phylo.per_protein.whole_polyprotein_tree`` is true, ALSO build
    one tree per spec on the **whole polyprotein CDS** (not the sliced
    peptides) — each rep's parent CDS is fed directly through the same
    engine, so the phyloXML ``<domain_architecture>`` shows every HMM hit
    across the entire polyprotein (the full nsp/peptide layout end-to-end).
    Outputs as ``{prefix}_<spec>_polyprotein_tree.{xml,pdf,png}`` etc. in
    the same directory; OFF by default (one extra large alignment + tree
    per spec). The unaligned ``{prefix}_<spec>_polyprotein.fasta`` is
    written independently by :func:`output.report.write_polyprotein_outputs`
    (always-on, parallel to the per-peptide FASTAs).

    Sparse peptides (fewer than ``phylo.per_protein.min_taxa`` ``ok``
    slices) are skipped with a log note; per-peptide build failures log
    and continue (the whole step never aborts a run).

    Peptide trees are **not** added to the per-protein incongruence
    table (that compares whole-genome markers); peptides within a
    polyprotein answer a different scientific question and live in
    their own subdirectory.

    Raises :class:`PhyloError` only when there's nothing to build (no
    specs configured, HMM tier didn't run, or every peptide is sparse).
    """
    from ..polyprotein import collect_polyprotein_specs, slice_polyprotein

    specs = collect_polyprotein_specs(cfg)
    if not specs:
        raise PhyloError(
            "no polyprotein specs configured "
            "(clustering.polyprotein / virus.polyprotein)"
        )

    hmm_active = _hmm_tier_ran(cfg, representatives)
    if not hmm_active:
        raise PhyloError(
            "the HMM tier did not run this session "
            "(need hmm.enabled plus configured hmms:), "
            "so polyprotein peptides can't be located"
        )

    min_taxa = _min_taxa(cfg)
    color_scheme = build_color_scheme(representatives, cfg)
    pp_mafft_args, pp_mafft_auto = _per_protein_mafft(cfg)
    pp_cfg = ((cfg or {}).get("phylo", {}) or {}).get("per_protein", {}) or {}
    emit_domains = bool(pp_cfg.get("domain_architecture", True))
    build_whole = bool(pp_cfg.get("whole_polyprotein_tree", False))
    tol = overlap_tolerance_from_cfg(cfg)

    sub_dir = out_dir / f"{prefix}_polyprotein"
    written: list[Path] = []
    built: list[str] = []
    sparse: list[str] = []

    from dataclasses import replace as _dc_replace

    def _proteins_for_segment(rep: Sequence, segment: Optional[str]) -> list[dict]:
        """Local copy of ``output.report._polyprotein_proteins`` to avoid
        reaching from phylo/ back into output/. Kept identical so the
        slicer sees the same CDS list the writer does.
        """
        if segment is None:
            return rep.proteins or []
        for seg in rep.concat_segments or []:
            if seg.segment == segment:
                return seg.proteins or []
        if rep.segment == segment:
            return rep.proteins or []
        return []

    for spec in specs:
        # Slice once per rep, then walk peptides — the slicer's parent
        # identification depends on the whole peptide set.
        per_rep_slices: dict[str, list] = {}
        per_rep_parent: dict[str, dict] = {}
        for rep in representatives:
            proteins = _proteins_for_segment(rep, spec.segment)
            parent_cds, sliced_list = slice_polyprotein(
                proteins, spec, overlap_tolerance=tol,
            )
            per_rep_slices[rep.id] = sliced_list
            if parent_cds is not None:
                per_rep_parent[rep.id] = parent_cds

        for pep_idx, pep in enumerate(spec.peptides):
            family_label = f"{_sanitize(spec.name)}_{_sanitize(pep.name)}"

            leaf_reps: list[Sequence] = []
            bodies: dict[str, str] = {}
            leaf_protein_ids: dict[str, set[str]] = {}
            for rep in representatives:
                sliced_list = per_rep_slices.get(rep.id) or []
                if pep_idx >= len(sliced_list):
                    continue
                sliced = sliced_list[pep_idx]
                if sliced.status not in ("ok", "overlap") or not sliced.sequence:
                    continue
                parent_cds = per_rep_parent.get(rep.id)
                if parent_cds is None:
                    continue
                synthetic = _build_peptide_synthetic_cds(
                    parent_cds, sliced, pep.name,
                )
                # Shallow-copy the rep with the synthetic peptide as its
                # only protein, so the phyloXML leaf shows the peptide
                # name + peptide-local domain architecture. Original
                # rep objects are NOT mutated.
                synth_rep = _dc_replace(
                    rep,
                    proteins=[synthetic],
                    concat_segments=None,
                    marker_protein_ids=[synthetic["protein_id"]],
                )
                leaf_reps.append(synth_rep)
                bodies[rep.id] = sliced.sequence
                leaf_protein_ids[rep.id] = {synthetic["protein_id"]}

            if len(leaf_reps) < min_taxa:
                sparse.append(f"{family_label} ({len(leaf_reps)})")
                logger.info(
                    "[polyprotein] peptide %s: %d rep(s) with ok slices "
                    "(< %d) — skipped", family_label, len(leaf_reps), min_taxa,
                )
                continue

            logger.info(
                "[polyprotein] building %s tree from %d rep(s)…",
                family_label, len(leaf_reps),
            )
            try:
                files = _build_tree(
                    leaf_reps,
                    bodies,
                    is_protein=True,
                    cfg=cfg,
                    out_dir=sub_dir,
                    file_prefix=family_label,
                    xml_name_prefix=f"{prefix}_{family_label}",
                    color_scheme=color_scheme,
                    leaf_protein_ids=leaf_protein_ids,
                    mafft_extra_args=pp_mafft_args,
                    mafft_use_auto=pp_mafft_auto,
                    domain_architecture=emit_domains,
                    trimal_settings=pp_cfg.get("trimal"),
                    basis_role="peptide",
                    basis_family=pep.name,
                    basis_parent=spec.name,
                    basis_segment=spec.segment,
                )
            except PhyloError as exc:
                logger.warning(
                    "[polyprotein] peptide %s failed: %s", family_label, exc,
                )
                continue
            written.extend(files)
            built.append(family_label)

        # Whole-polyprotein tree (opt-in, phylo.per_protein.whole_polyprotein_tree):
        # one tree per spec on the ENTIRE polyprotein CDS, not the sliced
        # peptides. Same per-spec loop minus the slicing — feed each rep's
        # parent CDS directly, so the phyloXML <domain_architecture> shows
        # every HMM hit across the whole polyprotein (the full nsp/peptide
        # layout end-to-end). Reuses the same _build_tree engine, colour
        # palette, MAFFT/trimAl settings, and min_taxa floor as the peptide
        # trees above.
        if build_whole:
            family_label = f"{_sanitize(spec.name)}_polyprotein"
            leaf_reps = []
            bodies = {}
            leaf_protein_ids = {}
            for rep in representatives:
                parent_cds = per_rep_parent.get(rep.id)
                if parent_cds is None:
                    continue
                seq_str = parent_cds.get("sequence")
                if not seq_str:
                    continue
                parent_pid = parent_cds.get("protein_id") or "polyprotein"
                # Shallow-copy the rep with the parent CDS as its only
                # protein, so the leaf shows the whole-polyprotein name +
                # full-length domain architecture. Originals untouched.
                synth_rep = _dc_replace(
                    rep,
                    proteins=[parent_cds],
                    concat_segments=None,
                    marker_protein_ids=[parent_pid],
                )
                leaf_reps.append(synth_rep)
                bodies[rep.id] = seq_str
                leaf_protein_ids[rep.id] = {parent_pid}

            if len(leaf_reps) < min_taxa:
                sparse.append(f"{family_label} ({len(leaf_reps)})")
                logger.info(
                    "[polyprotein] whole-polyprotein %s: %d rep(s) with a "
                    "parent CDS (< %d) — skipped",
                    family_label, len(leaf_reps), min_taxa,
                )
            else:
                logger.info(
                    "[polyprotein] building whole-polyprotein %s tree from "
                    "%d rep(s)…", family_label, len(leaf_reps),
                )
                try:
                    files = _build_tree(
                        leaf_reps,
                        bodies,
                        is_protein=True,
                        cfg=cfg,
                        out_dir=sub_dir,
                        file_prefix=family_label,
                        xml_name_prefix=f"{prefix}_{family_label}",
                        color_scheme=color_scheme,
                        leaf_protein_ids=leaf_protein_ids,
                        mafft_extra_args=pp_mafft_args,
                        mafft_use_auto=pp_mafft_auto,
                        domain_architecture=emit_domains,
                        trimal_settings=pp_cfg.get("trimal"),
                        basis_role="polyprotein",
                        basis_family=spec.name,
                        basis_segment=spec.segment,
                    )
                    written.extend(files)
                    built.append(family_label)
                except PhyloError as exc:
                    logger.warning(
                        "[polyprotein] whole-polyprotein %s failed: %s",
                        family_label, exc,
                    )

    if not built:
        detail = f" ({'; '.join(sparse)})" if sparse else ""
        raise PhyloError(
            f"no polyprotein peptide or whole-polyprotein CDS had >= "
            f"{min_taxa} representatives — nothing built{detail}"
        )
    return written


def run_per_segment_phylogeny(
    representatives: list[Sequence],
    cfg: dict[str, Any],
    out_dir: Path,
    prefix: str,
) -> list[Path]:
    """2K — build one **nucleotide** tree per segment, into
    ``{prefix}_per_segment/``.

    Segmented mode only. Complements 2F (per-marker protein trees) with a
    raw per-segment NT view: every representative isolate contributes its
    per-segment NT sequence (from ``concat_segments``) to one tree per
    segment. Useful because reassortment may not be visible in any single
    marker (which is one CDS per segment) but shows up as topological
    incongruence between the per-segment trees themselves.

    Soft-fail contract mirrors 2F: raises :class:`PhyloError` only when
    the run is non-segmented, no segments are declared, or no segment has
    >= ``phylo.per_protein.min_taxa`` representatives carrying it; per-
    segment build failures log + skip individual segments. Uses the
    shared :func:`_build_tree` engine so rooting / LCA / colouring /
    phyloXML are identical to 2E and 2F. The leaf colour palette is
    keyed on the same FULL representative set so a taxon keeps the same
    colour across 2E, 2F, and 2K.
    """
    segmented = bool((cfg.get("segmented", {}) or {}).get("enabled"))
    if not segmented:
        raise PhyloError(
            "per-segment trees only apply to segmented runs"
        )

    virus = get_virus_config(cfg) or {}
    segments = list(virus.get("segments") or [])
    if not segments:
        raise PhyloError(
            "no segments declared in virus config (need "
            "segmented.viruses.<v>.segments)"
        )

    min_taxa = _min_taxa(cfg)
    color_scheme = build_color_scheme(representatives, cfg)
    pp_cfg = ((cfg or {}).get("phylo", {}) or {}).get("per_protein", {}) or {}

    sub_dir = out_dir / f"{prefix}_per_segment"
    written: list[Path] = []
    built: list[str] = []
    sparse: list[str] = []

    for seg_name in segments:
        leaf_reps: list[Sequence] = []
        bodies: dict[str, str] = {}
        for rep in representatives:
            seg_seq = None
            for sub in rep.concat_segments or []:
                if sub.segment == seg_name:
                    seg_seq = sub
                    break
            if seg_seq is None or not seg_seq.sequence:
                continue
            leaf_reps.append(rep)
            bodies[rep.id] = seg_seq.sequence

        if len(leaf_reps) < min_taxa:
            sparse.append(f"{seg_name} ({len(leaf_reps)})")
            logger.info(
                "[per-segment] %s: %d rep(s) carry it (< %d) — skipped",
                seg_name, len(leaf_reps), min_taxa,
            )
            continue

        logger.info(
            "[per-segment] building %s NT tree from %d rep(s)…",
            seg_name, len(leaf_reps),
        )
        try:
            files = _build_tree(
                leaf_reps,
                bodies,
                is_protein=False,  # per-segment NT — always nucleotide
                cfg=cfg,
                out_dir=sub_dir,
                file_prefix=_sanitize(seg_name),
                xml_name_prefix=f"{prefix}_{seg_name}",
                color_scheme=color_scheme,
                leaf_protein_ids=None,
                mafft_extra_args=None,
                mafft_use_auto=True,
                domain_architecture=False,
                trimal_settings=(cfg.get("phylo", {}) or {}).get("trimal"),
                basis_role="segment_nt",
                basis_segment=seg_name,
            )
        except PhyloError as exc:
            logger.warning("[per-segment] %s failed: %s", seg_name, exc)
            continue
        written.extend(files)
        built.append(seg_name)

    if not built:
        detail = f" ({'; '.join(sparse)})" if sparse else ""
        raise PhyloError(
            f"no segment had >= {min_taxa} representatives carrying it{detail}"
        )
    return written


def _write_incongruence_table(
    built_labels: list[str],
    sub_dir: Path,
    out_dir: Path,
    prefix: str,
    cfg: dict[str, Any],
) -> Optional[Path]:
    """Pairwise unrooted-RF table across the built family trees (+ the 2E
    whole-genome tree when ``--phylo`` wrote one). Soft-fail: returns the
    path written, or None when disabled / too few trees / on any error.

    Reads each tree's on-disk Newick + ``_tree_id_map.tsv`` (the genome
    tree's land in ``out_dir`` from the ``--phylo`` step, which always
    runs before this one), so no tree objects need threading in.
    """
    pp_cfg = ((cfg or {}).get("phylo", {}) or {}).get("per_protein", {}) or {}
    if not pp_cfg.get("incongruence", True):
        return None
    try:
        trees: list[tuple[str, Path, Path]] = [
            (
                label,
                sub_dir / f"{label}_tree.nwk",
                sub_dir / f"{label}_tree_id_map.tsv",
            )
            for label in built_labels
        ]
        genome_nwk = out_dir / f"{prefix}_tree.nwk"
        genome_map = out_dir / f"{prefix}_tree_id_map.tsv"
        if genome_nwk.exists() and genome_map.exists():
            trees.append(("GENOME", genome_nwk, genome_map))
        if len(trees) < 2:
            return None  # need a pair to compare
        from .incongruence import compute_incongruence, write_incongruence_tsv

        rows = compute_incongruence(trees)
        if not rows:
            return None
        inc_path = sub_dir / f"{prefix}_incongruence.tsv"
        write_incongruence_tsv(rows, inc_path)
        logger.info(
            "[per-protein] wrote pairwise incongruence table (%d pair(s)): %s",
            len(rows), inc_path,
        )
        return inc_path
    except Exception as exc:  # never let a metric bug void the trees
        logger.warning("[per-protein] incongruence table failed: %s", exc)
        return None
