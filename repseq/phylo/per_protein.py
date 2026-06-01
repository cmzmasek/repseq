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
from pathlib import Path
from typing import Any, Optional

from ..config import get_virus_config
from ..hmm.runner import cds_satisfies_token, parse_hmm_token
from ..models import Sequence
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
    """Like :func:`_resolve_segment_marker` but also returns ``aliases``.

    Returns ``(name, hmm_tokens, aliases)``. ``segment_markers`` wins over
    ``cluster_protein`` (legacy); for the legacy form (a list of strings or
    dicts) plain-string entries are treated as alias strings, and each
    dict entry contributes its ``aliases:`` / ``hmms:`` lists. When
    ``segment_markers[seg]`` exists but its ``hmms:`` is empty (alias-only),
    we additionally pull ``hmms:`` from any ``cluster_protein[seg]`` dict
    entries so a split alias/HMM config still gets a usable HMM gate —
    ``segment_markers``'s own aliases are preserved and merged with any
    from ``cluster_protein``.
    """
    sm_name: Optional[str] = None
    sm_tokens: list[str] = []
    sm_aliases: list[str] = []
    if seg_name in segment_markers:
        spec = segment_markers[seg_name] or {}
        sm_name = spec.get("name")
        sm_tokens = list(spec.get("hmms") or [])
        sm_aliases = list(spec.get("aliases") or [])
        if sm_tokens:
            return sm_name, sm_tokens, sm_aliases
        # Fall through to cluster_protein for HMMs (and any extra aliases).
    if seg_name in cluster_protein_per_seg:
        tokens: list[str] = list(sm_tokens)
        aliases: list[str] = list(sm_aliases)
        name: Optional[str] = sm_name
        for entry in cluster_protein_per_seg[seg_name] or []:
            if isinstance(entry, str):
                aliases.append(entry)
            elif isinstance(entry, dict):
                name = name or entry.get("name")
                tokens.extend(entry.get("hmms") or [])
                aliases.extend(entry.get("aliases") or [])
        return name, tokens, aliases
    return sm_name, sm_tokens, sm_aliases


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
    """
    specs: list[tuple[str, list[str], list[str], Optional[str]]] = []

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
            specs.append(
                (_family_label(name, tokens, seg), tokens, aliases, seg)
            )
    else:
        for entry in (cfg.get("clustering", {}) or {}).get("cluster_protein", []) or []:
            if not isinstance(entry, dict):
                continue
            tokens = _dedup(list(entry.get("hmms") or []))
            aliases = _dedup(list(entry.get("aliases") or []))
            if not tokens and not aliases:
                continue
            specs.append(
                (_family_label(entry.get("name"), tokens, None), tokens, aliases, None)
            )
    return specs


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
    """
    specs: list[tuple[str, list[str], list[str], Optional[str]]] = []

    segmented = bool((cfg.get("segmented", {}) or {}).get("enabled"))
    if segmented:
        virus = get_virus_config(cfg) or {}
        extra = virus.get("extra_protein") or {}
        for seg_name, entries in extra.items():
            for entry in entries or []:
                if not isinstance(entry, dict):
                    continue
                name = (entry.get("name") or "").strip()
                if not name:
                    continue
                tokens = _dedup(list(entry.get("hmms") or []))
                aliases = _dedup(list(entry.get("aliases") or []))
                if not tokens and not aliases:
                    continue
                specs.append((_sanitize(name), tokens, aliases, seg_name))
    else:
        for entry in (cfg.get("clustering", {}) or {}).get("extra_protein", []) or []:
            if not isinstance(entry, dict):
                continue
            name = (entry.get("name") or "").strip()
            if not name:
                continue
            tokens = _dedup(list(entry.get("hmms") or []))
            aliases = _dedup(list(entry.get("aliases") or []))
            if not tokens and not aliases:
                continue
            specs.append((_sanitize(name), tokens, aliases, None))
    return specs


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
) -> Optional[dict]:
    """Pick one CDS to represent a marker spec on one sequence/segment.

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
    with_seq = [p for p in (proteins or []) if p.get("sequence")]
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
) -> tuple[list[Path], list[str], list[str]]:
    """Build one tree per spec into ``sub_dir``; shared by cluster + extra.

    ``basis_role`` is the :func:`repseq.phylo.basis.describe_tree_basis`
    role for these trees — ``"marker"`` for the cluster-driving markers
    (2F) and ``"extra_protein"`` for the accessory-protein trees.

    Spec form: ``(family_label, hmm_tokens, aliases, segment_or_None)`` —
    same shape ``collect_marker_specs`` / ``collect_extra_specs`` return.
    CDS picking honours :func:`pick_marker_cds` (HMM-first, alias-fallback)
    so a single tree leaf matches its FASTA-writer record.

    Returns ``(files_written, built_labels, sparse_summaries)``. ``label``
    appears in log messages so cluster vs. extra runs are distinguishable.
    """
    written: list[Path] = []
    built_labels: list[str] = []
    sparse: list[str] = []
    tol = overlap_tolerance_from_cfg(cfg)

    for family_label, tokens, aliases, segment in specs:
        leaf_reps: list[Sequence] = []
        bodies: dict[str, str] = {}
        leaf_protein_ids: dict[str, set[str]] = {}
        for rep in representatives:
            cds = pick_marker_cds(
                _segment_proteins(rep, segment),
                tokens, aliases, hmm_active,
                overlap_tolerance=tol,
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

    if not built:
        detail = f" ({'; '.join(sparse)})" if sparse else ""
        raise PhyloError(
            f"no polyprotein peptide had >= {min_taxa} ok-sliced "
            f"representatives — nothing built{detail}"
        )
    return written


def run_per_segment_phylogeny(
    representatives: list[Sequence],
    cfg: dict[str, Any],
    out_dir: Path,
    prefix: str,
) -> list[Path]:
    """2H — build one **nucleotide** tree per segment, into
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
    colour across 2E, 2F, and 2H.
    """
    segmented = bool((cfg.get("segmented", {}) or {}).get("enabled"))
    if not segmented:
        raise PhyloError(
            "per-segment trees (2H) only apply to segmented runs"
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
