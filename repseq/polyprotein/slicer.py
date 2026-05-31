"""Cut a polyprotein CDS into its mature peptides.

The pipeline:

1. :func:`identify_parent_cds` picks the CDS that carries hits from the
   most distinct peptide HMMs declared on the spec, requiring at least
   :attr:`PolyproteinSpec.min_peptides_hit` distinct matches.
2. :func:`compute_cuts` lays out the inter-peptide boundary positions
   (1-based amino-acid coordinates, inclusive) using one of three cut
   strategies (``boundary`` / ``bisect`` / ``motif``).
3. :func:`slice_polyprotein` ties it together: parent identification,
   ordering check, cut computation, peptide-string extraction, and a
   :class:`SlicedPeptide` audit record per declared peptide.

Hit dict shape mirrors what
:func:`repseq.hmm.hmmscan._parse_domtblout` produces — see
``feedback-test-fixtures-match-production``: ``target`` is the HMM
profile name, ``ali_from``/``ali_to`` are 1-based inclusive coords on
the protein, ``dom_evalue`` is the domain E-value used to break ties
when multiple hits to the same HMM exist on one CDS.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from ..hmm.runner import parse_hmm_token
from .specs import PeptideSpec, PolyproteinSpec


# Status codes the audit TSV (and the FASTA header `[cut_method=...]`
# tag) consume. ``ok`` is the only one that emits a FASTA record.
_OK = "ok"
_MISSING = "missing"
_OUT_OF_ORDER = "out_of_order"
_OVERLAP = "overlap"
_NO_PARENT = "no_parent_cds"


@dataclass
class SlicedPeptide:
    """One peptide's worth of audit info for a single representative.

    ``range_aa_from`` / ``range_aa_to`` are 1-based, inclusive (matching
    ``ali_from``/``ali_to`` convention). ``status`` says whether the
    peptide is usable: only ``ok`` peptides get a FASTA record; the
    others land in the audit TSV with the reason.

    ``cut_method_actual`` records what the slicer actually did at the
    peptide's *start* boundary (e.g. ``motif:LQ`` if a motif snap fired
    at the N-terminal cut, ``bisect`` if the motif window came up empty
    and we fell back, ``boundary`` if the spec used hit-boundary
    slicing, ``n-term`` for the first peptide whose start is just the
    polyprotein N-terminus). Mostly a debugging / scientific-audit
    field; the bench scientist uses it to spot peptides whose cut is
    purely heuristic vs. motif-supported.
    """

    peptide_name: str
    parent_protein_id: Optional[str]
    parent_accession: Optional[str]
    range_aa_from: int
    range_aa_to: int
    length_aa: int
    sequence: str  # empty when status != "ok"
    cut_method_actual: str
    status: str  # "ok" | "missing" | "out_of_order" | "overlap" | "no_parent_cds"
    note: str = ""  # free-text detail for the audit TSV
    # The peptide-token alternative that was satisfied for this slice
    # (one entry from ``PeptideSpec.hmms``). When a peptide has multiple
    # OR alternatives, this records which one actually fired — essential
    # for the alpha/beta-CoV NSP1 case where the architecture varies by
    # genus. Empty string when the peptide wasn't located (status
    # `missing` / `out_of_order` / `no_parent_cds`).
    matched_token: str = ""


def _satisfying_span_for_token(
    hits: list[dict],
    token: str,
    *,
    overlap_tolerance: int = 0,
) -> Optional[tuple[int, int]]:
    """Locate a peptide token's footprint on one CDS.

    Returns ``(span_from, span_to)`` (1-based, inclusive AA coords) when
    every HMM named in the token has a hit AND — for multidomain tokens
    — the hits appear in N-to-C order along the CDS with at most
    ``overlap_tolerance`` aa of overlap at each seam. Returns ``None``
    when the token can't be satisfied on this CDS.

    The span is the synthetic union: ``min(ali_from of chosen hits)``,
    ``max(ali_to of chosen hits)``. For single-HMM tokens this is just
    the hit's own ``ali_from..ali_to``; for multidomain tokens it covers
    the whole architecture from the N-terminal domain's start to the
    C-terminal domain's end.

    Filters by the per-hit ``passing`` flag set by ``_run_hmm_qc`` (only
    hits clearing the configured ``default_evalue`` / GA cutoff and
    ``relative_length_cutoff`` are considered). Hits without a
    ``passing`` key are treated as passing — a backward-compat shim for
    callers / tests that don't go through ``_run_hmm_qc``; the
    production pipeline always sets the flag (cli.py ``_run_hmm_qc``).

    The walk logic mirrors :func:`repseq.hmm.runner.cds_satisfies_token`
    so the slicer enforces the same architecture rules as the
    cluster_protein HMM gate.
    """
    try:
        token_hmms = parse_hmm_token(token)
    except ValueError:
        return None
    if not token_hmms:
        return None

    by_target: dict[str, list[dict]] = {}
    for h in hits or []:
        target = h.get("target") or h.get("hmm_name") or h.get("name")
        if not target:
            continue
        if not h.get("passing", True):
            continue
        by_target.setdefault(target, []).append(h)

    for name in token_hmms:
        if name not in by_target:
            return None

    if len(token_hmms) == 1:
        # Single-HMM token: best hit = lowest dom_evalue (tie-broken by
        # longer ali_span, since a longer hit gives a more informative
        # cut location).
        best = min(
            by_target[token_hmms[0]],
            key=lambda h: (
                float(h.get("dom_evalue", float("inf"))),
                -int(h.get("ali_span") or (
                    int(h.get("ali_to", 0)) - int(h.get("ali_from", 0)) + 1
                )),
            ),
        )
        return int(best["ali_from"]), int(best["ali_to"])

    # Multidomain token — greedy left-to-right walk, same shape as
    # cds_satisfies_token. Tiebreak: among forward-progressing
    # candidates for each named HMM, pick the leftmost by ``ali_to``,
    # then break ties by best (lowest) ``dom_evalue`` — so the walk
    # prefers higher-confidence hits when two candidates end at the
    # same position (v0.36.1+; pre-v0.36.1 the tiebreak ignored hit
    # quality entirely).
    chosen: list[dict] = []
    prev_from = -1
    prev_to = -1
    for name in token_hmms:
        candidates = [
            h for h in by_target[name]
            if int(h["ali_from"]) > prev_from
            and int(h["ali_to"]) > prev_to
            and (
                prev_to < 0
                or int(h["ali_from"]) >= prev_to - overlap_tolerance + 1
            )
        ]
        if not candidates:
            return None
        pick = min(
            candidates,
            key=lambda h: (
                int(h["ali_to"]),
                float(h.get("dom_evalue", float("inf"))),
            ),
        )
        chosen.append(pick)
        prev_from = int(pick["ali_from"])
        prev_to = int(pick["ali_to"])
    span_from = min(int(h["ali_from"]) for h in chosen)
    span_to = max(int(h["ali_to"]) for h in chosen)
    return span_from, span_to


def _best_satisfying_alternative(
    hits: list[dict],
    tokens: list[str],
    *,
    overlap_tolerance: int = 0,
) -> Optional[tuple[int, int, str]]:
    """OR across alternative peptide-token architectures.

    Tries each token in ``tokens`` against ``hits``. Among the
    alternatives that are satisfied, returns the one with the best
    (lowest) worst-domain E-value across its chosen hits — mirroring
    how :func:`repseq.phylo.per_protein._best_satisfying_cds_any` ranks
    OR-alternatives for ``cluster_protein``. Ties are broken by the
    declared order in ``tokens`` (first-declared wins).

    Returns ``(span_from, span_to, matched_token)`` or ``None`` when no
    alternative is satisfied. ``matched_token`` is the verbatim entry
    from ``tokens`` so the audit row can record which architecture
    fired (e.g. ``"aCoV_NSP1"`` vs. ``"bCoV_NSP1"``).
    """
    best: Optional[tuple[float, int, tuple[int, int], str]] = None
    for i, token in enumerate(tokens):
        span = _satisfying_span_for_token(
            hits, token, overlap_tolerance=overlap_tolerance,
        )
        if span is None:
            continue
        # Worst-domain E-value across the chosen hits for this token —
        # the same scalar the cluster_protein selector ranks by.
        worst_e = _token_worst_evalue(hits, token, overlap_tolerance)
        if worst_e is None:
            # Token satisfied but we can't score it (shouldn't happen);
            # treat as least-preferred so a scorable alternative wins.
            worst_e = float("inf")
        key = (worst_e, i)
        if best is None or key < (best[0], best[1]):
            best = (worst_e, i, span, token)
    if best is None:
        return None
    _, _, span, token = best
    return span[0], span[1], token


def _token_worst_evalue(
    hits: list[dict], token: str, overlap_tolerance: int,
) -> Optional[float]:
    """Return the worst (largest) dom_evalue across the hits that
    satisfy ``token`` on this CDS, or ``None`` if the token isn't
    satisfied. Replays the same walk as
    :func:`_satisfying_span_for_token` but reports the E-value scalar
    instead of the span. Same ``passing``-flag filter as the span
    helper, so the two stay in lockstep.
    """
    try:
        token_hmms = parse_hmm_token(token)
    except ValueError:
        return None
    by_target: dict[str, list[dict]] = {}
    for h in hits or []:
        target = h.get("target") or h.get("hmm_name") or h.get("name")
        if not target:
            continue
        if not h.get("passing", True):
            continue
        by_target.setdefault(target, []).append(h)
    for name in token_hmms:
        if name not in by_target:
            return None
    if len(token_hmms) == 1:
        # Single-HMM pick: same key as _satisfying_span_for_token so the
        # worst-E and span helpers stay in lockstep on tied dom_evalue
        # (best-E first, longer ali_span as tiebreak).
        best = min(
            by_target[token_hmms[0]],
            key=lambda h: (
                float(h.get("dom_evalue", float("inf"))),
                -int(h.get("ali_span") or (
                    int(h.get("ali_to", 0)) - int(h.get("ali_from", 0)) + 1
                )),
            ),
        )
        return float(best.get("dom_evalue", float("inf")))
    chosen: list[dict] = []
    prev_from = -1
    prev_to = -1
    for name in token_hmms:
        candidates = [
            h for h in by_target[name]
            if int(h["ali_from"]) > prev_from
            and int(h["ali_to"]) > prev_to
            and (
                prev_to < 0
                or int(h["ali_from"]) >= prev_to - overlap_tolerance + 1
            )
        ]
        if not candidates:
            return None
        pick = min(
            candidates,
            key=lambda h: (
                int(h["ali_to"]),
                float(h.get("dom_evalue", float("inf"))),
            ),
        )
        chosen.append(pick)
        prev_from = int(pick["ali_from"])
        prev_to = int(pick["ali_to"])
    return float(max(h.get("dom_evalue", float("inf")) for h in chosen))


def identify_parent_cds(
    proteins: list[dict],
    spec: PolyproteinSpec,
    *,
    overlap_tolerance: int = 0,
) -> Optional[dict]:
    """The CDS that best fits the polyprotein declaration.

    Counts the number of *satisfied peptide tokens* each CDS carries —
    a single-HMM peptide counts when its one HMM hits; a multidomain
    peptide counts only when the whole architecture is present in N→C
    order — AND requires those satisfied tokens to themselves appear in
    the declared N→C order along the CDS. The order check prevents a
    chimeric / contaminated CDS (peptides hit but scrambled) from being
    elected as parent and then producing a confusing
    ``parent_protein_id=...`` audit row for a spec that was always
    going to fail. The CDS satisfying the most tokens
    (≥ ``spec.min_peptides_hit``) wins; ties are broken by translation
    length (a polyprotein is usually the longest CDS on the segment,
    but we don't assume — the satisfied-token count is the decisive
    signal).

    Returns ``None`` when no CDS clears the threshold or the satisfied
    tokens are out of declared order; the caller emits one
    ``no_parent_cds`` audit row per declared peptide.
    """
    best: tuple[Optional[dict], int, int] = (None, 0, 0)
    for prot in proteins or []:
        hits = prot.get("hmm_hits") or []
        # Collect (declared_index, span_from) for satisfied peptides.
        satisfied_positions: list[tuple[int, int]] = []
        for pep_idx, pep in enumerate(spec.peptides):
            # A peptide counts as satisfied when ANY of its alternative
            # architectures (`pep.hmms`) hits — OR semantics. Each token
            # internally is AND across its named domains.
            chosen = _best_satisfying_alternative(
                hits, pep.hmms, overlap_tolerance=overlap_tolerance,
            )
            if chosen is not None:
                satisfied_positions.append((pep_idx, chosen[0]))
        satisfied = len(satisfied_positions)
        if satisfied < spec.min_peptides_hit:
            continue
        # Order check: satisfied peptides (which are already in
        # declared-N→C order — `satisfied_positions` is appended in
        # `pep_idx` order) must also appear in increasing span_from on
        # the CDS. Strict `>=` so equal span starts (two peptides
        # collocated) also disqualify — they can't both be the parent
        # CDS's mature peptides.
        starts = [sf for _, sf in satisfied_positions]
        if any(starts[i] >= starts[i + 1] for i in range(len(starts) - 1)):
            continue
        length = int(prot.get("length") or len(prot.get("sequence") or ""))
        if (satisfied, length) > (best[1], best[2]):
            best = (prot, satisfied, length)
    return best[0]


def _peptide_span_midpoint(span: tuple[int, int]) -> float:
    """Centre of the span (1-based AA coords, fractional)."""
    return (span[0] + span[1]) / 2.0


def _find_motif_snap(
    protein_seq: str,
    bisect_point: int,
    motif: str,
    window_aa: int,
) -> Optional[int]:
    """Snap a cut to the last occurrence of ``motif`` near ``bisect_point``.

    The motif is the residue(s) just N-terminal of the cut (3CL: ``"LQ"``,
    picornavirus 3C: ``"Q"``). We search the window
    ``[bisect_point - window_aa, bisect_point + window_aa]`` on the
    protein and snap to the rightmost match (the cleavage site after
    that motif). Returns the 1-based AA position of the *first* residue
    of the downstream peptide (i.e. the position just after the motif
    end), or ``None`` if no motif occurrence falls in the window.

    All coordinates are 1-based inclusive amino acids — the caller is
    responsible for slicing ``protein_seq`` with the Python convention.
    """
    if not motif or not protein_seq:
        return None
    n = len(protein_seq)
    lo = max(1, bisect_point - window_aa)
    hi = min(n, bisect_point + window_aa)
    # Search the window for occurrences of motif; pick the rightmost so
    # the downstream peptide starts as far C-terminal as the motif allows
    # (mirrors "last protease cut before the next domain").
    found: Optional[int] = None
    motif_len = len(motif)
    # Convert to 0-based python indices for substring search.
    start0 = lo - 1
    end0 = hi - 1
    idx = protein_seq.find(motif, start0, end0 + motif_len)
    while idx != -1 and idx <= end0:
        # Cut comes AFTER the motif; downstream peptide starts at idx+motif_len+1 (1-based).
        found = idx + motif_len + 1
        idx = protein_seq.find(motif, idx + 1, end0 + motif_len)
    return found


def compute_cuts(
    protein_seq: str,
    spec: PolyproteinSpec,
    peptide_spans: list[tuple[PeptideSpec, Optional[tuple[int, int]]]],
) -> tuple[list[Optional[tuple[int, int, str]]], list[str]]:
    """Place inter-peptide cuts according to the spec's cut strategy.

    Each element of ``peptide_spans`` is ``(peptide, span)`` where
    ``span`` is ``(from_aa, to_aa)`` 1-based inclusive — covering the
    full footprint of the peptide token on the parent CDS (a single
    HMM's hit, or the union of all named domains for a multidomain
    token). ``None`` means the token couldn't be satisfied on this CDS.

    Returns ``(ranges, notes)``:

    * ``ranges[i]`` is ``(from_aa, to_aa, cut_method_actual)`` for the
      ``i``-th peptide, or ``None`` if the peptide was missing /
      otherwise unusable.
    * ``notes`` is a list of free-text warnings for the audit TSV.

    The strategy:

    1. Drop peptides without a span. Note them as ``missing``.
    2. Verify the surviving peptides' spans are in N→C order. Return
       early (empty ranges, ``out_of_order`` note) if not — the spec
       fails for this rep.
    3. ``boundary``: each surviving peptide spans its token's footprint
       verbatim.
    4. ``bisect`` / ``motif``: each surviving peptide is extended into
       its **immediate** spec-neighbours only. The bisect (or motif-
       snapped) cut is placed between consecutive-in-spec surviving
       peptides; when the immediate neighbour is missing the peptide's
       boundary on that side stays at the HMM hit's ``ali_from`` /
       ``ali_to`` (``cut_method_actual`` = ``hit-boundary``). N- and
       C-termini extensions to AA 1 / ``n_aa`` likewise apply only when
       the first / last declared peptide is itself surviving — a
       leading or trailing missing peptide leaves the gap unassigned
       rather than letting its neighbour absorb the missing peptide's
       territory. Pre-v0.37.0 the bisect blindly split the gap evenly,
       which on SARS-CoV-2 (NSP2 missing) inflated NSP1 from 180 aa to
       505 aa and bled 313 aa of NSP2 into NSP3.
    """
    n_aa = len(protein_seq)
    notes: list[str] = []
    n_decl = len(peptide_spans)
    ranges: list[Optional[tuple[int, int, str]]] = [None] * n_decl

    if n_aa == 0:
        notes.append("parent CDS has no translation; cannot slice")
        return ranges, notes

    surviving_idx = [i for i in range(n_decl) if peptide_spans[i][1] is not None]
    if not surviving_idx:
        notes.append("no peptide token was satisfied on the parent CDS")
        return ranges, notes

    starts = [peptide_spans[i][1][0] for i in surviving_idx]
    # Strict `>=` so two peptides with identical span starts (collocated
    # tokens — biologically nonsensical but possible if a config declares
    # the same HMM for two peptides) also fail. Pre-v0.36.1 the check
    # was `>`, which let collocated tokens through and then produced
    # garbage bisect math downstream.
    if any(starts[k] >= starts[k + 1] for k in range(len(starts) - 1)):
        notes.append(
            "peptide tokens were satisfied out of N-to-C order on the "
            "parent CDS — spec fails for this representative"
        )
        return ranges, notes

    if spec.cut_strategy == "boundary":
        for i in surviving_idx:
            span = peptide_spans[i][1]
            f = max(1, min(span[0], n_aa))
            t = max(1, min(span[1], n_aa))
            if t < f:
                continue
            ranges[i] = (f, t, "boundary")
        return ranges, notes

    # bisect / motif — place a cut between each pair of CONSECUTIVE-IN-SPEC
    # surviving peptides. cuts_between[i] = (cut_aa, method_label) is the
    # cut sitting between declared peptide i and declared peptide i+1, set
    # only when both peptides are surviving. When a peptide in between is
    # missing, no cut is placed across the gap — the flanking peptides
    # keep their HMM-hit boundary on the missing side rather than splitting
    # the missing peptide's territory between them (the v0.37.0 fix).
    use_motif = spec.cut_strategy == "motif"
    overlap_note_emitted = False
    cuts_between: dict[int, tuple[int, str]] = {}

    for slot in range(len(surviving_idx) - 1):
        idx_a = surviving_idx[slot]
        idx_b = surviving_idx[slot + 1]
        if idx_b != idx_a + 1:
            # Missing peptide(s) between idx_a and idx_b — gap stays
            # unassigned; idx_a's C-side and idx_b's N-side both stay
            # at their HMM-hit edges.
            continue
        a_span = peptide_spans[idx_a][1]
        b_span = peptide_spans[idx_b][1]
        pep_b = peptide_spans[idx_b][0]
        a_to = a_span[1]
        b_from = b_span[0]
        overlap_here = b_from <= a_to
        if overlap_here:
            # Adjacent token footprints overlap on the parent — bisect
            # on the midpoint of their centres so we don't pin a sharp
            # boundary at a contested residue. Every overlap pair gets
            # this treatment (pre-v0.36.1 only the first did); the audit
            # note is gated separately so the user isn't spammed.
            if not overlap_note_emitted:
                notes.append(
                    "footprints for adjacent peptides overlap on the "
                    "parent; cut(s) placed at the midpoint of the "
                    "footprint centres"
                )
                overlap_note_emitted = True
            mid_centre = (
                _peptide_span_midpoint(a_span) + _peptide_span_midpoint(b_span)
            ) / 2.0
            bisect_point = int(round(mid_centre))
        else:
            bisect_point = int(round((a_to + b_from) / 2.0))
        bisect_point = max(2, min(n_aa, bisect_point))

        snapped: Optional[int] = None
        if use_motif and pep_b.cleavage_motif:
            snapped = _find_motif_snap(
                protein_seq, bisect_point, pep_b.cleavage_motif,
                spec.motif_window_aa,
            )
        if snapped is not None and pep_b.cleavage_motif:
            cuts_between[idx_a] = (snapped, f"motif:{pep_b.cleavage_motif}")
        else:
            cuts_between[idx_a] = (bisect_point, "bisect")

    # Audit-note when any declared peptide is missing — the surviving
    # peptides on either side use their HMM-hit boundary rather than
    # extending into the unassigned territory. Includes the leading-/
    # trailing-missing case (first/last declared not surviving).
    has_unassigned_gap = (
        any(surviving_idx[k + 1] - surviving_idx[k] > 1
            for k in range(len(surviving_idx) - 1))
        or surviving_idx[0] != 0
        or surviving_idx[-1] != n_decl - 1
    )
    if has_unassigned_gap:
        notes.append(
            "one or more declared peptides were not located on the "
            "parent CDS; flanking peptides do not extend into the "
            "unassigned gap (boundary kept at the HMM hit edge)"
        )

    for idx in surviving_idx:
        span = peptide_spans[idx][1]

        # N-side: extend to AA 1 only when this is the first declared
        # peptide; bisect with idx-1 only when idx-1 is itself surviving
        # (i.e. cuts_between[idx-1] exists); otherwise stay at the HMM
        # hit's N-edge.
        if idx == 0:
            range_from = 1
            method = "n-term"
        elif (idx - 1) in cuts_between:
            cut, method = cuts_between[idx - 1]
            range_from = cut
        else:
            range_from = span[0]
            method = "hit-boundary"

        # C-side: extend to n_aa only when this is the last declared
        # peptide; bisect with idx+1 only when idx+1 is itself surviving;
        # otherwise stay at the HMM hit's C-edge.
        if idx == n_decl - 1:
            range_to = n_aa
        elif idx in cuts_between:
            range_to = cuts_between[idx][0] - 1
        else:
            range_to = span[1]

        range_from = max(1, min(n_aa, range_from))
        range_to = max(1, min(n_aa, range_to))
        if range_to < range_from:
            continue
        ranges[idx] = (range_from, range_to, method)

    return ranges, notes


def slice_polyprotein(
    proteins: list[dict],
    spec: PolyproteinSpec,
    *,
    overlap_tolerance: int = 0,
) -> tuple[Optional[dict], list[SlicedPeptide]]:
    """Top-level entry: identify parent CDS, compute cuts, build records.

    Returns ``(parent_cds, sliced_peptides)``:

    * ``parent_cds`` is the CDS dict (from :attr:`Sequence.proteins`)
      identified as the polyprotein, or ``None`` if no CDS satisfies
      ≥ ``spec.min_peptides_hit`` peptide tokens. In the ``None``
      case ``sliced_peptides`` is a single ``no_parent_cds`` audit row
      per peptide so the user can see *why* the spec produced nothing
      on this representative.
    * ``sliced_peptides`` is one :class:`SlicedPeptide` per declared
      peptide of the spec (in N→C order). Peptides whose token isn't
      satisfied produce a ``missing`` row with no FASTA-eligible
      sequence; ``ok`` rows carry the spliced AA string.

    ``overlap_tolerance`` is forwarded to the multidomain-token walk so
    Pfam-boundary fuzz at adjacent named domains (e.g. ``A--B`` where
    A's HMM model overruns into B by a few residues) doesn't reject a
    biologically valid CDS.
    """
    parent = identify_parent_cds(
        proteins, spec, overlap_tolerance=overlap_tolerance,
    )
    if parent is None:
        return None, [
            SlicedPeptide(
                peptide_name=pep.name,
                parent_protein_id=None,
                parent_accession=None,
                range_aa_from=0,
                range_aa_to=0,
                length_aa=0,
                sequence="",
                cut_method_actual="",
                status=_NO_PARENT,
                note=(
                    f"no CDS on this representative satisfies "
                    f"≥ {spec.min_peptides_hit} of the declared peptide "
                    f"tokens"
                ),
            )
            for pep in spec.peptides
        ]

    protein_seq = parent.get("sequence") or ""
    parent_pid = parent.get("protein_id")
    parent_acc = parent.get("parent_accession") or parent.get("accession")
    parent_hits = parent.get("hmm_hits") or []

    # Resolve each declared peptide's footprint on the parent CDS.
    # Restrict the search to the chosen parent — a peptide token hit
    # elsewhere on the rep is irrelevant to this polyprotein's layout.
    # For OR-peptides (multiple alternative architectures), pick the
    # best-E alternative and remember which one fired.
    peptide_spans: list[
        tuple[PeptideSpec, Optional[tuple[int, int]], str]
    ] = []
    for pep in spec.peptides:
        chosen = _best_satisfying_alternative(
            parent_hits, pep.hmms, overlap_tolerance=overlap_tolerance,
        )
        if chosen is None:
            peptide_spans.append((pep, None, ""))
        else:
            f, t, matched = chosen
            peptide_spans.append((pep, (f, t), matched))

    # compute_cuts only cares about (peptide, span) pairs; the matched
    # token rides alongside and gets written into SlicedPeptide below.
    ranges, notes = compute_cuts(
        protein_seq, spec,
        [(pep, span) for pep, span, _ in peptide_spans],
    )

    # Was the global "out of order" or "no satisfied tokens" fatal?
    # Emit one row per declared peptide carrying the failure reason.
    if all(r is None for r in ranges) and notes:
        reason = notes[0]
        global_status = _OUT_OF_ORDER if "out of N-to-C order" in reason else _MISSING
        result: list[SlicedPeptide] = []
        for pep, span, matched in peptide_spans:
            if span is None:
                status = _MISSING
                rng = (0, 0)
            else:
                status = global_status
                rng = span
            result.append(SlicedPeptide(
                peptide_name=pep.name,
                parent_protein_id=parent_pid,
                parent_accession=parent_acc,
                range_aa_from=rng[0],
                range_aa_to=rng[1],
                length_aa=max(0, rng[1] - rng[0] + 1) if span else 0,
                sequence="",
                cut_method_actual="",
                status=status,
                note=(
                    reason if status != _MISSING
                    else "peptide token was not satisfied on the parent CDS"
                ),
                matched_token=matched,
            ))
        return parent, result

    overlap_seen = any("overlap" in n.lower() for n in notes)

    out: list[SlicedPeptide] = []
    for (pep, span, matched), rng in zip(peptide_spans, ranges):
        if span is None:
            out.append(SlicedPeptide(
                peptide_name=pep.name,
                parent_protein_id=parent_pid,
                parent_accession=parent_acc,
                range_aa_from=0,
                range_aa_to=0,
                length_aa=0,
                sequence="",
                cut_method_actual="",
                status=_MISSING,
                note="peptide token was not satisfied on the parent CDS",
                matched_token="",
            ))
            continue
        if rng is None:
            # Defensive: a satisfied token whose range resolved to None
            # means the parent has zero residues somehow.
            out.append(SlicedPeptide(
                peptide_name=pep.name,
                parent_protein_id=parent_pid,
                parent_accession=parent_acc,
                range_aa_from=span[0],
                range_aa_to=span[1],
                length_aa=0,
                sequence="",
                cut_method_actual="",
                status=_MISSING,
                note="cut math produced no slice for this peptide",
                matched_token=matched,
            ))
            continue
        f, t, method = rng
        body = protein_seq[f - 1: t]
        status = _OK
        note = ""
        if overlap_seen and "overlap" in (notes[0] if notes else "").lower():
            status = _OVERLAP
            note = "adjacent peptide token footprints overlap on the parent CDS"
        out.append(SlicedPeptide(
            peptide_name=pep.name,
            parent_protein_id=parent_pid,
            parent_accession=parent_acc,
            range_aa_from=f,
            range_aa_to=t,
            length_aa=len(body),
            sequence=body,
            cut_method_actual=method,
            status=status,
            note=note,
            matched_token=matched,
        ))

    return parent, out
