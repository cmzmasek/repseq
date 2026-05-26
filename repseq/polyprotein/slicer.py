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

    Deliberately ignores the ``passing`` flag from ``_run_hmm_qc``: by
    the time the slicer asks for a token's footprint, the parent CDS has
    already been chosen, and we want every domain placement the hmmscan
    found, not just the QC-clearing subset. The walk logic mirrors
    :func:`repseq.hmm.runner.cds_satisfies_token` so the slicer enforces
    the same architecture rules as the cluster_protein HMM gate.
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
    # cds_satisfies_token but without the passing-flag gate.
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
        pick = min(candidates, key=lambda h: int(h["ali_to"]))
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
    instead of the span.
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
        by_target.setdefault(target, []).append(h)
    for name in token_hmms:
        if name not in by_target:
            return None
    if len(token_hmms) == 1:
        best = min(
            by_target[token_hmms[0]],
            key=lambda h: float(h.get("dom_evalue", float("inf"))),
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
        pick = min(candidates, key=lambda h: int(h["ali_to"]))
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
    order. The CDS satisfying the most tokens
    (≥ ``spec.min_peptides_hit``) wins; ties are broken by translation
    length (a polyprotein is usually the longest CDS on the segment,
    but we don't assume — the satisfied-token count is the decisive
    signal).

    Returns ``None`` when no CDS clears the threshold; the caller emits
    one ``no_parent_cds`` audit row per declared peptide.
    """
    best: tuple[Optional[dict], int, int] = (None, 0, 0)
    for prot in proteins or []:
        hits = prot.get("hmm_hits") or []
        satisfied = 0
        for pep in spec.peptides:
            # A peptide counts as satisfied when ANY of its alternative
            # architectures (`pep.hmms`) hits — OR semantics. Each token
            # internally is AND across its named domains.
            if _best_satisfying_alternative(
                hits, pep.hmms, overlap_tolerance=overlap_tolerance,
            ) is not None:
                satisfied += 1
        if satisfied < spec.min_peptides_hit:
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
    4. ``bisect`` / ``motif``: each surviving peptide spans from the
       previous boundary to the bisect (or motif-snapped) point with
       the next surviving peptide. The first peptide starts at AA 1;
       the last ends at the protein C-term.
    """
    n_aa = len(protein_seq)
    notes: list[str] = []
    ranges: list[Optional[tuple[int, int, str]]] = [None] * len(peptide_spans)

    if n_aa == 0:
        notes.append("parent CDS has no translation; cannot slice")
        return ranges, notes

    surviving_idx = [i for i, (_, s) in enumerate(peptide_spans) if s is not None]
    if not surviving_idx:
        notes.append("no peptide token was satisfied on the parent CDS")
        return ranges, notes

    starts = [peptide_spans[i][1][0] for i in surviving_idx]
    if any(starts[k] > starts[k + 1] for k in range(len(starts) - 1)):
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

    # bisect / motif — chained cut placement across surviving peptides.
    use_motif = spec.cut_strategy == "motif"
    cuts: list[int] = []
    overlap_flagged = False
    for k in range(len(surviving_idx) - 1):
        a_span = peptide_spans[surviving_idx[k]][1]
        b_span = peptide_spans[surviving_idx[k + 1]][1]
        a_to = a_span[1]
        b_from = b_span[0]
        if b_from <= a_to and not overlap_flagged:
            # Adjacent token footprints overlap on the parent — bisect
            # on the midpoint of their centres so we don't pin a sharp
            # boundary at a contested residue.
            notes.append(
                f"footprints for adjacent peptides overlap on the parent "
                f"(a_to={a_to}, b_from={b_from}); cut placed at the "
                f"midpoint of the footprint centres"
            )
            overlap_flagged = True
            mid_centre = (
                _peptide_span_midpoint(a_span) + _peptide_span_midpoint(b_span)
            ) / 2.0
            bisect_point = int(round(mid_centre))
        else:
            bisect_point = int(round((a_to + b_from) / 2.0))
        bisect_point = max(2, min(n_aa, bisect_point))

        snapped: Optional[int] = None
        if use_motif:
            motif = peptide_spans[surviving_idx[k + 1]][0].cleavage_motif
            if motif:
                snapped = _find_motif_snap(
                    protein_seq, bisect_point, motif, spec.motif_window_aa,
                )
        cuts.append(snapped if snapped is not None else bisect_point)

    starts_iter = [1] + cuts
    ends_iter = [c - 1 for c in cuts] + [n_aa]

    methods = ["n-term"]
    for k in range(len(surviving_idx) - 1):
        pep = peptide_spans[surviving_idx[k + 1]][0]
        a_to = peptide_spans[surviving_idx[k]][1][1]
        b_from = peptide_spans[surviving_idx[k + 1]][1][0]
        bisect_point = int(round((a_to + b_from) / 2.0))
        snapped: Optional[int] = None
        if use_motif and pep.cleavage_motif:
            snapped = _find_motif_snap(
                protein_seq, bisect_point, pep.cleavage_motif,
                spec.motif_window_aa,
            )
        if snapped is not None and pep.cleavage_motif:
            methods.append(f"motif:{pep.cleavage_motif}")
        else:
            methods.append("bisect")

    for slot, idx in enumerate(surviving_idx):
        f, t = starts_iter[slot], ends_iter[slot]
        if t < f:
            continue
        ranges[idx] = (f, t, methods[slot])

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
