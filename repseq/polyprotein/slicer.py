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


def _best_hit_for_hmm(
    proteins: list[dict], hmm_name: str
) -> tuple[Optional[dict], Optional[dict]]:
    """Return ``(parent_cds, hit)`` for the best hit to ``hmm_name``.

    Best = lowest ``dom_evalue`` (tie-broken by longer ``ali_span``).
    Ignores the ``passing`` flag — at the slicing stage we WANT the
    locations of every hit, even ones that wouldn't have cleared
    QC. The QC-passing check has already happened upstream; here we
    are mapping out the polyprotein's domain layout.

    Returns ``(None, None)`` when no protein has any hit to this HMM.
    """
    best: tuple[Optional[dict], Optional[dict], float, int] = (None, None, float("inf"), 0)
    for prot in proteins or []:
        for hit in prot.get("hmm_hits") or []:
            name = hit.get("target") or hit.get("hmm_name") or hit.get("name")
            if not name or name != hmm_name:
                continue
            ev = float(hit.get("dom_evalue", hit.get("evalue", 1.0)))
            span = int(hit.get("ali_span") or (
                int(hit.get("ali_to", 0)) - int(hit.get("ali_from", 0)) + 1
            ))
            if (ev, -span) < (best[2], -best[3]):
                best = (prot, hit, ev, span)
    return best[0], best[1]


def identify_parent_cds(
    proteins: list[dict], spec: PolyproteinSpec,
) -> Optional[dict]:
    """The CDS that best fits the polyprotein declaration.

    Counts the number of *distinct* peptide HMMs each CDS carries hits
    for; the CDS with the most distinct hits (≥ ``spec.min_peptides_hit``)
    wins. Ties are broken by translation length (a polyprotein is usually
    the longest CDS on the segment, but we don't assume — the HMM-hit
    count is the decisive signal).

    Returns ``None`` when no CDS clears the threshold, which becomes a
    soft-fail at the caller (no peptides emitted for this rep × spec
    combination; a single ``no_parent_cds`` audit row records it).
    """
    target_hmms = {p.hmm for p in spec.peptides}
    best: tuple[Optional[dict], int, int] = (None, 0, 0)
    for prot in proteins or []:
        hits = prot.get("hmm_hits") or []
        present = {
            (h.get("target") or h.get("hmm_name") or h.get("name"))
            for h in hits
        }
        n_distinct = len(present & target_hmms)
        if n_distinct < spec.min_peptides_hit:
            continue
        length = int(prot.get("length") or len(prot.get("sequence") or ""))
        if (n_distinct, length) > (best[1], best[2]):
            best = (prot, n_distinct, length)
    return best[0]


def _peptide_hit_midpoint(hit: dict) -> float:
    """Centre of the hit (1-based AA coords, fractional)."""
    return (int(hit["ali_from"]) + int(hit["ali_to"])) / 2.0


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
    peptide_hits: list[tuple[PeptideSpec, Optional[dict]]],
) -> tuple[list[Optional[tuple[int, int, str]]], list[str]]:
    """Place inter-peptide cuts according to the spec's cut strategy.

    Returns ``(ranges, notes)``:

    * ``ranges[i]`` is ``(from_aa, to_aa, cut_method_actual)`` for the
      ``i``-th peptide (1-based inclusive AA coords), or ``None`` if the
      peptide's HMM didn't hit on the parent CDS (peptide skipped; the
      hole is closed by the neighbouring peptides extending across it).
    * ``notes`` is a list of free-text warnings for the audit TSV
      (e.g. ``"NSP3 hit overlaps NSP5 hit; bisecting from averaged starts"``).

    The strategy:

    1. Drop peptides whose HMM didn't hit. Note them as ``missing``.
    2. Verify the surviving peptides' hits are in N→C order. Return
       early (empty ranges, ``out_of_order`` note) if not — the spec
       fails for this rep.
    3. ``boundary``: each surviving peptide spans its hit's
       ``ali_from..ali_to`` verbatim.
    4. ``bisect`` / ``motif``: each surviving peptide spans from the
       previous boundary to the bisect (or motif-snapped) point with the
       next surviving peptide. The first peptide starts at AA 1; the
       last ends at the protein C-term.
    """
    n_aa = len(protein_seq)
    notes: list[str] = []
    ranges: list[Optional[tuple[int, int, str]]] = [None] * len(peptide_hits)

    if n_aa == 0:
        notes.append("parent CDS has no translation; cannot slice")
        return ranges, notes

    # Index of each surviving (peptide_hits[i] has a real hit).
    surviving_idx = [i for i, (_, h) in enumerate(peptide_hits) if h is not None]
    if not surviving_idx:
        notes.append("no peptide HMM hit the parent CDS")
        return ranges, notes

    # Order check: hits must be N→C on the parent.
    starts = [int(peptide_hits[i][1]["ali_from"]) for i in surviving_idx]
    if any(starts[k] > starts[k + 1] for k in range(len(starts) - 1)):
        notes.append(
            "peptide HMMs hit out of N-to-C order on the parent CDS — "
            "spec fails for this representative"
        )
        return ranges, notes

    if spec.cut_strategy == "boundary":
        for i in surviving_idx:
            hit = peptide_hits[i][1]
            f, t = int(hit["ali_from"]), int(hit["ali_to"])
            f = max(1, min(f, n_aa))
            t = max(1, min(t, n_aa))
            if t < f:
                continue
            ranges[i] = (f, t, "boundary")
        return ranges, notes

    # bisect / motif — chained cut placement across surviving peptides.
    # Compute the bisect point between each adjacent surviving pair.
    use_motif = spec.cut_strategy == "motif"
    cuts: list[int] = []  # 1-based start position of each surviving peptide after the first
    overlap_flagged = False
    for k in range(len(surviving_idx) - 1):
        a_hit = peptide_hits[surviving_idx[k]][1]
        b_hit = peptide_hits[surviving_idx[k + 1]][1]
        a_to = int(a_hit["ali_to"])
        b_from = int(b_hit["ali_from"])
        if b_from <= a_to and not overlap_flagged:
            # Hits overlap; bisect on the midpoints of their *centres*
            # so we don't pin a sharp boundary at a contested residue.
            notes.append(
                f"hits for adjacent peptides overlap on the parent "
                f"(a_to={a_to}, b_from={b_from}); cut placed at the "
                f"midpoint of the hit centres"
            )
            overlap_flagged = True
            mid_centre = (
                _peptide_hit_midpoint(a_hit) + _peptide_hit_midpoint(b_hit)
            ) / 2.0
            bisect_point = int(round(mid_centre))
        else:
            bisect_point = int(round((a_to + b_from) / 2.0))
        bisect_point = max(2, min(n_aa, bisect_point))

        snapped: Optional[int] = None
        if use_motif:
            motif = peptide_hits[surviving_idx[k + 1]][0].cleavage_motif
            if motif:
                snapped = _find_motif_snap(
                    protein_seq, bisect_point, motif, spec.motif_window_aa,
                )
        if snapped is not None:
            cuts.append(snapped)
        else:
            cuts.append(bisect_point)

    # Assemble ranges: each surviving peptide owns AA positions
    # [start, end] where start = previous-cut (or 1 for the first) and
    # end = this-cut - 1 (or n_aa for the last). The method tag records
    # what produced the START boundary of the peptide (the more
    # informative one — the C-term cut becomes the next peptide's N-term).
    starts_iter = [1] + cuts
    ends_iter = [c - 1 for c in cuts] + [n_aa]

    methods = ["n-term"]  # first surviving peptide's start is the protein N-term
    for k in range(len(surviving_idx) - 1):
        pep = peptide_hits[surviving_idx[k + 1]][0]
        bisect_point = int(round(
            (int(peptide_hits[surviving_idx[k]][1]["ali_to"])
             + int(peptide_hits[surviving_idx[k + 1]][1]["ali_from"])) / 2.0
        ))
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
) -> tuple[Optional[dict], list[SlicedPeptide]]:
    """Top-level entry: identify parent CDS, compute cuts, build records.

    Returns ``(parent_cds, sliced_peptides)``:

    * ``parent_cds`` is the CDS dict (from :attr:`Sequence.proteins`)
      identified as the polyprotein, or ``None`` if no CDS cleared
      ``min_peptides_hit`` distinct peptide-HMM hits. In the ``None``
      case ``sliced_peptides`` is a single ``no_parent_cds`` audit row
      so the user can see *why* the spec produced nothing on this
      representative.
    * ``sliced_peptides`` is one :class:`SlicedPeptide` per declared
      peptide of the spec (in N→C order). Peptides whose HMM didn't
      hit produce a ``missing`` row with no FASTA-eligible sequence;
      ``ok`` rows carry the spliced AA string.
    """
    parent = identify_parent_cds(proteins, spec)
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
                    f"no CDS on this representative carries hits from "
                    f"≥ {spec.min_peptides_hit} of the declared peptide HMMs"
                ),
            )
            for pep in spec.peptides
        ]

    protein_seq = parent.get("sequence") or ""
    parent_pid = parent.get("protein_id")
    parent_acc = parent.get("parent_accession") or parent.get("accession")

    # Resolve each declared peptide's best hit on the parent CDS.
    peptide_hits: list[tuple[PeptideSpec, Optional[dict]]] = []
    for pep in spec.peptides:
        # Restrict the hit search to the chosen parent CDS — a peptide
        # HMM that hits another CDS on the rep is not relevant to this
        # polyprotein's layout.
        _, hit = _best_hit_for_hmm([parent], pep.hmm)
        peptide_hits.append((pep, hit))

    ranges, notes = compute_cuts(protein_seq, spec, peptide_hits)

    # Was the global "out of order" or "no hits" fatal? If so emit one
    # status row per declared peptide carrying the failure reason.
    if all(r is None for r in ranges) and notes:
        # Distinguish out-of-order from "no hits at all": a single
        # combined audit reason per peptide.
        reason = notes[0]
        global_status = _OUT_OF_ORDER if "out of N-to-C order" in reason else _MISSING
        result: list[SlicedPeptide] = []
        for slot, (pep, hit) in enumerate(peptide_hits):
            if hit is None:
                status = _MISSING
                rng = (0, 0)
            else:
                status = global_status
                rng = (int(hit["ali_from"]), int(hit["ali_to"]))
            result.append(SlicedPeptide(
                peptide_name=pep.name,
                parent_protein_id=parent_pid,
                parent_accession=parent_acc,
                range_aa_from=rng[0],
                range_aa_to=rng[1],
                length_aa=max(0, rng[1] - rng[0] + 1) if hit else 0,
                sequence="",
                cut_method_actual="",
                status=status,
                note=reason if status != _MISSING else "peptide HMM did not hit on the parent CDS",
            ))
        return parent, result

    overlap_seen = any("overlap" in n.lower() for n in notes)

    out: list[SlicedPeptide] = []
    for (pep, hit), rng in zip(peptide_hits, ranges):
        if hit is None:
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
                note="peptide HMM did not hit on the parent CDS",
            ))
            continue
        if rng is None:
            # Should be unreachable now that we drop None-hit slots
            # above, but defensive: a non-None hit that still resolved
            # to no range means we ran out of room (e.g. parent has
            # zero residues somehow).
            out.append(SlicedPeptide(
                peptide_name=pep.name,
                parent_protein_id=parent_pid,
                parent_accession=parent_acc,
                range_aa_from=int(hit["ali_from"]),
                range_aa_to=int(hit["ali_to"]),
                length_aa=0,
                sequence="",
                cut_method_actual="",
                status=_MISSING,
                note="cut math produced no slice for this peptide",
            ))
            continue
        f, t, method = rng
        body = protein_seq[f - 1: t]
        status = _OK
        note = ""
        if overlap_seen and "overlap" in (notes[0] if notes else "").lower():
            # Tag every peptide of this slicing as overlap-affected, but
            # still keep the sequences (the cut math fell back to
            # midpoint-of-centres for the overlap — usable, just flagged).
            status = _OVERLAP
            note = "adjacent peptide HMM hits overlap on the parent CDS"
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
        ))

    return parent, out
