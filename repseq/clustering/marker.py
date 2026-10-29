"""Marker-protein selection for amino-acid clustering.

When ``clustering.alphabet_for_clustering`` is ``protein``, every
sequence (non-segmented) or segment (segmented) must contribute one
amino-acid string to the clustering input. This module picks that
string from the sequence's CDS list (``seq.proteins``).

Marker selection follows the configured marker specs. Each spec is
either an alias string (legacy form) or a dict ``{name, aliases?,
hmms?}`` where ``hmms`` is a list of **token** strings:

    - ``"Name"``         — single-HMM token
    - ``"A--B--C"``      — multidomain token, HMMs listed in N-to-C
                           order (A most N-terminal, C most C-terminal)

The semantic:

    - A CDS satisfies a single-HMM token when that HMM has a passing
      hit on it.
    - A CDS satisfies a multidomain token when every named HMM has a
      passing hit AND the hits appear in N-to-C order on the protein
      (strict non-overlap). Extra domains on the same CDS are fine.
    - Within one spec the tokens in ``hmms:`` are **alternative
      architectures (OR)**: a CDS satisfies the spec when it satisfies
      *any one* of the tokens. This lets a single marker cover divergent
      architectures of the same protein — e.g. coronavirus Spike as
      ``["CoV_S1--CoV_S2", "bCoV_S1_N--bCoV_S1_RBD--CoV_S2"]`` so alpha-
      and beta-CoV Spikes are both accepted.
    - Per-spec: ``aliases`` and ``hmms`` are both consulted when both
      are set, BUT once an HMM tier is active and any spec defines
      ``hmms``, the HMM tier is AUTHORITATIVE for that spec — if no
      CDS satisfies any token in the spec, aliases are NOT consulted
      as a fallback and the spec is recorded as a hard failure.
    - Marker selection: from the set of CDSes that satisfy at least
      one token in the spec, pick the **longest**. (Alias-only specs
      keep the legacy "first alias substring match, longest CDS"
      behaviour.)

Per-spec QC (drop vs keep) is enforced upstream by ``cli._run_hmm_qc``;
this module only picks the marker among already-surviving CDSes. The
distinction matters because QC checks "every token satisfied by some
CDS in the segment", while marker selection picks one CDS — different
uses of the same hit cache.

Length bounds (``length_aa: [min, max]``, v0.76.0): a dict spec may carry
an inclusive amino-acid range, and only a CDS inside it can be picked for
that spec. The same bound is applied by every picker in the pipeline —
this one (clustering), the marker-length QC stage, and the per-protein
FASTA / tree / supermatrix pickers in :mod:`repseq.phylo` — so a tree never
shows a CDS that clustering did not use. A bounded spec is strict: when the
protein is present but outside the range the result is a
``length_out_of_bounds`` failure, and when a bounded alias-only spec matches
no product at all it is a ``no_alias_match`` failure (missing, not partial) —
never a silent fall-through to the longest CDS.

Marker completeness (``qc.marker_completeness.require_met_start``,
v0.77.0): a global switch, not a per-spec one. When on, only a CDS whose
translation begins with methionine can be picked for ANY spec — a complete
viral CDS always starts with Met, so a translation that does not is a
5'-partial record, which no length range can detect (a short N-terminal
truncation stays inside the range). Composed with ``length_aa`` as a second
candidate filter in every picker, with length taking precedence in the
reported failure: a candidate outside the range is ``length_out_of_bounds``,
one inside the range but lacking the Met start is ``no_met_start``.

When no specs are configured at all, the longest CDS with a translation
is returned (legacy fallback). When specs ARE configured but none
yield a satisfying CDS, the function returns ``(None, MarkerFailure)``
so the caller can drop the sequence/isolate with a structured reason.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Optional

from ..hmm.runner import cds_satisfies_token, coverage_of, parse_hmm_token
from ..models import QCReport, Sequence


# ---------------------------------------------------------------------------
# length_aa bounds — shared by every marker / extra-protein picker
# ---------------------------------------------------------------------------
#
# ``length_aa: [min, max]`` (amino acids, inclusive; either end may be null)
# restricts which CDS a spec may pick. The helpers live here, in the
# clustering picker's module, because the SAME bound must hold wherever a
# marker is picked (clustering, the marker-length QC stage, per-protein
# FASTAs and trees, the partitioned supermatrix); ``repseq.phylo.per_protein``
# re-exports them for the output-side pickers.

LengthBounds = tuple[Optional[int], Optional[int]]


def parse_length_bounds(raw: Any) -> Optional[LengthBounds]:
    """``[min, max]`` → ``(min, max)``; ``None`` when absent or malformed.

    Validation rejects malformed values at load; this stays defensive so a
    hand-built cfg (tests, replay) can't crash a picker.
    """
    if not isinstance(raw, (list, tuple)) or len(raw) != 2:
        return None
    lo, hi = raw
    lo = lo if isinstance(lo, int) and not isinstance(lo, bool) else None
    hi = hi if isinstance(hi, int) and not isinstance(hi, bool) else None
    if lo is None and hi is None:
        return None
    return (lo, hi)


def format_length_bounds(bounds: Optional[LengthBounds]) -> str:
    """``(250, 400)`` → ``"[250, 400]"``; an open end renders as ``null``."""
    if not bounds:
        return "[null, null]"
    lo, hi = bounds
    return f"[{'null' if lo is None else lo}, {'null' if hi is None else hi}]"


def within_length_bounds(prot: dict, bounds: Optional[LengthBounds]) -> bool:
    """True when the CDS translation length lies inside ``bounds``."""
    if not bounds:
        return True
    n = len(prot.get("sequence") or "")
    lo, hi = bounds
    return (lo is None or n >= lo) and (hi is None or n <= hi)


def starts_with_met(prot: dict) -> bool:
    """True when the CDS translation begins with methionine.

    The decisive 5'-completeness test: a complete viral CDS translation
    always starts with Met, so one that does not is a 5'-partial record.
    A CDS with no translation is not judged here (callers filter those out
    first) and counts as passing, so this never becomes a second
    "no translation" drop.
    """
    seq = prot.get("sequence") or ""
    return not seq or seq[0].upper() == "M"


def marker_candidate_ok(
    prot: dict,
    bounds: Optional[LengthBounds],
    *,
    require_met_start: bool = False,
) -> bool:
    """Whether ``prot`` may be picked for a spec with these restrictions.

    Composes the two candidate filters for pickers that only need a
    yes/no — ``phylo.per_protein.pick_marker_cds`` and everything routed
    through it (per-protein FASTAs, 2F/2K trees, the supermatrix, the
    coverage report). ``select_marker_protein`` and ``pick_extra_cds``
    deliberately apply the filters in SEQUENCE instead, because they must
    report WHICH restriction withheld a CDS and length takes precedence
    over completeness; a single boolean cannot carry that.
    """
    if not within_length_bounds(prot, bounds):
        return False
    return not require_met_start or starts_with_met(prot)


def require_met_start_from_cfg(cfg: Any) -> bool:
    """``qc.marker_completeness.require_met_start`` (default False).

    One global switch rather than a per-spec flag (designer's call,
    2026-10-29): the requirement is identical on every marker, and a
    per-spec boolean would let a config gate seven of eight segments and
    look entirely normal.
    """
    qc = (cfg or {}).get("qc", {}) or {}
    block = qc.get("marker_completeness", {}) or {}
    return bool(block.get("require_met_start", False))


MARKER_QUALITY_RUNTIME = "_marker_quality_runtime"


def require_met_start_active(cfg: Any) -> bool:
    """``require_met_start``, but False once the QC stage has skipped it.

    The stage skips under ``--no-resolve`` (no GenBank translations to
    judge) and records that on ``cfg[MARKER_QUALITY_RUNTIME]``. Every
    PICKER asks this, so the FASTAs, trees, supermatrix and coverage
    tables cannot enforce a requirement the Methods summary reports as
    not applied. ``require_met_start_from_cfg`` stays the CONFIGURED
    value, which is what that summary prose needs.
    """
    if not require_met_start_from_cfg(cfg):
        return False
    rt = (cfg or {}).get(MARKER_QUALITY_RUNTIME) or {}
    return not rt.get("require_met_start_skipped", False)


def describe_length_miss(
    length: Optional[int], bounds: Optional[LengthBounds],
) -> str:
    """Reason fragment for a length failure, e.g. ``481aa_below_485``.

    ``length`` is the CDS the picker would have taken without the bound;
    ``None`` means no candidate matched at all (``no_candidate``).
    """
    if length is None:
        return "no_candidate"
    lo, hi = bounds or (None, None)
    if lo is not None and length < lo:
        return f"{length}aa_below_{lo}"
    if hi is not None and length > hi:
        return f"{length}aa_above_{hi}"
    return f"{length}aa"


@dataclass
class MarkerFailure:
    """Structured reason why ``select_marker_protein`` returned no marker.

    ``reason`` is one of:
        - ``no_proteins``     — proteins list is None or empty
        - ``no_translated``   — every CDS lacks a /translation
        - ``hmm_failed``      — an HMM-gated spec had no CDS satisfy
                                any token; ``failed_tokens`` lists the
                                token strings that no CDS satisfied,
                                and ``best_evalue`` is the best E-value
                                seen across passing hits for HMMs that
                                appear in those failed tokens (best-of
                                the partial evidence, for triage).
        - ``no_alias_match``  — a bounded alias-only spec matched no
                                ``/product`` at all (the protein was not
                                found; unbounded alias-only specs fall back
                                to the longest CDS instead)
        - ``length_out_of_bounds`` — a spec with ``length_aa`` found its
                                protein only outside the range. ``length_aa``
                                is the length of the CDS that would have been
                                picked without the bound; ``bounds`` is the
                                spec's range.
        - ``no_met_start``    — the spec's protein was found INSIDE any
                                ``length_aa`` range but its translation does
                                not begin with Met, i.e. the CDS is
                                5'-partial (only when
                                ``qc.marker_completeness.require_met_start``
                                is on). ``length_aa`` carries the length of
                                the CDS that would have been picked without
                                the requirement. Reported only when the
                                length bound itself is satisfied, so length
                                failures never read as truncations.
    """

    reason: str
    marker_name: Optional[str] = None
    failed_tokens: list[str] = field(default_factory=list)
    best_evalue: Optional[float] = None
    length_aa: Optional[int] = None
    bounds: Optional[LengthBounds] = None


def _normalize_marker_specs(specs: Any) -> list[dict]:
    """Normalize a mixed list of strings / dicts into uniform spec dicts.

    Each output dict has keys ``name``, ``aliases``, ``hmms`` (lists are
    always present even if empty) and ``bounds`` (the parsed ``length_aa``,
    or ``None``). A bare string ``"foo"`` normalises to
    ``{name: "foo", aliases: ["foo"], hmms: [], bounds: None}``.
    """
    out: list[dict] = []
    if not specs:
        return out
    for entry in specs:
        if isinstance(entry, str):
            s = entry.strip()
            if s:
                out.append({"name": s, "aliases": [s], "hmms": [], "bounds": None})
        elif isinstance(entry, dict):
            name = entry.get("name") or ""
            aliases = list(entry.get("aliases") or [])
            hmms = list(entry.get("hmms") or [])
            if name or aliases or hmms:
                out.append({
                    "name": name, "aliases": aliases, "hmms": hmms,
                    "bounds": (
                        entry["bounds"] if "bounds" in entry
                        else parse_length_bounds(entry.get("length_aa"))
                    ),
                })
    return out


def _passing_hits(p: dict) -> list[dict]:
    """Return the CDS's passing HMM hits (pre-annotated by _run_hmm_qc)."""
    return [h for h in (p.get("hmm_hits") or []) if h.get("passing")]


def _cds_satisfies_any_token(
    p: dict, tokens: list[str], overlap_tolerance: int = 0,
) -> Optional[str]:
    """Return the token string this CDS satisfies (any token), or None.

    Used by marker selection to pick the longest CDS satisfying at least
    one token in the spec. ``overlap_tolerance`` is forwarded to
    :func:`cds_satisfies_token`.
    """
    hits = p.get("hmm_hits") or []
    for token in tokens:
        try:
            parsed = parse_hmm_token(token)
        except ValueError:
            continue
        if cds_satisfies_token(
            hits, parsed, overlap_tolerance=overlap_tolerance,
        ) is not None:
            return token
    return None


def select_marker_protein(
    proteins: Optional[list[dict]],
    marker_specs: Any = None,
    *,
    hmm_active: bool = False,
    ga_cutoffs: Optional[dict[str, Optional[float]]] = None,
    hmm_cfg: Optional[dict[str, Any]] = None,
    ignore_restrictions: bool = False,
    require_met_start: bool = False,
) -> tuple[Optional[dict], Optional[MarkerFailure]]:
    """Pick the marker CDS for one sequence/segment.

    Returns ``(marker_dict, None)`` on success or
    ``(None, MarkerFailure)`` on failure. See module docstring for the
    gating semantics. ``marker_specs`` accepts either the legacy alias
    string list (back-compat with pre-v0.13 callers) or the new mixed
    list of strings / dicts.

    ``length_aa`` bounds on a spec restrict its candidates; a spec whose
    protein exists only outside the range yields a ``length_out_of_bounds``
    failure (HMM failures take precedence, since the marker-length QC stage
    runs after the HMM gate). ``require_met_start`` adds the second
    candidate filter: a CDS inside the range whose translation does not
    begin with Met is 5'-partial, and a spec left with only such candidates
    yields a ``no_met_start`` failure.

    When several specs fail, the failure KIND decides what is returned
    (hmm → length → met → missing, the order the QC stages attribute in),
    and within one kind the FIRST declared spec that failed that way names
    the marker. One rule for every branch: the reported ``marker_name``
    reaches ``_qc_removed.tsv`` and the per-marker breakdown, so it must not
    depend on which code path produced the failure.

    ``ignore_restrictions=True`` waives BOTH restrictions — the length
    bound and the Met requirement — which is exactly what the
    ``marker_quality`` protect stage buys (``overrides.protect_qc``), so
    clustering uses a protected isolate's real protein exactly as the FASTA
    / tree pickers do. (Named for what it does rather than for the length
    half alone: a caller that waived one and kept the other would silently
    half-protect a reference strain.)

    Note: the HMM hits have already been scanned and annotated with
    ``passing`` flags by ``cli._run_hmm_qc`` (v0.14.0+). The
    ``hmm_active`` / ``ga_cutoffs`` / ``hmm_cfg`` kwargs are kept on
    the signature for back-compat with unit tests that construct hits
    by hand without the pre-annotation; in normal pipeline use they
    are not consulted.
    """
    if not proteins:
        return None, MarkerFailure("no_proteins")
    with_seq = [p for p in proteins if p.get("sequence")]
    if not with_seq:
        return None, MarkerFailure("no_translated")

    # ONE place decides what protection waives: the `marker_quality` protect
    # stage covers both restrictions, so a caller passes a single flag and
    # cannot accidentally waive one and keep the other.
    met_required = require_met_start and not ignore_restrictions

    # Test back-compat: if hits lack the pre-computed `passing` flag,
    # recompute it on the fly from the supplied cutoffs.
    if hmm_active and hmm_cfg is not None:
        _ensure_passing_annotation(with_seq, ga_cutoffs, hmm_cfg)

    specs = _normalize_marker_specs(marker_specs)

    if not specs:
        # No markers configured → legacy longest-CDS fallback, still
        # restricted to a complete CDS when the run asks for one (the
        # picked protein is what clustering sees either way).
        pool = [p for p in with_seq if starts_with_met(p)] if met_required else with_seq
        if not pool:
            return None, MarkerFailure(
                reason="no_met_start",
                length_aa=len(max(with_seq, key=lambda p: len(p["sequence"]))["sequence"]),
            )
        return max(pool, key=lambda p: len(p["sequence"])), None

    def _longest(cands: list[dict]) -> dict:
        return max(cands, key=lambda p: len(p["sequence"]))

    last_hmm_failure: Optional[MarkerFailure] = None
    last_length_failure: Optional[MarkerFailure] = None
    last_met_failure: Optional[MarkerFailure] = None
    last_missing_failure: Optional[MarkerFailure] = None
    # Bounds of the first spec that cannot identify its protein here (no
    # aliases, and its hmms unusable because the HMM tier is off). The legacy
    # longest-CDS fallback below must still respect that bound.
    fallback_bounds: Optional[LengthBounds] = None
    fallback_name: Optional[str] = None
    overlap_tol = int((hmm_cfg or {}).get("multidomain_overlap_tolerance", 30))

    def _met_ok(p: dict) -> bool:
        return not met_required or starts_with_met(p)

    for spec in specs:
        bounds = None if ignore_restrictions else spec.get("bounds")
        hmms = spec["hmms"]
        if hmms and hmm_active:
            # HMM tier is authoritative for this spec.
            satisfying: list[dict] = []
            for p in with_seq:
                if _cds_satisfies_any_token(
                    p, hmms, overlap_tolerance=overlap_tol,
                ) is not None:
                    satisfying.append(p)
            in_range = [p for p in satisfying if within_length_bounds(p, bounds)]
            pickable = [p for p in in_range if _met_ok(p)]
            if pickable:
                return _longest(pickable), None
            if in_range:
                # Length is fine; every candidate is 5'-partial. Reported
                # separately so a truncation never reads as a length miss.
                if last_met_failure is None:
                    last_met_failure = MarkerFailure(
                        reason="no_met_start",
                        marker_name=spec["name"] or ",".join(hmms),
                        length_aa=len(_longest(in_range)["sequence"]),
                        bounds=bounds,
                    )
                continue
            if satisfying:
                # HMM-verified protein present, but only outside length_aa.
                if last_length_failure is None:
                    last_length_failure = MarkerFailure(
                        reason="length_out_of_bounds",
                        marker_name=spec["name"] or ",".join(hmms),
                        length_aa=len(_longest(satisfying)["sequence"]),
                        bounds=bounds,
                    )
                continue
            # Record the failure and move on — do NOT consult aliases.
            if last_hmm_failure is None:
                last_hmm_failure = MarkerFailure(
                    reason="hmm_failed",
                    marker_name=spec["name"] or ",".join(hmms),
                    failed_tokens=list(hmms),
                    best_evalue=_best_evalue_across_token_hmms(with_seq, hmms),
                )
            continue

        # Alias-only spec (or HMM tier inactive). Try the aliases.
        if not spec["aliases"]:
            # Nothing can identify this spec's protein here; it reaches the
            # legacy longest-CDS fallback, which must honour its bound.
            if bounds and fallback_bounds is None:
                fallback_bounds = bounds
                fallback_name = spec["name"] or ",".join(hmms)
            continue
        matched: list[dict] = []
        matched_in_range: list[dict] = []
        for alias in spec["aliases"]:
            needle = alias.lower().strip()
            if not needle:
                continue
            matches = [
                p for p in with_seq
                if needle in (p.get("product") or "").lower()
            ]
            if not matches:
                continue
            in_range = [p for p in matches if within_length_bounds(p, bounds)]
            pickable = [p for p in in_range if _met_ok(p)]
            if pickable:
                return _longest(pickable), None
            matched.extend(matches)
            matched_in_range.extend(in_range)
        if met_required and matched_in_range and last_met_failure is None:
            # An alias matched an in-range CDS, but every match is
            # 5'-partial. Strict for the same reason a bounded alias spec
            # is: falling through to "the longest CDS" could substitute an
            # unrelated protein that happens to be complete.
            last_met_failure = MarkerFailure(
                reason="no_met_start",
                marker_name=spec["name"] or ",".join(spec["aliases"]),
                length_aa=len(_longest(matched_in_range)["sequence"]),
                bounds=bounds,
            )
            continue
        if bounds:
            # A bounded alias-only spec is strict: with no in-range alias
            # match it never falls back to "the longest CDS", which could be
            # an unrelated protein that happens to fit the range. An alias
            # that matched only out-of-range CDSes is a length failure; one
            # that matched nothing at all means the protein was not found,
            # which is a missing marker, not a partial one.
            name = spec["name"] or ",".join(spec["aliases"])
            if matched:
                if last_length_failure is None:
                    last_length_failure = MarkerFailure(
                        reason="length_out_of_bounds",
                        marker_name=name,
                        length_aa=len(_longest(matched)["sequence"]),
                        bounds=bounds,
                    )
            elif last_missing_failure is None:
                last_missing_failure = MarkerFailure(
                    reason="no_alias_match", marker_name=name,
                )

    if last_hmm_failure is not None:
        # Any HMM-spec failure is a HARD drop — never fall through to
        # longest CDS, since the user has declared the HMM authoritative
        # for that marker.
        return None, last_hmm_failure
    if last_length_failure is not None:
        # A bound is a hard requirement too: no fall-through. Length is
        # reported ahead of a completeness miss so the two stages' counters
        # agree with the reason written to _qc_removed.tsv.
        return None, last_length_failure
    if last_met_failure is not None:
        # The marker was found at an acceptable length but is 5'-partial.
        return None, last_met_failure
    if last_missing_failure is not None:
        # A bounded alias-only spec found no protein: strict, no fall-through.
        return None, last_missing_failure
    if fallback_bounds is not None:
        # Identity-less spec with a bound (e.g. HMM tier soft-failed):
        # longest CDS, but only inside the range and complete.
        in_range = [p for p in with_seq if within_length_bounds(p, fallback_bounds)]
        pickable = [p for p in in_range if _met_ok(p)]
        if pickable:
            return _longest(pickable), None
        if in_range:
            return None, MarkerFailure(
                reason="no_met_start",
                marker_name=fallback_name,
                length_aa=len(_longest(in_range)["sequence"]),
                bounds=fallback_bounds,
            )
        return None, MarkerFailure(
            reason="length_out_of_bounds",
            marker_name=fallback_name,
            length_aa=len(_longest(with_seq)["sequence"]),
            bounds=fallback_bounds,
        )
    # Legacy behaviour: an unbounded alias-only spec that matched no product
    # is advisory — fall through to the longest CDS rather than dropping.
    # Same for a run with only empty specs (no aliases AND no hmms). Both
    # still honour require_met_start, so clustering never takes a 5'-partial
    # CDS once the run asks for complete markers. A spec whose alias DID
    # match, but only truncated CDSes, is strict: it set last_met_failure
    # above, which is returned before this point.
    pool = [p for p in with_seq if _met_ok(p)]
    if not pool:
        return None, MarkerFailure(
            reason="no_met_start",
            length_aa=len(_longest(with_seq)["sequence"]),
        )
    return _longest(pool), None

def _ensure_passing_annotation(
    proteins: list[dict],
    ga_cutoffs: Optional[dict[str, Optional[float]]],
    hmm_cfg: dict[str, Any],
) -> None:
    """Annotate each hit with ``passing`` if it isn't already.

    The pipeline path goes through ``cli._run_hmm_qc`` which sets the
    flag once per hit. Direct unit tests of ``select_marker_protein``
    construct hits by hand without the flag — this helper covers that
    case so the test surface stays small.
    """
    use_ga = hmm_cfg.get("use_ga_when_available", True)
    default_evalue = hmm_cfg.get("default_evalue", 1.0e-5)
    rel_len = hmm_cfg.get("relative_length_cutoff", 0.5)
    for p in proteins:
        for hit in (p.get("hmm_hits") or []):
            if "passing" in hit:
                continue
            ga = None
            if use_ga and ga_cutoffs is not None:
                ga = ga_cutoffs.get(hit["target"])
            if ga is not None:
                similarity_pass = hit["dom_score"] >= ga
            else:
                similarity_pass = hit["dom_evalue"] <= default_evalue
            hit["passing"] = similarity_pass and (coverage_of(hit) >= rel_len)


def _best_evalue_across_token_hmms(
    proteins: list[dict], tokens: list[str]
) -> Optional[float]:
    """Best (smallest) passing-hit E-value seen for any HMM appearing in any
    of the failed tokens. Used to surface partial evidence in MarkerFailure
    so the _qc_removed.tsv reason is informative for triage.
    """
    names: set[str] = set()
    for token in tokens:
        try:
            names.update(parse_hmm_token(token))
        except ValueError:
            continue
    if not names:
        return None
    best: Optional[float] = None
    for p in proteins:
        for hit in _passing_hits(p):
            if hit["target"] not in names:
                continue
            ev = hit.get("dom_evalue")
            if ev is None:
                continue
            if best is None or ev < best:
                best = ev
    return best


def select_concatenated_markers(
    proteins: Optional[list[dict]],
    marker_specs: Any = None,
    *,
    hmm_active: bool = False,
    ga_cutoffs: Optional[dict[str, Optional[float]]] = None,
    hmm_cfg: Optional[dict[str, Any]] = None,
    ignore_restrictions: bool = False,
    require_met_start: bool = False,
) -> tuple[Optional[list[dict]], Optional[MarkerFailure]]:
    """Pick the marker CDS from EVERY spec, for non-segmented concat clustering.

    Returns ``(markers, None)`` — one CDS per declared spec, in declared
    order — or ``(None, MarkerFailure)`` on the first spec that yields no
    marker. Each spec is selected independently via
    :func:`select_marker_protein` (so the per-spec HMM/alias gating is
    identical to single-marker mode); the caller concatenates the AA
    bodies in the returned order.

    When no specs are configured the single longest CDS is returned as a
    one-element list (parity with the no-specs branch of
    ``select_marker_protein``). A spec failure short-circuits with that
    spec's ``MarkerFailure`` so the drop reason is informative.
    """
    if not proteins:
        return None, MarkerFailure("no_proteins")
    specs = _normalize_marker_specs(marker_specs)
    if not specs:
        marker, failure = select_marker_protein(
            proteins, marker_specs,
            hmm_active=hmm_active, ga_cutoffs=ga_cutoffs, hmm_cfg=hmm_cfg,
            ignore_restrictions=ignore_restrictions, require_met_start=require_met_start,
        )
        return ([marker], None) if marker is not None else (None, failure)

    markers: list[dict] = []
    for spec in specs:
        marker, failure = select_marker_protein(
            proteins, [spec],
            hmm_active=hmm_active, ga_cutoffs=ga_cutoffs, hmm_cfg=hmm_cfg,
            ignore_restrictions=ignore_restrictions, require_met_start=require_met_start,
        )
        if marker is None:
            return None, failure
        markers.append(marker)
    return markers, None


def populate_protein_sequences(
    sequences: list[Sequence],
    marker_specs: Any = None,
    report: Optional[QCReport] = None,
    *,
    hmm_active: bool = False,
    ga_cutoffs: Optional[dict[str, Optional[float]]] = None,
    hmm_cfg: Optional[dict[str, Any]] = None,
    concatenate: bool = False,
    policy=None,
    require_met_start: bool = False,
) -> list[Sequence]:
    """Set ``seq.protein_sequence`` on each sequence to its marker protein.

    Used by non-segmented inputs only — segmented isolates get their
    concatenated marker via ``build_concatenated_sequences``.

    When ``concatenate`` is true (``clustering.concatenate_markers``), the
    marker CDS from EVERY spec is selected and the AA bodies are joined in
    declared spec order (e.g. Spike+Nucleocapsid); a sequence missing any
    required marker is dropped. When false (default) the single marker
    from the first satisfying spec is used (legacy behaviour).

    Sequences that fail marker selection are dropped. Drops with an
    HMM-tier failure are counted under ``report.removed_hmm_failed``
    and broken out by marker in ``report.removed_hmm_by_marker``; a
    ``length_aa`` failure under ``report.removed_marker_length`` /
    ``removed_marker_length_by_marker``. All other failure modes are
    counted under ``report.removed_proteins``. Returns the surviving list.

    ``require_met_start`` (``qc.marker_completeness``) additionally requires
    the picked CDS to start with Met; such a failure is counted under
    ``report.removed_marker_completeness`` / ``_by_marker``.

    ``policy`` (an ``overrides.ProtectionPolicy``): a sequence protected
    against the ``marker_quality`` stage is picked as if neither the length
    bound nor the Met requirement existed, matching what the output pickers
    show for it.

    Note: in v0.14.0+ this should rarely fire under HMM control, because
    ``cli._run_hmm_qc`` has already dropped sequences whose HMM gate
    failed. The HMM-failure branch here remains as a defence-in-depth /
    direct-test entry point.
    """
    def _select(seq: Sequence, ignore_restrictions: bool):
        if concatenate:
            return select_concatenated_markers(
                seq.proteins,
                marker_specs,
                hmm_active=hmm_active,
                ga_cutoffs=ga_cutoffs,
                hmm_cfg=hmm_cfg,
                ignore_restrictions=ignore_restrictions,
                require_met_start=require_met_start,
            )
        marker, failure = select_marker_protein(
            seq.proteins,
            marker_specs,
            hmm_active=hmm_active,
            ga_cutoffs=ga_cutoffs,
            hmm_cfg=hmm_cfg,
            ignore_restrictions=ignore_restrictions,
            require_met_start=require_met_start,
        )
        return ([marker] if marker is not None else None), failure

    kept: list[Sequence] = []
    for seq in sequences:
        markers, failure = _select(seq, False)
        if (
            markers is None
            and failure is not None
            and failure.reason in ("length_out_of_bounds", "no_met_start")
            and policy is not None
            and policy.protects(seq, "marker_quality")
        ):
            markers, failure = _select(seq, True)
        if markers is None:
            if report is not None:
                reason = _format_failure_reason(failure)
                seq.qc_passed = False
                seq.qc_fail_reason = reason
                if failure is not None and failure.reason == "hmm_failed":
                    report.removed_hmm_failed += 1
                    key = failure.marker_name or "?"
                    report.removed_hmm_by_marker[key] = (
                        report.removed_hmm_by_marker.get(key, 0) + 1
                    )
                elif failure is not None and failure.reason == "length_out_of_bounds":
                    report.removed_marker_length += 1
                    key = failure.marker_name or "?"
                    report.removed_marker_length_by_marker[key] = (
                        report.removed_marker_length_by_marker.get(key, 0) + 1
                    )
                elif failure is not None and failure.reason == "no_met_start":
                    report.removed_marker_completeness += 1
                    key = failure.marker_name or "?"
                    report.removed_marker_completeness_by_marker[key] = (
                        report.removed_marker_completeness_by_marker.get(key, 0) + 1
                    )
                else:
                    report.removed_proteins += 1
                report.add_removed(seq.id, reason)
            continue
        seq.protein_sequence = "".join(m["sequence"] for m in markers)
        ids = [m["protein_id"] for m in markers if m.get("protein_id")]
        if ids:
            seq.marker_protein_ids = ids
        kept.append(seq)
    return kept


def _format_failure_reason(failure: Optional[MarkerFailure]) -> str:
    """Render a MarkerFailure as the string written to _qc_removed.tsv.

    Non-HMM failures keep the legacy ``no_marker_protein_for_clustering``
    string for back-compat with downstream parsers; HMM failures get the
    structured ``hmm_failed:<marker>:<tokens>(E=...)`` form and length
    failures ``marker_length:<marker>:<len>aa_below_<min>`` (or
    ``_above_<max>`` / ``no_candidate``) and 5'-partial markers
    ``marker_completeness:<marker>:no_met_start``, so they're visually
    distinct in the output and grep-able for triage.
    """
    if failure is None:
        return "no_marker_protein_for_clustering"
    if failure.reason == "length_out_of_bounds":
        return (
            f"marker_length:{failure.marker_name or '?'}:"
            f"{describe_length_miss(failure.length_aa, failure.bounds)}"
        )
    if failure.reason == "no_met_start":
        return f"marker_completeness:{failure.marker_name or '?'}:no_met_start"
    if failure.reason == "hmm_failed":
        parts = ["hmm_failed", failure.marker_name or "?"]
        if failure.failed_tokens:
            parts.append(",".join(failure.failed_tokens))
        suffix = ":".join(parts)
        if failure.best_evalue is not None:
            suffix = f"{suffix}(E={failure.best_evalue:.2g})"
        return suffix
    return "no_marker_protein_for_clustering"
