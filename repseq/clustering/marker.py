"""Marker-protein selection for amino-acid clustering.

When ``clustering.alphabet_for_clustering`` is ``protein``, every
sequence (non-segmented) or segment (segmented) must contribute one
amino-acid string to the clustering input. This module picks that
string from the sequence's CDS list (``seq.proteins``).

Marker selection follows the configured marker specs. Each spec is
either an alias string (legacy form) or a dict ``{name, aliases?,
hmms?}`` where ``hmms`` is a list of **token** strings:

    - ``"Name"``         — single-HMM token
    - ``"A--B--C"``      — multidomain token, HMMs listed in C-to-N
                           order (A most C-terminal, C most N-terminal)

The semantic for v0.14.0 (hard cutover from v0.13.0's list-AND):

    - A CDS satisfies a single-HMM token when that HMM has a passing
      hit on it.
    - A CDS satisfies a multidomain token when every named HMM has a
      passing hit AND the hits appear in C-to-N order on the protein
      (strict non-overlap). Extra domains on the same CDS are fine.
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
        - ``no_alias_match``  — only alias-only specs were tried and
                                none matched
    """

    reason: str
    marker_name: Optional[str] = None
    failed_tokens: list[str] = field(default_factory=list)
    best_evalue: Optional[float] = None


def _normalize_marker_specs(specs: Any) -> list[dict]:
    """Normalize a mixed list of strings / dicts into uniform spec dicts.

    Each output dict has keys ``name``, ``aliases``, ``hmms`` (lists are
    always present even if empty). A bare string ``"foo"`` normalises to
    ``{name: "foo", aliases: ["foo"], hmms: []}``.
    """
    out: list[dict] = []
    if not specs:
        return out
    for entry in specs:
        if isinstance(entry, str):
            s = entry.strip()
            if s:
                out.append({"name": s, "aliases": [s], "hmms": []})
        elif isinstance(entry, dict):
            name = entry.get("name") or ""
            aliases = list(entry.get("aliases") or [])
            hmms = list(entry.get("hmms") or [])
            if name or aliases or hmms:
                out.append({"name": name, "aliases": aliases, "hmms": hmms})
    return out


def _passing_hits(p: dict) -> list[dict]:
    """Return the CDS's passing HMM hits (pre-annotated by _run_hmm_qc)."""
    return [h for h in (p.get("hmm_hits") or []) if h.get("passing")]


def _cds_satisfies_any_token(p: dict, tokens: list[str]) -> Optional[str]:
    """Return the token string this CDS satisfies (any token), or None.

    Used by marker selection to pick the longest CDS satisfying at least
    one token in the spec.
    """
    hits = p.get("hmm_hits") or []
    for token in tokens:
        try:
            parsed = parse_hmm_token(token)
        except ValueError:
            continue
        if cds_satisfies_token(hits, parsed) is not None:
            return token
    return None


def select_marker_protein(
    proteins: Optional[list[dict]],
    marker_specs: Any = None,
    *,
    hmm_active: bool = False,
    ga_cutoffs: Optional[dict[str, Optional[float]]] = None,
    hmm_cfg: Optional[dict[str, Any]] = None,
) -> tuple[Optional[dict], Optional[MarkerFailure]]:
    """Pick the marker CDS for one sequence/segment.

    Returns ``(marker_dict, None)`` on success or
    ``(None, MarkerFailure)`` on failure. See module docstring for the
    gating semantics. ``marker_specs`` accepts either the legacy alias
    string list (back-compat with pre-v0.13 callers) or the new mixed
    list of strings / dicts.

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

    # Test back-compat: if hits lack the pre-computed `passing` flag,
    # recompute it on the fly from the supplied cutoffs.
    if hmm_active and hmm_cfg is not None:
        _ensure_passing_annotation(with_seq, ga_cutoffs, hmm_cfg)

    specs = _normalize_marker_specs(marker_specs)

    if not specs:
        # No markers configured → legacy longest-CDS fallback.
        return max(with_seq, key=lambda p: len(p["sequence"])), None

    last_hmm_failure: Optional[MarkerFailure] = None
    any_alias_only_attempted = False

    for spec in specs:
        hmms = spec["hmms"]
        if hmms and hmm_active:
            # HMM tier is authoritative for this spec.
            satisfying: list[dict] = []
            for p in with_seq:
                if _cds_satisfies_any_token(p, hmms) is not None:
                    satisfying.append(p)
            if satisfying:
                return max(satisfying, key=lambda p: len(p["sequence"])), None
            # Record the failure and move on — do NOT consult aliases.
            last_hmm_failure = MarkerFailure(
                reason="hmm_failed",
                marker_name=spec["name"] or ",".join(hmms),
                failed_tokens=list(hmms),
                best_evalue=_best_evalue_across_token_hmms(with_seq, hmms),
            )
            continue

        # Alias-only spec (or HMM tier inactive). Try the aliases.
        if spec["aliases"]:
            any_alias_only_attempted = True
        for alias in spec["aliases"]:
            needle = alias.lower().strip()
            if not needle:
                continue
            matches = [
                p for p in with_seq
                if needle in (p.get("product") or "").lower()
            ]
            if matches:
                return max(matches, key=lambda p: len(p["sequence"])), None

    if last_hmm_failure is not None:
        # Any HMM-spec failure is a HARD drop — never fall through to
        # longest CDS, since the user has declared the HMM authoritative
        # for that marker.
        return None, last_hmm_failure
    if any_alias_only_attempted:
        # Legacy behaviour: alias-only specs are advisory; if none matched,
        # fall through to longest CDS rather than dropping.
        return max(with_seq, key=lambda p: len(p["sequence"])), None
    # All specs were empty (no aliases AND no hmms). Treat as no specs.
    return max(with_seq, key=lambda p: len(p["sequence"])), None


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


def populate_protein_sequences(
    sequences: list[Sequence],
    marker_specs: Any = None,
    report: Optional[QCReport] = None,
    *,
    hmm_active: bool = False,
    ga_cutoffs: Optional[dict[str, Optional[float]]] = None,
    hmm_cfg: Optional[dict[str, Any]] = None,
) -> list[Sequence]:
    """Set ``seq.protein_sequence`` on each sequence to its marker protein.

    Used by non-segmented inputs only — segmented isolates get their
    concatenated marker via ``build_concatenated_sequences``.

    Sequences that fail marker selection are dropped. Drops with an
    HMM-tier failure are counted under ``report.removed_hmm_failed``
    and broken out by marker in ``report.removed_hmm_by_marker``. All
    other failure modes are counted under ``report.removed_proteins``.
    Returns the surviving list.

    Note: in v0.14.0+ this should rarely fire under HMM control, because
    ``cli._run_hmm_qc`` has already dropped sequences whose HMM gate
    failed. The HMM-failure branch here remains as a defence-in-depth /
    direct-test entry point.
    """
    kept: list[Sequence] = []
    for seq in sequences:
        marker, failure = select_marker_protein(
            seq.proteins,
            marker_specs,
            hmm_active=hmm_active,
            ga_cutoffs=ga_cutoffs,
            hmm_cfg=hmm_cfg,
        )
        if marker is None:
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
                else:
                    report.removed_proteins += 1
                report.add_removed(seq.id, reason)
            continue
        seq.protein_sequence = marker["sequence"]
        if marker.get("protein_id"):
            seq.marker_protein_ids = [marker["protein_id"]]
        kept.append(seq)
    return kept


def _format_failure_reason(failure: Optional[MarkerFailure]) -> str:
    """Render a MarkerFailure as the string written to _qc_removed.tsv.

    Non-HMM failures keep the legacy ``no_marker_protein_for_clustering``
    string for back-compat with downstream parsers; only HMM failures
    get the structured ``hmm_failed:<marker>:<tokens>(E=...)`` form so
    they're visually distinct in the output and grep-able for triage.
    """
    if failure is None:
        return "no_marker_protein_for_clustering"
    if failure.reason == "hmm_failed":
        parts = ["hmm_failed", failure.marker_name or "?"]
        if failure.failed_tokens:
            parts.append(",".join(failure.failed_tokens))
        suffix = ":".join(parts)
        if failure.best_evalue is not None:
            suffix = f"{suffix}(E={failure.best_evalue:.2g})"
        return suffix
    return "no_marker_protein_for_clustering"
