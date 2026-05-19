"""Marker-protein selection for amino-acid clustering.

When ``clustering.alphabet_for_clustering`` is ``protein``, every
sequence (non-segmented) or segment (segmented) must contribute one
amino-acid string to the clustering input. This module picks that
string from the sequence's CDS list (``seq.proteins``).

Marker selection follows the configured marker specs. Each spec is
either an alias string (legacy form) or a dict ``{name, aliases?,
hmms?}``. Specs are tried in order; the first to yield a passing CDS
wins.

Per-spec gating:
    * If the spec defines ``hmms`` AND the HMM tier is active (hmmscan
      available, DB indexed), HMM hits are AUTHORITATIVE. A CDS passes
      this spec only if every HMM in the list has at least one hit on
      it that clears both the E-value/GA threshold and the relative-
      length threshold. Aliases for this spec are NOT consulted —
      that's the whole point of the strict gate.
    * Otherwise, the spec's aliases are matched case-insensitively as
      substrings against /product. First spec with any alias match
      wins; among multiple CDSes matching the same alias, the longest
      is returned.

When no specs are configured at all, the longest CDS with a translation
is returned (legacy fallback). When specs ARE configured but none
yield a passing CDS, the function returns ``(None, MarkerFailure)`` so
the caller can drop the sequence/isolate with a structured reason.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Optional

from ..models import QCReport, Sequence


@dataclass
class MarkerFailure:
    """Structured reason why ``select_marker_protein`` returned no marker.

    ``reason`` is one of:
        - ``no_proteins``     — proteins list is None or empty
        - ``no_translated``   — every CDS lacks a /translation
        - ``hmm_failed``      — an HMM-gated spec had no CDS pass; the
                                most-specific HMM failure is reported
                                in ``failed_hmms`` + ``best_evalue``
        - ``no_alias_match``  — only alias-only specs were tried and
                                none matched
    """

    reason: str
    marker_name: Optional[str] = None
    failed_hmms: list[str] = field(default_factory=list)
    best_evalue: Optional[float] = None


def _normalize_marker_specs(specs: Any) -> list[dict]:
    """Normalize a mixed list of strings / dicts into uniform spec dicts.

    Each output dict has keys ``name``, ``aliases``, ``hmms`` (lists are
    always present even if empty). A bare string ``"foo"`` normalises
    to ``{name: "foo", aliases: ["foo"], hmms: []}``.
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


def _hit_passes(
    hit: dict,
    ga_cutoffs: Optional[dict[str, Optional[float]]],
    hmm_cfg: Optional[dict[str, Any]],
) -> bool:
    """Did this hit pass the configured E-value/GA + coverage cutoffs?

    The HMM scan step in ``cli._run_hmm_scan`` pre-annotates each hit
    with a ``passing`` bool, so the fast path just reads that. Fall-
    back recompute (used by direct unit tests of ``select_marker_protein``
    that build hits by hand) mirrors ``hmm.runner.passes_cutoffs`` —
    inlined to avoid importing the hmm package from clustering, which
    must work without HMMER installed when alphabet=nucleotide.
    """
    if "passing" in hit:
        return bool(hit["passing"])
    if hmm_cfg is None:
        return False
    use_ga = hmm_cfg.get("use_ga_when_available", True)
    default_evalue = hmm_cfg.get("default_evalue", 1.0e-5)
    rel_len = hmm_cfg.get("relative_length_cutoff", 0.5)
    ga = None
    if use_ga and ga_cutoffs is not None:
        ga = ga_cutoffs.get(hit["target"])
    if ga is not None:
        similarity_pass = hit["dom_score"] >= ga
    else:
        similarity_pass = hit["dom_evalue"] <= default_evalue
    hmm_len = max(int(hit.get("hmm_len", 0)), 1)
    coverage = hit["ali_span"] / hmm_len
    return similarity_pass and coverage >= rel_len


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
    """
    if not proteins:
        return None, MarkerFailure("no_proteins")
    with_seq = [p for p in proteins if p.get("sequence")]
    if not with_seq:
        return None, MarkerFailure("no_translated")

    specs = _normalize_marker_specs(marker_specs)

    if not specs:
        # No markers configured → legacy longest-CDS fallback.
        return max(with_seq, key=lambda p: len(p["sequence"])), None

    # Track the most informative HMM failure we encountered. HMM
    # failures are strict: when a spec defines hmms and none of the
    # configured HMMs hit on any CDS, the function returns no marker.
    # Alias-only specs preserve the legacy fall-through to longest CDS
    # when nothing matches — aliases were always advisory, and changing
    # that silently would surprise existing users.
    any_hmm_attempted = False
    last_hmm_failure: Optional[MarkerFailure] = None
    any_alias_only_attempted = False

    for spec in specs:
        hmms = spec["hmms"]
        if hmms and hmm_active:
            any_hmm_attempted = True
            passing: list[dict] = []
            best_evalues: dict[str, float] = {}
            for p in with_seq:
                pass_targets: set[str] = set()
                for hit in (p.get("hmm_hits") or []):
                    if _hit_passes(hit, ga_cutoffs, hmm_cfg):
                        pass_targets.add(hit["target"])
                    ev = hit.get("dom_evalue")
                    if ev is not None:
                        t = hit["target"]
                        if t not in best_evalues or ev < best_evalues[t]:
                            best_evalues[t] = ev
                if all(h in pass_targets for h in hmms):
                    passing.append(p)
            if passing:
                return max(passing, key=lambda p: len(p["sequence"])), None
            failed = [
                h for h in hmms
                if all(
                    h not in {
                        hit["target"]
                        for hit in (p.get("hmm_hits") or [])
                        if _hit_passes(hit, ga_cutoffs, hmm_cfg)
                    }
                    for p in with_seq
                )
            ]
            best_e = None
            for h in failed:
                if h in best_evalues:
                    if best_e is None or best_evalues[h] < best_e:
                        best_e = best_evalues[h]
            last_hmm_failure = MarkerFailure(
                reason="hmm_failed",
                marker_name=spec["name"] or ",".join(hmms),
                failed_hmms=failed,
                best_evalue=best_e,
            )
            # Strict-gate: do NOT consult this spec's aliases as a
            # fallback for an HMM failure. Move on to the next spec.
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
        # Legacy behaviour: alias-only specs are advisory; if none
        # matched, fall through to longest CDS rather than dropping.
        return max(with_seq, key=lambda p: len(p["sequence"])), None
    # All specs were empty (no aliases AND no hmms). Treat as no specs.
    return max(with_seq, key=lambda p: len(p["sequence"])), None


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
    get the structured ``hmm_failed:<marker>:<hmm>(E=...)`` form so
    they're visually distinct in the output and grep-able for triage.
    """
    if failure is None:
        return "no_marker_protein_for_clustering"
    if failure.reason == "hmm_failed":
        parts = ["hmm_failed", failure.marker_name or "?"]
        if failure.failed_hmms:
            parts.append(",".join(failure.failed_hmms))
        suffix = ":".join(parts)
        if failure.best_evalue is not None:
            suffix = f"{suffix}(E={failure.best_evalue:.2g})"
        return suffix
    return "no_marker_protein_for_clustering"
