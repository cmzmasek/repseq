"""Overrides for sequences of special importance.

A user can name sequences that should be treated specially regardless of
what the automated pipeline would do. This module implements the
**QC-protection** half: named sequences bypass the QC removal stages they
would otherwise fail (a force-keep / whitelist), so a curated reference
strain that is slightly noisy, missing an expected protein, or just below a
length bound is never silently dropped.

Two capabilities over one id list, independently toggled:

* ``protect_qc`` — **force-keep**: named sequences bypass the QC removal
  stages they would otherwise fail (this module's stage helpers).
* ``force_select`` — **force-select**: named sequences are guaranteed to
  appear among the representatives, regardless of how clustering would have
  collapsed them (:func:`apply_force_select`).

A third capability, **exclude** (:func:`apply_exclusions`), is the mirror
image: named sequences are dropped from the input the moment it is read —
before metadata resolution, QC, or clustering — exactly as if they had been
deleted from the FASTA. It carries its own dedicated id list
(``overrides.exclude.ids`` / ``ids_file``), kept separate from the shared
keep/pin list above; ``validate_config`` rejects a config that lists the
same id under both (an id cannot be both removed and guaranteed).

Config (``overrides:`` block)::

    overrides:
      ids: ["NC_045512.2"]      # accession (non-seg) / isolate_id (seg)
      ids_file: vip.txt         # one id per line; '#' comments allowed
      protect_qc: true
      protect_stages: all       # or a subset of QC_PROTECT_STAGES
      force_select: true

**Matching** is by ``seq.accession`` OR ``seq.id`` OR ``seq.isolate_id``,
case- and whitespace-normalised, and **version-insensitive** (``NC_045512``
matches ``NC_045512.2`` and vice versa). In segmented mode you typically
list ``isolate_id``s; note that the early QC stages (``ambiguous``,
``annotation``) run before ``isolate_id`` is populated, so those only match
on accession.

**Transparency**: a protected sequence that *would* have been removed is
recorded on :class:`~repseq.models.QCReport` via ``add_protected`` and
surfaces in ``{prefix}_overrides.tsv`` and the run summary — protection is
never silent.
"""

from __future__ import annotations

import re
from typing import Any, Iterable, Optional

# QC removal stages a sequence can be protected against, in pipeline order.
# `completeness` is deliberately absent: protection cannot synthesize a
# missing segment, so an incomplete isolate still drops (with a distinct
# message) — it can't be concatenated or treed.
QC_PROTECT_STAGES: tuple[str, ...] = (
    "duplicates",
    "length",
    "ambiguous",
    "annotation",
    "protein_count",
    "taxonomy_consistency",
    "protein_quality",
    "hmm",
)

_WS = re.compile(r"\s+")
_VERSION = re.compile(r"\.\d+$")


def _norm_id(value: Optional[str]) -> Optional[str]:
    """Normalise an id for matching: lower-case, then replace every
    whitespace run AND pipe with ``_``. Returns ``None`` for empty input.

    This is the **single source of truth** for id normalisation:
    ``segmented.completeness._normalise_isolate_id`` (which builds the
    concatenated isolate's ``isolate_id`` / ``seq.id``) delegates to it, so
    the two cannot drift. They MUST agree, because a force_select /
    protect_qc / exclude id is matched against BOTH a per-segment
    ``seq.isolate_id`` (the raw strain, e.g. ``"A/Foo Bar/1"``) and the
    concatenated isolate's normalised id (``"a/foo_bar/1"``). The
    pre-2026-06 behaviour collapsed whitespace to a single space here while
    the isolate-id builder used ``_``; the natural space-form id then
    silently failed to bind any segmented isolate whose name contained
    whitespace — i.e. most influenza / segmented strain names.
    """
    if value is None:
        return None
    s = _WS.sub("_", str(value).strip().lower()).replace("|", "_")
    return s or None


def _strip_version(norm: str) -> str:
    """Return the version-insensitive match key for a normalised id.

    Drops a trailing ``.<digits>`` accession version (so ``nc_045512.2``
    matches ``nc_045512``), and first trims any leading/trailing ``_`` — the
    residue a delimiter pipe or edge whitespace leaves behind once
    :func:`_norm_id` has mapped it to ``_``. A pipe-bearing header id from
    the UNKNOWN-source FASTA fallback (e.g. ``AB123456.1|``) normalises to
    ``ab123456.1_``; without trimming the trailing ``_`` the ``\\.\\d+$``
    anchor can't reach the version, and a bare/versioned override id would
    silently fail to bind it. Trimming is applied symmetrically to both the
    configured ids and the sequence fields (every caller routes through this
    helper), so matching stays consistent and version-insensitive in either
    direction. Internal ``_`` (real accession underscores, ``nc_045512``) is
    untouched — only the ends are trimmed.
    """
    return _VERSION.sub("", norm.strip("_"))


def resolve_stages(spec: Any) -> frozenset[str]:
    """Resolve ``protect_stages`` (``"all"`` or a list) to a token set."""
    if spec is None or spec == "all":
        return frozenset(QC_PROTECT_STAGES)
    if isinstance(spec, str):
        spec = [spec]
    return frozenset(str(s).strip() for s in spec if str(s).strip())


def load_ids_file(path: str) -> list[str]:
    """Read one id per line from ``path`` (blank lines and ``#`` comments skipped).

    Raises ``FileNotFoundError`` if the file is missing — the caller turns
    that into a friendly, actionable error rather than a stack trace.
    """
    ids: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.split("#", 1)[0].strip()
            if line:
                ids.append(line)
    return ids


def resolve_raw_ids(ov: dict[str, Any]) -> tuple[str, ...]:
    """Return the original (un-normalised) configured ids, order-preserving.

    Inline ``ids`` + ``ids_file`` (best-effort, missing file ignored),
    de-duplicated. Used for the force-select "unavailable" audit, which
    must echo the id the user actually typed.
    """
    raw: list[str] = list(ov.get("ids") or [])
    path = ov.get("ids_file")
    if path:
        try:
            raw.extend(load_ids_file(path))
        except OSError:
            pass
    seen: set[str] = set()
    out: list[str] = []
    for r in raw:
        if r not in seen:
            seen.add(r)
            out.append(r)
    return tuple(out)


def resolve_ids(ov: dict[str, Any]) -> frozenset[str]:
    """Build the normalised, version-augmented id set from an ``overrides`` dict.

    Both the full normalised id and its version-stripped form are stored so
    matching is version-insensitive in either direction. ``ids_file`` (if
    present and readable) is unioned with the inline ``ids`` list; a missing
    file is the caller's responsibility to validate — here it is ignored so
    a library call (e.g. a test passing only inline ids) never crashes.
    """
    raw: list[str] = list(ov.get("ids") or [])
    path = ov.get("ids_file")
    if path:
        try:
            raw.extend(load_ids_file(path))
        except OSError:
            pass
    out: set[str] = set()
    for item in raw:
        n = _norm_id(item)
        if not n:
            continue
        out.add(n)
        out.add(_strip_version(n))
    return frozenset(out)


class ProtectionPolicy:
    """Decides whether a sequence is protected (QC) and/or pinned (selection)."""

    def __init__(self, ids: frozenset[str], stages: frozenset[str], enabled: bool,
                 *, force_select: bool = False, raw_ids: tuple[str, ...] = ()):
        self.ids = ids
        self.stages = stages
        self.raw_ids = raw_ids
        # Each capability is active only when its flag is on AND at least one
        # id is listed — an empty list protects / pins nothing.
        self.enabled = bool(enabled and ids)
        self.force_select = bool(force_select and ids)

    @classmethod
    def from_cfg(cls, cfg: dict[str, Any]) -> "ProtectionPolicy":
        """Build a policy from ``cfg``.

        Reads the pre-resolved ``cfg["_overrides_runtime"]`` cache when the
        CLI populated it (one id-file read per run); otherwise resolves the
        ``overrides`` block inline so direct callers / tests work without the
        CLI driver.
        """
        rt = cfg.get("_overrides_runtime")
        if rt is not None:
            return cls(
                ids=rt.get("ids", frozenset()),
                stages=rt.get("stages", frozenset()),
                enabled=rt.get("protect_qc", False),
                force_select=rt.get("force_select", False),
                raw_ids=rt.get("raw_ids", ()),
            )
        ov = cfg.get("overrides") or {}
        return cls(
            ids=resolve_ids(ov),
            stages=resolve_stages(ov.get("protect_stages", "all")),
            enabled=bool(ov.get("protect_qc", False)),
            force_select=bool(ov.get("force_select", False)),
            raw_ids=resolve_raw_ids(ov),
        )

    def _matches(self, seq) -> bool:
        for key in (getattr(seq, "accession", None), getattr(seq, "id", None),
                    getattr(seq, "isolate_id", None)):
            n = _norm_id(key)
            if n is not None and (n in self.ids or _strip_version(n) in self.ids):
                return True
        return False

    def protects(self, seq, stage: str) -> bool:
        """True iff ``seq`` is protected against removal at ``stage``."""
        if not self.enabled or stage not in self.stages:
            return False
        return self._matches(seq)

    def protects_any(self, seqs: Iterable, stage: str) -> bool:
        """True iff any sequence in ``seqs`` is protected against ``stage``.

        Used by the isolate-grouped stages (taxonomy-consistency,
        protein-quality, HMM) where naming any one segment protects the
        whole isolate from being dropped.
        """
        if not self.enabled or stage not in self.stages:
            return False
        return any(self._matches(s) for s in seqs)

    def pins(self, seq) -> bool:
        """True iff ``seq`` is force-selected (guaranteed a representative)."""
        return self.force_select and self._matches(seq)


def protected_keep(seq, stage: str, reason: str, policy: Optional[ProtectionPolicy],
                   report) -> bool:
    """Record-and-keep helper for a single-sequence QC stage.

    If ``policy`` protects ``seq`` against ``stage``, log the would-be
    removal ``reason`` on ``report`` (for ``_overrides.tsv`` / the summary)
    and return ``True`` so the caller keeps the sequence and skips its
    removal bookkeeping. Returns ``False`` (and records nothing) otherwise.
    """
    if policy is not None and policy.protects(seq, stage):
        report.add_protected(seq.id, stage, reason)
        return True
    return False


# ---------------------------------------------------------------------------
# Force-select (guaranteed representative)
# ---------------------------------------------------------------------------

def _audit_id(seq) -> str:
    """The most user-recognisable id for the force-select audit."""
    return (getattr(seq, "isolate_id", None)
            or getattr(seq, "accession", None)
            or getattr(seq, "id", "") or "")


def apply_force_select(result, pool, cfg) -> None:
    """Guarantee every ``force_select``-pinned sequence is a representative.

    **Hybrid policy** (the chosen design): within a cluster a pinned member
    *wins* the representative slot — beating the configured
    refseq/reviewed/longest priority (the demoted old representative drops
    back into the cluster's members). When several pinned members land in
    the same cluster, the best one (by that same priority) wins and the
    rest are *split* into their own singleton clusters, so all pinned
    members survive. Pinned sequences that were diversity-deselected (and
    so belong to no cluster — the ``global -n`` path) are *added* as new
    singleton representatives. Pinned ids that no surviving sequence
    matched (dropped in QC and not protected, or never in the input) are
    recorded as ``unavailable`` — the visible half of the
    force-keep/force-select coupling.

    Mutates ``result.representatives`` / ``result.clusters`` in place and
    writes an audit list onto ``result.force_selected`` (entries are
    ``{id, action, detail}``). No-op unless ``overrides.force_select`` is
    on and at least one id is listed. ``pool`` is the post-QC sequence list
    fed to the mode (the same objects that populate the clusters, so
    identity comparison is valid).
    """
    from .models import Cluster
    from .representative.selector import select_representative

    policy = ProtectionPolicy.from_cfg(cfg)
    if not policy.force_select:
        return
    priority = cfg.get("representative", {}).get(
        "priority", ["refseq", "reviewed_uniprot", "longest"]
    )

    pool = pool or []
    pinned_pool = [s for s in pool if policy.pins(s)]
    audit: list[dict] = []

    # Unavailable: configured ids that no surviving sequence matched.
    matched_forms: set[str] = set()
    for s in pinned_pool:
        for key in (getattr(s, "accession", None), getattr(s, "id", None),
                    getattr(s, "isolate_id", None)):
            n = _norm_id(key)
            if n is not None:
                matched_forms.add(n)
                matched_forms.add(_strip_version(n))
    for raw in policy.raw_ids:
        n = _norm_id(raw)
        if n is not None and n not in matched_forms and _strip_version(n) not in matched_forms:
            audit.append({
                "id": raw, "action": "unavailable",
                "detail": "no surviving sequence matched (dropped in QC, or "
                          "not in input)",
            })

    if pinned_pool:
        rep_ids = {id(r) for r in result.representatives}
        cluster_of: dict[int, Cluster] = {}
        for c in result.clusters:
            for m in [c.representative, *c.members]:
                cluster_of[id(m)] = c

        pins_by_cluster: dict[int, list] = {}
        cluster_by_key: dict[int, Cluster] = {}
        orphans: list = []
        for p in pinned_pool:
            if id(p) in rep_ids:
                audit.append({"id": _audit_id(p),
                              "action": "already_representative", "detail": ""})
                continue
            c = cluster_of.get(id(p))
            if c is None:
                orphans.append(p)
                continue
            pins_by_cluster.setdefault(id(c), []).append(p)
            cluster_by_key[id(c)] = c

        for ckey, pins in pins_by_cluster.items():
            cluster = cluster_by_key[ckey]
            winner = None
            if not policy._matches(cluster.representative):
                # Election: a pinned member beats the current representative.
                winner = select_representative(pins, priority)
                old = cluster.representative
                cluster.members = [old] + [m for m in cluster.members if m is not winner]
                cluster.representative = winner
                result.representatives = [
                    winner if r is old else r for r in result.representatives
                ]
                audit.append({"id": _audit_id(winner),
                              "action": "elected_representative",
                              "detail": f"cluster={cluster.cluster_id}"})
            # Split the remaining pinned members into singleton clusters.
            n = 0
            for p in pins:
                if p is winner:
                    continue
                n += 1
                cluster.members = [m for m in cluster.members if m is not p]
                result.clusters.append(
                    Cluster(cluster_id=f"{cluster.cluster_id}_pin{n}",
                            representative=p)
                )
                result.representatives.append(p)
                audit.append({"id": _audit_id(p), "action": "split_singleton",
                              "detail": f"from_cluster={cluster.cluster_id}"})

        for i, p in enumerate(orphans, 1):
            result.clusters.append(
                Cluster(cluster_id=f"pin_{i:06d}", representative=p)
            )
            result.representatives.append(p)
            audit.append({"id": _audit_id(p), "action": "added_representative",
                          "detail": "diversity-deselected; added as singleton"})

    result.force_selected = audit


# ---------------------------------------------------------------------------
# Exclude (input blocklist)
# ---------------------------------------------------------------------------

def apply_exclusions(sequences, cfg) -> tuple[list, list[dict]]:
    """Drop blocklisted sequences before any processing.

    Reads the resolved exclude id set from ``cfg["_overrides_runtime"]``
    (``exclude_ids`` / ``exclude_raw_ids``, populated by the CLI's
    ``_resolve_overrides``) and removes every input sequence whose
    ``accession`` or ``id`` matches — case- and version-insensitively, the
    same normalisation the keep/pin path uses. ``isolate_id`` is **not**
    consulted: this runs before the GenBank lookup that populates it, so it
    would always be empty here; matching is deliberately header-id-only,
    mirroring "delete this record from the FASTA".

    Returns ``(kept_sequences, audit)``. ``audit`` is a list of
    ``{id, action, detail}`` entries (parallel to the force-select audit):
    one ``action="excluded"`` row per dropped sequence, plus one
    ``action="unavailable"`` row per configured id that matched nothing (a
    typo-catcher). No-op (returns the input list and an empty audit) when no
    exclude id is configured.
    """
    rt = cfg.get("_overrides_runtime") or {}
    exclude_ids: frozenset[str] = rt.get("exclude_ids", frozenset())
    if not exclude_ids:
        return sequences, []

    kept: list = []
    audit: list[dict] = []
    matched_forms: set[str] = set()
    for seq in sequences:
        hit_on: Optional[str] = None
        for field in ("accession", "id"):
            n = _norm_id(getattr(seq, field, None))
            if n is not None and (n in exclude_ids or _strip_version(n) in exclude_ids):
                hit_on = field
                matched_forms.add(n)
                matched_forms.add(_strip_version(n))
                break
        if hit_on is not None:
            audit.append({
                "id": getattr(seq, "accession", None) or getattr(seq, "id", "") or "",
                "action": "excluded",
                "detail": f"matched on {hit_on}",
            })
        else:
            kept.append(seq)

    # Unavailable: configured ids that matched no input sequence.
    for raw in rt.get("exclude_raw_ids", ()):
        n = _norm_id(raw)
        if n is not None and n not in matched_forms and _strip_version(n) not in matched_forms:
            audit.append({
                "id": raw, "action": "unavailable",
                "detail": "no input sequence matched (typo, or already absent)",
            })

    return kept, audit
