"""Overrides for sequences of special importance.

A user can name sequences that should be treated specially regardless of
what the automated pipeline would do. This module implements the
**QC-protection** half: named sequences bypass the QC removal stages they
would otherwise fail (a force-keep / whitelist), so a curated reference
strain that is slightly noisy, missing an expected protein, or just below a
length bound is never silently dropped.

Two capabilities are planned over one id list; this module currently
implements `protect_qc` (force-keep). `force_select` (guaranteed
representative) is a separate, later step.

Config (``overrides:`` block)::

    overrides:
      ids: ["NC_045512.2"]      # accession (non-seg) / isolate_id (seg)
      ids_file: vip.txt         # one id per line; '#' comments allowed
      protect_qc: true
      protect_stages: all       # or a subset of QC_PROTECT_STAGES

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
    """Lower-case, whitespace-collapse, and pipe-strip an id for matching.

    Mirrors the normalisation used for ``complete_isolates`` keys so an
    isolate id round-trips identically. Returns ``None`` for empty input.
    """
    if value is None:
        return None
    s = _WS.sub(" ", str(value).strip().strip("|")).lower()
    return s or None


def _strip_version(norm: str) -> str:
    """Drop a trailing ``.<digits>`` accession version from a normalised id."""
    return _VERSION.sub("", norm)


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
    """Decides whether a sequence is protected against a given QC stage."""

    def __init__(self, ids: frozenset[str], stages: frozenset[str], enabled: bool):
        self.ids = ids
        self.stages = stages
        # Only active when force-keep is on AND at least one id is listed —
        # an empty list with protect_qc=true protects nothing.
        self.enabled = bool(enabled and ids)

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
            )
        ov = cfg.get("overrides") or {}
        return cls(
            ids=resolve_ids(ov),
            stages=resolve_stages(ov.get("protect_stages", "all")),
            enabled=bool(ov.get("protect_qc", False)),
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
