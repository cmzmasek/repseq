"""Metadata resolver: DB-first, header-fallback, strain-label-fallback."""

from __future__ import annotations

import logging
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Any, Optional

from ..models import Sequence, SequenceSource, TaxonomyInfo
from . import http as _http
from .cache import TaxonomyCache
from .ncbi import NCBITaxonomy
from .uniprot import UniProtAPI

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Strain-label parsers (last-resort fallback)
# ---------------------------------------------------------------------------

# Influenza: A/duck/Hong_Kong/1/1997(H5N1)
_INFLUENZA_RE = re.compile(
    r"(?P<type>[AB])/(?P<host>[^/]+)/(?P<location>[^/]+)/(?P<number>[^/]+)/"
    r"(?P<year>\d{4})(?:\((?P<subtype>[^)]+)\))?"
)

# Strict influenza HA/NA subtype: "H" + 1-2 digits, optionally "N" + 1-2
# digits (H5N1, H3N2, H5). The strain regex above captures EVERYTHING inside
# the trailing parens (``[^)]+``), so the captured text is validated against
# this with ``fullmatch`` before it is accepted — otherwise arbitrary
# parenthetical notes after the year (e.g. "(mixed, partial sequence)") would
# be mistaken for a subtype. Deliberately strict: variant tags like
# "H1N1pdm09" are not accepted here (the authoritative subtype still comes
# from the GenBank /serotype path); this fallback only fills the clean case.
_SUBTYPE_RE = re.compile(r"H\d{1,2}(?:N\d{1,2})?", re.IGNORECASE)

# Generic year extraction from strain or description
_YEAR_RE = re.compile(r"\b(19\d{2}|20\d{2})\b")


def _parse_strain_label(seq: Sequence) -> dict[str, Any]:
    """Extract metadata from viral strain nomenclature."""
    target = seq.strain or seq.description or seq.header
    if not target:
        return {}

    result: dict[str, Any] = {}

    m = _INFLUENZA_RE.search(target)
    if m:
        if not seq.host:
            host = m.group("host")
            if host.lower() not in ("human", "unknown"):
                result["host"] = host
            else:
                result["host"] = "Homo sapiens" if host.lower() == "human" else host
        if not seq.country:
            result["country"] = m.group("location").replace("_", " ")
        if not seq.collection_date:
            result["collection_date"] = m.group("year")
        if not seq.subtype:
            sub = (m.group("subtype") or "").strip()
            # Accept only a strict HA/NA subtype — the capture group grabbed
            # everything in the parens, so fullmatch ensures arbitrary text
            # (e.g. "mixed, partial sequence") is never taken as a subtype.
            if sub and _SUBTYPE_RE.fullmatch(sub):
                result["subtype"] = sub.upper()
        return result

    # Generic year
    if not seq.collection_date:
        ym = _YEAR_RE.search(target)
        if ym:
            result["collection_date"] = ym.group(1)

    return result


# ---------------------------------------------------------------------------
# Resolver
# ---------------------------------------------------------------------------

class MetadataResolver:
    def __init__(
        self,
        cache: TaxonomyCache,
        ncbi: NCBITaxonomy,
        uniprot: UniProtAPI,
        threads: int = 4,
    ) -> None:
        self._cache = cache
        self._ncbi = ncbi
        self._uniprot = uniprot
        self._threads = threads
        self.failures: list[tuple[str, str]] = []  # (seq_id, error_message)

    # ------------------------------------------------------------------
    # Resolve a single sequence
    # ------------------------------------------------------------------

    def resolve(self, seq: Sequence) -> None:
        """Populate taxonomy and missing metadata fields on seq in-place."""
        # 1. Database query
        db_meta = self._db_query(seq)
        # A lookup that failed for a network reason (after the client's
        # retries) returns nothing, exactly like "no such record" — so ask
        # the failure log which it was, and count it. Without this the
        # record silently fell through to header parsing with no taxonomy.
        reason = self._lookup_failure(seq, db_meta)
        if reason:
            self.failures.append((seq.id, reason))

        # 2. Fill fields from the DB result.
        #    The database is authoritative: a DB-provided value overrides
        #    the heuristic header parse (which, for NCBI Virus headers, is
        #    fragile bracket-guessing). A field is only left to the header
        #    value when the DB has nothing for it.
        if db_meta:
            if db_meta.get("organism"):
                seq.organism = db_meta["organism"]
            if db_meta.get("host"):
                seq.host = db_meta["host"]
            if db_meta.get("collection_date"):
                seq.collection_date = db_meta["collection_date"]
            if db_meta.get("country"):
                seq.country = db_meta["country"]
            if db_meta.get("strain"):
                seq.strain = db_meta["strain"]
            if db_meta.get("subtype"):
                seq.subtype = db_meta["subtype"]
            if db_meta.get("is_reviewed") is not None:
                seq.is_reviewed = db_meta["is_reviewed"]

            # Build TaxonomyInfo from DB result
            lineage_data = db_meta.get("lineage") or db_meta
            seq.taxonomy = _build_taxonomy(lineage_data)
            if seq.taxonomy and not seq.taxonomy.taxid and db_meta.get("taxid"):
                seq.taxonomy.taxid = db_meta["taxid"]

        # 3. Fallback: strain label parsing
        if _needs_metadata(seq):
            extra = _parse_strain_label(seq)
            for field, value in extra.items():
                if not getattr(seq, field, None):
                    setattr(seq, field, value)

    @staticmethod
    def _lookup_failure(
        seq: Sequence, db_meta: Optional[dict[str, Any]],
    ) -> Optional[str]:
        """Why this record's database lookup came back short, if a network
        failure (rather than a genuine "not found") is the reason."""
        acc = seq.accession
        if not acc:
            return None
        if not db_meta:
            reason = _http.failure_for(_http.STAGE_METADATA, acc)
            return f"metadata lookup failed: {reason}" if reason else None
        taxid = db_meta.get("taxid")
        if taxid and not db_meta.get("lineage"):
            reason = _http.failure_for(_http.STAGE_LINEAGE_TAXID, taxid)
            if reason:
                msg = f"taxonomy lookup for taxid {taxid} failed: {reason}"
                # Pin it to THIS record. The taxid-keyed entry is shared
                # and is cleared as soon as any other record's fetch of the
                # same taxid succeeds, but this record has already been
                # resolved without its lineage and must stay on the report.
                _http.record_failure(_http.STAGE_LINEAGE, acc, msg)
                return msg
        return None

    def _ncbi_metadata(self, acc: str) -> Optional[dict[str, Any]]:
        """NCBI metadata for ``acc``, with the lineage back-filled.

        ``fetch_accession_metadata`` caches its result even when the lineage
        fetch inside it failed, so a cache hit can come back without one.
        Asking for the lineage again here is what lets a re-run recover it —
        on EVERY path that uses NCBI, the unknown-source one included (it
        used to skip this, leaving such a record without taxonomy until the
        cache entry expired).
        """
        meta = self._ncbi.fetch_accession_metadata(acc)
        if meta and meta.get("taxid") and not meta.get("lineage"):
            lineage = self._ncbi.fetch_lineage(meta["taxid"])
            if lineage:
                meta["lineage"] = lineage
        return meta

    def _db_query(self, seq: Sequence) -> Optional[dict[str, Any]]:
        """Query the appropriate database based on sequence source."""
        acc = seq.accession
        if not acc:
            return None

        if seq.source == SequenceSource.UNIPROT:
            return self._uniprot.fetch_entry(acc)

        if seq.source in (SequenceSource.NCBI, SequenceSource.NCBI_VIRUS):
            return self._ncbi_metadata(acc)

        # Unknown source: try NCBI first, then UniProt
        meta = self._ncbi_metadata(acc)
        if meta:
            return meta
        entry = self._uniprot.fetch_entry(acc)
        if entry:
            # Found at UniProt, so a failed NCBI attempt for the same
            # accession is not a failure of this record.
            _http.clear_failure(_http.STAGE_METADATA, acc)
        return entry

    # ------------------------------------------------------------------
    # Batch resolve with thread pool
    # ------------------------------------------------------------------

    def resolve_batch(self, sequences: list[Sequence], progress=None) -> None:
        """Resolve metadata for a list of sequences using a thread pool.

        Per-sequence resolution errors are non-fatal: they are logged and
        recorded on `self.failures` so the caller can summarise them, but
        the batch continues so one bad accession doesn't halt the run.
        """
        with ThreadPoolExecutor(max_workers=self._threads) as executor:
            futures = {executor.submit(self.resolve, seq): seq for seq in sequences}
            for future in as_completed(futures):
                seq = futures[future]
                try:
                    future.result()
                except Exception as exc:
                    msg = f"{type(exc).__name__}: {exc}"
                    self.failures.append((seq.id, msg))
                    _http.record_failure(
                        _http.STAGE_METADATA, seq.accession or seq.id, msg,
                    )
                    logger.warning("Resolution failed for %s: %s", seq.id, msg)
                if progress:
                    progress.update(1)
        if self.failures:
            logger.warning(
                "Metadata resolution: %d of %d sequences failed (see warnings above).",
                len(self.failures), len(sequences),
            )

    # ------------------------------------------------------------------
    # Taxonomy by taxid
    # ------------------------------------------------------------------

    def fetch_taxonomy(self, taxid: int, source: SequenceSource) -> Optional[TaxonomyInfo]:
        """Fetch and return TaxonomyInfo for a taxid."""
        if source == SequenceSource.UNIPROT:
            data = self._uniprot.fetch_lineage(taxid)
        else:
            data = self._ncbi.fetch_lineage(taxid)
        if data:
            return _build_taxonomy(data)
        return None


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _build_taxonomy(data: dict[str, Any]) -> TaxonomyInfo:
    lineage = data.get("lineage", {})
    if not isinstance(lineage, dict):
        lineage = {}
    return TaxonomyInfo(
        taxid=data.get("taxid"),
        species=data.get("species") or lineage.get("species"),
        genus=data.get("genus") or lineage.get("genus"),
        family=data.get("family") or lineage.get("family"),
        order=data.get("order") or lineage.get("order"),
        class_=data.get("class") or lineage.get("class"),
        phylum=data.get("phylum") or lineage.get("phylum"),
        kingdom=data.get("kingdom") or lineage.get("kingdom"),
        superkingdom=data.get("superkingdom") or lineage.get("superkingdom"),
        lineage=lineage,
    )


def _needs_metadata(seq: Sequence) -> bool:
    return not all([seq.host, seq.country, seq.collection_date])
