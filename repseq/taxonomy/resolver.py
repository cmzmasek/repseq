"""Metadata resolver: DB-first, header-fallback, strain-label-fallback."""

from __future__ import annotations

import logging
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Any, Optional

from ..models import Sequence, SequenceSource, TaxonomyInfo
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

    def _db_query(self, seq: Sequence) -> Optional[dict[str, Any]]:
        """Query the appropriate database based on sequence source."""
        acc = seq.accession
        if not acc:
            return None

        if seq.source == SequenceSource.UNIPROT:
            return self._uniprot.fetch_entry(acc)

        if seq.source in (SequenceSource.NCBI, SequenceSource.NCBI_VIRUS):
            meta = self._ncbi.fetch_accession_metadata(acc)
            if meta and meta.get("taxid") and not (meta.get("lineage")):
                lineage = self._ncbi.fetch_lineage(meta["taxid"])
                if lineage:
                    meta["lineage"] = lineage
            return meta

        # Unknown source: try NCBI first, then UniProt
        meta = self._ncbi.fetch_accession_metadata(acc)
        if meta:
            return meta
        return self._uniprot.fetch_entry(acc)

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
