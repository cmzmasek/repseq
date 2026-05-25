"""NCBI Entrez taxonomy and metadata queries."""

from __future__ import annotations

import logging
import threading
import time
from io import StringIO
from typing import Any, Optional

import requests

from .cache import TaxonomyCache

logger = logging.getLogger(__name__)

_ENTREZ_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
_RATE_LIMIT_DELAY = 0.34  # ~3 req/s without API key; 0.1 with key
_SOURCE = "ncbi_taxonomy"
_SOURCE_NUCCORE = "ncbi_nuccore"
_SOURCE_PROTEINS = "ncbi_proteins"
_GENBANK_BATCH_SIZE = 200
_NOT_FOUND: dict = {"_not_found": True}


class NCBITaxonomy:
    def __init__(
        self,
        cache: TaxonomyCache,
        email: Optional[str] = None,
        api_key: Optional[str] = None,
    ) -> None:
        self._cache = cache
        self._email = email
        self._api_key = api_key
        self._delay = 0.11 if api_key else _RATE_LIMIT_DELAY
        self._last_request: float = 0.0
        self._throttle_lock = threading.Lock()

    @property
    def cache(self) -> TaxonomyCache:
        """Read-only handle to the underlying TaxonomyCache.

        Exposed so siblings that share the cache (e.g. the HMM scan
        cache) don't have to reach into ``_cache`` privately.
        """
        return self._cache

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _auth_params(self, extra: dict) -> dict:
        """Attach email/api_key auth params without forcing a retmode."""
        p = dict(extra)
        if self._email:
            p["email"] = self._email
        if self._api_key:
            p["api_key"] = self._api_key
        return p

    def _params(self, extra: dict) -> dict:
        return self._auth_params({"retmode": "json", **extra})

    def _throttle(self) -> None:
        """Space successive Entrez requests at least ``self._delay`` apart.

        The sleep happens inside the lock so concurrent resolver threads
        queue up and each one waits for its slot — without the lock the
        shared ``_last_request`` timestamp races and the NCBI rate limit
        (3 req/s, or 10 with an API key) is not actually enforced.
        """
        with self._throttle_lock:
            now = time.time()
            wait = self._last_request + self._delay - now
            if wait > 0:
                time.sleep(wait)
                now = time.time()
            self._last_request = now

    def _get(self, endpoint: str, params: dict) -> dict:
        self._throttle()
        resp = requests.get(f"{_ENTREZ_BASE}/{endpoint}", params=params, timeout=30)
        resp.raise_for_status()
        return resp.json()

    def _get_text(self, endpoint: str, params: dict) -> str:
        """Like ``_get`` but returns the raw response body — for the Entrez
        endpoints that only speak XML (taxonomy ``efetch``)."""
        self._throttle()
        resp = requests.get(f"{_ENTREZ_BASE}/{endpoint}", params=params, timeout=30)
        resp.raise_for_status()
        return resp.text

    # ------------------------------------------------------------------
    # Taxonomy by taxid
    # ------------------------------------------------------------------

    def fetch_lineage(self, taxid: int) -> Optional[dict[str, Any]]:
        """Return the full taxonomic lineage for a taxid (cached).

        Uses Entrez ``efetch`` on the taxonomy DB (XML). The taxonomy
        *esummary* response — which a previous implementation parsed — does
        not carry the lineage at all: its ``genus``/``species`` fields are
        blank for viruses and there is no ``LineageEx``. That silently gave
        every viral sequence an empty lineage, so taxonomic modes grouped
        everything under "Unknown". ``efetch`` is the only Entrez endpoint
        that returns the ranked lineage.

        Returns a dict with the standard rank keys plus a ``lineage``
        rank→name map, or ``None`` if the taxid has no usable lineage.
        """
        key = str(taxid)
        cached = self._cache.get(_SOURCE, key)
        if cached is not None:
            return None if cached.get("_not_found") else cached

        try:
            xml_text = self._get_text(
                "efetch.fcgi",
                self._auth_params(
                    {"db": "taxonomy", "id": str(taxid), "retmode": "xml"}
                ),
            )
        except Exception:
            return None

        rank_map = _parse_taxonomy_xml(xml_text)
        if not rank_map:
            self._cache.set(_SOURCE, key, _NOT_FOUND)
            return None

        result = {
            "taxid": taxid,
            "species": rank_map.get("species"),
            "genus": rank_map.get("genus"),
            "family": rank_map.get("family"),
            "order": rank_map.get("order"),
            "class": rank_map.get("class"),
            "phylum": rank_map.get("phylum"),
            "kingdom": rank_map.get("kingdom"),
            "superkingdom": rank_map.get("superkingdom"),
            "lineage": rank_map,
        }
        self._cache.set(_SOURCE, key, result)
        return result

    # ------------------------------------------------------------------
    # Metadata from accession (nuccore / protein)
    # ------------------------------------------------------------------

    def fetch_accession_metadata(self, accession: str) -> Optional[dict[str, Any]]:
        """Fetch organism, taxid, host, collection_date, country for an accession."""
        cached = self._cache.get(_SOURCE_NUCCORE, accession)
        if cached is not None:
            return None if cached.get("_not_found") else cached

        db = "protein" if _looks_like_protein_acc(accession) else "nuccore"
        try:
            # First esearch to get uid
            search = self._get(
                "esearch.fcgi",
                self._params({"db": db, "term": accession}),
            )
            ids = search.get("esearchresult", {}).get("idlist", [])
            if not ids:
                self._cache.set(_SOURCE_NUCCORE, accession, _NOT_FOUND)
                return None

            # Then esummary
            summary = self._get(
                "esummary.fcgi",
                self._params({"db": db, "id": ids[0]}),
            )
            rec = summary.get("result", {}).get(ids[0], {})
            if not rec:
                self._cache.set(_SOURCE_NUCCORE, accession, _NOT_FOUND)
                return None

            taxid = rec.get("taxid")
            result: dict[str, Any] = {
                "accession": accession,
                "organism": rec.get("organism"),
                "taxid": int(taxid) if taxid else None,
                "title": rec.get("title"),
            }

            # Source-feature qualifiers (host, country, collection_date,
            # strain/isolate) are exposed by esummary as two parallel
            # pipe-delimited fields: 'subtype' holds the qualifier names,
            # 'subname' the corresponding values. Without this, NCBI
            # sequences get host/country/date only from fragile header
            # parsing, never from the authoritative database.
            subtype = rec.get("subtype") or ""
            subname = rec.get("subname") or ""
            if subtype and subname:
                quals = dict(zip(subtype.split("|"), subname.split("|")))
                if quals.get("host"):
                    result["host"] = quals["host"]
                # NCBI is migrating 'country' → 'geo_loc_name'; accept both.
                location = quals.get("geo_loc_name") or quals.get("country")
                if location:
                    result["country"] = location
                if quals.get("collection_date"):
                    result["collection_date"] = quals["collection_date"]
                strain = quals.get("strain") or quals.get("isolate")
                if strain:
                    result["strain"] = strain

            # Fetch lineage if we have a taxid
            if result["taxid"]:
                lineage = self.fetch_lineage(result["taxid"])
                if lineage:
                    result["lineage"] = lineage

            self._cache.set(_SOURCE_NUCCORE, accession, result)
            return result
        except Exception:
            return None

    # ------------------------------------------------------------------
    # Batched GenBank protein fetching
    # ------------------------------------------------------------------

    def fetch_proteins_batch(
        self,
        accessions: list[str],
        batch_size: int = _GENBANK_BATCH_SIZE,
        progress: Optional[Any] = None,
    ) -> dict[str, list[dict]]:
        """Fetch CDS protein annotations for nucleotide accessions.

        Returns a mapping ``accession → [{"protein_id", "product", "length"}, ...]``.
        Accessions with no record (or no CDS features) map to ``[]``.

        Network calls are batched (up to ``batch_size`` per ``efetch``).
        Per-accession results are cached under the ``ncbi_proteins`` source,
        so subsequent runs incur no network cost for the same inputs. The
        cache entry also carries source-feature qualifiers (isolate, strain,
        segment) extracted from the same GenBank record — read via
        ``fetch_source_metadata_batch`` to avoid a second network round trip.
        """
        records = self._fetch_genbank_batch(accessions, batch_size, progress=progress)
        return {acc: rec.get("proteins", []) for acc, rec in records.items()}

    def fetch_source_metadata_batch(
        self,
        accessions: list[str],
        batch_size: int = _GENBANK_BATCH_SIZE,
        progress: Optional[Any] = None,
    ) -> dict[str, dict[str, Optional[str]]]:
        """Fetch source-feature qualifiers (isolate, strain, segment).

        Returns ``accession → {"isolate": ..., "strain": ..., "segment": ...}``.
        Values are ``None`` when the qualifier is absent on the GenBank source
        feature. Accessions with no record map to all-None.

        Shares the ``ncbi_proteins`` cache with ``fetch_proteins_batch`` — a
        single efetch populates both the protein list and the source metadata,
        so running protein QC and segmented metadata extraction together costs
        one round trip. Cache entries written by older repseq versions did not
        store source metadata; for those, all fields come back ``None`` (the
        caller falls back to header parsing).
        """
        records = self._fetch_genbank_batch(accessions, batch_size, progress=progress)
        empty: dict[str, Optional[str]] = {
            "isolate": None, "strain": None, "segment": None,
        }
        return {acc: dict(rec.get("source") or empty) for acc, rec in records.items()}

    def _fetch_genbank_batch(
        self,
        accessions: list[str],
        batch_size: int,
        progress: Optional[Any] = None,
    ) -> dict[str, dict[str, Any]]:
        """Cached-batched GenBank fetch shared by protein and source-metadata APIs.

        Returns ``accession → {"proteins": [...], "source": {...}}``. The
        ``source`` value may be ``None`` for accessions cached by earlier
        repseq versions that did not capture source qualifiers — callers
        should treat that as "fall back to other means".

        When ``progress`` is callable it is invoked once per batch with the
        signature ``progress(done_batches, total_batches, batch_size_actual,
        cached_count)`` so the caller can emit a heartbeat for what is
        otherwise a many-minute silent network round-trip. The first call
        carries ``done_batches=0`` and reports the cache hit rate.
        """
        results: dict[str, dict[str, Any]] = {}
        to_fetch: list[str] = []
        for acc in accessions:
            cached = self._cache.get(_SOURCE_PROTEINS, acc)
            if cached is not None:
                results[acc] = {
                    "proteins": cached.get("proteins", []),
                    "source": cached.get("source"),
                }
            else:
                to_fetch.append(acc)

        total_batches = (len(to_fetch) + batch_size - 1) // batch_size
        cached_count = len(accessions) - len(to_fetch)
        if callable(progress):
            progress(0, total_batches, 0, cached_count)

        for i in range(0, len(to_fetch), batch_size):
            chunk = to_fetch[i : i + batch_size]
            try:
                fetched = self._fetch_genbank_chunk(chunk)
            except Exception as exc:
                logger.warning(
                    "GenBank batch fetch failed for %d accessions: %s",
                    len(chunk), exc,
                )
                fetched = {}

            for acc in chunk:
                rec = fetched.get(acc, {"proteins": [], "source": None})
                self._cache.set(_SOURCE_PROTEINS, acc, rec)
                results[acc] = rec

            if callable(progress):
                progress(i // batch_size + 1, total_batches, len(chunk), cached_count)

        return results

    def _fetch_genbank_chunk(
        self,
        accessions: list[str],
    ) -> dict[str, dict[str, Any]]:
        """One efetch call → parsed CDS features and source qualifiers per accession.

        Records from NCBI carry version suffixes (e.g. ``MW626064.1``); we
        match each parsed record back to the caller's accession with and
        without the version.
        """
        from Bio import SeqIO

        self._throttle()

        params: dict[str, Any] = {
            "db": "nuccore",
            "id": ",".join(accessions),
            "rettype": "gb",
            "retmode": "text",
        }
        if self._email:
            params["email"] = self._email
        if self._api_key:
            params["api_key"] = self._api_key

        resp = requests.get(
            f"{_ENTREZ_BASE}/efetch.fcgi", params=params, timeout=120,
        )
        resp.raise_for_status()

        by_acc_full: dict[str, dict[str, Any]] = {}
        by_acc_no_version: dict[str, dict[str, Any]] = {}
        for record in SeqIO.parse(StringIO(resp.text), "genbank"):
            proteins: list[dict] = []
            source: dict[str, Optional[str]] = {
                "isolate": None, "strain": None, "segment": None,
            }
            for feat in record.features:
                if feat.type == "source":
                    q = feat.qualifiers
                    source["isolate"] = q.get("isolate", [None])[0]
                    source["strain"] = q.get("strain", [None])[0]
                    source["segment"] = q.get("segment", [None])[0]
                elif feat.type == "CDS":
                    q = feat.qualifiers
                    translation = q.get("translation", [None])[0]
                    proteins.append({
                        "protein_id": q.get("protein_id", [None])[0],
                        "product": q.get("product", [None])[0],
                        "length": len(translation) if translation else None,
                        # Amino-acid sequence from the GenBank /translation=
                        # qualifier. Stored so we can later emit a proteins.fasta
                        # without a second network call.
                        "sequence": translation,
                    })
            rec = {"proteins": proteins, "source": source}
            by_acc_full[record.id] = rec
            by_acc_no_version[record.id.split(".")[0]] = rec

        results: dict[str, dict[str, Any]] = {}
        for acc in accessions:
            if acc in by_acc_full:
                results[acc] = by_acc_full[acc]
            else:
                results[acc] = by_acc_no_version.get(
                    acc.split(".")[0],
                    {"proteins": [], "source": None},
                )
        return results


def _parse_taxonomy_xml(xml_text: str) -> dict[str, str]:
    """Parse an Entrez taxonomy ``efetch`` XML payload into a rank→name map.

    Collects every ranked entry in the queried taxon's ``<LineageEx>`` plus
    the queried taxon's own rank. Entries with rank ``no rank`` or ``clade``
    are skipped — only named ranks (genus, family, order, …) are kept.
    """
    import xml.etree.ElementTree as ET

    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return {}

    taxon = root.find("Taxon")
    if taxon is None:
        return {}

    rank_map: dict[str, str] = {}

    def _add(node) -> None:
        rank = (node.findtext("Rank") or "").strip().lower()
        name = (node.findtext("ScientificName") or "").strip()
        if rank and name and rank not in ("no rank", "clade"):
            rank_map.setdefault(rank, name)

    lineage_ex = taxon.find("LineageEx")
    if lineage_ex is not None:
        for node in lineage_ex.findall("Taxon"):
            _add(node)
    # The queried taxon itself may carry a useful rank (e.g. it *is* the
    # species); add it last so a LineageEx entry of the same rank wins.
    _add(taxon)
    return rank_map


def _looks_like_protein_acc(acc: str) -> bool:
    """Return True if the accession looks like a protein record.

    RefSeq protein accessions are NP_/XP_/YP_/WP_/AP_ — the discriminator
    is a 'P' in the second position followed by '_'. RefSeq *nucleotide*
    accessions (NC_, NM_, NG_, NR_, NT_, NW_, NZ_, XM_, XR_) must NOT match,
    or their metadata would be looked up in the wrong Entrez database.
    GenBank protein accessions are three letters followed by digits.
    """
    import re
    return bool(re.match(r"^[A-Z]P_\d|^[A-Z]{3}\d", acc))
