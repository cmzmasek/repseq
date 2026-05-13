"""NCBI Entrez taxonomy and metadata queries."""

from __future__ import annotations

import time
from typing import Any, Optional

import requests

from .cache import TaxonomyCache

_ENTREZ_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
_RATE_LIMIT_DELAY = 0.34  # ~3 req/s without API key; 0.1 with key
_SOURCE = "ncbi_taxonomy"
_SOURCE_NUCCORE = "ncbi_nuccore"


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

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _params(self, extra: dict) -> dict:
        p: dict = {"retmode": "json", **extra}
        if self._email:
            p["email"] = self._email
        if self._api_key:
            p["api_key"] = self._api_key
        return p

    def _get(self, endpoint: str, params: dict) -> dict:
        elapsed = time.time() - self._last_request
        if elapsed < self._delay:
            time.sleep(self._delay - elapsed)
        resp = requests.get(f"{_ENTREZ_BASE}/{endpoint}", params=params, timeout=30)
        self._last_request = time.time()
        resp.raise_for_status()
        return resp.json()

    # ------------------------------------------------------------------
    # Taxonomy by taxid
    # ------------------------------------------------------------------

    def fetch_lineage(self, taxid: int) -> Optional[dict[str, Any]]:
        """Return lineage dict for a taxid, using cache."""
        key = str(taxid)
        cached = self._cache.get(_SOURCE, key)
        if cached is not None:
            return cached

        try:
            data = self._get(
                "efetch.fcgi",
                self._params({"db": "taxonomy", "id": str(taxid), "retmode": "xml"}),
            )
        except Exception:
            return None

        # efetch with retmode=json for taxonomy doesn't work; use XML via requests
        # Fallback: use esummary
        return self._fetch_lineage_esummary(taxid)

    def _fetch_lineage_esummary(self, taxid: int) -> Optional[dict[str, Any]]:
        key = str(taxid)
        try:
            data = self._get(
                "esummary.fcgi",
                self._params({"db": "taxonomy", "id": str(taxid)}),
            )
            result_dict = data.get("result", {})
            rec = result_dict.get(str(taxid), {})
            if not rec:
                return None

            lineage_str = rec.get("lineage", "")
            lineage_names = [x.strip() for x in lineage_str.split(";") if x.strip()]

            lineage_ex = rec.get("lineageex", [])
            rank_map: dict[str, str] = {}
            for entry in lineage_ex:
                rank = entry.get("rank", "").lower()
                name = entry.get("scientificname", "")
                if rank and name and rank != "no rank":
                    rank_map[rank] = name

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
        except Exception:
            return None

    # ------------------------------------------------------------------
    # Metadata from accession (nuccore / protein)
    # ------------------------------------------------------------------

    def fetch_accession_metadata(self, accession: str) -> Optional[dict[str, Any]]:
        """Fetch organism, taxid, host, collection_date, country for an accession."""
        cached = self._cache.get(_SOURCE_NUCCORE, accession)
        if cached is not None:
            return cached

        db = "protein" if _looks_like_protein_acc(accession) else "nuccore"
        try:
            # First esearch to get uid
            search = self._get(
                "esearch.fcgi",
                self._params({"db": db, "term": accession}),
            )
            ids = search.get("esearchresult", {}).get("idlist", [])
            if not ids:
                return None

            # Then esummary
            summary = self._get(
                "esummary.fcgi",
                self._params({"db": db, "id": ids[0]}),
            )
            rec = summary.get("result", {}).get(ids[0], {})
            if not rec:
                return None

            taxid = rec.get("taxid")
            result: dict[str, Any] = {
                "accession": accession,
                "organism": rec.get("organism"),
                "taxid": int(taxid) if taxid else None,
                "title": rec.get("title"),
            }

            # Fetch lineage if we have a taxid
            if result["taxid"]:
                lineage = self._fetch_lineage_esummary(result["taxid"])
                if lineage:
                    result["lineage"] = lineage

            self._cache.set(_SOURCE_NUCCORE, accession, result)
            return result
        except Exception:
            return None


def _looks_like_protein_acc(acc: str) -> bool:
    # Protein accessions: [A-Z]{2}_\d+ or [A-Z]{3}\d+
    import re
    return bool(re.match(r"^[A-NR-Z][A-Z]_\d|^[A-Z]{3}\d", acc))
