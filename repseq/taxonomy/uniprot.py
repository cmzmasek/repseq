"""UniProt REST API queries for taxonomy and metadata."""

from __future__ import annotations

import time
from typing import Any, Optional

import requests

from .cache import TaxonomyCache

_UNIPROT_API = "https://rest.uniprot.org/uniprotkb"
_TAXONOMY_API = "https://rest.uniprot.org/taxonomy"
_RATE_LIMIT_DELAY = 0.2
_SOURCE = "uniprot"
_SOURCE_TAX = "uniprot_taxonomy"


class UniProtAPI:
    def __init__(self, cache: TaxonomyCache) -> None:
        self._cache = cache
        self._last_request: float = 0.0

    def _get(self, url: str, params: Optional[dict] = None) -> dict:
        elapsed = time.time() - self._last_request
        if elapsed < _RATE_LIMIT_DELAY:
            time.sleep(_RATE_LIMIT_DELAY - elapsed)
        resp = requests.get(url, params=params, timeout=30)
        self._last_request = time.time()
        resp.raise_for_status()
        return resp.json()

    # ------------------------------------------------------------------
    # Entry metadata by accession
    # ------------------------------------------------------------------

    def fetch_entry(self, accession: str) -> Optional[dict[str, Any]]:
        """Fetch organism, taxonomy, host, and review status for a UniProt accession."""
        cached = self._cache.get(_SOURCE, accession)
        if cached is not None:
            return cached

        try:
            data = self._get(
                f"{_UNIPROT_API}/{accession}",
                params={"format": "json"},
            )
        except Exception:
            return None

        result = _parse_uniprot_entry(accession, data)
        if result:
            self._cache.set(_SOURCE, accession, result)
        return result

    # ------------------------------------------------------------------
    # Taxonomy by taxid
    # ------------------------------------------------------------------

    def fetch_lineage(self, taxid: int) -> Optional[dict[str, Any]]:
        key = str(taxid)
        cached = self._cache.get(_SOURCE_TAX, key)
        if cached is not None:
            return cached

        try:
            data = self._get(f"{_TAXONOMY_API}/{taxid}")
        except Exception:
            return None

        lineage_entries = data.get("lineage", [])
        rank_map: dict[str, str] = {}
        for entry in lineage_entries:
            rank = entry.get("rank", "").lower()
            name = entry.get("scientificName", "")
            if rank and name and rank != "no rank":
                rank_map[rank] = name

        result = {
            "taxid": taxid,
            "species": rank_map.get("species") or data.get("scientificName"),
            "genus": rank_map.get("genus"),
            "family": rank_map.get("family"),
            "order": rank_map.get("order"),
            "class": rank_map.get("class"),
            "phylum": rank_map.get("phylum"),
            "kingdom": rank_map.get("kingdom"),
            "superkingdom": rank_map.get("superkingdom"),
            "lineage": rank_map,
        }
        self._cache.set(_SOURCE_TAX, key, result)
        return result


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------

def _parse_uniprot_entry(accession: str, data: dict) -> dict[str, Any]:
    result: dict[str, Any] = {"accession": accession}

    # Review status
    result["is_reviewed"] = data.get("entryType") == "UniProtKB reviewed (Swiss-Prot)"

    # Organism
    organism = data.get("organism", {})
    result["organism"] = organism.get("scientificName")
    result["taxid"] = organism.get("taxonId")

    # Lineage from organism.lineage
    lineage_names = organism.get("lineage", [])
    result["lineage_names"] = lineage_names

    # Host
    hosts = data.get("organismHosts", [])
    if hosts:
        result["host"] = hosts[0].get("scientificName")

    # Description
    protein = data.get("proteinDescription", {})
    recommended = protein.get("recommendedName", {})
    full_name = recommended.get("fullName", {})
    result["description"] = full_name.get("value")

    return result
