"""NCBI Entrez taxonomy and metadata queries."""

from __future__ import annotations

import logging
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

    # ------------------------------------------------------------------
    # Batched GenBank protein fetching
    # ------------------------------------------------------------------

    def fetch_proteins_batch(
        self,
        accessions: list[str],
        batch_size: int = _GENBANK_BATCH_SIZE,
    ) -> dict[str, list[dict]]:
        """Fetch CDS protein annotations for nucleotide accessions.

        Returns a mapping ``accession → [{"protein_id", "product", "length"}, ...]``.
        Accessions with no record (or no CDS features) map to ``[]``.

        Network calls are batched (up to ``batch_size`` per ``efetch``).
        Per-accession results are cached under the ``ncbi_proteins`` source,
        so subsequent runs incur no network cost for the same inputs.
        """
        results: dict[str, list[dict]] = {}
        to_fetch: list[str] = []
        for acc in accessions:
            cached = self._cache.get(_SOURCE_PROTEINS, acc)
            if cached is not None:
                results[acc] = cached.get("proteins", [])
            else:
                to_fetch.append(acc)

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
                proteins = fetched.get(acc, [])
                self._cache.set(_SOURCE_PROTEINS, acc, {"proteins": proteins})
                results[acc] = proteins

        return results

    def _fetch_genbank_chunk(self, accessions: list[str]) -> dict[str, list[dict]]:
        """One efetch call → parsed CDS features per accession.

        Records from NCBI carry version suffixes (e.g. ``MW626064.1``); we
        match each parsed record back to the caller's accession with and
        without the version.
        """
        from Bio import SeqIO

        elapsed = time.time() - self._last_request
        if elapsed < self._delay:
            time.sleep(self._delay - elapsed)

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
        self._last_request = time.time()
        resp.raise_for_status()

        by_acc_full: dict[str, list[dict]] = {}
        by_acc_no_version: dict[str, list[dict]] = {}
        for record in SeqIO.parse(StringIO(resp.text), "genbank"):
            proteins: list[dict] = []
            for feat in record.features:
                if feat.type != "CDS":
                    continue
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
            by_acc_full[record.id] = proteins
            by_acc_no_version[record.id.split(".")[0]] = proteins

        results: dict[str, list[dict]] = {}
        for acc in accessions:
            if acc in by_acc_full:
                results[acc] = by_acc_full[acc]
            else:
                results[acc] = by_acc_no_version.get(acc.split(".")[0], [])
        return results


def _looks_like_protein_acc(acc: str) -> bool:
    # Protein accessions: [A-Z]{2}_\d+ or [A-Z]{3}\d+
    import re
    return bool(re.match(r"^[A-NR-Z][A-Z]_\d|^[A-Z]{3}\d", acc))
