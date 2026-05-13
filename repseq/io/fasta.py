"""FASTA parsing with header field extraction for UniProt, NCBI, and NCBI Virus."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Iterator, Optional

from ..models import Sequence, SequenceSource, SequenceType


# ---------------------------------------------------------------------------
# Sequence type detection
# ---------------------------------------------------------------------------

_NUCLEOTIDE_CHARS = set("ATCGUNRYWSKMBDHVatcgunrywskmbdhv")
_PROTEIN_ONLY_CHARS = set("EFILPQZefilpqz")


def detect_seq_type(sequence: str) -> SequenceType:
    if not sequence:
        return SequenceType.UNKNOWN
    sample = sequence[:200].upper()
    if any(c in _PROTEIN_ONLY_CHARS for c in sample):
        return SequenceType.PROTEIN
    na_count = sum(1 for c in sample if c in _NUCLEOTIDE_CHARS)
    if na_count / len(sample) >= 0.90:
        return SequenceType.NUCLEOTIDE
    return SequenceType.PROTEIN


# ---------------------------------------------------------------------------
# Header parsers
# ---------------------------------------------------------------------------

# UniProt: >sp|P12345|PROT_HUMAN Description OS=Homo sapiens OX=9606 GN=gene PE=1 SV=1
_UNIPROT_RE = re.compile(
    r"^(?P<db>sp|tr)\|(?P<accession>[A-Z0-9]+)\|(?P<entry>\S+)"
    r"\s+(?P<description>.*?)"
    r"(?:\s+OS=(?P<organism>[^=]+?)(?=\s+\w+=|$))?"
    r"(?:\s+OX=(?P<taxid>\d+))?"
    r"(?:\s+GN=(?P<gene>[^=]+?)(?=\s+\w+=|$))?"
    r"(?:\s+PE=\d+)?(?:\s+SV=\d+)?$"
)

# NCBI standard: >accession.version description [organism]
_NCBI_RE = re.compile(
    r"^(?P<accession>[A-Z_]+[\d_]+(?:\.\d+)?)\s+(?P<description>.*?)"
    r"(?:\s+\[(?P<organism>[^\]]+)\])?$"
)

# NCBI Virus: >accession.version strain [host] [organism] [country] [collection_date] ...
# Example: >MW626064.1 Influenza A virus (A/duck/Philippines/14-0056/2014(H5N1)) segment 6
_NCBI_VIRUS_RE = re.compile(
    r"^(?P<accession>[A-Z_]+[\d_]+(?:\.\d+)?)\s+(?P<description>.*)"
)

# Influenza strain nomenclature: A/host/location/number/year(subtype)
_INFLUENZA_STRAIN_RE = re.compile(
    r"\((?P<strain>[AB]/[^/()]+/[^()]+)\)"
    r"|(?P<strain2>[AB]/[^/\s(]+/[^/\s(]+/[^/\s(]+/[^/\s(]+)"
)

# Collection date in header: [YYYY] or [YYYY-MM] or [YYYY-MM-DD]
_DATE_RE = re.compile(r"\[(\d{4}(?:-\d{2}(?:-\d{2})?)?)\]")

# Country in header: [country_name]
_COUNTRY_RE = re.compile(r"\[([A-Z][a-zA-Z\s]+(?::[A-Za-z\s]+)?)\]")

# Host in header
_HOST_RE = re.compile(r"\[(?:host=)?([^\]]+)\]", re.IGNORECASE)

# Segment in header
_SEGMENT_RE = re.compile(r"\bsegment\s+(\d+|\w+)\b", re.IGNORECASE)

# RefSeq accession prefixes
_REFSEQ_PREFIXES = (
    "NC_", "NM_", "NR_", "NP_", "XM_", "XR_", "XP_", "NG_", "NT_", "NW_", "NZ_",
)


def _is_refseq(accession: Optional[str]) -> bool:
    if not accession:
        return False
    return any(accession.startswith(p) for p in _REFSEQ_PREFIXES)


def _parse_uniprot_header(header: str) -> dict:
    m = _UNIPROT_RE.match(header)
    if not m:
        return {}
    db = m.group("db")
    return {
        "accession": m.group("accession"),
        "description": (m.group("description") or "").strip(),
        "organism": (m.group("organism") or "").strip() or None,
        "taxid": int(m.group("taxid")) if m.group("taxid") else None,
        "is_reviewed": db == "sp",
        "source": SequenceSource.UNIPROT,
    }


def _parse_ncbi_virus_header(header: str) -> dict:
    """Parse NCBI Virus FASTA headers, extracting all embedded metadata."""
    m = _NCBI_VIRUS_RE.match(header)
    if not m:
        return {}

    accession = m.group("accession")
    description = m.group("description")

    result: dict = {
        "accession": accession,
        "description": description,
        "is_refseq": _is_refseq(accession),
        "source": SequenceSource.NCBI_VIRUS,
    }

    # Organism: last [...] block that looks like a species name
    brackets = re.findall(r"\[([^\]]+)\]", description)
    for b in reversed(brackets):
        if re.match(r"[A-Z][a-z]", b) and len(b.split()) >= 2:
            result["organism"] = b
            break

    # Strain from influenza-style nomenclature
    sm = _INFLUENZA_STRAIN_RE.search(description)
    if sm:
        result["strain"] = sm.group("strain") or sm.group("strain2")

    # Collection date
    dm = _DATE_RE.search(description)
    if dm:
        result["collection_date"] = dm.group(1)

    # Segment
    segm = _SEGMENT_RE.search(description)
    if segm:
        result["segment"] = segm.group(1)

    # Country — look for [Country] or [Country: Region] not already captured
    remaining = description
    for b in brackets:
        if b == result.get("organism"):
            continue
        if re.match(r"[A-Z][a-zA-Z\s]+(?::[A-Za-z\s]+)?$", b):
            if not re.match(r"\d{4}", b):
                result.setdefault("country", b)
                break

    # Host — look for [host] bracket early in the description
    host_candidates = [b for b in brackets if b not in (result.get("organism"), result.get("country"))]
    if host_candidates:
        result.setdefault("host", host_candidates[0])

    return result


def _parse_ncbi_header(header: str) -> dict:
    m = _NCBI_RE.match(header)
    if not m:
        return {}
    accession = m.group("accession")
    return {
        "accession": accession,
        "description": (m.group("description") or "").strip(),
        "organism": (m.group("organism") or "").strip() or None,
        "is_refseq": _is_refseq(accession),
        "source": SequenceSource.NCBI,
    }


def parse_header(header: str) -> dict:
    """Auto-detect header format and return extracted fields."""
    h = header.lstrip(">").strip()

    # UniProt
    if h.startswith("sp|") or h.startswith("tr|"):
        result = _parse_uniprot_header(h)
        if result:
            return result

    # NCBI Virus — detect by presence of virus-like bracket patterns or accession pattern
    accession_match = re.match(r"^([A-Z_]+[\d_]+(?:\.\d+)?)\s", h)
    if accession_match:
        acc = accession_match.group(1)
        if _is_refseq(acc) or re.search(r"\[[^\]]+\].*\[[^\]]+\]", h):
            result = _parse_ncbi_virus_header(h)
            if result:
                return result
        result = _parse_ncbi_header(h)
        if result:
            return result

    # Fallback: use full header as description
    return {
        "accession": h.split()[0] if h.split() else None,
        "description": h,
        "source": SequenceSource.UNKNOWN,
    }


# ---------------------------------------------------------------------------
# FASTA reader
# ---------------------------------------------------------------------------

def read_fasta(
    path: str | Path,
    source_override: Optional[SequenceSource] = None,
) -> Iterator[Sequence]:
    """Yield Sequence objects from a FASTA file.

    Args:
        path: Path to the FASTA file.
        source_override: If set, force this source on every sequence instead
                         of auto-detecting from the header.
    """
    path = Path(path)
    current_header: Optional[str] = None
    current_lines: list[str] = []

    def _emit(header: str, lines: list[str]) -> Sequence:
        seq_str = "".join(lines).upper().replace(" ", "").replace("\r", "")
        fields = parse_header(header)
        seq_type = detect_seq_type(seq_str)
        seq = Sequence(
            id=fields.get("accession") or header.split()[0].lstrip(">"),
            header=header.lstrip(">").strip(),
            sequence=seq_str,
            seq_type=seq_type,
            source=source_override if source_override is not None else fields.get("source", SequenceSource.UNKNOWN),
            accession=fields.get("accession"),
            organism=fields.get("organism"),
            description=fields.get("description"),
            strain=fields.get("strain"),
            host=fields.get("host"),
            collection_date=fields.get("collection_date"),
            country=fields.get("country"),
            segment=fields.get("segment"),
            is_refseq=fields.get("is_refseq", False),
            is_reviewed=fields.get("is_reviewed", False),
        )
        if fields.get("taxid"):
            from ..models import TaxonomyInfo
            seq.taxonomy = TaxonomyInfo(taxid=fields["taxid"])
        return seq

    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_header is not None:
                    yield _emit(current_header, current_lines)
                current_header = line[1:].strip()
                current_lines = []
            else:
                current_lines.append(line.strip())

    if current_header is not None:
        yield _emit(current_header, current_lines)


def write_fasta(sequences: list[Sequence], path: str | Path, line_width: int = 70) -> None:
    """Write a list of Sequence objects to a FASTA file."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        for seq in sequences:
            fh.write(f">{seq.header}\n")
            for i in range(0, len(seq.sequence), line_width):
                fh.write(seq.sequence[i : i + line_width] + "\n")
