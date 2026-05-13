"""Core data models for repseq."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional


class SequenceType(Enum):
    PROTEIN = "protein"
    NUCLEOTIDE = "nucleotide"
    UNKNOWN = "unknown"


class SequenceSource(Enum):
    UNIPROT = "uniprot"
    NCBI = "ncbi"
    NCBI_VIRUS = "ncbi_virus"
    UNKNOWN = "unknown"


@dataclass
class TaxonomyInfo:
    taxid: Optional[int] = None
    species: Optional[str] = None
    genus: Optional[str] = None
    family: Optional[str] = None
    order: Optional[str] = None
    class_: Optional[str] = None
    phylum: Optional[str] = None
    kingdom: Optional[str] = None
    superkingdom: Optional[str] = None
    # Full rank → name mapping for arbitrary rank lookups
    lineage: dict[str, str] = field(default_factory=dict)

    def get_rank(self, rank: str) -> Optional[str]:
        """Return taxon name at the given rank, checking both standard fields and lineage."""
        rank_lower = rank.lower()
        standard = {
            "species": self.species,
            "genus": self.genus,
            "family": self.family,
            "order": self.order,
            "class": self.class_,
            "phylum": self.phylum,
            "kingdom": self.kingdom,
            "superkingdom": self.superkingdom,
        }
        if rank_lower in standard:
            return standard[rank_lower]
        return self.lineage.get(rank_lower)


@dataclass
class Sequence:
    id: str
    header: str
    sequence: str
    seq_type: SequenceType = SequenceType.UNKNOWN
    source: SequenceSource = SequenceSource.UNKNOWN

    # Parsed metadata
    accession: Optional[str] = None
    organism: Optional[str] = None
    description: Optional[str] = None
    strain: Optional[str] = None
    host: Optional[str] = None
    collection_date: Optional[str] = None
    country: Optional[str] = None
    segment: Optional[str] = None
    isolate_id: Optional[str] = None

    # Quality flags
    is_refseq: bool = False
    is_reviewed: bool = False  # UniProt Swiss-Prot reviewed

    # Taxonomy (resolved after DB lookup)
    taxonomy: Optional[TaxonomyInfo] = None

    # Protein annotations (None = not fetched, [] = fetched but none found).
    # Each dict has keys: protein_id, product, length.
    proteins: Optional[list[dict]] = None

    # QC state
    qc_passed: bool = True
    qc_fail_reason: Optional[str] = None

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def ambiguous_fraction(self) -> float:
        if not self.sequence:
            return 0.0
        seq = self.sequence.upper()
        if self.seq_type == SequenceType.PROTEIN:
            ambiguous = set("XBZJU")
        else:
            ambiguous = set("NRYWSKMBDHV")
        return sum(1 for c in seq if c in ambiguous) / len(seq)

    def to_fasta(self) -> str:
        return f">{self.header}\n{self.sequence}\n"


@dataclass
class Cluster:
    cluster_id: str
    representative: Sequence
    members: list[Sequence] = field(default_factory=list)
    identity_threshold: Optional[float] = None

    @property
    def size(self) -> int:
        return len(self.members) + 1  # +1 for representative


@dataclass
class QCReport:
    total_input: int = 0
    passed: int = 0
    removed_duplicates: int = 0
    removed_length: int = 0
    removed_ambiguous: int = 0
    removed_annotation: int = 0
    removed_proteins: int = 0
    removed_incomplete_isolates: int = 0
    details: list[dict] = field(default_factory=list)

    def add_removed(self, seq_id: str, reason: str) -> None:
        self.details.append({"id": seq_id, "reason": reason})

    def summary(self) -> str:
        lines = [
            f"QC Summary",
            f"  Input sequences     : {self.total_input}",
            f"  Passed QC           : {self.passed}",
            f"  Removed (duplicates): {self.removed_duplicates}",
            f"  Removed (length)    : {self.removed_length}",
            f"  Removed (ambiguous) : {self.removed_ambiguous}",
            f"  Removed (annotation): {self.removed_annotation}",
            f"  Removed (proteins)  : {self.removed_proteins}",
            f"  Removed (incomplete): {self.removed_incomplete_isolates}",
        ]
        return "\n".join(lines)


@dataclass
class RunResult:
    mode: str
    representatives: list[Sequence] = field(default_factory=list)
    clusters: list[Cluster] = field(default_factory=list)
    qc_report: Optional[QCReport] = None
    config_snapshot: dict = field(default_factory=dict)
