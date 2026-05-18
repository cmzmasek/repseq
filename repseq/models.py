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
    # Provenance of ``isolate_id``: ``"isolate"`` if the value came from
    # the GenBank source feature's ``/isolate`` qualifier,
    # ``"strain"`` if the populator fell back to ``/strain`` (no
    # ``/isolate`` qualifier was present), ``"regex"`` if the
    # header-regex fallback fired (``--no-resolve`` / UniProt input /
    # missing accession). ``None`` when ``isolate_id`` is unset. Strain
    # provenance is the over-merge risk: a single named strain is often
    # shared across distinct biological samples, so two accessions can
    # legitimately share ``isolate_id`` while being different isolates —
    # see the strain-collision detector for the active check.
    isolate_id_source: Optional[str] = None

    # Quality flags
    is_refseq: bool = False
    is_reviewed: bool = False  # UniProt Swiss-Prot reviewed

    # Taxonomy (resolved after DB lookup)
    taxonomy: Optional[TaxonomyInfo] = None

    # Protein annotations (None = not fetched, [] = fetched but none found).
    # Each dict has keys: protein_id, product, length, sequence.
    proteins: Optional[list[dict]] = None

    # Amino-acid string fed to the clustering backend when
    # clustering.alphabet=protein. For non-segmented input this is the
    # chosen marker protein (longest CDS, or first cluster_protein alias
    # that matches). For segmented input the concat sequence holds the
    # in-segment-order concatenation of each segment's marker protein.
    # Always None when alphabet=nucleotide.
    protein_sequence: Optional[str] = None

    # ``protein_id`` of each marker CDS used for clustering / phylogeny.
    # Non-segmented: a single-element list (the one marker chosen on
    # this sequence). Segmented CONCAT: one per segment, in segment
    # order — so a phyloXML writer can mark those proteins as
    # "this is what fed the tree" and list them first.
    marker_protein_ids: Optional[list[str]] = None

    # For segmented CONCAT records: the per-segment Sequence objects
    # that were concatenated, in segment order. Each carries its own
    # accession + .proteins list, so downstream output (phyloXML
    # multi-<sequence> emission, isolate_proteins.tsv) can list every
    # underlying nuc accession and protein without re-fetching. ``None``
    # for non-segmented input — in that case the leaf IS its own
    # single "segment", and the writer reads from ``seq`` directly.
    concat_segments: Optional[list["Sequence"]] = None

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
            # X=any, B=Asx, Z=Glx, J=Xle. U (selenocysteine) and O
            # (pyrrolysine) are definite residues, not ambiguity codes.
            ambiguous = set("XBZJ")
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
class GroupStat:
    """Per-group before/after counts for one stratified-selection group.

    Recorded by every grouping mode (taxonomic1/2, host, time, geographic,
    custom, hybrid) and by global mode. ``n_before`` is the number of
    sequences entering the group — in segmented runs these are the
    concatenated per-isolate sequences, since segmentation happens before
    the mode runs. ``cutoff`` is the MMseqs2 identity threshold the binary
    search settled on, or ``None`` when the group was small enough to keep
    without clustering.
    """
    grouping: str            # the grouping dimension (rank name, "host", field, …)
    group: str               # the group label/value within that dimension
    n_before: int            # sequences entering the group
    n_after: int             # representatives selected from the group
    clustered: bool          # whether MMseqs2 clustering ran (vs. kept whole)
    cutoff: Optional[float] = None   # MMseqs2 identity threshold used, if clustered


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
    # Segments dropped because their isolate's segments disagreed on
    # the configured taxonomy rank (segmented.taxonomy_consistency).
    # Counted in *segments removed*, not isolates, so the units stay
    # consistent with the other ``removed_*`` counters.
    removed_taxonomy_mismatch: int = 0
    # Segments dropped by the strain-collision detector — accessions
    # sharing a strain-derived ``isolate_id`` AND a segment, which is
    # the over-merge signature of the /strain → isolate_id fallback.
    # Always zero when ``segmented.strain_collision_action`` is
    # ``"warn"`` (the default); only increments under ``"drop"``.
    # Counted in segments, consistent with the other ``removed_*``
    # counters.
    removed_strain_collisions: int = 0
    # Whole-pool length filter is skipped in segmented mode (the input is a
    # mix of segments of very different lengths; per-segment bounds are
    # applied later via segmented.viruses.<v>.segment_lengths instead).
    length_filter_skipped: bool = False
    # Per-segment isolate drops from the segmented length filter. Shape:
    # {"L": {"too_short": 257, "too_long": 0}, "M": {...}, ...}. Counted
    # in *isolates*, not segments — the filter is isolate-level (one bad
    # segment drops the whole isolate), so "257 isolates lost their L"
    # is the actionable number. This is the only ``removed_*`` field in
    # isolate units; all the others count segments.
    removed_length_by_segment: dict = field(default_factory=dict)
    # Exact-duplicate removal is likewise skipped on the segment pool in
    # segmented mode: a segment can be byte-identical between two otherwise
    # distinct isolates, and dropping it would leave one isolate incomplete.
    # Dedup is instead applied to the concatenated per-isolate sequences
    # after build_concatenated_sequences.
    dedup_skipped: bool = False
    details: list[dict] = field(default_factory=list)

    def add_removed(self, seq_id: str, reason: str) -> None:
        self.details.append({"id": seq_id, "reason": reason})

    def summary(self) -> str:
        length_lines: list[str]
        if self.removed_length_by_segment:
            total = sum(
                c["too_short"] + c["too_long"]
                for c in self.removed_length_by_segment.values()
            )
            length_lines = [
                f"  Removed (length)    : {total} isolate(s) "
                "(per-segment, isolate-level)"
            ]
            for seg_name in sorted(self.removed_length_by_segment.keys()):
                counts = self.removed_length_by_segment[seg_name]
                if counts["too_short"]:
                    length_lines.append(
                        f"    {seg_name} too short  : {counts['too_short']}"
                    )
                if counts["too_long"]:
                    length_lines.append(
                        f"    {seg_name} too long   : {counts['too_long']}"
                    )
        elif self.length_filter_skipped:
            length_lines = ["  Removed (length)    : skipped (segmented mode)"]
        else:
            length_lines = [f"  Removed (length)    : {self.removed_length}"]
        dup_line = (
            f"  Removed (duplicates): {self.removed_duplicates} "
            "(applied to concatenated isolates)"
            if self.dedup_skipped
            else f"  Removed (duplicates): {self.removed_duplicates}"
        )
        lines = [
            f"QC Summary",
            f"  Input sequences     : {self.total_input}",
            f"  Passed QC           : {self.passed}",
            dup_line,
            *length_lines,
            f"  Removed (ambiguous) : {self.removed_ambiguous}",
            f"  Removed (annotation): {self.removed_annotation}",
            f"  Removed (proteins)  : {self.removed_proteins}",
            f"  Removed (incomplete): {self.removed_incomplete_isolates}",
            f"  Removed (tax-mismatch): {self.removed_taxonomy_mismatch}",
            f"  Removed (strain-collision): {self.removed_strain_collisions}",
        ]
        return "\n".join(lines)


@dataclass
class RunResult:
    mode: str
    representatives: list[Sequence] = field(default_factory=list)
    clusters: list[Cluster] = field(default_factory=list)
    group_stats: list[GroupStat] = field(default_factory=list)
    qc_report: Optional[QCReport] = None
    config_snapshot: dict = field(default_factory=dict)
