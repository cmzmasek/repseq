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
    # Viral subtype / serotype from the GenBank ``/serotype`` source-feature
    # qualifier — e.g. influenza A "H5N1". Sample metadata (like host /
    # country), NOT an NCBI taxonomy rank, so it never enters the 9-rank
    # ``_TAX_RANKS`` ladder. Blank when the qualifier is absent.
    subtype: Optional[str] = None
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
    # Diversity curve: cluster count at each configured "standard" identity
    # threshold (e.g. {0.99: 87, 0.95: 42, 0.90: 18, 0.80: 5, 0.70: None}).
    # ``None`` for a cutoff means it is below the active backend's identity
    # floor (cd-hit-est refuses <0.80) and was not run. Populated only when
    # ``clustered=True`` and ``clustering.diversity_curve_cutoffs`` is
    # non-empty; reporting only — does not influence selection.
    cutoff_counts: Optional[dict[float, Optional[int]]] = None


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
    # Isolates dropped because their seg_map carried segment names outside
    # the configured ``segments`` list (e.g. an extra fourth segment, or
    # a non-canonical identifier that ``identify_segment`` returned
    # unchanged). Always zero when ``segmented.extra_segments_action`` is
    # ``"warn"`` (the default); only increments under ``"drop"``. Counted
    # in **isolates** (the actionable number for a bench scientist), not
    # segments — every segment of a dropped isolate still lands in
    # ``_qc_removed.tsv`` with reason ``extra_segments:<extras>``.
    removed_extra_segments: int = 0
    # Sequences (non-segmented) or isolates (segmented) dropped because a
    # CDS protein translation carried too high a fraction of ambiguous
    # residues (X/B/Z/J) — the protein-quality QC step, the amino-acid
    # analogue of the nucleotide ambiguous-character filter. A bad protein
    # fails its segment, which drops the whole isolate (segmented) or the
    # sequence (non-segmented). Counted in **isolates** (segmented) /
    # **sequences** (non-segmented) — the unit a bench scientist acts on —
    # while each dropped segment still lands in ``_qc_removed.tsv``. Always
    # zero unless ``qc.protein_quality.enabled`` is set.
    removed_protein_quality: int = 0
    # Sequences (non-segmented) or isolates (segmented) dropped because
    # an HMM-gated marker spec had no CDS pass the HMM check (E-value or
    # coverage). Distinct from ``removed_proteins`` /
    # ``removed_incomplete_isolates`` so the actionable HMM failures
    # don't blur into generic "marker missing" drops — bench scientists
    # need to see "12 isolates lost their L because RdRp_4 failed" as a
    # separate signal from "12 isolates were missing an L segment".
    removed_hmm_failed: int = 0
    # Per-marker breakdown of the same drops, keyed by the marker spec's
    # ``name`` (e.g. ``{"L": 8, "M": 4}``). Used by the live stderr
    # line and the QC Summary breakdown.
    removed_hmm_by_marker: dict = field(default_factory=dict)
    # Whole-pool length filter is skipped in segmented mode (the input is a
    # mix of segments of very different lengths; per-segment bounds are
    # applied later via segmented.viruses.<v>.segment_lengths instead).
    length_filter_skipped: bool = False
    # Why the whole-genome length filter did not run, so the summary line
    # reads honestly: "segmented" (segmented mode — per-segment bounds apply
    # instead) or "disabled" (non-segmented run with the filter off, the
    # default). ``None`` defaults to the segmented wording for backward
    # compatibility with callers that set length_filter_skipped directly.
    length_filter_skip_reason: Optional[str] = None
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
    # Number of records that survived the *entire* QC pipeline (basic
    # QC plus protein-annotation QC, segmented completeness, taxonomy
    # consistency, strain-collision, and per-segment length). Set by
    # the CLI driver at the end of QC, right before mode selection
    # runs. Units match what the mode actually consumes:
    # ``"sequences"`` for non-segmented runs (one per input record that
    # survived) and ``"isolates"`` for segmented runs (one per
    # concatenated complete isolate). Distinct from ``passed``, which
    # only counts what made it past the basic-QC stages and is set
    # inside ``run_qc`` before the segmented/protein steps run.
    final_survivors: Optional[int] = None
    final_survivors_unit: str = "sequences"
    details: list[dict] = field(default_factory=list)

    # Sequences of special importance that bypassed a QC removal stage via
    # the `overrides.protect_qc` whitelist. Each entry records the id, the
    # stage it was protected against, and the reason it *would* have been
    # removed (so `_overrides.tsv` and the run summary stay transparent —
    # protection is never silent). Populated by overrides.protected_keep.
    protected: list[dict] = field(default_factory=list)

    def add_removed(self, seq_id: str, reason: str) -> None:
        self.details.append({"id": seq_id, "reason": reason})

    def add_protected(self, seq_id: str, stage: str, would_be_reason: str) -> None:
        """Record a sequence kept despite failing ``stage`` (override list)."""
        self.protected.append(
            {"id": seq_id, "stage": stage, "reason": would_be_reason}
        )

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
            why = (
                "filter disabled"
                if self.length_filter_skip_reason == "disabled"
                else "segmented mode"
            )
            length_lines = [f"  Removed (length)    : skipped ({why})"]
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
            f"  Passed basic QC     : {self.passed}",
            dup_line,
            *length_lines,
            f"  Removed (ambiguous) : {self.removed_ambiguous}",
            f"  Removed (annotation): {self.removed_annotation}",
            f"  Removed (proteins)  : {self.removed_proteins}",
            f"  Removed (protein-quality): {self.removed_protein_quality}",
            f"  Removed (incomplete): {self.removed_incomplete_isolates}",
            f"  Removed (tax-mismatch): {self.removed_taxonomy_mismatch}",
            f"  Removed (strain-collision): {self.removed_strain_collisions}",
            f"  Removed (extra-segments): {self.removed_extra_segments}",
        ]
        if self.removed_hmm_failed or self.removed_hmm_by_marker:
            hmm_line = f"  Removed (HMM QC)    : {self.removed_hmm_failed}"
            if self.removed_hmm_by_marker:
                parts = [
                    f"{k}={v}"
                    for k, v in sorted(self.removed_hmm_by_marker.items())
                ]
                hmm_line += f" ({', '.join(parts)})"
            lines.append(hmm_line)
        if self.protected:
            n_ids = len({p["id"] for p in self.protected})
            lines.append(
                f"  Protected (overrides): {n_ids} record(s) kept despite "
                f"failing QC ({len(self.protected)} stage-hit(s))"
            )
        # Add a Final survivors line when the CLI driver populated it.
        # Older callers / tests that build a QCReport directly (without
        # going through _handle_segmented) will see the field as None
        # and the line is omitted, preserving prior behaviour.
        if self.final_survivors is not None:
            lines.append(
                f"  Final survivors     : {self.final_survivors} "
                f"{self.final_survivors_unit} "
                f"(after every QC stage)"
            )
        return "\n".join(lines)


@dataclass
class RunResult:
    mode: str
    representatives: list[Sequence] = field(default_factory=list)
    clusters: list[Cluster] = field(default_factory=list)
    group_stats: list[GroupStat] = field(default_factory=list)
    qc_report: Optional[QCReport] = None
    config_snapshot: dict = field(default_factory=dict)
    # Audit of the force-select override (overrides.force_select): one entry
    # per pinned sequence, {id, action, detail}, where action is one of
    # elected_representative / split_singleton / added_representative /
    # already_representative / unavailable. Written to
    # {prefix}_force_selected.tsv and summarised in {prefix}_summary.md.
    force_selected: list[dict] = field(default_factory=list)
