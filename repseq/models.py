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

    # Display-only context for the console ``summary()`` block. Not part of
    # the QC logic — purely what the rendered tally needs to show live
    # thresholds in labels (``ambiguous_pct``, ``protein_quality_pct``,
    # ``annotation_keywords``, ``length_bounds``), pick the segmented vs
    # non-segmented layout (``segmented``, ``n_segments``), and reconcile
    # the segment→isolate unit change (``segments_entering_assembly`` and
    # the per-cause Phase-2 magnitudes ``incomplete_segments`` /
    # ``concat_marker_isolates`` / ``concat_hmm_isolates`` /
    # ``hmm_isolates_pre_assembly`` / ``duplicate_isolates``). Populated by
    # the CLI driver (``_run_qc`` for the thresholds, ``_handle_segmented``
    # for the assembly magnitudes); empty for directly-constructed reports,
    # in which case ``summary()`` infers the layout and omits the live
    # threshold/magnitude annotations.
    display: dict = field(default_factory=dict)

    # Pre-QC taxonomic snapshot: ``{rank: Counter(taxon -> count)}`` over the
    # input pool *before any QC removal* (after exclusions + metadata
    # resolution). Taken once at the QC chokepoint so the taxonomic report can
    # show a leading pre-QC column and the "eliminated taxa" alarm can compare
    # it against the post-QC pool — without holding the raw input list to end
    # of run. Counts raw input records (segments, in segmented mode). Empty
    # for directly-constructed reports and ``--no-resolve`` runs (no taxonomy).
    pre_qc_taxa: dict = field(default_factory=dict)
    # Taxa (genus and higher) present pre-QC but with zero survivors after QC —
    # the silent-drop alarm. Each entry is ``{rank, taxon, pre_qc_count}``.
    # Computed by the CLI driver from ``pre_qc_taxa`` vs the post-QC pool
    # (output.report.find_eliminated_taxa); surfaced in the console final
    # summary and, via the taxonomic-report TSV, in ``_flags.txt``.
    eliminated_taxa: list = field(default_factory=list)
    # Polyprotein peptide-coverage "walls of zeros": a clade whose home
    # polyprotein spec leaves most of its peptides at 0 % coverage (mistuned
    # profiles), or that no spec slices at all. Each entry is a dict from
    # ``output.report.find_polyprotein_coverage_walls``; surfaced in the console
    # final summary and, via the polyprotein-report TSV, in ``_flags.txt``.
    polyprotein_walls: list = field(default_factory=list)

    def add_removed(self, seq_id: str, reason: str) -> None:
        self.details.append({"id": seq_id, "reason": reason})

    def add_protected(self, seq_id: str, stage: str, would_be_reason: str) -> None:
        """Record a sequence kept despite failing ``stage`` (override list)."""
        self.protected.append(
            {"id": seq_id, "stage": stage, "reason": would_be_reason}
        )

    # ------------------------------------------------------------------
    # Console QC summary
    #
    # The block is rendered as a *running survivor tally* so a bench
    # scientist can see, after each step, both how many records were
    # removed and how many remain. The one genuine subtlety it makes
    # explicit (rather than hiding, as the old flat list did) is the
    # **unit change**: per-record screening counts segments/sequences,
    # but segmented isolate assembly groups segments into isolates
    # (N segments → 1 isolate), so the count both drops *and* changes
    # units at that boundary. That seam — not any single removal step —
    # is what makes "1.6M segments → 12K isolates" look paradoxical.
    # ------------------------------------------------------------------
    @staticmethod
    def _n(v) -> str:
        try:
            return f"{int(v):,}"
        except (TypeError, ValueError):
            return str(v)

    @staticmethod
    def _pct(frac) -> str:
        try:
            return f"{float(frac) * 100:g}%"
        except (TypeError, ValueError):
            return "?%"

    @staticmethod
    def _kw_hint(keywords) -> str:
        kws = [k for k in (keywords or []) if k]
        if not kws:
            return ""
        shown = ", ".join(kws[:3])
        more = ", …" if len(kws) > 3 else ""
        return f" ({shown}{more})"

    @staticmethod
    def _block(rows: list[tuple]) -> list[str]:
        """Render ``(label, removed, running)`` rows with the label left-
        justified and the two value columns right-justified to their own
        max widths. A blank value column collapses (zero width) so blocks
        that only use one column don't carry empty padding."""
        if not rows:
            return []
        lw = max(len(r[0]) for r in rows)
        rw = max(len(r[1]) for r in rows)
        nw = max(len(r[2]) for r in rows)
        out: list[str] = []
        for label, removed, running in rows:
            line = "  " + label.ljust(lw)
            if rw:
                line += "  " + removed.rjust(rw)
            if nw:
                line += "  " + running.rjust(nw)
            out.append(line.rstrip())
        return out

    def _length_label(self, d: dict) -> str:
        bounds = d.get("length_bounds") or {}
        mn, mx = bounds.get("min"), bounds.get("max")
        bits = []
        if mn is not None:
            bits.append(f"<{self._n(mn)}")
        if mx is not None:
            bits.append(f">{self._n(mx)}")
        rng = f" ({' or '.join(bits)} nt)" if bits else ""
        return "genome length out of bounds" + rng

    def summary(self) -> str:
        d = self.display or {}
        segmented = d.get("segmented")
        if segmented is None:
            # Infer for directly-constructed reports (tests, legacy callers).
            segmented = (
                self.final_survivors_unit == "isolates"
                or self.length_filter_skip_reason == "segmented"
                or self.dedup_skipped
                or bool(self.removed_length_by_segment)
            )
        if segmented:
            return self._summary_segmented(d)
        return self._summary_nonsegmented(d)

    def _summary_nonsegmented(self, d: dict) -> str:
        n = self._n
        lines = ["QC summary", "", "Per-record screening  (unit: sequences)"]
        run = self.total_input
        rows: list[tuple] = [("input sequences", "", n(run))]
        stages = [
            (self.removed_duplicates, "exact duplicates"),
            (0 if self.length_filter_skipped else self.removed_length,
             self._length_label(d)),
            (self.removed_annotation,
             "description keywords" + self._kw_hint(d.get("annotation_keywords"))),
            (self.removed_ambiguous,
             f"ambiguous nucleotides (>{self._pct(d.get('ambiguous_pct', 0.05))} N)"),
            (self.removed_proteins, "wrong / missing CDS count"),
            (self.removed_protein_quality,
             f"ambiguous protein residues "
             f"(>{self._pct(d.get('protein_quality_pct', 0.05))} X)"),
            (self.removed_hmm_failed, "marker fails HMM identity"),
        ]
        for removed, label in stages:
            if not removed:
                continue
            run -= removed
            rows.append((label, "−" + n(removed), n(run)))
        final = self.final_survivors if self.final_survivors is not None else run
        rows.append(("sequences passing QC", "", n(final)))
        lines += self._block(rows)
        lines += self._protected_lines()
        lines.append("")
        lines.append(f"{n(self.total_input)} input sequences → {n(final)} passing QC")
        return "\n".join(lines)

    def _summary_segmented(self, d: dict) -> str:
        n = self._n
        lines = ["QC summary", "", "Phase 1 · per-record screening  (unit: segments)"]

        # Pure per-segment removal stages — these advance the segment tally.
        run = self.total_input
        rows: list[tuple] = [("input records", "", n(run))]
        seg_stages = [
            (self.removed_annotation,
             "description keywords" + self._kw_hint(d.get("annotation_keywords"))),
            (self.removed_ambiguous,
             f"ambiguous nucleotides (>{self._pct(d.get('ambiguous_pct', 0.05))} N)"),
            (self.removed_proteins, "wrong / missing CDS count"),
            (self.removed_taxonomy_mismatch, "cross-segment species mismatch"),
            (self.removed_strain_collisions, "strain-id collision"),
        ]
        for removed, label in seg_stages:
            if not removed:
                continue
            run -= removed
            rows.append((label, "−" + n(removed), n(run)))
        rows.append(("segments passing per-record screening", "", n(run)))
        lines += self._block(rows)

        # Isolate-level identity gates run on the segment pool but drop the
        # whole isolate, so they are counted (and shown) in isolates — a
        # separate sub-tally, reconciled back to the pool by the
        # "segments entering assembly" ground-truth line below.
        pq = self.removed_protein_quality
        hmm_pre = d.get("hmm_isolates_pre_assembly", self.removed_hmm_failed)
        if pq or hmm_pre:
            lines.append("")
            lines.append(
                "  isolate-level identity gates — a failure drops the whole "
                "isolate  (unit: isolates)"
            )
            gate_rows: list[tuple] = []
            if pq:
                gate_rows.append((
                    f"ambiguous protein residues "
                    f"(>{self._pct(d.get('protein_quality_pct', 0.05))} X)",
                    "", "−" + n(pq) + " isolates"))
            if hmm_pre:
                gate_rows.append(
                    ("marker fails HMM identity", "", "−" + n(hmm_pre) + " isolates"))
            # Reconcile the isolate-gate removals back to the segment pool.
            # Only meaningful here — with no gates the Phase-1 subtotal already
            # is the pool entering assembly, so the line would be redundant.
            sea = d.get("segments_entering_assembly")
            if sea is not None:
                gate_rows.append(("segments entering assembly", "", n(sea)))
            lines += self._block(gate_rows)

        # Phase 2 — isolate assembly. Units are now isolates.
        lines.append("")
        nseg = d.get("n_segments")
        seg_word = f"{nseg} segments → 1 isolate" if nseg else "segments → 1 isolate"
        lines.append(f"Phase 2 · isolate assembly  ({seg_word}; unit: isolates)")
        p2: list[tuple] = []
        inc = d.get("incomplete_segments", self.removed_incomplete_isolates)
        if inc:
            miss = (f"missing ≥1 of {nseg} segments / unidentifiable"
                    if nseg else "missing segments / unidentifiable")
            p2.append((miss, "", "≈" + n(inc) + " segments"))
        if self.removed_length_by_segment:
            tot = sum(c["too_short"] + c["too_long"]
                      for c in self.removed_length_by_segment.values())
            if tot:
                p2.append(("segment length out of bounds", "",
                           "−" + n(tot) + " isolates"))
        if self.removed_extra_segments:
            p2.append(("extra / non-canonical segments", "",
                       "−" + n(self.removed_extra_segments) + " isolates"))
        cm = d.get("concat_marker_isolates", 0)
        if cm:
            hmm_c = d.get("concat_hmm_isolates", 0)
            note = f" ({n(hmm_c)} via HMM)" if hmm_c else ""
            p2.append(("marker lost at concatenation", "",
                       "−" + n(cm) + " isolates" + note))
        dup = d.get("duplicate_isolates", 0)
        if dup:
            p2.append(("identical duplicate genomes", "",
                       "−" + n(dup) + " isolates"))
        final = self.final_survivors
        if final is not None:
            p2.append(("complete isolates assembled", "", n(final)))
        lines += self._block(p2)

        # Optional per-segment length breakdown (which segment's bound bit).
        if self.removed_length_by_segment:
            for seg_name in sorted(self.removed_length_by_segment.keys()):
                counts = self.removed_length_by_segment[seg_name]
                if counts.get("too_short"):
                    lines.append(f"    {seg_name} too short: {n(counts['too_short'])}")
                if counts.get("too_long"):
                    lines.append(f"    {seg_name} too long : {n(counts['too_long'])}")

        lines += self._protected_lines()
        lines.append("")
        if final is not None:
            lines.append(
                f"{n(self.total_input)} input segments → {n(final)} complete isolates"
            )
        lines.append(
            "not applied in segmented mode: exact-duplicate & whole-genome "
            "length filters"
        )
        return "\n".join(lines)

    def _protected_lines(self) -> list[str]:
        if not self.protected:
            return []
        n_ids = len({p["id"] for p in self.protected})
        return [
            "",
            f"  Protected (overrides): {n_ids} record(s) kept despite failing "
            f"QC ({len(self.protected)} stage-hit(s))",
        ]


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
