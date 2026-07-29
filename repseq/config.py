"""YAML configuration loading, validation, and defaults."""

from __future__ import annotations

import copy
import difflib
import os
import re
from pathlib import Path
from typing import Any, Optional

import yaml

from .errors import ConfigError


# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

DEFAULTS: dict[str, Any] = {
    "cache_dir": "~/.repseq/cache",
    "temp_dir": "/tmp/repseq",
    "threads": 4,
    "seed": 42,
    "qc": {
        "remove_duplicates": True,
        "genome_length_filter": {
            # Absolute nucleotide-length bounds on whole input sequences.
            # NON-SEGMENTED MODE ONLY — enabling this with segmented.enabled
            # is a config error (segmented runs use per-segment bounds via
            # segmented.viruses.<v>.segment_lengths instead). Default off:
            # the user must opt in and supply at least one of min/max.
            "enabled": False,
            "min": None,   # drop sequences shorter than this many nt
            "max": None,   # drop sequences longer than this many nt
        },
        "ambiguous_threshold": 0.05,
        "protein_annotation": {
            # Drop sequences whose NCBI GenBank record has fewer than
            # min_proteins annotated CDS features. Requires network access
            # (skipped automatically with --no-resolve).
            "enabled": False,
            "min_proteins": 1,
        },
        "protein_quality": {
            # Amino-acid analogue of ambiguous_threshold. Drop a CDS
            # protein when the fraction of ambiguous residues (X/B/Z/J) in
            # its translation exceeds max_bad_fraction; a bad protein fails
            # its segment, which drops the whole isolate (segmented mode) or
            # the sequence (non-segmented). An empty/absent translation
            # counts as fully bad. Network-dependent: when enabled it
            # force-fetches GenBank CDS translations if no earlier step
            # already did. Skipped under --no-resolve.
            "enabled": False,
            "max_bad_fraction": 0.05,
        },
        "annotation_filter": {
            "enabled": True,
            "keywords": [
                "MAG:",
                "metagenome-assembled",
                "synthetic",
                "artificial",
                "fragment",
                "partial",
                "environmental sample",
                "uncultured",
                "unclassified",
                "unidentified",
                "hypothetical",
            ],
        },
    },
    "segmented": {
        "enabled": False,
        "virus": None,
        "viruses": {},
        # When true (the default), repseq fetches the GenBank source feature
        # for each NCBI-sourced sequence and uses its /isolate, /strain, and
        # /segment qualifiers in preference to the header-regex parse. The
        # regex still runs as a fallback for sequences without an accession,
        # for UniProt input, when --no-resolve is set, or when the GenBank
        # record lacks the qualifier. Set to false to bypass the GenBank
        # lookup entirely (header-regex only).
        "use_genbank_metadata": True,
        # Drop any segmented isolate whose segments disagree on the
        # taxonomic rank named in ``rank``. Reassortment between
        # different parent species is real biology for many segmented
        # viruses (peribunyaviruses, orthomyxoviruses, …), so the user
        # who turns this on is asking for monophyletic-at-this-rank
        # isolates only. Off-by-default would let those isolates
        # through unflagged, so the default is on. Missing values are
        # ignored — an isolate is only dropped when *populated* labels
        # disagree.
        "taxonomy_consistency": {
            "enabled": True,
            "rank": "species",
        },
        # What to do when the strain-collision detector finds two or more
        # distinct accessions sharing the same strain-derived isolate_id
        # AND the same segment — the over-merge signature of the
        # /strain → isolate_id fallback in _populate_genbank_isolate_segment.
        # "warn" (default) prints one line per collision to stderr and
        # lets the pipeline continue (the over-merged isolate keeps the
        # longest sequence per segment, the rest get dedup-dropped
        # downstream). "drop" removes every accession involved in any
        # collision before the completeness filter runs, adds them to
        # _qc_removed.tsv with reason "strain_collision:<segment>", and
        # increments QCReport.removed_strain_collisions.
        "strain_collision_action": "warn",
        # What to do when an isolate's segments include names outside
        # the configured ``segments`` list (e.g. an L/M/S virus that
        # has a fourth segment, or a non-canonical segment identifier
        # that ``identify_segment`` couldn't map). "warn" (default)
        # prints one line per affected isolate to stderr and pipes the
        # isolate through with just its expected segments — today's
        # silent prune becomes visible without changing analysis.
        # "drop" removes the entire isolate (every segment lands in
        # _qc_removed.tsv with reason "extra_segments:<extras>") and
        # increments QCReport.removed_extra_segments (units: isolates).
        "extra_segments_action": "warn",
    },
    "clustering": {
        "backend": "mmseqs2",              # "mmseqs2" | "cdhit"
        # Alphabet fed to the clustering backend.
        #
        # IMPORTANT: this setting only chooses what the clustering backend
        # (mmseqs2 / cd-hit) sees. It does NOT disable GenBank CDS download,
        # protein-count QC (qc.protein_annotation), or the
        # virus.expected_proteins_per_segment check — those run on every
        # isolate regardless of this value.
        #
        #   "protein"    — cluster on amino acid sequences (recommended for
        #                  diverged virus families). Non-segmented: the
        #                  marker protein (longest CDS, or the first
        #                  matching alias in `cluster_protein`).
        #                  Segmented: in-order concat of each segment's
        #                  marker protein. Triggers a one-shot GenBank CDS
        #                  fetch when proteins aren't already cached.
        #   "nucleotide" — cluster on the raw nucleotide sequence.
        #                  Non-segmented: the input FASTA sequence as-is.
        #                  Segmented: concat of all segments in
        #                  `segments` order.
        "alphabet_for_clustering": "protein",
        # Non-segmented marker-protein override. Alias list, matched
        # case-insensitively as substrings against /product. First alias
        # that matches a CDS wins; if no aliases match (or the list is
        # empty), the longest CDS on the sequence is used.
        "cluster_protein": [],
        # Non-segmented multi-marker clustering (v0.38.0+). When false
        # (default) the clustering string is the SINGLE marker from the
        # first cluster_protein spec that a CDS satisfies (legacy
        # behaviour — e.g. Spike alone for coronaviruses). When true, the
        # marker CDS from EVERY cluster_protein spec is selected and the
        # AA strings are concatenated in declared spec order
        # (e.g. Spike+Nucleocapsid) into seq.protein_sequence, mirroring
        # the segmented per-segment concat. A sequence missing any
        # required marker is dropped (the HMM AND-gate already enforces
        # this upstream when every spec declares `hmms`). No effect in
        # segmented mode (segmented always concatenates its per-segment
        # markers) or when alphabet_for_clustering="nucleotide". The
        # whole-genome tree (2E) is unaffected by this flag — for a
        # multi-gene tree use phylo.tool: auto/iqtree so the partitioned
        # supermatrix path aligns each marker separately.
        "concatenate_markers": False,
        # Additional proteins to extract and (optionally) phylogenize but
        # NOT use for clustering or the whole-genome tree. Same schema as
        # cluster_protein (list of {name, aliases?, hmms?} dicts). For
        # every entry, a single CDS per representative (HMM-gated when
        # `hmms` is declared and the HMM tier ran; alias fallback against
        # /product otherwise) is written to
        # {prefix}_extra_protein_fasta/{prefix}_<name>.fasta. With
        # --per-protein-phylo, one tree per entry is also emitted into
        # {prefix}_extra_protein/. Intended for proteins that are sparse
        # across taxa (e.g. coronavirus ORF7) — clustering them would
        # systematically drop representatives that lack the protein.
        "extra_protein": [],
        # Polyprotein cutting (v0.33.0+): slice each representative's
        # polyprotein CDS into its mature peptides using one HMM per
        # peptide. PURELY ADDITIVE — the polyprotein still drives
        # clustering and the whole-genome tree; mature peptides are
        # emitted as accessory artifacts (per-peptide FASTAs in
        # {prefix}_polyprotein/, plus an audit TSV per spec). Each list
        # entry is a dict:
        #   { name: <spec_name, unique>,
        #     peptides: [ {name, hmm, cleavage_motif?}, ... ],  # N→C order
        #     cut_strategy: "boundary" | "bisect" | "motif",
        #     motif_window_aa: 10,
        #     min_peptides_hit: 2 }     # parent-CDS identification floor
        # cut_strategy defaults to "motif" if any peptide carries
        # cleavage_motif, else "bisect". Requires the HMM tier active and
        # each peptide HMM present in the HMM database. Soft-fails with a
        # stderr line when the HMM tier didn't run.
        "polyprotein": [],
        # Diagnostic: for each stratum where binary-search clustering ran,
        # also run the clustering backend at each of these identity
        # thresholds and report the cluster count as additional columns
        # (n_clusters_0.99, n_clusters_0.95, ...) in
        # {prefix}_group_counts.tsv. Reporting-only — does NOT influence
        # representative selection. Cutoffs below the backend's identity
        # floor (cd-hit-est: 0.80, cd-hit protein: 0.40, MMseqs2: 0) are
        # reported as NA. Set to [] to disable (e.g. on huge runs where
        # the extra clustering work would dominate runtime).
        "diversity_curve_cutoffs": [0.99, 0.95, 0.9, 0.8, 0.7],
        "mmseqs2_mode": "easy-linclust",   # "easy-linclust" | "easy-cluster"
        "coverage": 0.8,
        "coverage_mode": 0,
        "extra_args": [],
        "cdhit": {
            # Binary auto-selected from input alphabet:
            #   protein  -> cd-hit
            #   nucleic  -> cd-hit-est
            # Override to pin a specific path or variant.
            "binary": None,
            # Word size (-n). None = auto-pick from threshold per the cd-hit
            # user guide; cd-hit refuses out-of-range -n for a given -c.
            "word_size": None,
            # Shorter-sequence coverage (-aS); only honoured when
            # global_alignment is False (cd-hit's -G 0).
            "coverage": 0.8,
            # -G: True = global identity (cd-hit default), False = local.
            "global_alignment": True,
            # -g: False = greedy (fast, default), True = accurate
            # (slower; compares each input against every existing cluster).
            "accurate": False,
            # -M: memory cap in MB. 0 = unlimited (cd-hit default).
            "memory_mb": 0,
            # Raw cd-hit flags appended verbatim (e.g. ["-s", "0.8"]).
            "extra_args": [],
        },
    },
    "hmm": {
        # HMM-based marker-protein selection. When enabled AND a marker
        # has `hmms: [...]` configured AND hmmscan (HMMER) is on PATH,
        # HMM hits become the AUTHORITATIVE gate for that marker — a CDS
        # whose /product matches an alias but FAILS the HMM check is
        # rejected, and the segment / sequence is dropped. Markers
        # without `hmms` fall through to the legacy alias → longest
        # chain. Soft-fails (warns and falls back) when hmmscan is
        # missing or the database is unavailable.
        "enabled": True,
        # null = use the bundled viral-core set (`repseq/data/hmms/
        # repseq_viral_core.hmm`, 19 Pfam-A profiles for RdRp,
        # nucleocapsid, glycoprotein, protease, etc.). Absolute path = user-
        # supplied .hmm file; auto-`hmmpress`-ed on first use if the
        # .h3* index files are missing.
        "database": None,
        # E-value cutoff used when a profile has no curated Pfam GA
        # (gathering threshold). Hits with E ≤ this value pass.
        "default_evalue": 1.0e-5,
        # When true, use each profile's GA cutoff when available; fall
        # back to default_evalue otherwise. When false, always use
        # default_evalue regardless of GA availability.
        "use_ga_when_available": True,
        # Length cutoff: ali_span / hmm_model_length must be ≥ this.
        # Guards against tiny single-domain hits being accepted as full
        # marker matches. Range (0, 1].
        "relative_length_cutoff": 0.5,
        # Multidomain-token overlap tolerance, in amino acids. Consecutive
        # named domains in an N-to-C token (`A--B--C`) may overlap by up
        # to this many residues at their boundary AND still satisfy the
        # token. The walk also requires forward progression — the next
        # hit must start strictly C-terminal to the prior hit's start AND
        # end strictly C-terminal to the prior hit's end — so a single
        # fused region can never masquerade as two domains, and full
        # containment is always rejected. Set 0 for strict non-overlap
        # (pre-v0.22 behaviour). 30 (default) is comfortable for typical
        # Pfam viral-glycoprotein profiles whose boundaries don't align
        # exactly to the molecular cleavage site (e.g. coronavirus
        # S1/S2 across the furin site, bunya G1/G2).
        "multidomain_overlap_tolerance": 30,
        # null = use cfg.threads.
        "threads": None,
    },
    "representative": {
        "priority": ["refseq", "reviewed_uniprot", "longest"],
    },
    "phylo": {
        # Optional MSA + phylogeny step. Triggered with --phylo on any
        # mode subcommand; skipped automatically if fewer than 3
        # representatives survive selection.
        #
        # Tree-builder selection:
        #   "auto"     — IQ-TREE for protein alignments, FastTree for
        #                nucleotide. Default; matches what each tool is
        #                best at.
        #   "iqtree"   — always IQ-TREE (slower; ModelFinder + UFBoot)
        #   "fasttree" — always FastTree (faster; approximate-ML)
        "tool": "auto",
        # Retain the plain-text Newick (*_tree.nwk) for every tree built
        # this run (whole-genome 2E, per-protein/extra/segment 2F,
        # pre-cluster, partition). OFF by default to reduce output clutter:
        # the annotated phyloXML (*_tree.xml) is a topological superset, and
        # the short-id *_tree_id_map.tsv (which decodes the retained
        # *_msa.fasta leaves) is kept regardless. The Newick is still
        # generated internally during the run — the phyloXML is re-parsed
        # from it and the incongruence RF table reads it — then dropped at
        # the end unless this is true. CLI: --newick / --no-newick.
        "newick": False,
        # Render a graphical PDF + PNG of every phyloXML tree built this run
        # (whole-genome 2E, per-protein/extra/segment 2F, pre-cluster,
        # partition). ON by default. Each {prefix}..._tree.xml gets a sibling
        # {prefix}..._tree.pdf and {prefix}..._tree.png — a ladderized
        # rectangular phylogram (matplotlib + Bio.Phylo) with taxonomy-coloured
        # leaf labels, a genus/subfamily legend, and branch-support labels,
        # reconstructed entirely from the phyloXML. Soft-fails with a single
        # stderr line (and emits no figures) when matplotlib is unavailable —
        # same posture as --plot; matplotlib is the [viz] extra. CLI:
        # --pdf / --no-pdf.
        "pdf": True,
        "mafft": {
            # Raw mafft flags appended to "mafft --auto --thread N <input>".
            # Examples: ["--maxiterate", "1000"] for L-INS-i; or
            # ["--retree", "1"] for a faster pass on a very large input.
            "extra_args": [],
            # Pass --auto to MAFFT. Set false to drop it and rely on
            # extra_args alone (used by the --fast CLI flag to force
            # single-pass FFT-NS-1).
            "use_auto": True,
        },
        # Optional alignment trimming (trimAl) between MAFFT and the
        # tree-builder, for the whole-genome tree (2E). OFF by default.
        # In partitioned mode each per-family alignment is trimmed BEFORE
        # concatenation (so the partition charset ranges stay valid).
        # Soft-fails (loud warning + builds on the UNTRIMMED alignment)
        # when trimal is missing, errors, or strips the alignment to
        # nothing. The per-protein trees have their own knob at
        # phylo.per_protein.trimal.
        "trimal": {
            "enabled": False,
            # trimAl method → the `-<mode>` flag. "automated1" is trimAl's
            # heuristic best-for-ML-trees default; other column-trimming
            # methods: gappyout, strict, strictplus, nogaps, noallgaps.
            # Threshold trimming (-gt / -st / -cons) goes in extra_args.
            "mode": "automated1",
            # Raw trimal flags appended verbatim, e.g. ["-gt", "0.8"].
            "extra_args": [],
        },
        "fasttree": {
            # Raw FastTree flags appended to its argv. The protein /
            # nucleotide model is picked automatically from the rep
            # alphabet (default JTT for protein, -nt -gtr for nucleotide).
            "extra_args": [],
        },
        "iqtree": {
            # Binary auto-detected: tries iqtree2 first, then iqtree.
            # Set to a name or absolute path to pin a specific build.
            "binary": None,
            # Substitution model. "MFP" runs ModelFinder Plus (recommended
            # for protein — JTT/WAG/LG/etc. tested by BIC). Pin a model
            # like "LG+G4" or "JTT+G4" for a faster fixed-model pass.
            "model": "MFP",
            # Ultrafast bootstrap replicates. 0 disables bootstrap (faster).
            # IQ-TREE recommends >= 1000 when reporting branch support.
            "ultrafast_bootstrap": 1000,
            # Raw flags appended verbatim, e.g. ["-alrt", "1000"] for SH-aLRT.
            "extra_args": [],
        },
        # Partitioned-supermatrix tree (the principled multi-marker
        # analysis). When enabled AND the run resolves to protein + IQ-TREE
        # AND the HMM tier resolved >= 2 marker families, the whole-genome
        # tree (2E) is built by aligning each marker family separately, then
        # concatenating the per-family MSAs column-wise into a supermatrix
        # and letting IQ-TREE fit a model PER partition (-p/-q/-Q). This
        # replaces gluing the markers into one string + one model + one
        # MAFFT (which can align unrelated proteins across segment seams).
        # Soft-falls back to that concat-then-align path when the run can't
        # be partitioned (FastTree, no HMM families, < 2 families).
        "partition": {
            # Master switch. Default ON for protein + IQ-TREE runs; set
            # false to force the legacy concat-then-align behaviour.
            "enabled": True,
            # IQ-TREE partition linkage (Chernomor et al. 2016):
            #   "proportional" — -p, edge-linked proportional: one shared
            #                    branch-length set + a per-partition rate
            #                    multiplier. The standard default.
            #   "equal"        — -q, edge-equal: all partitions share branch
            #                    lengths (fewest parameters).
            #   "unlinked"     — -Q, edge-unlinked: each partition gets its
            #                    own branch lengths (most flexible; many more
            #                    parameters — useful when segments may have
            #                    different histories, e.g. reassortment).
            "linkage": "proportional",
            # Optional per-family substitution-model pin, keyed by the family
            # LABEL (segment_token, e.g. "L_RdRP_4" or "M_Bunya_G1--Bunya_G2";
            # the segment prefix is dropped in non-segmented runs). A family
            # absent here gets "MFP" → ModelFinder picks its model. Example:
            #   models: {L_RdRP_4: "LG+G4", S_Bunya_nucleocap: "WAG+G4"}
            "models": {},
        },
        # Per-leaf display labels on the phyloXML tree.
        # Supported placeholders: {species}, {genus}, {subgenus},
        # {subfamily}, {family}, {order}, {class}, {phylum}, {id},
        # {accession}, {host}, {strain}, {isolate_id}, {country},
        # {date}, {year}, {organism}.
        "labeling": {
            # Default for non-segmented runs.
            "format": "{species}|{id}|{host}",
            # Used when segmented.enabled is true. Falls back to {format}
            # if null. When {strain} is requested but the GenBank record
            # has no /strain qualifier, the writer substitutes
            # {isolate_id} (which segmented mode always has) so the
            # label never collapses to ``...||host``.
            "segmented_format": "{species}|{strain}|{subtype}|{host}",
            # Replace internal whitespace runs in each placeholder value
            # with underscores. Keeps the label round-trippable through
            # tree viewers that treat whitespace as a token boundary.
            "replace_whitespace": True,
            # When a placeholder resolves to empty (and isn't a {strain}
            # that can fall back to {isolate_id}), drop the placeholder
            # AND the single separator character immediately before it,
            # so the rendered label never contains ``||`` or trailing
            # ``|``. Set true to keep all separators verbatim.
            "keep_separator_on_empty": False,
        },
        # Tree rooting. Tree-building tools (FastTree, IQ-TREE) produce
        # unrooted trees; the post-processing step picks a root before
        # the writer ladderizes + serialises.
        #
        # method:
        #   auto       — try taxonomy-guided → MAD → midpoint, first
        #                success wins. The most-likely-correct default.
        #   taxonomy   — root at the branch that maximises mean LCA
        #                specificity of internal clades against the
        #                resolved NCBI lineages. Falls through to
        #                midpoint if no leaves carry lineage data.
        #   mad        — Minimal Ancestor Deviation (Tria et al. 2017).
        #                Pure-Python implementation; robust when
        #                taxonomy is sparse.
        #   midpoint   — Bio.Phylo's root_at_midpoint. Last-resort
        #                fallback; always succeeds.
        #   none       — leave the tree as parsed (use when the input
        #                is already rooted, e.g. by an outgroup).
        "rooting": {
            # method: "auto" | "taxonomy" | "mad" | "midpoint" |
            #         "outgroup" | "none"
            #
            # "outgroup" uses a user-specified accession or clade as the
            # outgroup. Configure via either or both:
            #   outgroup:      "AB123456"  OR  ["AB123456", "CD789012"]
            #   outgroup_rank: {family: "Hantaviridae"}
            # Accessions match against rep.accession / rep.id / rep.isolate_id
            # (case-insensitive). outgroup_rank pulls every rep whose
            # taxonomy carries that taxon at that rank. The two sources are
            # unioned; multi-leaf outgroups root at the leaves' MRCA. When
            # the spec matches no representative the rooter falls through
            # to midpoint with a stderr note (so a typo never voids the tree).
            "method": "auto",
            "outgroup": None,
            "outgroup_rank": None,
        },
        # Internal-node LCA labels. After rooting, every internal
        # clade is labelled with the lowest common ancestor of its
        # terminals (read from the resolved NCBI lineage) and the
        # label's rank attached as a PhyloXML <rank>.
        "lca": {
            # Master switch.
            "enabled": True,
            # Leaves whose lineage doesn't reach this rank are
            # excluded from the LCA *vote* — they stay on the tree
            # but don't pull internal labels toward an over-coarse
            # taxon. "none" disables the gate. Default "genus" suits
            # viral data, where lots of leaves lack species-level
            # classification.
            "min_rank": "genus",
            # An internal node is annotated only if at least this
            # fraction of its terminals carry usable lineage data.
            # Guards against a handful of well-annotated leaves
            # dictating the label of a much larger bare clade.
            "coverage_threshold": 0.5,
        },
        # Phylogeny-based taxonomy review (v0.39.0). OFF by default — a
        # new scientific inference step. When enabled (and --phylo built a
        # tree), each representative leaf is checked against the smallest
        # well-supported, taxonomically-pure clade enclosing it: a blank
        # rank is imputed from the clade's majority value, a populated rank
        # that disagrees is flagged (never auto-changed). Writes
        # {prefix}_taxonomy_review.tsv; with write_corrected, also emits
        # *_corrected copies of the rep TSV + protein FASTA with
        # high-confidence imputed blanks filled (each tagged "imputed").
        # The tree, its colouring, and its LCA labels are NOT modified.
        "taxonomy_review": {
            "enabled": False,
            # Ranks to evaluate, coarse→fine; evaluated in that order so
            # imputations stay hierarchy-consistent (a finer rank is only
            # imputed from neighbours agreeing on the coarser ones).
            # Species is deliberately omitted — viral species monophyly is
            # too often violated to trust this way.
            "ranks": ["family", "genus", "subgenus"],
            # "high"-confidence bar (only high-confidence imputations are
            # written into the corrected copies; everything down to the
            # relaxed "medium" bar is still listed in the review TSV).
            "min_support": 90,        # enclosing-clade branch support, 0-100
            "min_purity": 0.9,        # fraction of labelled neighbours agreeing
            "min_agreeing": 3,        # min labelled neighbours backing the call
            "require_refseq_anchor": True,  # a RefSeq/reviewed leaf must agree
            # Also write *_corrected copies of the rep TSV + protein FASTA
            # with high-confidence imputed blanks filled (originals kept).
            "write_corrected": True,
        },
        # MSA conservation scoring (always-on when MSAs are produced).
        # After the phylo steps, every alignment written this run (the
        # genome tree 2E, the partition per-family alignments +
        # supermatrix, the per-protein 2F trees, the extra-protein and
        # polyprotein-peptide trees, the per-segment NT trees) is scored
        # for overall conservation and the results collected into one
        # file, {prefix}_msa_conservation.tsv. The metric is the mean
        # per-column Jensen-Shannon divergence to a residue background
        # (Capra & Singh 2007), with Henikoff & Henikoff (1994) sequence
        # weighting and a (1 - gap-fraction) gap penalty. Bounded [0,1]:
        # ~0 = unrelated/background columns, ~0.85-0.95 = a perfectly
        # conserved column (the JSD ceiling depends on the conserved
        # residue's background frequency, so it never reaches exactly 1).
        "conservation": {
            "enabled": True,
        },
        # Per-taxon monophyly report ({prefix}_monophyly.tsv), an always-on
        # sweep over every tree built this run. ``min_support`` makes it
        # support-aware: internal branches with support below this value are
        # collapsed into polytomies before assessing, and a taxon is called
        # monophyletic when no *well-supported* node contradicts it (it could
        # be a clade under some resolution of the collapsed polytomies). So
        # only confident non-monophyly is flagged — a taxon broken solely by
        # weakly-supported branches reads as monophyletic. 0 disables the
        # collapse (topology-only, every branch trusted). Range [0, 100].
        # ``include_species`` adds species-rank rows to the report. Off by
        # default because viral species labels are annotation-noisy (the same
        # reason the taxonomic reports skip species); turn it on for
        # reassortment / misannotation analysis, where species is the rank at
        # which a taxon broken on a marker tree but clean on the genome tree
        # localises the signal.
        "monophyly": {
            "min_support": 70,
            "include_species": False,
        },
        # PhyloXML writer knobs.
        "phyloxml": {
            # Override the <confidence type="..."> attribute. ``auto``
            # picks ``sh_like`` for FastTree and ``ufboot`` for IQ-TREE,
            # which matches what each tool actually produces by default.
            # Set explicitly if you pass non-default tree args (e.g.
            # IQ-TREE ``-b`` for classical bootstrap, in which case use
            # ``bootstrap``; ``-alrt`` only, use ``sh_alrt``).
            "confidence_type": "auto",
        },
        # Taxonomy-driven leaf colouring of the phyloXML output. Each
        # external node gets a <property ref="style:font_color"
        # datatype="xsd:token" applies_to="node">#RRGGBB</property> so a
        # viewer (Archaeopteryx) tints the leaf label by taxonomy. The
        # palette is shared across the whole-genome tree (2E) and every
        # per-protein tree (2F), so a taxon is the same colour everywhere.
        "coloring": {
            # Master switch. Default ON, colouring by genus.
            "enabled": True,
            # One or two taxonomy ranks (from the 9-rank ladder or
            # phylum/kingdom/superkingdom). ONE rank → each value gets a
            # distinct hue. TWO ranks [parent, child] → the child fans
            # out across hues of its parent's base colour (e.g.
            # [genus, subgenus]: subgenera as shades-of-genus). List the
            # coarser rank first.
            "ranks": ["genus"],
            # HSV saturation / value for the generated palette, in [0,1].
            # Defaults give legible mid-tone colours on a white canvas.
            "saturation": 0.65,
            "value": 0.90,
            # Colour for leaves whose rank is empty / unknown / na / ? /
            # etc. — a medium grey by default. Must be #RRGGBB.
            "missing_color": "#808080",
        },
        # Per-protein trees (2F), triggered with --per-protein-phylo.
        # One tree per declared HMM domain-architecture token (the same
        # hmms: used for QC). Requires the HMM tier to have run.
        "per_protein": {
            # A family is built only when at least this many
            # representatives carry the architecture. Never drops below
            # 3 (the tree-builder floor) regardless of what's set here.
            "min_taxa": 3,
            # After building the per-protein trees, write a pairwise
            # unrooted Robinson-Foulds distance table
            # ({prefix}_incongruence.tsv) quantifying topological
            # incongruence between the marker trees (and the
            # whole-genome tree from --phylo, when present). High RF
            # between marker trees is the reassortment/recombination
            # signal. Needs >= 2 trees to compare; soft-fails otherwise.
            "incongruence": True,
            # MAFFT args for the per-protein (single-gene) alignments.
            # Default is empty → MAFFT --auto (fast; size-adaptive), the
            # same as the whole-genome tree. For a high-accuracy
            # publication run set this to ["--maxiterate", "1000",
            # "--localpair"] (L-INS-i) — these single-gene alignments are
            # small enough to afford it, but it's noticeably slower, so
            # it's opt-in. A non-empty list is passed to MAFFT WITHOUT
            # --auto (so the explicit strategy takes effect).
            "mafft": {
                "extra_args": [],
            },
            # Emit a phyloXML <domain_architecture> on each per-protein
            # tree leaf, built from that CDS's HMM hits (every hit, with
            # its E-value as the domain confidence — so Archaeopteryx's
            # interactive E-value slider can filter them). Default on.
            "domain_architecture": True,
            # In addition to the per-mature-peptide trees, build ONE tree
            # per declared polyprotein spec on the WHOLE polyprotein CDS
            # (e.g. the entire ORF1ab), one leaf per representative carrying
            # it. The headline benefit is the phyloXML <domain_architecture>:
            # every HMM hit across the whole polyprotein is drawn as a domain
            # box, so a leaf shows the full nsp/peptide layout end-to-end (vs
            # the single-peptide box on a peptide tree). OFF by default —
            # it's one extra (large) alignment + tree per spec. Only fires
            # when --per-protein-phylo runs with a polyprotein spec and the
            # HMM tier is active. Outputs to {prefix}_polyprotein/ as
            # {prefix}_<spec>_polyprotein_tree.{xml,pdf,png}. The unaligned
            # whole-polyprotein FASTA ({prefix}_<spec>_polyprotein.fasta) is
            # written independently of this knob (always-on, parallel to the
            # per-peptide FASTAs).
            "whole_polyprotein_tree": False,
            # Optional trimAl trimming for the per-protein (single-gene)
            # alignments — independent of the whole-genome phylo.trimal.
            # OFF by default; same shape (mode → -<mode>, default
            # automated1; extra_args raw passthrough). Soft-fails to the
            # untrimmed alignment exactly like phylo.trimal.
            "trimal": {
                "enabled": False,
                "mode": "automated1",
                "extra_args": [],
            },
        },
        # Pre-cluster overview tree (2H, v0.32.0+), triggered with
        # --pre-cluster-tree or by setting `enabled: true` here. Builds a
        # single rough tree of EVERY post-QC sequence (one leaf per
        # CONCAT isolate in segmented mode) so the bench scientist can
        # see at a glance where the elected representatives land in the
        # broader diversity. Capped at max_leaves leaves (default 5000;
        # reps always kept, background randomly sampled) so a huge pool
        # builds + stays legible. Pipeline is hard-coded for speed
        # regardless of the rest of phylo: MAFFT (--retree 1, or --parttree
        # above parttree_threshold — see below), FastTree (no IQ-TREE /
        # ModelFinder / UFBoot), midpoint root only, no LCA, no trimAl.
        # Representative leaves are prefixed with "[repr] " in the
        # phyloXML <name> for visual identification.
        "pre_cluster_tree": {
            "enabled": False,
            # Cap on the number of leaves (post-QC sequences) in the
            # pre-cluster tree. FastTree memory scales ~linearly with leaf
            # count, so an uncapped pool of tens of thousands of leaves
            # needs hundreds of GB and gets OOM-killed — and a tree that
            # large isn't legible anyway. Above this many sequences the
            # tree is built on ALL representatives + a random sample of
            # the non-representative background, up to max_leaves (reps are
            # never dropped — they're the point of the overview). At/below
            # it the full pool is used. 0 = no cap (build everything; only
            # do this on a big-memory machine).
            "max_leaves": 5000,
            # Above this many post-QC sequences (AFTER the max_leaves cap),
            # the pre-cluster MSA switches from MAFFT --retree 1 (single-
            # pass FFT-NS-1, which builds an O(N^2) distance matrix and
            # runs out of memory on very large pools) to MAFFT --retree 2
            # --parttree (the PartTree guide — no full matrix, scales to
            # 10^5+ sequences). Below it the standard --retree 1 pass is
            # used. 0 = always PartTree; a very large value effectively
            # never. PartTree is rougher, but the pre-cluster tree is a
            # rough overview anyway. (With the default max_leaves: 5000
            # cap this never triggers; it only matters when max_leaves: 0
            # uncaps the pool.)
            "parttree_threshold": 10000,
        },
    },
    "taxonomy": {
        "ncbi_email": None,
        "ncbi_api_key": None,
        "cache_ttl_days": 90,
    },
    "output": {
        "dir": "./repseq_output",
        "prefix": "repseq",
        # {prefix}_protein_taxonomic_report.txt — protein coverage and length
        # statistics per taxonomic rank (subgenus → class). max_breakdown caps
        # the rows per rank table; extra taxa are summarised on a "+N more"
        # line so column widths stay readable on a terminal. Set higher for
        # publication runs where the full per-rank table is expected.
        "protein_report": {
            "max_breakdown": 20,
        },
        # {prefix}_polyprotein_taxonomic_report.* — the peptide-coverage
        # "wall of zeros" alarm. When a declared polyprotein spec's peptide
        # profiles fail to cover a whole clade, that clade's peptide columns
        # go 0% (the pre-fix Orthohepacivirus/Pegivirus/Orthopestivirus case
        # under the Orthoflavivirus-only slicing). This warns loudly (console
        # + _flags.txt + HTML) when it sees such a wall. Each taxon is judged
        # against the ONE spec that covers the most of its peptides (its "home"
        # spec), so the expected all-zero rows a taxon has under OTHER clades'
        # specs never false-fire.
        "polyprotein_report": {
            "wall_warning": {
                "enabled": True,
                # Rank to assess (one rank avoids redundant multi-rank alarms).
                "rank": "genus",
                # Warn when >= this fraction of the home spec's peptides are at
                # 0% coverage for the taxon.
                "wall_fraction": 0.6,
                # Ignore taxa with fewer than this many representatives (noise
                # guard — a 1-2 rep taxon isn't a systemic coverage wall).
                "min_reps": 3,
            },
        },
    },
    # Sequences of special importance. `protect_qc` force-keeps the named
    # sequences through the QC removal stages they would otherwise fail (a
    # whitelist), so a curated reference strain is never silently dropped.
    # Matching is by accession (non-segmented) / isolate_id (segmented),
    # case- and version-insensitive. `protect_stages` is "all" or a subset
    # of the QC stage tokens (duplicates, length, ambiguous, annotation,
    # protein_count, taxonomy_consistency, protein_quality, hmm). Kept
    # sequences are recorded in {prefix}_overrides.tsv and the run summary —
    # protection is never silent. Default off (empty id list = no-op).
    "overrides": {
        "ids": [],
        "ids_file": None,
        "protect_qc": False,
        "protect_stages": "all",
        # force_select: guarantee the named sequences appear among the
        # representatives (independent of protect_qc). Hybrid policy: a
        # pinned sequence wins its cluster's representative slot; pins that
        # collide in one cluster are split into singletons; diversity-
        # deselected pins are added as singletons. Force-select needs the
        # sequence to survive QC first — pair with protect_qc for a
        # "present no matter what" guarantee. Audit in
        # {prefix}_force_selected.tsv.
        "force_select": False,
        # exclude: a blocklist. Named sequences are dropped from the input
        # the moment it is read — before metadata resolution, QC, or
        # clustering — exactly as if the records had been deleted from the
        # FASTA. For sequences known to be "bad" (chimeric, mislabelled,
        # contaminated). Matching is by accession / id only, case- and
        # version-insensitive (isolate_id is not yet populated this early).
        # Has its own dedicated id list (NOT the shared `ids` above); a
        # config that lists the same id under both exclude and a keep/pin
        # capability is a hard error. Dropped records are audited in
        # {prefix}_excluded.tsv and the run summary. Default off.
        "exclude": {
            "enabled": False,
            "ids": [],
            "ids_file": None,
        },
    },
}


# ---------------------------------------------------------------------------
# Required fields for specific scenarios
# ---------------------------------------------------------------------------

SEGMENTED_VIRUS_REQUIRED_FIELDS = ["expected_segments", "segments", "isolate_regex"]


# ---------------------------------------------------------------------------
# Unknown-key audit
# ---------------------------------------------------------------------------
#
# repseq rejects config keys it doesn't recognise: a typo'd or misplaced key
# is a HARD ERROR, not a silent no-op (a bench scientist who sets a key that's
# quietly ignored has no way to discover the run didn't honour it). Two
# cooperating mechanisms:
#
#   * ``_audit_unknown_keys`` walks the config against DEFAULTS, auditing the
#     keys of every CLOSED-SCHEMA node — a node that is a NON-EMPTY dict in
#     DEFAULTS. Empty dicts ({}) / None / lists in DEFAULTS mark user-keyed
#     maps (``segmented.viruses``, ``phylo.partition.models``,
#     ``phylo.rooting.outgroup_rank``) and spec collections
#     (``cluster_protein``, ``polyprotein``, ...) whose contents are DATA, not
#     schema — those are left to the spec validators.
#   * the spec validators (``_validate_marker_entry``,
#     ``_validate_segment_markers``, ``_validate_polyprotein_list``, and the
#     segmented virus block) carry their own ``_ALLOWED_*_KEYS`` so a stray key
#     INSIDE a list-entry / virus dict (e.g. ``whole_polyprotein_tree`` on a
#     polyprotein spec — the bug this audit was built for) is caught too.
#
# Keys starting with ``_`` are always skipped: both pipeline-injected runtime
# state (``_hmm_runtime``, ...) and the sanctioned channel for a user who
# wants to stash an annotation key in their YAML.
#
# WHEN ADDING A CONFIG KEY: a key added to DEFAULTS is picked up by
# ``_audit_unknown_keys`` automatically. A key added to a SPEC dict (marker /
# segment-marker / polyprotein / peptide / virus) MUST also be added to the
# matching ``_ALLOWED_*_KEYS`` set below, or valid configs using it are
# rejected. ``tests/test_config.py`` asserts ``repseq/data/default_config.yaml``
# (the bundled reference config that ``repseq init-config`` emits) still
# validates clean as a drift guard.

# Top-level keys the CLI injects into ``cfg`` before validation that are not
# in DEFAULTS (so the audit must not flag them as unknown).
_AUDIT_IGNORE_KEYS: frozenset[str] = frozenset({"verbose"})

# Renamed/removed keys that already carry a tailored migration message; the
# generic audit skips them so the user sees the specific guidance rather than
# a duplicate "unknown key" line.
_RENAMED_CONFIG_PATHS: frozenset[str] = frozenset({"qc.length_filter"})

_ALLOWED_MARKER_KEYS: frozenset[str] = frozenset({"name", "aliases", "hmms"})
_ALLOWED_SEGMENT_MARKER_KEYS: frozenset[str] = frozenset({"aliases", "hmms"})
_ALLOWED_POLYPROTEIN_SPEC_KEYS: frozenset[str] = frozenset(
    {"name", "peptides", "cut_strategy", "motif_window_aa", "min_peptides_hit"}
)
_ALLOWED_PEPTIDE_KEYS: frozenset[str] = frozenset(
    {"name", "hmm", "hmms", "cleavage_motif"}
)
_ALLOWED_VIRUS_KEYS: frozenset[str] = frozenset({
    "expected_segments", "segments", "isolate_regex", "segment_regex",
    "segment_aliases", "cluster_protein", "segment_markers", "extra_protein",
    "expected_proteins_per_segment", "segment_lengths", "polyprotein",
})


def _did_you_mean(key: Any, allowed) -> str:
    """Return a ``' — did you mean 'X'?'`` suffix for the nearest allowed key."""
    hint = difflib.get_close_matches(
        str(key), sorted(map(str, allowed)), n=1, cutoff=0.6
    )
    return f" — did you mean '{hint[0]}'?" if hint else ""


def _spec_unknown_key_errors(
    entry: Any, allowed: frozenset[str], path: str
) -> list[str]:
    """Flag keys on a spec dict that aren't in ``allowed`` (``_`` keys skipped).

    Non-dict ``entry`` returns no errors — the caller's type check owns that.
    """
    if not isinstance(entry, dict):
        return []
    errs: list[str] = []
    for key in entry:
        if isinstance(key, str) and key.startswith("_"):
            continue
        if key not in allowed:
            errs.append(
                f"{path}: unknown key '{key}'{_did_you_mean(key, allowed)} "
                f"(allowed: {', '.join(sorted(allowed))})"
            )
    return errs


def _audit_unknown_keys(
    node: Any, default_node: Any, path: str, errors: list[str]
) -> None:
    """Append an error for every key in ``node`` absent from ``default_node``.

    Recurses only into sub-nodes that are NON-EMPTY dicts in DEFAULTS (closed
    schemas); stops at empty dicts / None / lists (user-keyed maps + spec
    collections, handled by the spec validators). ``_``-prefixed keys and the
    CLI-injected ``_AUDIT_IGNORE_KEYS`` are skipped.
    """
    if not isinstance(node, dict) or not isinstance(default_node, dict):
        return
    allowed = set(default_node)
    for key, value in node.items():
        if isinstance(key, str) and (
            key.startswith("_") or key in _AUDIT_IGNORE_KEYS
        ):
            continue
        dotted = f"{path}.{key}" if path else str(key)
        if key not in allowed:
            if dotted in _RENAMED_CONFIG_PATHS:
                continue
            errors.append(
                f"unknown config key '{dotted}'{_did_you_mean(key, allowed)}"
            )
            continue
        dv = default_node[key]
        if isinstance(dv, dict) and dv and isinstance(value, dict):
            _audit_unknown_keys(value, dv, dotted, errors)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _deep_merge(base: dict, override: dict) -> dict:
    """Recursively merge override into base, returning a new dict."""
    result = copy.deepcopy(base)
    for key, value in override.items():
        if key in result and isinstance(result[key], dict) and isinstance(value, dict):
            result[key] = _deep_merge(result[key], value)
        else:
            result[key] = copy.deepcopy(value)
    return result


def _expand_paths(cfg: dict) -> dict:
    """Expand ~ in path fields."""
    for key in ("cache_dir", "temp_dir"):
        if cfg.get(key):
            cfg[key] = str(Path(cfg[key]).expanduser())
    if cfg.get("output", {}).get("dir"):
        cfg["output"]["dir"] = str(Path(cfg["output"]["dir"]).expanduser())
    return cfg


# Config keys whose VALUES are secrets and must never be persisted into
# an output artifact — blanked to ``None`` wherever they appear (any
# depth), by name. Currently the NCBI credentials (also injectable via
# the REPSEQ_NCBI_EMAIL / REPSEQ_NCBI_API_KEY env vars).
SECRET_CONFIG_KEYS: frozenset[str] = frozenset({"ncbi_email", "ncbi_api_key"})


def sanitize_config(cfg: dict[str, Any], *, drop_private: bool = True) -> dict[str, Any]:
    """Return a copy of ``cfg`` safe to write to disk.

    * **Secrets blanked** — any key in :data:`SECRET_CONFIG_KEYS`, at any
      depth, is set to ``None`` (the key stays so the structure is
      complete and the file remains a valid, re-runnable config; the
      user fills in their own credential).
    * **Private runtime keys dropped** (when ``drop_private``) — any
      top-level-or-nested key starting with ``_`` (e.g. ``_hmm_runtime``,
      ``_taxonomy_review``) is removed. These are pipeline-injected
      runtime state, not configuration, and are recomputed on every run.

    Pure / non-mutating: builds fresh dicts and lists, leaves scalars
    (immutable) as-is, and coerces tuples to lists so the result is
    YAML/JSON-native. ``drop_private=False`` keeps the ``_`` keys (used by
    the lockfile, which records the full post-mutation runtime config for
    replay but still must not persist credentials)."""
    def _clean(obj: Any) -> Any:
        if isinstance(obj, dict):
            out: dict[Any, Any] = {}
            for key, value in obj.items():
                if drop_private and isinstance(key, str) and key.startswith("_"):
                    continue
                if isinstance(key, str) and key in SECRET_CONFIG_KEYS:
                    out[key] = None
                else:
                    out[key] = _clean(value)
            return out
        if isinstance(obj, (list, tuple)):
            return [_clean(x) for x in obj]
        return obj

    return _clean(cfg)


def write_effective_config(cfg: dict[str, Any], path: Path) -> Path:
    """Write the sanitized, fully-resolved effective config to ``path``.

    The dumped config is the run-time ``cfg`` — i.e. the user's YAML
    already deep-merged over :data:`DEFAULTS`, so **every setting is
    present at the value it actually ran with** (defaults filled in for
    anything the user didn't set). Secrets are blanked and runtime ``_``
    keys dropped (see :func:`sanitize_config`); ``yaml.safe_dump`` emits
    no comments. Block style + insertion (DEFAULTS) ordering are kept so
    the file is human-readable and re-loadable verbatim via
    ``repseq <mode> -c <this-file>``.
    """
    sanitized = sanitize_config(cfg, drop_private=True)
    body = yaml.safe_dump(
        sanitized,
        sort_keys=False,
        default_flow_style=False,
        allow_unicode=True,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(body)
    return path


def effective_config_filename(prefix: str) -> str:
    """Filename for the effective-config snapshot:
    ``{prefix}_config_repseq{version}.yaml`` with the repseq version's
    periods replaced by underscores (e.g. ``cov_config_repseq0_42_0.yaml``)."""
    from . import __version__
    return f"{prefix}_config_repseq{__version__.replace('.', '_')}.yaml"


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def _format_yaml_error(path: Path, err: yaml.YAMLError) -> str:
    """Turn a PyYAML parse error into a one-glance, line-located message."""
    mark = getattr(err, "problem_mark", None)
    loc = f" at line {mark.line + 1}, column {mark.column + 1}" if mark else ""
    problem = getattr(err, "problem", None) or str(err).splitlines()[0]
    return (
        f"Config file {path} is not valid YAML{loc}: {problem}.\n"
        f"       Check the indentation (use spaces, not tabs) and quoting."
    )


def load_config(path: Optional[str | Path] = None) -> dict[str, Any]:
    """Load config from YAML file, merging over defaults.

    Raises :class:`~repseq.errors.ConfigError` (rendered without a traceback
    at the CLI boundary) for the mistakes a user is likely to make: a path
    that doesn't exist, a file that isn't valid YAML, or a top level that
    isn't a mapping.
    """
    cfg = copy.deepcopy(DEFAULTS)
    if path is not None:
        path = Path(path)
        if not path.exists():
            raise ConfigError(
                f"Config file not found: {path}\n"
                f"       Check the path you passed to -c/--config."
            )
        try:
            with open(path) as fh:
                user_cfg = yaml.safe_load(fh)
        except yaml.YAMLError as e:
            raise ConfigError(_format_yaml_error(path, e)) from e
        except OSError as e:
            raise ConfigError(
                f"Could not read config file {path}: {e.strerror or e}."
            ) from e
        if user_cfg is None:
            user_cfg = {}
        if not isinstance(user_cfg, dict):
            raise ConfigError(
                f"Config file {path} must contain a mapping (key: value pairs) "
                f"at the top level, but found a {type(user_cfg).__name__}.\n"
                f"       Check the file's structure — a repseq config is a set "
                f"of named sections, not a bare list or value."
            )
        cfg = _deep_merge(cfg, user_cfg)

    # Environment variable overrides
    if os.environ.get("REPSEQ_NCBI_EMAIL"):
        cfg["taxonomy"]["ncbi_email"] = os.environ["REPSEQ_NCBI_EMAIL"]
    if os.environ.get("REPSEQ_NCBI_API_KEY"):
        cfg["taxonomy"]["ncbi_api_key"] = os.environ["REPSEQ_NCBI_API_KEY"]

    cfg = _expand_paths(cfg)
    return cfg


def _validate_hmm_tokens(hmms: Any, path: str) -> tuple[list[str], list[str]]:
    """Validate an ``hmms:`` list of token strings.

    Each token is either a single HMM name (``"Name"``) or a multidomain
    spec joined with ``--`` (``"A--B--C"``, HMMs listed in N-to-C order).
    Multiple tokens in one list are **alternative architectures (OR)** — a
    CDS satisfying any one of them satisfies the marker. Returns
    ``(errors, validated_tokens)``. Invalid tokens are dropped from
    ``validated_tokens`` so the caller can still check the "at least one of
    aliases / hmms" invariant.
    """
    from .hmm.runner import parse_hmm_token

    errs: list[str] = []
    validated: list[str] = []
    if not isinstance(hmms, list):
        errs.append(
            f"{path} must be a list of HMM token strings "
            "(single 'Name' or multidomain 'A--B--C')"
        )
        return errs, validated
    for i, token in enumerate(hmms):
        if not isinstance(token, str):
            errs.append(f"{path}[{i}] must be a string, got {type(token).__name__}")
            continue
        try:
            parse_hmm_token(token)
        except ValueError as e:
            errs.append(f"{path}[{i}]: {e}")
            continue
        validated.append(token)
    return errs, validated


def _validate_marker_entry(entry: Any, path: str) -> list[str]:
    """Validate one cluster_protein marker entry.

    Each entry is either a non-empty alias string (legacy, alias-only)
    or a dict with required ``name`` and at least one of ``aliases``
    (list of non-empty strings) or ``hmms`` (list of HMM token strings,
    where a token is either a single HMM name like ``"Name"`` or a
    multidomain spec like ``"A--B"`` in N-to-C order; multiple tokens are
    alternative architectures, OR).
    """
    errs: list[str] = []
    if isinstance(entry, str):
        if not entry.strip():
            errs.append(f"{path}: alias string must be non-empty")
        return errs
    if not isinstance(entry, dict):
        errs.append(
            f"{path} must be an alias string or a dict "
            "{name, aliases?, hmms?}"
        )
        return errs
    errs.extend(_spec_unknown_key_errors(entry, _ALLOWED_MARKER_KEYS, path))
    name = entry.get("name")
    if not isinstance(name, str) or not name.strip():
        errs.append(f"{path}: dict-form entry must include a non-empty 'name'")
    aliases = entry.get("aliases", [])
    if not isinstance(aliases, list) or not all(
        isinstance(a, str) and a.strip() for a in aliases
    ):
        errs.append(f"{path}: 'aliases' must be a list of non-empty strings")
        aliases = []
    token_errs, hmms = _validate_hmm_tokens(entry.get("hmms", []), f"{path}.hmms")
    errs.extend(token_errs)
    if not aliases and not hmms:
        errs.append(
            f"{path}: dict-form entry must define at least one of "
            "'aliases' or 'hmms' (otherwise the marker can't be matched)"
        )
    return errs


def _validate_segment_markers(
    sm: Any, virus_name: str, seg_names: set[str]
) -> list[str]:
    """Validate per-virus segment_markers block.

    Shape: ``{segment_name: {aliases: [...], hmms: [...]}}``. Each
    segment-spec must define at least one of aliases / hmms. Coexists
    with the legacy per-segment ``cluster_protein`` block; when both
    define a marker for the same segment, ``segment_markers`` wins.
    """
    errs: list[str] = []
    if not isinstance(sm, dict):
        errs.append(
            f"segmented.viruses.{virus_name}.segment_markers must be a "
            "mapping of segment-name → {aliases: [...], hmms: [...]}"
        )
        return errs
    for seg_name, spec in sm.items():
        prefix = f"segmented.viruses.{virus_name}.segment_markers.{seg_name}"
        if seg_name not in seg_names:
            errs.append(
                f"segmented.viruses.{virus_name}.segment_markers: "
                f"unknown segment '{seg_name}'"
            )
        if not isinstance(spec, dict):
            errs.append(
                f"{prefix} must be a dict with 'aliases' and/or 'hmms' keys"
            )
            continue
        errs.extend(
            _spec_unknown_key_errors(spec, _ALLOWED_SEGMENT_MARKER_KEYS, prefix)
        )
        aliases = spec.get("aliases", [])
        if not isinstance(aliases, list) or not all(
            isinstance(a, str) and a.strip() for a in aliases
        ):
            errs.append(f"{prefix}.aliases must be a list of non-empty strings")
            aliases = []
        token_errs, hmms = _validate_hmm_tokens(spec.get("hmms", []), f"{prefix}.hmms")
        errs.extend(token_errs)
        if not aliases and not hmms:
            errs.append(
                f"{prefix}: must define at least one of 'aliases' or 'hmms'"
            )
    return errs


_POLYPROTEIN_CUT_STRATEGIES = ("boundary", "bisect", "motif")


def _validate_polyprotein_list(entries: Any, path: str) -> list[str]:
    """Validate a list of polyprotein specs (clustering or per-segment).

    Each entry must be a dict with:
      * ``name`` (non-empty string; used as the output basename component)
      * ``peptides`` (list, length ≥ 2, ordered N→C). Each peptide is a dict
        with ``name`` and ``hmm`` (both non-empty strings); optional
        ``cleavage_motif`` (non-empty string, residues just N-terminal of
        the cut that liberates this peptide).
      * Optional ``cut_strategy`` in {"boundary", "bisect", "motif"};
        defaults to "motif" if any peptide declares ``cleavage_motif``,
        else "bisect".
      * Optional ``motif_window_aa`` (positive int; default 10).
      * Optional ``min_peptides_hit`` (positive int; default 2; clamped at
        ≥ 1 since 0 would accept any CDS).
    Spec names must be unique within the same list (the caller is
    responsible for cross-list uniqueness — handled separately in the
    segmented validator since names key the output filenames globally).
    """
    errs: list[str] = []
    if not isinstance(entries, list):
        errs.append(
            f"{path} must be a list of {{name, peptides, ...}} dicts"
        )
        return errs

    seen_names: set[str] = set()
    for i, entry in enumerate(entries):
        ipath = f"{path}[{i}]"
        if not isinstance(entry, dict):
            errs.append(
                f"{ipath} must be a dict "
                f"{{name, peptides, cut_strategy?, ...}}"
            )
            continue
        errs.extend(_spec_unknown_key_errors(
            entry, _ALLOWED_POLYPROTEIN_SPEC_KEYS, ipath
        ))

        name = entry.get("name")
        if not isinstance(name, str) or not name.strip():
            errs.append(f"{ipath} must include a non-empty 'name'")
        else:
            key = name.strip().lower()
            if key in seen_names:
                errs.append(
                    f"{ipath}: duplicate spec name '{name}' "
                    f"(names key the output filenames)"
                )
            seen_names.add(key)

        peptides = entry.get("peptides")
        if not isinstance(peptides, list) or len(peptides) < 2:
            errs.append(
                f"{ipath}.peptides must be a list of at least 2 peptide "
                f"dicts ordered N-to-C"
            )
        else:
            seen_pep: set[str] = set()
            for j, pep in enumerate(peptides):
                ppath = f"{ipath}.peptides[{j}]"
                if not isinstance(pep, dict):
                    errs.append(
                        f"{ppath} must be a dict "
                        f"{{name, hmm, cleavage_motif?}}"
                    )
                    continue
                errs.extend(_spec_unknown_key_errors(
                    pep, _ALLOWED_PEPTIDE_KEYS, ppath
                ))
                pname = pep.get("name")
                if not isinstance(pname, str) or not pname.strip():
                    errs.append(f"{ppath}.name must be a non-empty string")
                else:
                    pkey = pname.strip().lower()
                    if pkey in seen_pep:
                        errs.append(
                            f"{ppath}: duplicate peptide name '{pname}' "
                            f"within this polyprotein"
                        )
                    seen_pep.add(pkey)

                # Accept either ``hmm: <token>`` (singular, legacy) or
                # ``hmms: [<token>, ...]`` (the v0.34.0 OR form).
                # Exactly one must be set; both is an error
                # (ambiguous), neither is an error (peptide has no
                # locator).
                phmm = pep.get("hmm")
                phmms = pep.get("hmms")
                has_singular = phmm is not None
                has_plural = phmms is not None
                if has_singular and has_plural:
                    errs.append(
                        f"{ppath}: set exactly one of 'hmm' (single "
                        f"token) or 'hmms' (list of alternative "
                        f"tokens), not both"
                    )
                elif not has_singular and not has_plural:
                    errs.append(
                        f"{ppath}: must set either 'hmm' (single HMM "
                        f"token like \"CoV_NSP8\") or 'hmms' (a list of "
                        f"alternative tokens, OR-joined — e.g. "
                        f"[aCoV_NSP1, bCoV_NSP1] for a peptide whose "
                        f"architecture differs across genera)"
                    )
                else:
                    from .hmm.runner import parse_hmm_token
                    if has_singular:
                        if not isinstance(phmm, str) or not phmm.strip():
                            errs.append(
                                f"{ppath}.hmm must be a non-empty HMM "
                                f"token (single profile name or "
                                f"multidomain \"A--B--C\")"
                            )
                        else:
                            try:
                                parse_hmm_token(phmm)
                            except ValueError as exc:
                                errs.append(f"{ppath}.hmm: {exc}")
                    else:
                        if (
                            not isinstance(phmms, list)
                            or not phmms
                            or not all(
                                isinstance(t, str) and t.strip()
                                for t in phmms
                            )
                        ):
                            errs.append(
                                f"{ppath}.hmms must be a non-empty list "
                                f"of HMM token strings (each either a "
                                f"single profile name or a multidomain "
                                f"\"A--B--C\" token)"
                            )
                        else:
                            for ti, tok in enumerate(phmms):
                                try:
                                    parse_hmm_token(tok)
                                except ValueError as exc:
                                    errs.append(
                                        f"{ppath}.hmms[{ti}]: {exc}"
                                    )
                motif = pep.get("cleavage_motif")
                if motif is not None and (
                    not isinstance(motif, str) or not motif.strip()
                ):
                    errs.append(
                        f"{ppath}.cleavage_motif must be a non-empty "
                        f"string (residues N-terminal of the cut), or omitted"
                    )

        strat = entry.get("cut_strategy")
        if strat is not None and strat not in _POLYPROTEIN_CUT_STRATEGIES:
            errs.append(
                f"{ipath}.cut_strategy '{strat}' is not supported "
                f"(use one of {list(_POLYPROTEIN_CUT_STRATEGIES)})"
            )

        mwin = entry.get("motif_window_aa", 10)
        if (
            not isinstance(mwin, int)
            or isinstance(mwin, bool)
            or mwin < 1
        ):
            errs.append(
                f"{ipath}.motif_window_aa must be a positive integer "
                f"(amino acids; default 10)"
            )

        mph = entry.get("min_peptides_hit", 2)
        if (
            not isinstance(mph, int)
            or isinstance(mph, bool)
            or mph < 1
        ):
            errs.append(
                f"{ipath}.min_peptides_hit must be a positive integer "
                f"(parent-CDS identification threshold; default 2)"
            )

    return errs


def validate_config(cfg: dict[str, Any]) -> list[str]:
    """Return a list of validation error messages (empty = valid)."""
    errors: list[str] = []

    # Reject unrecognised keys up front: a typo'd or misplaced key is a hard
    # error, not a silent no-op. Closed-schema sections (non-empty dicts in
    # DEFAULTS) are audited here; spec dicts (markers, peptides, polyprotein,
    # the segmented virus block) are audited by the spec validators below.
    _audit_unknown_keys(cfg, DEFAULTS, "", errors)

    # Genome length filter (non-segmented only). The old median-percent
    # filter (qc.length_filter) was removed in favour of explicit absolute
    # bounds — reject the renamed key so configs migrate consciously rather
    # than silently losing their length filtering.
    qc = cfg.get("qc", {}) or {}
    if "length_filter" in qc:
        errors.append(
            "qc.length_filter was renamed to qc.genome_length_filter and the "
            "median-percent mode was removed. Use "
            "qc.genome_length_filter: {enabled: true, min: <nt>, max: <nt>} "
            "with absolute nucleotide bounds (non-segmented mode only)."
        )
    glf = qc.get("genome_length_filter", {}) or {}
    if "enabled" in glf and not isinstance(glf.get("enabled"), bool):
        errors.append("qc.genome_length_filter.enabled must be a boolean")
    if glf.get("enabled"):
        if cfg.get("segmented", {}).get("enabled"):
            errors.append(
                "qc.genome_length_filter.enabled cannot be true when "
                "segmented.enabled is true — the whole-genome length filter "
                "only applies to non-segmented runs. For segmented viruses, "
                "set per-segment bounds via "
                "segmented.viruses.<virus>.segment_lengths instead."
            )
        mn = glf.get("min")
        mx = glf.get("max")
        if mn is None and mx is None:
            errors.append(
                "qc.genome_length_filter.enabled is true but neither min nor "
                "max is set — supply at least one absolute nucleotide bound."
            )
        for label, val in (("min", mn), ("max", mx)):
            if val is not None and (not isinstance(val, int) or isinstance(val, bool) or val <= 0):
                errors.append(
                    f"qc.genome_length_filter.{label} must be a positive integer "
                    "(nucleotides) or null"
                )
        if (
            isinstance(mn, int) and not isinstance(mn, bool)
            and isinstance(mx, int) and not isinstance(mx, bool)
            and mn > mx
        ):
            errors.append(
                "qc.genome_length_filter.min must be <= qc.genome_length_filter.max"
            )

    # Ambiguous threshold
    thresh = cfg.get("qc", {}).get("ambiguous_threshold")
    if not isinstance(thresh, (int, float)) or not (0 <= thresh <= 1):
        errors.append("qc.ambiguous_threshold must be a number between 0 and 1")

    # Protein annotation QC
    pa = cfg.get("qc", {}).get("protein_annotation", {})
    if pa.get("enabled"):
        mp = pa.get("min_proteins")
        if not isinstance(mp, int) or mp < 0:
            errors.append("qc.protein_annotation.min_proteins must be a non-negative integer")

    # Protein quality QC
    pq = cfg.get("qc", {}).get("protein_quality", {})
    if pq.get("enabled"):
        mbf = pq.get("max_bad_fraction")
        if not isinstance(mbf, (int, float)) or not (0 <= mbf <= 1):
            errors.append(
                "qc.protein_quality.max_bad_fraction must be a number between 0 and 1"
            )

    # Segmented virus
    seg = cfg.get("segmented", {})
    if "use_genbank_metadata" in seg and not isinstance(
        seg["use_genbank_metadata"], bool
    ):
        errors.append("segmented.use_genbank_metadata must be a boolean")

    sca = seg.get("strain_collision_action", "warn")
    if sca not in ("warn", "drop"):
        errors.append(
            f"segmented.strain_collision_action '{sca}' is not supported "
            "(use 'warn' or 'drop')"
        )

    esa = seg.get("extra_segments_action", "warn")
    if esa not in ("warn", "drop"):
        errors.append(
            f"segmented.extra_segments_action '{esa}' is not supported "
            "(use 'warn' or 'drop')"
        )

    tc = seg.get("taxonomy_consistency", {}) or {}
    if "enabled" in tc and not isinstance(tc["enabled"], bool):
        errors.append("segmented.taxonomy_consistency.enabled must be a boolean")
    # The rank must be one a resolved TaxonomyInfo can actually answer.
    # Standard fields on TaxonomyInfo plus arbitrary entries from
    # ``lineage`` are both accepted by ``get_rank``, so we only
    # validate against the named-attribute set — anything else still
    # works, just relies on the lineage map being populated.
    valid_consistency_ranks = {
        "species", "subgenus", "genus", "subfamily", "family",
        "suborder", "order", "subclass", "class", "phylum",
        "kingdom", "superkingdom",
    }
    rank = tc.get("rank", "species")
    if not isinstance(rank, str) or rank.lower() not in valid_consistency_ranks:
        errors.append(
            f"segmented.taxonomy_consistency.rank '{rank}' is not supported "
            f"(use one of {sorted(valid_consistency_ranks)})"
        )
    if seg.get("enabled"):
        virus_name = seg.get("virus")
        if not virus_name:
            errors.append("segmented.virus must be set when segmented.enabled is true")
        else:
            viruses = seg.get("viruses", {})
            if virus_name not in viruses:
                errors.append(
                    f"segmented.virus '{virus_name}' not found in segmented.viruses"
                )
            else:
                vdef = viruses[virus_name]
                errors.extend(_spec_unknown_key_errors(
                    vdef, _ALLOWED_VIRUS_KEYS,
                    f"segmented.viruses.{virus_name}",
                ))
                for field in SEGMENTED_VIRUS_REQUIRED_FIELDS:
                    if field not in vdef:
                        errors.append(
                            f"segmented.viruses.{virus_name} missing required field '{field}'"
                        )
                # Segment aliases: optional dict[canonical → list[str]]
                aliases = vdef.get("segment_aliases")
                if aliases is not None:
                    if not isinstance(aliases, dict):
                        errors.append(
                            f"segmented.viruses.{virus_name}.segment_aliases "
                            f"must be a mapping of canonical-segment-name → list of strings"
                        )
                    else:
                        seg_names = set(vdef.get("segments", []))
                        for canonical, syns in aliases.items():
                            if canonical not in seg_names:
                                errors.append(
                                    f"segmented.viruses.{virus_name}."
                                    f"segment_aliases: unknown segment '{canonical}'"
                                )
                            if not isinstance(syns, list) or not all(
                                isinstance(s, str) and s.strip() for s in syns
                            ):
                                errors.append(
                                    f"segmented.viruses.{virus_name}."
                                    f"segment_aliases.{canonical} "
                                    f"must be a list of non-empty strings"
                                )

                cp = vdef.get("cluster_protein")
                if cp is not None:
                    if not isinstance(cp, dict):
                        errors.append(
                            f"segmented.viruses.{virus_name}.cluster_protein "
                            f"must be a mapping of segment-name → list of "
                            f"alias strings and/or {{name, aliases?, hmms?}} dicts"
                        )
                    else:
                        seg_names = set(vdef.get("segments", []))
                        for seg_name, entries in cp.items():
                            if seg_name not in seg_names:
                                errors.append(
                                    f"segmented.viruses.{virus_name}."
                                    f"cluster_protein: unknown segment '{seg_name}'"
                                )
                            if not isinstance(entries, list):
                                errors.append(
                                    f"segmented.viruses.{virus_name}."
                                    f"cluster_protein.{seg_name} "
                                    f"must be a list of alias strings and/or "
                                    f"{{name, aliases?, hmms?}} dicts"
                                )
                            else:
                                for j, entry in enumerate(entries):
                                    errors.extend(_validate_marker_entry(
                                        entry,
                                        f"segmented.viruses.{virus_name}."
                                        f"cluster_protein.{seg_name}[{j}]",
                                    ))

                sm = vdef.get("segment_markers")
                if sm is not None:
                    errors.extend(_validate_segment_markers(
                        sm, virus_name, set(vdef.get("segments", []))
                    ))

                ep = vdef.get("extra_protein")
                if ep is not None:
                    if not isinstance(ep, dict):
                        errors.append(
                            f"segmented.viruses.{virus_name}.extra_protein "
                            f"must be a mapping of segment-name → list of "
                            f"{{name, aliases?, hmms?}} dicts"
                        )
                    else:
                        seg_names = set(vdef.get("segments", []))
                        seen_extra: set[str] = set()
                        for seg_name, entries in ep.items():
                            if seg_name not in seg_names:
                                errors.append(
                                    f"segmented.viruses.{virus_name}."
                                    f"extra_protein: unknown segment "
                                    f"'{seg_name}'"
                                )
                            if not isinstance(entries, list):
                                errors.append(
                                    f"segmented.viruses.{virus_name}."
                                    f"extra_protein.{seg_name} must be a "
                                    f"list of {{name, aliases?, hmms?}} dicts"
                                )
                                continue
                            for j, entry in enumerate(entries):
                                ipath = (
                                    f"segmented.viruses.{virus_name}."
                                    f"extra_protein.{seg_name}[{j}]"
                                )
                                if isinstance(entry, str):
                                    errors.append(
                                        f"{ipath} must be a dict "
                                        f"{{name, aliases?, hmms?}}, not a "
                                        f"bare alias string"
                                    )
                                    continue
                                errors.extend(
                                    _validate_marker_entry(entry, ipath)
                                )
                                if isinstance(entry, dict):
                                    nm = (entry.get("name") or "").strip()
                                    if nm:
                                        # Names must be unique across the
                                        # whole virus (not just within one
                                        # segment) — the output basename is
                                        # the name alone, no segment prefix.
                                        key = nm.lower()
                                        if key in seen_extra:
                                            errors.append(
                                                f"{ipath}: duplicate name "
                                                f"'{nm}' — extra_protein "
                                                f"names must be unique "
                                                f"across all segments "
                                                f"(they key the output "
                                                f"filenames)"
                                            )
                                        seen_extra.add(key)

                epps = vdef.get("expected_proteins_per_segment")
                if epps is not None:
                    if not isinstance(epps, dict):
                        errors.append(
                            f"segmented.viruses.{virus_name}.expected_proteins_per_segment "
                            f"must be a mapping of segment-name → int or list[int]"
                        )
                    else:
                        seg_names = set(vdef.get("segments", []))
                        for seg_name, count in epps.items():
                            if seg_name not in seg_names:
                                errors.append(
                                    f"segmented.viruses.{virus_name}."
                                    f"expected_proteins_per_segment: unknown segment '{seg_name}'"
                                )
                            # Allow int (exact) or list[int] (any-of).
                            # bool is a subclass of int in Python — reject it explicitly.
                            valid = (
                                (isinstance(count, int) and not isinstance(count, bool) and count >= 0)
                                or (
                                    isinstance(count, list)
                                    and len(count) > 0
                                    and all(
                                        isinstance(x, int)
                                        and not isinstance(x, bool)
                                        and x >= 0
                                        for x in count
                                    )
                                )
                            )
                            if not valid:
                                errors.append(
                                    f"segmented.viruses.{virus_name}."
                                    f"expected_proteins_per_segment.{seg_name} "
                                    f"must be a non-negative integer or a non-empty "
                                    f"list of non-negative integers"
                                )

                sl = vdef.get("segment_lengths")
                if sl is not None:
                    if not isinstance(sl, dict):
                        errors.append(
                            f"segmented.viruses.{virus_name}.segment_lengths "
                            f"must be a mapping of segment-name → {{min: N, max: M}}"
                        )
                    else:
                        seg_names = set(vdef.get("segments", []))
                        for seg_name, bounds in sl.items():
                            if seg_name not in seg_names:
                                errors.append(
                                    f"segmented.viruses.{virus_name}."
                                    f"segment_lengths: unknown segment '{seg_name}'"
                                )
                            if not isinstance(bounds, dict):
                                errors.append(
                                    f"segmented.viruses.{virus_name}."
                                    f"segment_lengths.{seg_name} must be a dict "
                                    f"with optional 'min' and/or 'max' integer keys"
                                )
                            else:
                                mn = bounds.get("min")
                                mx = bounds.get("max")
                                if mn is not None and (
                                    not isinstance(mn, int) or isinstance(mn, bool) or mn < 0
                                ):
                                    errors.append(
                                        f"segmented.viruses.{virus_name}."
                                        f"segment_lengths.{seg_name}.min must be a "
                                        f"non-negative integer"
                                    )
                                if mx is not None and (
                                    not isinstance(mx, int) or isinstance(mx, bool) or mx < 0
                                ):
                                    errors.append(
                                        f"segmented.viruses.{virus_name}."
                                        f"segment_lengths.{seg_name}.max must be a "
                                        f"non-negative integer"
                                    )
                                if (
                                    mn is not None and mx is not None
                                    and isinstance(mn, int) and isinstance(mx, int)
                                    and mn >= mx
                                ):
                                    errors.append(
                                        f"segmented.viruses.{virus_name}."
                                        f"segment_lengths.{seg_name}: min must be "
                                        f"less than max"
                                    )

    # Clustering backend
    backend = cfg.get("clustering", {}).get("backend")
    if backend not in ("mmseqs2", "cdhit"):
        errors.append(
            f"clustering.backend '{backend}' is not supported "
            f"(use 'mmseqs2' or 'cdhit')"
        )

    alphabet = cfg.get("clustering", {}).get("alphabet_for_clustering", "protein")
    if alphabet not in ("protein", "nucleotide"):
        errors.append(
            f"clustering.alphabet_for_clustering '{alphabet}' is not supported "
            f"(use 'protein' or 'nucleotide')"
        )

    cluster_protein_global = cfg.get("clustering", {}).get("cluster_protein", [])
    if not isinstance(cluster_protein_global, list):
        errors.append(
            "clustering.cluster_protein must be a list of alias strings "
            "and/or {name, aliases?, hmms?} dicts"
        )
    else:
        for i, entry in enumerate(cluster_protein_global):
            errors.extend(
                _validate_marker_entry(entry, f"clustering.cluster_protein[{i}]")
            )

    extra_protein_global = cfg.get("clustering", {}).get("extra_protein", [])
    if not isinstance(extra_protein_global, list):
        errors.append(
            "clustering.extra_protein must be a list of "
            "{name, aliases?, hmms?} dicts"
        )
    else:
        seen_extra_names: set[str] = set()
        for i, entry in enumerate(extra_protein_global):
            # extra_protein entries MUST be dicts with a name — they drive
            # output filenames, so a bare alias string (which has no name)
            # has nothing to key the artifact on.
            if isinstance(entry, str):
                errors.append(
                    f"clustering.extra_protein[{i}] must be a dict "
                    f"{{name, aliases?, hmms?}}, not a bare alias string"
                )
                continue
            errors.extend(
                _validate_marker_entry(entry, f"clustering.extra_protein[{i}]")
            )
            if isinstance(entry, dict):
                name = (entry.get("name") or "").strip()
                if name:
                    key = name.lower()
                    if key in seen_extra_names:
                        errors.append(
                            f"clustering.extra_protein[{i}]: duplicate name "
                            f"'{name}' — extra_protein names must be unique "
                            f"(they key the output filenames)"
                        )
                    seen_extra_names.add(key)

    # Polyprotein cutting specs (non-segmented at clustering.polyprotein;
    # segmented at virus.polyprotein, validated alongside the virus block
    # above).
    poly_global = cfg.get("clustering", {}).get("polyprotein", [])
    if not isinstance(poly_global, list):
        errors.append(
            "clustering.polyprotein must be a list of "
            "{name, peptides, cut_strategy?, motif_window_aa?, "
            "min_peptides_hit?} dicts"
        )
    else:
        errors.extend(_validate_polyprotein_list(
            poly_global, "clustering.polyprotein",
        ))

    # Segmented per-segment polyprotein blocks (parity with extra_protein).
    if seg.get("enabled"):
        virus_name = seg.get("virus")
        viruses = seg.get("viruses", {}) or {}
        if virus_name and virus_name in viruses:
            vdef = viruses[virus_name]
            pp = vdef.get("polyprotein")
            if pp is not None:
                if not isinstance(pp, dict):
                    errors.append(
                        f"segmented.viruses.{virus_name}.polyprotein "
                        f"must be a mapping of segment-name → list of "
                        f"{{name, peptides, ...}} dicts"
                    )
                else:
                    seg_names = set(vdef.get("segments", []))
                    flat: list[Any] = []
                    for seg_name, entries in pp.items():
                        if seg_name not in seg_names:
                            errors.append(
                                f"segmented.viruses.{virus_name}."
                                f"polyprotein: unknown segment '{seg_name}'"
                            )
                        if not isinstance(entries, list):
                            errors.append(
                                f"segmented.viruses.{virus_name}."
                                f"polyprotein.{seg_name} must be a list of "
                                f"{{name, peptides, ...}} dicts"
                            )
                            continue
                        errors.extend(_validate_polyprotein_list(
                            entries,
                            f"segmented.viruses.{virus_name}."
                            f"polyprotein.{seg_name}",
                        ))
                        flat.extend(entries)
                    # Spec names must be unique across all segments — they
                    # drive the output filenames.
                    seen: set[str] = set()
                    for entry in flat:
                        if isinstance(entry, dict):
                            nm = (entry.get("name") or "").strip().lower()
                            if nm and nm in seen:
                                errors.append(
                                    f"segmented.viruses.{virus_name}."
                                    f"polyprotein: duplicate spec name "
                                    f"'{nm}' across segments — names must be "
                                    f"unique (they key the output filenames)"
                                )
                            if nm:
                                seen.add(nm)

    diversity_cutoffs = cfg.get("clustering", {}).get("diversity_curve_cutoffs", [])
    if not isinstance(diversity_cutoffs, list):
        errors.append(
            "clustering.diversity_curve_cutoffs must be a list of floats in (0, 1]"
        )
    else:
        for v in diversity_cutoffs:
            if not isinstance(v, (int, float)) or not (0.0 < float(v) <= 1.0):
                errors.append(
                    f"clustering.diversity_curve_cutoffs entry {v!r} is not a "
                    f"number in (0, 1]"
                )
                break

    mmseqs2_mode = cfg.get("clustering", {}).get("mmseqs2_mode")
    if mmseqs2_mode not in ("easy-linclust", "easy-cluster"):
        errors.append(
            f"clustering.mmseqs2_mode must be 'easy-linclust' or 'easy-cluster', got '{mmseqs2_mode}'"
        )

    cdhit_cfg = cfg.get("clustering", {}).get("cdhit", {}) or {}
    if cdhit_cfg:
        ws = cdhit_cfg.get("word_size")
        if ws is not None and (
            not isinstance(ws, int) or isinstance(ws, bool) or not (2 <= ws <= 11)
        ):
            errors.append(
                "clustering.cdhit.word_size must be an integer in [2, 11] or null "
                "(auto-pick from threshold)"
            )
        cov = cdhit_cfg.get("coverage", 0.8)
        if not isinstance(cov, (int, float)) or not (0 <= cov <= 1):
            errors.append("clustering.cdhit.coverage must be a number between 0 and 1")
        mem = cdhit_cfg.get("memory_mb", 0)
        if not isinstance(mem, int) or isinstance(mem, bool) or mem < 0:
            errors.append(
                "clustering.cdhit.memory_mb must be a non-negative integer "
                "(0 = unlimited)"
            )
        for flag in ("global_alignment", "accurate"):
            if flag in cdhit_cfg and not isinstance(cdhit_cfg[flag], bool):
                errors.append(f"clustering.cdhit.{flag} must be a boolean")
        extra = cdhit_cfg.get("extra_args", [])
        if not isinstance(extra, list) or not all(isinstance(x, str) for x in extra):
            errors.append("clustering.cdhit.extra_args must be a list of strings")

    # Phylo
    phylo_cfg = cfg.get("phylo", {}) or {}
    tool = phylo_cfg.get("tool", "auto")
    if tool not in ("auto", "iqtree", "fasttree"):
        errors.append(
            f"phylo.tool '{tool}' is not supported "
            f"(use 'auto', 'iqtree', or 'fasttree')"
        )
    for tool_name in ("mafft", "fasttree", "iqtree"):
        tool_cfg = phylo_cfg.get(tool_name, {}) or {}
        extra = tool_cfg.get("extra_args", [])
        if not isinstance(extra, list) or not all(isinstance(x, str) for x in extra):
            errors.append(f"phylo.{tool_name}.extra_args must be a list of strings")

    mafft_cfg = phylo_cfg.get("mafft", {}) or {}
    if "use_auto" in mafft_cfg and not isinstance(mafft_cfg["use_auto"], bool):
        errors.append("phylo.mafft.use_auto must be a boolean")

    iq_cfg = phylo_cfg.get("iqtree", {}) or {}
    if "model" in iq_cfg and not isinstance(iq_cfg["model"], str):
        errors.append("phylo.iqtree.model must be a string (e.g. 'MFP', 'LG+G4')")
    if "binary" in iq_cfg and iq_cfg["binary"] is not None and not isinstance(
        iq_cfg["binary"], str
    ):
        errors.append("phylo.iqtree.binary must be a string or null")
    ufb = iq_cfg.get("ultrafast_bootstrap", 0)
    if not isinstance(ufb, int) or isinstance(ufb, bool) or ufb < 0:
        errors.append(
            "phylo.iqtree.ultrafast_bootstrap must be a non-negative integer "
            "(0 disables; IQ-TREE recommends >= 1000 for interpretable support)"
        )

    labeling_cfg = phylo_cfg.get("labeling", {}) or {}
    for key in ("format", "segmented_format"):
        if key in labeling_cfg and labeling_cfg[key] is not None and not isinstance(
            labeling_cfg[key], str
        ):
            errors.append(f"phylo.labeling.{key} must be a string or null")
    for key in ("replace_whitespace", "keep_separator_on_empty"):
        if key in labeling_cfg and not isinstance(labeling_cfg[key], bool):
            errors.append(f"phylo.labeling.{key} must be a boolean")

    phyloxml_cfg = phylo_cfg.get("phyloxml", {}) or {}
    ct = phyloxml_cfg.get("confidence_type", "auto")
    if ct not in ("auto", "sh_like", "sh_alrt", "ufboot", "bootstrap"):
        errors.append(
            f"phylo.phyloxml.confidence_type '{ct}' is not supported "
            "(use 'auto', 'sh_like', 'sh_alrt', 'ufboot', or 'bootstrap')"
        )

    coloring_cfg = phylo_cfg.get("coloring", {}) or {}
    if "enabled" in coloring_cfg and not isinstance(coloring_cfg["enabled"], bool):
        errors.append("phylo.coloring.enabled must be a boolean")
    # Ranks usable for colouring — the 9-rank ladder plus the coarse
    # ranks TaxonomyInfo.get_rank answers from its named attributes.
    valid_color_ranks = {
        "species", "subgenus", "genus", "subfamily", "family",
        "suborder", "order", "subclass", "class", "phylum",
        "kingdom", "superkingdom",
    }
    cranks = coloring_cfg.get("ranks", ["genus"])
    if not isinstance(cranks, list) or not (1 <= len(cranks) <= 2):
        errors.append(
            "phylo.coloring.ranks must be a list of one or two rank names "
            "(one rank = distinct hue per value; two = [parent, child] with "
            "the child shaded within its parent's hue)"
        )
    else:
        for r in cranks:
            if not isinstance(r, str) or r.lower() not in valid_color_ranks:
                errors.append(
                    f"phylo.coloring.ranks contains unsupported rank '{r}' "
                    f"(use one of {sorted(valid_color_ranks)})"
                )
    for key in ("saturation", "value"):
        v = coloring_cfg.get(key, 0.65 if key == "saturation" else 0.90)
        if not isinstance(v, (int, float)) or isinstance(v, bool) or not (0 <= v <= 1):
            errors.append(f"phylo.coloring.{key} must be a number between 0 and 1")
    mc = coloring_cfg.get("missing_color", "#808080")
    if not (isinstance(mc, str) and re.fullmatch(r"#[0-9A-Fa-f]{6}", mc)):
        errors.append(
            "phylo.coloring.missing_color must be a '#RRGGBB' hex colour string"
        )

    pp_cfg = phylo_cfg.get("per_protein", {}) or {}
    mt = pp_cfg.get("min_taxa", 3)
    if not isinstance(mt, int) or isinstance(mt, bool) or mt < 1:
        errors.append(
            "phylo.per_protein.min_taxa must be a positive integer "
            "(values below 3 are clamped to 3 at runtime)"
        )
    if "incongruence" in pp_cfg and not isinstance(pp_cfg["incongruence"], bool):
        errors.append("phylo.per_protein.incongruence must be a boolean")
    pp_mafft = pp_cfg.get("mafft", {}) or {}
    pp_mafft_extra = pp_mafft.get("extra_args", [])
    if not isinstance(pp_mafft_extra, list) or not all(
        isinstance(x, str) for x in pp_mafft_extra
    ):
        errors.append("phylo.per_protein.mafft.extra_args must be a list of strings")
    if "domain_architecture" in pp_cfg and not isinstance(
        pp_cfg["domain_architecture"], bool
    ):
        errors.append("phylo.per_protein.domain_architecture must be a boolean")
    if "whole_polyprotein_tree" in pp_cfg and not isinstance(
        pp_cfg["whole_polyprotein_tree"], bool
    ):
        errors.append(
            "phylo.per_protein.whole_polyprotein_tree must be a boolean"
        )

    # trimAl blocks (whole-genome phylo.trimal + per-protein
    # phylo.per_protein.trimal share the same shape).
    _TRIMAL_MODES = (
        "automated1", "gappyout", "strict", "strictplus", "nogaps", "noallgaps",
    )

    def _check_trimal_block(block: dict, where: str) -> None:
        if "enabled" in block and not isinstance(block["enabled"], bool):
            errors.append(f"{where}.enabled must be a boolean")
        mode = block.get("mode", "automated1")
        if mode not in _TRIMAL_MODES:
            errors.append(
                f"{where}.mode '{mode}' is not supported "
                f"(use one of {list(_TRIMAL_MODES)}; threshold trimming goes "
                f"in {where}.extra_args)"
            )
        ta_extra = block.get("extra_args", [])
        if not isinstance(ta_extra, list) or not all(
            isinstance(x, str) for x in ta_extra
        ):
            errors.append(f"{where}.extra_args must be a list of strings")

    _check_trimal_block(phylo_cfg.get("trimal", {}) or {}, "phylo.trimal")
    _check_trimal_block(
        pp_cfg.get("trimal", {}) or {}, "phylo.per_protein.trimal",
    )

    part_cfg = phylo_cfg.get("partition", {}) or {}
    if "enabled" in part_cfg and not isinstance(part_cfg["enabled"], bool):
        errors.append("phylo.partition.enabled must be a boolean")
    plink = part_cfg.get("linkage", "proportional")
    if plink not in ("proportional", "equal", "unlinked"):
        errors.append(
            f"phylo.partition.linkage '{plink}' is not supported "
            "(use 'proportional', 'equal', or 'unlinked')"
        )
    pmodels = part_cfg.get("models", {})
    if not isinstance(pmodels, dict) or not all(
        isinstance(k, str) and isinstance(v, str) for k, v in pmodels.items()
    ):
        errors.append(
            "phylo.partition.models must be a mapping of family-label → "
            "model string"
        )

    rooting_cfg = phylo_cfg.get("rooting", {}) or {}
    rmethod = rooting_cfg.get("method", "auto")
    if rmethod not in ("auto", "taxonomy", "mad", "midpoint", "outgroup", "none"):
        errors.append(
            f"phylo.rooting.method '{rmethod}' is not supported "
            "(use 'auto', 'taxonomy', 'mad', 'midpoint', 'outgroup', "
            "or 'none')"
        )
    og = rooting_cfg.get("outgroup")
    if og is not None and not (
        isinstance(og, str)
        or (isinstance(og, list) and all(isinstance(x, str) for x in og))
    ):
        errors.append(
            "phylo.rooting.outgroup must be a string (single accession) "
            "or a list of accession strings"
        )
    ogr = rooting_cfg.get("outgroup_rank")
    if ogr is not None and not (
        isinstance(ogr, dict)
        and all(
            isinstance(k, str) and isinstance(v, str) and v.strip()
            for k, v in ogr.items()
        )
    ):
        errors.append(
            "phylo.rooting.outgroup_rank must be a mapping of "
            "rank-name → taxon name (strings)"
        )
    if rmethod == "outgroup" and not og and not ogr:
        errors.append(
            "phylo.rooting.method='outgroup' requires either "
            "phylo.rooting.outgroup (an accession or list) or "
            "phylo.rooting.outgroup_rank (a rank → taxon mapping)"
        )

    lca_cfg = phylo_cfg.get("lca", {}) or {}
    if "enabled" in lca_cfg and not isinstance(lca_cfg["enabled"], bool):
        errors.append("phylo.lca.enabled must be a boolean")
    valid_min_ranks = {
        "none", "superkingdom", "realm", "kingdom", "subkingdom",
        "phylum", "subphylum", "class", "subclass",
        "order", "suborder", "family", "subfamily",
        "genus", "subgenus", "species",
    }
    mr = lca_cfg.get("min_rank", "genus")
    if mr not in valid_min_ranks:
        errors.append(
            f"phylo.lca.min_rank '{mr}' is not supported "
            f"(use one of {sorted(valid_min_ranks)})"
        )
    ctv = lca_cfg.get("coverage_threshold", 0.5)
    if not isinstance(ctv, (int, float)) or isinstance(ctv, bool) or not (0 <= ctv <= 1):
        errors.append(
            "phylo.lca.coverage_threshold must be a number between 0 and 1"
        )

    mono_cfg = phylo_cfg.get("monophyly", {}) or {}
    ms = mono_cfg.get("min_support", 70)
    if not isinstance(ms, (int, float)) or isinstance(ms, bool) or not (0 <= ms <= 100):
        errors.append(
            "phylo.monophyly.min_support must be a number in [0, 100] "
            "(0 disables the support-aware collapse)"
        )
    if "include_species" in mono_cfg and not isinstance(
        mono_cfg["include_species"], bool
    ):
        errors.append("phylo.monophyly.include_species must be a boolean")

    pc_cfg = phylo_cfg.get("pre_cluster_tree", {}) or {}
    if "enabled" in pc_cfg and not isinstance(pc_cfg["enabled"], bool):
        errors.append("phylo.pre_cluster_tree.enabled must be a boolean")
    pml = pc_cfg.get("max_leaves", 5000)
    if not isinstance(pml, int) or isinstance(pml, bool) or pml < 0:
        errors.append(
            "phylo.pre_cluster_tree.max_leaves must be a non-negative "
            "integer (0 = no cap, build every post-QC sequence as a leaf "
            "— only on a big-memory machine)"
        )
    ptt = pc_cfg.get("parttree_threshold", 10000)
    if not isinstance(ptt, int) or isinstance(ptt, bool) or ptt < 0:
        errors.append(
            "phylo.pre_cluster_tree.parttree_threshold must be a "
            "non-negative integer (0 = always use MAFFT --parttree; a "
            "very large value effectively never)"
        )

    # HMM block
    hmm = cfg.get("hmm", {}) or {}
    if "enabled" in hmm and not isinstance(hmm["enabled"], bool):
        errors.append("hmm.enabled must be a boolean")
    if "use_ga_when_available" in hmm and not isinstance(
        hmm["use_ga_when_available"], bool
    ):
        errors.append("hmm.use_ga_when_available must be a boolean")
    db = hmm.get("database")
    if db is not None and not isinstance(db, str):
        errors.append(
            "hmm.database must be null (use bundled set) or a path string"
        )
    ev = hmm.get("default_evalue", 1.0e-5)
    if not isinstance(ev, (int, float)) or isinstance(ev, bool) or ev <= 0:
        errors.append(
            "hmm.default_evalue must be a positive number (e.g. 1.0e-5)"
        )
    rc = hmm.get("relative_length_cutoff", 0.5)
    if not isinstance(rc, (int, float)) or isinstance(rc, bool) or not (0 < rc <= 1):
        errors.append(
            "hmm.relative_length_cutoff must be a number in (0, 1] "
            "(fraction of HMM model length the alignment must span)"
        )
    ht = hmm.get("threads")
    if ht is not None and (
        not isinstance(ht, int) or isinstance(ht, bool) or ht < 1
    ):
        errors.append(
            "hmm.threads must be null (use cfg.threads) or a positive integer"
        )
    mot = hmm.get("multidomain_overlap_tolerance", 30)
    if not isinstance(mot, int) or isinstance(mot, bool) or mot < 0:
        errors.append(
            "hmm.multidomain_overlap_tolerance must be a non-negative integer "
            "(amino acids; 0 = strict non-overlap, 30 = default)"
        )

    # Representative priority
    priority = cfg.get("representative", {}).get("priority", [])
    valid_priorities = {"refseq", "reviewed_uniprot", "longest"}
    for p in priority:
        if p not in valid_priorities:
            errors.append(
                f"representative.priority contains unknown value '{p}'. "
                f"Valid values: {sorted(valid_priorities)}"
            )
    if "longest" not in priority:
        errors.append("representative.priority must include 'longest' as a fallback")

    pr_cfg = cfg.get("output", {}).get("protein_report", {}) or {}
    max_breakdown = pr_cfg.get("max_breakdown", 20)
    if (
        not isinstance(max_breakdown, int)
        or isinstance(max_breakdown, bool)
        or max_breakdown < 1
    ):
        errors.append(
            "output.protein_report.max_breakdown must be a positive integer"
        )

    # Overrides (force-keep whitelist).
    ov = cfg.get("overrides", {}) or {}
    ids = ov.get("ids", [])
    if ids is not None and not isinstance(ids, list):
        errors.append("overrides.ids must be a list of accession/isolate-id strings")
    ids_file = ov.get("ids_file")
    if ids_file is not None and not isinstance(ids_file, str):
        errors.append("overrides.ids_file must be a path string (or null)")
    if "protect_qc" in ov and not isinstance(ov.get("protect_qc"), bool):
        errors.append("overrides.protect_qc must be a boolean")
    if "force_select" in ov and not isinstance(ov.get("force_select"), bool):
        errors.append("overrides.force_select must be a boolean")
    stages = ov.get("protect_stages", "all")
    if stages != "all":
        from .overrides import QC_PROTECT_STAGES

        if isinstance(stages, str):
            stages_list = [stages]
        elif isinstance(stages, list):
            stages_list = stages
        else:
            stages_list = None
            errors.append(
                "overrides.protect_stages must be \"all\" or a list of stage tokens"
            )
        if stages_list is not None:
            bad = [s for s in stages_list if s not in QC_PROTECT_STAGES]
            if bad:
                errors.append(
                    f"overrides.protect_stages has unknown stage(s): "
                    f"{', '.join(map(str, bad))}. Valid tokens: "
                    f"{', '.join(QC_PROTECT_STAGES)} (or \"all\")."
                )

    # Exclude blocklist (overrides.exclude): own dedicated id list.
    ex = ov.get("exclude", {}) or {}
    if "enabled" in ex and not isinstance(ex.get("enabled"), bool):
        errors.append("overrides.exclude.enabled must be a boolean")
    ex_ids = ex.get("ids", [])
    if ex_ids is not None and not isinstance(ex_ids, list):
        errors.append(
            "overrides.exclude.ids must be a list of accession/id strings"
        )
    ex_file = ex.get("ids_file")
    if ex_file is not None and not isinstance(ex_file, str):
        errors.append("overrides.exclude.ids_file must be a path string (or null)")
    # Contradiction guard: the same id cannot be both blocklisted and
    # force-kept / force-selected. Only meaningful when a keep/pin capability
    # is actually active (the shared `ids` list is inert otherwise), and only
    # attempted once the id fields have valid types (else the type errors
    # above already explain the problem).
    types_ok = isinstance(ex_ids, list) and (ex_file is None or isinstance(ex_file, str))
    if (ov.get("protect_qc") or ov.get("force_select")) and types_ok:
        from .overrides import _norm_id, _strip_version, resolve_ids, resolve_raw_ids

        keep_set = resolve_ids(ov)
        exclude_set = resolve_ids({"ids": ex_ids, "ids_file": ex_file})
        clash = keep_set & exclude_set
        if clash:
            # Echo the ids as the user typed them (original case), not the
            # normalised forms, so the message is actionable.
            raw = resolve_raw_ids({"ids": ex_ids, "ids_file": ex_file})
            shown = [
                r for r in raw
                if (_norm_id(r) in clash
                    or _strip_version(_norm_id(r) or "") in clash)
            ] or sorted(clash)
            errors.append(
                "overrides.exclude lists id(s) that are also force-kept "
                "(overrides.protect_qc) or force-selected "
                "(overrides.force_select): "
                f"{', '.join(shown[:10])}"
                f"{' …' if len(shown) > 10 else ''}. An id cannot be both "
                "removed and guaranteed — remove it from one list."
            )

    # output.polyprotein_report.wall_warning — the peptide-coverage alarm knobs.
    ww = (
        cfg.get("output", {}).get("polyprotein_report", {}) or {}
    ).get("wall_warning", {}) or {}
    if not isinstance(ww.get("enabled", True), bool):
        errors.append("output.polyprotein_report.wall_warning.enabled must be true or false")
    wall_ranks = (
        "species", "subgenus", "genus", "subfamily", "family",
        "suborder", "order", "subclass", "class",
    )
    wrank = ww.get("rank", "genus")
    if not isinstance(wrank, str) or wrank not in wall_ranks:
        errors.append(
            "output.polyprotein_report.wall_warning.rank must be one of "
            f"{', '.join(wall_ranks)}"
        )
    wfrac = ww.get("wall_fraction", 0.6)
    if not isinstance(wfrac, (int, float)) or isinstance(wfrac, bool) or not (0.0 <= wfrac <= 1.0):
        errors.append(
            "output.polyprotein_report.wall_warning.wall_fraction must be a "
            "number in [0, 1]"
        )
    wmin = ww.get("min_reps", 3)
    if not isinstance(wmin, int) or isinstance(wmin, bool) or wmin < 1:
        errors.append(
            "output.polyprotein_report.wall_warning.min_reps must be a "
            "positive integer"
        )

    return errors


def get_virus_config(cfg: dict[str, Any]) -> Optional[dict[str, Any]]:
    """Return the active virus config dict, or None if segmented mode is off."""
    seg = cfg.get("segmented", {})
    if not seg.get("enabled"):
        return None
    virus_name = seg.get("virus")
    return seg.get("viruses", {}).get(virus_name)
