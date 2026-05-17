# CLAUDE.md

Guidance for Claude Code (or any future agent) working on this repo.

## What this project is

`repseq` is a Python CLI for **selecting representative sequences** from
large FASTA datasets — protein or nucleotide, with a strong viral-sequence
focus. The user is a JCVI bioinformatics researcher; treat them as fluent in
the domain (taxonomy, clustering, FASTA conventions) but not necessarily in
every implementation detail, since most of the code was written by Claude.

Entry point: `repseq` (defined in `repseq/cli.py`).

## Purpose & scope (designer's mental model)

The points below are the load-bearing scope statement from the project
designer. Treat them as durable requirements: when changes would conflict
with one of them, raise it explicitly rather than silently breaking
intent.

1. **Audience: bench scientists, not developers.** End users are
   inexperienced in computer science, math, and statistics. Ease of use is
   paramount: `--help` screens, README/CLAUDE.md-level documentation, log
   messages during long-running steps (so the user knows the pipeline is
   alive), and **clear, actionable error messages** (plain English, named
   next step — never a raw stack trace) are first-class deliverables, not
   polish.

2. **Main outputs.** Every run is expected to produce, where applicable:
   - **2A — Representative genome nucleotide sequences (FASTA):** the
     selected representatives as NT sequences. For segmented viruses,
     emit *both* the per-segment files and the per-isolate concatenated
     file.
   - **2B — Proteins of the representatives (FASTA):** the amino-acid
     translations of every CDS of every representative.
   - **2C — Detailed TSV** that ties AA sequences to NT sequences: one
     row per protein, with isolate, segment, accession, length, and full
     taxonomy / metadata columns (see `{prefix}_isolate_proteins.tsv`).
   - **2D — Clustering summary files:** per-group cluster counts and
     thresholds (`{prefix}_group_counts.tsv`), plus the textual run
     report.
   - **2E — Phylogenetic tree of the representatives (phyloXML), with
     MSA kept:** built from the same concatenated protein sequences used
     for clustering when `alphabet=protein`, otherwise from genome /
     clustered segment NT sequences. The MSA must be retained alongside
     the tree.
   - **2F — Per-protein phylogenetic trees (phyloXML), with MSAs kept**
     — *not implemented yet.* One tree per protein found across the
     clustered representatives. The designer will request this in a
     future session.
   - **2G — Graphical clustering visualisation (UMAP scatter):**
     implemented under `[viz]` extras but currently struggles with
     `umap-learn` / `matplotlib` dependency installation; treat
     dependency robustness as an open issue.

3. **Clustering alphabet is selectable.** Clustering may operate on
   either nucleotide sequences (concatenated for segmented viruses) or
   amino-acid sequences (concatenated for segmented viruses) — see
   `clustering.alphabet` and the v0.6.0 status block.

4. **Segmented viruses are first-class.** The pipeline must process
   them seamlessly: completeness filtering, per-segment metadata,
   concatenation, and per-segment output are core behaviour, not
   add-ons.

5. **Phylogenetic inference stack:** MAFFT (MSA), **Trimal (not
   implemented yet — alignment trimming will be added)**, IQ-TREE
   (protein default), FastTree (NT default / fallback). Trimal would
   sit between MAFFT and the tree-builder.

6. **Robustness and caching are core.** All external lookups go through
   the SQLite-backed taxonomy cache; optional steps (cd-hit, `--phylo`,
   `--plot`) fail soft with a clear stderr message rather than aborting
   the run; `repseq doctor` exists so a bench scientist can self-diagnose
   install state.

7. **Database fetches drive metadata enrichment.** The pipeline reaches
   out to NCBI Entrez and UniProt to obtain taxonomy, isolate ID,
   strain, segment, host, collection date, and country — see the
   `taxonomy/` and segmented-metadata code paths. Network usage is
   batched, cached, and rate-limited; `--no-resolve` exists for fully
   offline runs.

## Top-level layout

```
repseq/
├── cli.py                  ← click commands, dispatch into modes
├── config.py               ← YAML defaults + validation
├── models.py               ← Sequence, Cluster, QCReport, RunResult dataclasses
├── io/fasta.py             ← parse 3 header formats, read/write FASTA
├── qc/
│   ├── pipeline.py         ← duplicates → length → ambiguous → annotation
│   └── protein_qc.py       ← NCBI-backed protein-count filter (opt-in)
├── segmented/completeness.py ← isolate grouping, completeness, concat
├── taxonomy/
│   ├── cache.py            ← SQLite TTL cache
│   ├── ncbi.py             ← Entrez API client
│   ├── uniprot.py          ← UniProt REST client
│   └── resolver.py         ← DB → header → strain-label fallback chain
├── clustering/
│   ├── __init__.py         ← `run_clustering` dispatcher + `min_threshold`
│   ├── mmseqs2.py          ← shells out to `mmseqs`
│   ├── cdhit.py            ← shells out to `cd-hit` / `cd-hit-est`
│   ├── marker.py           ← `select_marker_protein` + `populate_protein_sequences`
│   └── diversity.py        ← MaxMin selection (alignment-free, k-mer Jaccard)
├── representative/selector.py ← RefSeq > reviewed > longest priority
├── modes/                  ← one file per selection mode, all extend BaseMode
├── output/{writer,report}.py ← FASTA + TSV/text report writers
├── viz/clustering_plot.py  ← UMAP scatter (optional, [viz] extras)
├── phylo/                  ← optional MSA + tree step (--phylo)
│   ├── mafft.py            ← shells out to `mafft --auto`
│   ├── fasttree.py         ← shells out to `FastTree` (auto-picks -nt/-gtr)
│   ├── iqtree.py           ← shells out to `iqtree2` (ModelFinder + UFBoot)
│   └── pipeline.py         ← short-id remap → MAFFT → IQ-TREE/FastTree → phyloXML
└── doctor.py               ← `repseq doctor` self-test (deps, tools, net, config)
```

`config/default_config.yaml` is the documented config schema; the same dict
shape is hard-coded in `repseq/config.py:DEFAULTS`.

## Data flow

1. `_load_sequences` (cli.py) → `read_fasta` yields `Sequence` objects with
   parsed `accession`, `organism`, `is_refseq`, `is_reviewed`, etc.
2. `_resolve_metadata` builds a `MetadataResolver` and runs
   `resolve_batch` in a thread pool. Failures are non-fatal but **tracked**
   on `resolver.failures` and surfaced via stderr. Returns the `NCBITaxonomy`
   client so downstream steps can reuse it.
3. `_run_qc` produces a passing-list + `QCReport`. In segmented mode it
   skips both the whole-pool length filter and exact-duplicate removal
   (`QCReport.length_filter_skipped` / `dedup_skipped`) — see step 5.
4. `_run_protein_qc` (optional) — if `qc.protein_annotation.enabled` or
   `virus.expected_proteins_per_segment` is set, fetches CDS features via
   `NCBITaxonomy.fetch_proteins_batch` (200 accessions per request, cached
   in SQLite under source `ncbi_proteins`), populates `seq.proteins`, and
   drops sequences failing the count check.
5. If `segmented.enabled` *and* `segmented.use_genbank_metadata` (default
   true) *and* an NCBI client is available, `_populate_genbank_isolate_segment`
   pulls `/isolate`, `/strain`, `/segment` from the GenBank source feature
   for every NCBI-sourced sequence and writes them onto `seq.isolate_id`,
   `seq.strain`, `seq.segment`. Shares the `ncbi_proteins` cache with
   protein QC via `fetch_source_metadata_batch`, so a run that also does
   protein QC pays one efetch round trip, not two. UniProt sequences and
   sequences without an accession are skipped — they fall through to the
   header-regex parse in step 7. The toggle exists so a user who would
   rather not depend on NCBI lookups can revert to the pre-v0.6 behaviour.
6. `_setup_protein_alphabet` reads `clustering.alphabet` (default `protein`):
   for `protein` it auto-triggers `attach_proteins` if any sequence's
   `proteins` is `None` (one-shot GenBank CDS fetch, same `ncbi_proteins`
   cache as protein-QC), then — for *non-segmented* runs — selects each
   sequence's marker via `clustering.marker.select_marker_protein`
   (longest CDS, or first `clustering.cluster_protein` alias that matches
   `/product` as a case-insensitive substring) and stores its AA on
   `seq.protein_sequence`. Sequences with no viable marker are dropped
   under `removed_proteins`. For segmented runs this step only fetches
   proteins — the per-isolate AA concat is built in step 7. `auto` falls
   back to `nucleotide` if no proteins are fetchable; pure `protein` with
   `--no-resolve` aborts at startup. No-op when `alphabet=nucleotide`.
7. If `segmented.enabled`, `_handle_segmented` runs `filter_complete_isolates`
   then `build_concatenated_sequences`, replacing the sequence list with one
   concatenated sequence per complete isolate. `filter_complete_isolates`
   prefers `seq.isolate_id` / `seq.segment` (populated in step 5) and only
   falls back to `isolate_regex` and the segment alias/header scan when
   those fields are unset. When `alphabet=protein`,
   `build_concatenated_sequences` also picks each segment's marker protein
   (per-segment `cluster_protein` aliases, else longest CDS) and stores the
   in-segments-order concat on the concat's `protein_sequence`; isolates
   whose marker is missing on any segment are dropped and counted under
   `removed_incomplete_isolates`. Exact-duplicate removal is then applied
   on the *concatenated* sequences (not the segment pool — a conserved
   segment shared between two distinct isolates must not knock either out
   as incomplete); de-duplicated isolates are also pruned from
   `complete_isolates` so the per-segment output files stay consistent.
8. The chosen mode runs (`modes/<mode>.py`), returning a `RunResult` with
   `representatives` + `clusters`. Clustering backends honour
   `clustering.alphabet`: with `protein`, `_write_id_fasta` writes
   `seq.protein_sequence` (the marker / AA concat) as the FASTA body and
   the cd-hit dispatcher forces `cd-hit` (not `cd-hit-est`) regardless of
   `seq_type`. Cluster objects still carry the original NT-bearing
   `Sequence`, so rep selection and downstream output are unchanged.
9. `write_results` writes the FASTA(s); `write_all_reports` writes the
   plain-text/TSV reports. Column names are deliberately harmonised
   across files so the same concept (accession, organism, length, the
   9-rank taxonomy ladder, boolean cells as `TRUE`/`FALSE`) carries the
   same name everywhere — see `repseq/output/report.py` constants
   `_TAX_RANKS` and `_tsv_bool`. The full set:

   - `{prefix}_qc_removed.tsv`: `accession, reason`.
   - `{prefix}_representative_sequences.tsv` (non-segmented mode):
     `accession, organism, description, strain, host, collection_date,
     country, segment, isolate_id, molecule_type, length_nt,
     is_refseq, is_reviewed, ncbi_taxon_id`, then the 9-rank
     `_TAX_RANKS` ladder (`species, subgenus, genus, subfamily,
     family, suborder, order, subclass, class`). One row per
     representative sequence.
   - `{prefix}_representative_isolates.tsv` (segmented mode): one row
     per representative isolate (synthetic `CONCAT|<isolate_id>`
     Sequence with `concat_segments` populated). Columns: `isolate_id,
     organism, strain, host, collection_date, country, n_segments,
     segments` (comma-joined segment names in concat order),
     `accessions` (comma-joined per-segment GenBank accessions in
     concat order), `total_length_nt` (sum of per-segment NT lengths),
     `is_refseq, is_reviewed, ncbi_taxon_id`, then the same
     `_TAX_RANKS` ladder. The per-sequence columns (`accession`,
     `segment`, `description`, `molecule_type`, `length_nt`) are
     deliberately *absent* — they have no isolate-level meaning. Two
     distinct writers (`write_representative_sequences_tsv` /
     `write_representative_isolates_tsv`) so each schema stays clean.
     Sub-ranks populate only from the resolver lineage map and are
     commonly blank for viruses.
   - `{prefix}_isolate_proteins.tsv` (segmented mode, one row per CDS
     per passing isolate): `protein_id, product, length_aa, isolate_id,
     segment, segment_length_nt, accession, representative`, then the
     same `_TAX_RANKS` ladder. `representative` is `TRUE` if the
     isolate survived clustering, `FALSE` otherwise — derived from
     `result.representatives` (handles both `isolate_id` reps and
     synthetic `CONCAT|<isolate_id>` reps).
   - `{prefix}_representative_isolate_proteins.tsv` (segmented mode):
     row-filtered companion to `_isolate_proteins.tsv`. Same exact
     column schema (including the `representative` column, which is
     always `TRUE` here) but only rows belonging to representative
     isolates. Written by the same `write_isolate_proteins_tsv` writer
     with a pre-filtered `complete_isolates` dict.
   - `{prefix}_clusters.tsv`: `cluster_id, accession, organism,
     cluster_size, is_refseq, is_reviewed`. `accession` is the
     representative's accession.
   - `{prefix}_group_counts.tsv` (one row per stratum):
     `stratified_by, stratum, stratum_size_before, stratum_size_after,
     clustered, cutoff` — populated from `RunResult.group_stats`,
     which every mode fills in. Internal `GroupStat` field names are
     still `grouping`/`group`/`n_before`/`n_after` for historical
     reasons; only the column labels were harmonised.
   - `{prefix}_tree_id_map.tsv` (only when `--phylo` ran):
     `short_id, accession`.
   - `{prefix}_proteins.fasta` (amino-acid sequences for all proteins
     of the selected representatives). Reconstructed from the same
     cached GenBank records — no extra network calls.
   - `{prefix}_representatives_protein.fasta` (only when any
     representative carries a populated `protein_sequence`, i.e.
     `alphabet=protein` actually fired): the AA strings that were fed
     into the clustering step — per-isolate marker concat in segmented
     mode, per-rep marker in non-segmented mode.
10. If `--plot` is passed, `viz.clustering_plot.write_clustering_plot`
    embeds the clustered sequences with UMAP on k-mer Jaccard distance
    and writes `{prefix}_clustering.png`. Cost-bounded by a default
    2000-point subsample (representatives always kept); skipped when
    the run produced no clusters. Requires the `[viz]` extras
    (`matplotlib` + `umap-learn`) — `ImportError` is surfaced gracefully.
11. If `--phylo` is passed, `phylo.run_phylogeny` builds an MSA (MAFFT
    `--auto`) and an ML tree over the representatives. The tree-builder
    is selected by `phylo.tool` (`auto` | `iqtree` | `fasttree`); `auto`
    (the default) picks **IQ-TREE for protein alignments, FastTree for
    nucleotide** — IQ-TREE's ModelFinder is the killer feature for
    protein (JTT/WAG/LG/etc. can change topology), while FastTree's
    `-nt -gtr` is well-understood and orders of magnitude faster on NT.
    When `clustering.alphabet=protein` and every rep carries a
    `protein_sequence`, the MSA/tree are built on the AA strings;
    otherwise the alphabet comes from the rep `seq_type`. Every rep
    gets a deterministic short id `S0001…SNNNN` for the MAFFT pipeline
    (long names, whitespace, and pipes break many phylo tools); the
    intermediate Newick keeps the short ids. The orchestrator then
    **roots** the tree (`repseq/phylo/rooting.py`: taxonomy-guided →
    MAD → midpoint chain by default; pinnable via `phylo.rooting.method`)
    and labels every internal node with the LCA of its terminals
    (`repseq/phylo/lca.py`: ≥50% lineage coverage gate, `min_rank=genus`
    leaf-vote filter, `keep_deepest_labels` cleanup so each
    monophyletic group's label lands on the crown,
    `suppress_same_species_pairs` to avoid duplicating species names on
    a 2-leaf internal). The final phyloXML is emitted by
    `repseq/phylo/phyloxml_writer.py` — every leaf gets a formatted
    `<name>` (driven by `phylo.labeling.format` or `segmented_format`),
    a `<taxonomy>` block with NCBI taxon id, a `<sequence>` block with
    the GenBank accession + title, and repseq-namespaced `<property>`
    elements for host, collection_date, country, strain, isolate_id,
    year, species, genus, subfamily, and family (empty values omitted);
    every annotated internal clade gets a `<name>` and a `<taxonomy>`
    with `<scientific_name>` + `<rank>` (validated against PhyloXML's
    enum, falling back to `"other"`). The `<phylogeny>` element carries
    a `<name>` and a `<description>` recording MAFFT/IQ-TREE/FastTree
    versions, the selected model, bootstrap settings, and which rooting
    method actually fired. The tree is ladderized (`reverse=True`) and
    confidence values are normalised to 0-100 integers
    (`type="sh_like"` for FastTree, `"ufboot"` for IQ-TREE). IQ-TREE's `--prefix` lands all its scratch files in a
    temp dir that's wiped at end of run; only the canonical `.treefile`
    is copied to `{prefix}_tree.nwk` and the `.iqtree` model-selection
    report is copied to `{prefix}_iqtree_summary.txt`. IQ-TREE refuses
    UFBoot with `<4` sequences, so the wrapper auto-drops bootstrap
    (keeps the tree) when the MSA has 3 reps. Outputs:
    `{prefix}_msa.fasta`, `{prefix}_tree.nwk`, `{prefix}_tree.xml`,
    `{prefix}_tree_id_map.tsv`, and (IQ-TREE only)
    `{prefix}_iqtree_summary.txt`. Skipped with a stderr warning if
    `<3` reps or if any of the binaries are missing or fail (mirrors
    `--plot` failure handling).

## Invariants worth knowing

- **Representative priority** (`representative.priority`) must include
  `"longest"` — `validate_config` enforces this. The selector always
  uses `length` as final tiebreaker even when not listed.
- **Cluster representative swapping**: MMseqs2 picks the longest as rep by
  default; `apply_representative_selection` re-elects across all
  members so RefSeq/reviewed wins. When swapping, the old rep is appended
  to `cluster.members`.
- **Segmented mode** requires three fields per virus in config:
  `expected_segments`, `segments`, `isolate_regex`. The regex must capture
  either a named group `isolate` or group 1. The regex is now a *fallback*:
  with `segmented.use_genbank_metadata: true` (default), `seq.isolate_id`
  and `seq.segment` come from the GenBank source feature
  (`/isolate`/`/strain`/`/segment`) first; the regex only fires for sequences
  the GenBank lookup can't fill (UniProt input, missing accession, missing
  qualifier, or `--no-resolve` / toggle off).
- **Taxonomy cache** keys are `(source, key)`; TTL eviction happens lazily
  on `get` (expired entries are deleted in place).
- **Threads**: only `MetadataResolver.resolve_batch` and the clustering
  binary (MMseqs2 or cd-hit) consume `cfg["threads"]`. Internal Python
  work is single-threaded.
- **Clustering backend**: `cfg["clustering"]["backend"]` selects
  `"mmseqs2"` (default) or `"cdhit"`. Modes import
  `from ..clustering import run_clustering` (the dispatcher), never from
  a backend module directly. The cd-hit wrapper auto-picks `cd-hit` for
  protein input and `cd-hit-est` for nucleotide. Identity floors differ
  by backend: mmseqs2 = 0.0, cd-hit protein = 0.40, cd-hit-est = 0.80.
  `clustering.min_threshold(cfg, sequences)` returns the active floor;
  `_binary_search_threshold` clamps `lo` to it so the search never asks
  a backend for a value it would refuse.
- **Clustering alphabet**: `cfg["clustering"]["alphabet"]` selects what's
  fed to the backend — `"protein"` (default, since v0.6.0), `"nucleotide"`,
  or `"auto"`. Backends read `seq.protein_sequence` instead of
  `seq.sequence` when alphabet=`protein`; the cd-hit dispatcher likewise
  forces `cd-hit` (with the 0.40 floor) over `cd-hit-est` regardless of
  `seq_type` because the concat sequence is NT-typed even when its
  `protein_sequence` is set. `seq.protein_sequence` is populated by
  `_setup_protein_alphabet` (non-segmented) or `build_concatenated_sequences`
  (segmented); both use `clustering.marker.select_marker_protein` (longest
  CDS by default, or first matching `cluster_protein` alias against
  `/product`, case-insensitive substring). Sequences/isolates without a
  viable marker are dropped, not chimerised with NT. The clustering output
  filenames are unchanged; an extra `{prefix}_representatives_protein.fasta`
  carries the AA strings when alphabet=protein actually fired.

## Modes

Every mode subclasses `BaseMode` and returns a `RunResult`. Most modes
share the same shape:

1. Group sequences by some key (rank / host / decade / country / custom field).
2. For each group, either keep all (if `<= n_per_group`) or run a
   **binary-search over clustering thresholds** to land at `n_per_group`
   clusters (see `taxonomic1._binary_search_threshold` for the canonical
   implementation; other modes import or mirror this). The search lower
   bound is clamped to the active backend's identity floor.
3. `--overflow keep` returns whatever cluster count the binary search
   produced; `--overflow trim` follows up with `select_diverse` to enforce
   an exact count.

`GlobalMode` is the exception: it runs *either* a single threshold pass
(`-T`) or one shot of `select_diverse` (`-n`).

## How to run things

```bash
# Tests (offline, mocked)
pytest tests/

# Smoke a single mode
repseq global -c config/default_config.yaml -i some.fasta -T 0.95 --no-resolve

# Dry-run config validation
repseq taxonomic1 -c my.yaml -i x.fasta --rank genus -n 5 --dry-run
```

`--no-resolve` skips all network lookups — useful for fast local iteration.

## Common gotchas

- **MMseqs2 not in PATH** ⇒ `MMseqs2Error` from `clustering/mmseqs2.py`.
  Cluster-based modes fail; `--no-resolve` does not help — only `global -n`
  (diversity-only) and modes where every group is small enough to skip
  clustering will run without it. Same failure mode for `CDHitError` when
  `clustering.backend: cdhit` and `cd-hit` / `cd-hit-est` are missing.
- **No taxonomy info** ⇒ taxonomic modes fall back to grouping under
  `"Unknown"`. Always run with metadata resolution unless you've
  pre-resolved.
- **Source auto-detection on weird headers**: prefer `--source <name>` for
  guaranteed parsing. NCBI Virus headers can resemble plain NCBI; the
  auto-detector checks for multiple `[...]` brackets as a tell.
- **Concatenated sequences in output**: in segmented mode, the
  representatives list contains synthetic `CONCAT|<isolate_id>` entries.
  `output/writer._write_segmented` expands these back into per-segment
  files. The `--phylo` step also operates on those concatenated reps
  directly, so the tree groups isolates (one leaf per isolate) rather
  than segments.
- **Output directory must be empty**: `_load_and_validate` calls
  `_check_output_dir`, which aborts (exit 1) if `output.dir` already
  exists and is non-empty — runs never overwrite or mix into prior
  results. Tests that exercise the writers point at a fresh `tmp_path`.

## When making changes

- Keep the `Sequence` dataclass surface stable — it's threaded through
  every layer.
- Don't add network calls outside `taxonomy/`. Network-touching tests
  should `unittest.mock.MagicMock` the NCBI/UniProt clients (see
  `tests/test_resolver.py`).
- New modes: subclass `BaseMode`, register a new click subcommand in
  `cli.py`, and follow the existing `_load_sequences → _resolve_metadata
  → _run_qc → _run_protein_qc → _setup_protein_alphabet →
  _populate_genbank_isolate_segment → _handle_segmented → mode.run →
  _write_output` flow.
- Config additions: extend `DEFAULTS` and (if validation matters) add a
  check in `validate_config`. Document in `config/default_config.yaml`.

## Status

`v0.8.0` is a TSV output-schema overhaul. Column names are now
harmonised across the six TSVs the program writes — the canonical
sequence identifier is `accession` everywhere (was a mix of
`sequence_id` / `representative_id` / `original_id`), booleans are
`TRUE` / `FALSE` everywhere (was a mix of cases), length columns
carry their alphabet (`length_nt` / `length_aa` / `segment_length_nt`
/ `total_length_nt`), the taxonomic ladder is the same 9-rank
sequence everywhere (`_TAX_RANKS`), and `_group_counts.tsv` adopts
statistical vocabulary (`stratified_by`, `stratum`,
`stratum_size_before/after`). `seq_type` became `molecule_type`,
`taxid` became `ncbi_taxon_id`. The rep TSV split by mode:
`_representative_sequences.tsv` (non-segmented, per-sequence schema)
or `_representative_isolates.tsv` (segmented, per-isolate schema with
`n_segments`, `segments`, `accessions`, `total_length_nt` instead of
the per-sequence `accession`/`segment`/`description`/`molecule_type`/
`length_nt`). `_isolate_proteins.tsv` gained a `representative`
column (TRUE/FALSE) marking whether the isolate survived clustering.
A new `_representative_isolate_proteins.tsv` is a row-filtered
companion to `_isolate_proteins.tsv` carrying only proteins of
representative isolates (same schema). Breaking change for any
downstream script that parses output by column name or filename.

