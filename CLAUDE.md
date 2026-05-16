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
   plain-text/TSV reports — including `{prefix}_group_counts.tsv`
   (one row per stratification group: `grouping, group, n_before,
   n_after, clustered, cutoff` — populated from `RunResult.group_stats`,
   which every mode fills in), `{prefix}_isolate_proteins.tsv`
   (one row per CDS per passing isolate, segmented mode — columns:
   `protein_id, product, length, isolate_id, segment, segment_length,
   accession, species, subgenus, genus, subfamily, family, suborder,
   order, subclass, class`; the four sub-ranks are populated only from
   the resolver lineage map and commonly blank for viruses) and
   `{prefix}_proteins.fasta` (amino-acid sequences for all proteins of
   the selected representatives). The protein FASTA is reconstructed
   from the same cached GenBank records — no extra network calls. When
   any representative carries a populated `protein_sequence` (i.e.
   `alphabet=protein` actually fired), an additional
   `{prefix}_representatives_protein.fasta` is written, holding the AA
   strings that were fed into the clustering step (per-isolate marker
   concat in segmented mode, per-rep marker in non-segmented mode).
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

`v0.7.0` overhauls phylogenetic-tree annotation in two passes, both
released together.

**Pass B — rooting + internal-node LCA annotation** —
builds on Pass A's rich phyloXML by post-processing the tree before
serialisation. Two new modules:

* `repseq/phylo/rooting.py` — picks a root via the chain
  **taxonomy-guided → MAD → midpoint** (first success wins), or by
  user-pinned method (`phylo.rooting.method`: `auto` default,
  `taxonomy`, `mad`, `midpoint`, or `none`). Taxonomy-guided scoring
  is mean LCA specificity across internal clades, weighted by clade
  size, gated by ≥50% lineage coverage. MAD is a pure-Python port
  (no external dep) that minimises Σρ² across all leaf pairs over
  per-branch split points; analytical solution per branch.
* `repseq/phylo/lca.py` — labels every internal clade with the
  lowest common ancestor of its terminals (read from
  `TaxonomyInfo.lineage`). Coverage gate `phylo.lca.coverage_threshold`
  (default 0.5) skips bare clades; `phylo.lca.min_rank` (default
  `"genus"` — viral data rarely reaches species universally)
  excludes thinly-annotated leaves from the LCA vote without
  removing them from the tree. `keep_deepest_labels` walks
  largest-first and clears every nested duplicate so each
  monophyletic group's label lands on the crown, not the
  intermediate internals. `suppress_same_species_pairs` clears the
  LCA name from any internal whose only children are two leaves of
  the same species (the species name is already on the leaves).

ICTV rank inference (`_infer_rank_from_name`) handles the common
suffix → rank mapping (`-viridae` → family, `-virales` → order,
`-viricetes` → class, single-word `…virus` → genus, etc.) for when
the lineage map doesn't carry an explicit rank. PhyloXML's `<rank>`
enum is validated via `phyloxml_rank` — unknown ranks (`"no rank"`,
`"clade"`, NCBI's odd custom ranks) fall back to `"other"` so the
file always validates.

Writer (`phyloxml_writer.py`) now accepts a pre-parsed `tree=` (so
the pipeline can root + LCA-annotate before serialisation). Internal
clades emit `<name>` and `<taxonomy>` (`<scientific_name>` + `<rank>`)
when `_lca_name` is set; otherwise stay bare. The
`<phylogeny><description>` records which rooting method actually
fired (`rooting=taxonomy` / `mad` / `midpoint`), so an `auto`-chain
result is auditable from the file alone.

New tests:
- `tests/test_rooting.py` (10): LCA-prefix correctness, method=none
  passthrough, midpoint fallback when no lineage, taxonomy method
  fall-through when nothing to root by, taxonomy success on a tree
  with clean family groups, MAD success, invalid method falls back to
  auto, mean-LCA-specificity rewards consistent grouping, zero-score
  with no lineage data.
- `tests/test_lca.py` (15): ICTV-suffix inference (4: family, order,
  phylum, class, subfamily — and the multi-word "X virus" → no rank
  case), min_rank gate (species, too-coarse, "none" disables),
  LCA-prefix on (rank, name) tuples; annotate_internal_nodes
  family-rollup, coverage-gate skip, min_rank-filter excludes-but-keeps
  leaves; keep_deepest clears nested duplicates and preserves
  distinct labels; same-species pair suppression on 2-leaf internal,
  no-op when species differ, no-op when >2 children; phyloxml_rank
  accepts standard ranks and falls back to "other" for unknown.
- `tests/test_phyloxml_writer.py` (+4): internal `<name>` +
  `<taxonomy>(<scientific_name>,<rank>)` when LCA set; rank fallback
  to "other"; no `<name>`/`<taxonomy>` on bare internals; rooting
  method recorded in description.
- `tests/test_config.py` (+9): rooting defaults, all-five-methods
  accepted, unknown rejected, LCA defaults, enabled bool check,
  min_rank accept-list, unknown min_rank rejected,
  coverage_threshold range check.

596 offline tests total (Pass A + Pass B combined).

**Pass A — phyloXML annotation overhaul** — replaces the
one-line `Bio.Phylo.write(...)` phyloXML emission with a hand-rolled
writer (`repseq/phylo/phyloxml_writer.py`) that produces a
**richly-annotated** tree. Each leaf now carries:

* A formatted `<name>` driven by `phylo.labeling.format` (defaults
  `"{species}|{id}|{host}"`) or `phylo.labeling.segmented_format`
  (defaults `"{species}|{strain}|{host}"`) when segmented mode is on.
  Empty placeholders drop the preceding separator so labels never
  read `||` or end with `|`; `{strain}` falls back to `{isolate_id}`
  when the GenBank `/strain` qualifier is missing — one template
  works for both segmented and non-segmented runs.
* A `<taxonomy>` block with `<id provider="ncbi">` (from
  `TaxonomyInfo.taxid`) and `<scientific_name>` (species).
* A `<sequence type="dna|protein">` block with
  `<accession source="ncbi">` and `<name>` = original GenBank
  header / description.
* Repseq-namespaced `<property>` elements (under `xmlns:repseq=
  https://github.com/cmzmasek/repseq`) for `host`, `collection_date`,
  `country`, `strain`, `isolate_id`, `year` (parsed from
  `collection_date`), `species`, `genus`, `subfamily`, `family`.
  Empty values are **omitted**, not emitted as empty stubs.

Tree-level: `<phylogeny>` gains `<name>` (e.g.
`peribunyaviridae_genus5 [protein|MAFFT|IQ-TREE LG+G4]`) and
`<description>` (run timestamp + MAFFT/IQ-TREE/FastTree versions
captured via new `tool_version()` helpers + selected model + bootstrap
count + extra args). Tree is **ladderized** (`reverse=True`, larger
clades top) before write. Confidence values are normalised to **0-100
integers**: FastTree's SH-like `[0,1]` is rescaled, IQ-TREE's UFBoot
`[0,100]` passes through; the `<confidence type="...">` attribute
records the metric (`sh_like` / `ufboot` / `sh_alrt` / `bootstrap`),
overridable via `phylo.phyloxml.confidence_type`. Optional
`phylo.phyloxml.embed_alignment: true` inlines the per-leaf aligned
residues as `<mol_seq is_aligned="true">` (off by default —
`{prefix}_msa.fasta` is always written separately).

Per-isolate metadata inheritance in `concatenate_isolate` now uses
**first non-empty** across segments instead of always segment 0, so a
single blank field on the first segment no longer wipes the value
from the concat record (e.g. host present on segment M but blank on
segment L).

New files: `repseq/phylo/labels.py` (placeholder substitution +
separator-drop), `repseq/phylo/phyloxml_writer.py` (stdlib
`xml.etree.ElementTree`, namespace-aware, schema-ordered children).
The orchestrator (`repseq/phylo/pipeline.py`) now collects MAFFT /
tree-tool versions and the resolved model + bootstrap count and
passes them to the writer.

New tests:
- `tests/test_labels.py` (21): placeholder substitution, taxonomy
  rank resolution, year parsing, empty-token handling, separator-drop
  on empty, strain→isolate_id fallback, literal-text preservation,
  config helpers (`pick_format_string` segmented/non-segmented,
  `labeling_options` defaults + overrides).
- `tests/test_phyloxml_writer.py` (20): confidence normalisation
  (SH-like rescale, UFBoot passthrough, clamp, None → skip),
  confidence-type mapping (defaults, override, auto, unknown-tool),
  end-to-end write (XML well-formed, taxonomy block, sequence block,
  property namespace + datatype + applies_to, omitted-empty
  properties, year-from-date, phylogeny name + description, internal
  confidence rescale, type attribute per tool, schema element order,
  embed_alignment on/off, configured label format, segmented label
  uses strain, ladderize).
- `tests/test_segmented.py` (+2): first-non-empty inheritance picks
  up later-segment metadata; segment-0-wins when all set.
- `tests/test_config.py` (+8): labeling defaults, format type,
  null-segmented-format accepted, replace_whitespace bool check;
  phyloxml defaults, embed_alignment bool, all confidence_type
  values accepted, unknown confidence_type rejected.

Removed tests: two `_newick_to_phyloxml` direct tests in
`tests/test_phylo.py` — the helper is gone (the orchestrator now
calls `write_phyloxml` directly), and the dedicated writer test file
covers the same ground in more detail.

`v0.6.1` adds **IQ-TREE as the protein tree-builder**. New
`phylo.tool` knob (`auto` | `iqtree` | `fasttree`); `auto` (default)
picks IQ-TREE for protein alignments and FastTree for nucleotide.
IQ-TREE is invoked with **ModelFinder Plus** (`-m MFP`, scans
substitution models and picks by BIC — JTT vs WAG vs LG can change
topology on viral proteins) and **1000 ultrafast-bootstrap replicates**
(`-B 1000`, IQ-TREE's recommended minimum for interpretable support).
Binary auto-detected: tries `iqtree2` first, then `iqtree`; overridable
via `phylo.iqtree.binary`. The wrapper runs IQ-TREE under a temp
`--prefix` so its ~8 auxiliary files (`.iqtree`, `.log`, `.bionj`,
`.mldist`, `.contree`, `.splits.nex`, `.ckp.gz`) land in scratch and
get deleted; only the canonical `.treefile` (→ `{prefix}_tree.nwk`)
and the human-readable `.iqtree` model-selection report (→
`{prefix}_iqtree_summary.txt`) are kept. UFBoot refuses `<4`
sequences, so the wrapper auto-drops bootstrap with a stderr note when
the MSA has 3 reps (keeps the tree). `doctor` gains `iqtree2` /
`iqtree` to the external-binaries check (WARN-only). 27 new tests
cover IQ-TREE binary auto-detect (4: prefer iqtree2, fall back to
iqtree, override honoured, all-missing raises), argv construction
(model / threads / seed / UFBoot / extra_args, 3), the UFBoot
auto-skip below 4 seqs (1), config-disable of UFBoot (1), missing
summary path (1), subprocess failure → IQTreeError (1), missing
treefile detection (1), MSA record counting (1); dispatcher
auto-by-alphabet (2: iqtree-for-protein, fasttree-for-NT), explicit
tool override (1), IQ-TREE summary appended to output list (2:
present and absent), IQ-TREE-error → PhyloError wrapping (1); config
validation (8: defaults, tool override accepted/rejected, ufboot
non-negative + integer, model + binary string typing); doctor
binary check (2: dual-name accepted, missing-is-warn). 501 offline
tests total.

`v0.6.0` switches **clustering to amino-acid sequences by default**.
The new `clustering.alphabet` knob (`protein` default, `nucleotide`,
or `auto`) picks what's fed to MMseqs2 / cd-hit. The motivation is the
v0.5.x Orthobunyavirus run that drove `--min-seq-id` to 0.30 looking
for a target cluster count: synonymous substitutions inflate NT
divergence by 30–40% with no biological signal, the binary search
drifted into mmseqs2's sensitivity dead zone, and the reported cutoff
ceased to mean anything. Protein clustering avoids both problems —
reliable homology to ~25–30% identity, biologically meaningful
thresholds. The new `clustering.marker.select_marker_protein` picks
the marker per sequence (longest CDS by default; first matching
`cluster_protein` alias against `/product` as a case-insensitive
substring otherwise — alias order encodes preference). For segmented
viruses, each segment contributes its own marker; the per-isolate
concat lives on `concat.protein_sequence` and is what the clustering
backend sees. Isolates whose marker is missing on any segment are
dropped under `removed_incomplete_isolates`; non-segmented sequences
without a viable marker are dropped under `removed_proteins`. The
backends thread the alphabet through `_write_id_fasta`; the cd-hit
dispatcher's `_is_protein` consults `cfg["clustering"]["alphabet"]`
first, so the protein binary (`cd-hit`, 0.40 floor) is used even on
NT-typed CONCAT records carrying a `protein_sequence`. The new
`_setup_protein_alphabet` step auto-triggers `attach_proteins` if QC
didn't already (one-shot GenBank CDS fetch, same `ncbi_proteins` cache);
`--no-resolve` + `alphabet=protein` aborts at startup, `alphabet=auto`
silently falls back to `nucleotide`. Output gains
`{prefix}_representatives_protein.fasta` (AA strings actually fed to
clustering) alongside the existing NT outputs; `--phylo` builds the
MSA / tree on the AA strings when alphabet=protein actually fired
(single `_msa.fasta`, FastTree JTT). 37 new tests cover the marker
selector (8 cases: empty inputs, longest-CDS default, alias-order
preference, case-insensitive substring match, alias fallback,
translation-missing skip, all-translations-missing → None, alias-tie
length break), non-segmented `populate_protein_sequences` (3:
populates field, alias override, drops-no-marker sequences), segmented
`build_concatenated_sequences` protein path (4: longest-CDS-per-segment
fallback, per-segment alias selection, missing-marker drop, no-alias-
match fallback), `_setup_protein_alphabet` (6: nucleotide no-op,
auto-fetch when proteins missing, `--no-resolve` abort, auto fallback
to nucleotide, segmented skips per-sequence populate, `_resolve_alphabet`
auto branches), backend alphabet threading (`_write_id_fasta`
protein-body and missing-protein error; cd-hit `_is_protein` honours
alphabet override; cd-hit `min_threshold` floor follows override),
config validation (5: default alphabet=protein, nucleotide/auto
accepted, unknown rejected, non-list `cluster_protein` rejected,
per-segment `cluster_protein` accepted + unknown-segment + empty-alias
rejected), the writer's AA output (2: emitted when reps carry it,
skipped otherwise), and the phylo AA path (1: AA bodies + protein
model when alphabet=protein on NT-typed reps). 474 offline tests total.

Cache compatibility: v0.5.9+ `ncbi_proteins` cache entries already
carry `proteins[i]["sequence"]` (the GenBank `/translation`), so
existing runs benefit immediately with no refetch. Older entries
without translations return `None` from the marker selector and the
isolate is dropped — clear the `ncbi_proteins` cache to refresh.

`v0.5.10` makes segmented mode prefer **GenBank source-feature
qualifiers** (`/isolate`, `/strain`, `/segment`) over the
header-regex parse. Behind the new `segmented.use_genbank_metadata`
toggle (default `true`). The new
`NCBITaxonomy.fetch_source_metadata_batch` reuses the existing
`ncbi_proteins` SQLite cache — a run with protein QC and segmented
metadata extraction pays one efetch round trip, not two. Fetch happens
in a new CLI helper `_populate_genbank_isolate_segment` that runs after
`_run_protein_qc` and before `_handle_segmented`; UniProt sequences,
sequences without an accession, and `--no-resolve` runs all fall back
to the regex transparently (no warning, no error). `extract_isolate_id`
and `identify_segment` in `segmented/completeness.py` were already
written to short-circuit when `seq.isolate_id` / `seq.segment` are set
— no changes needed there. Cache entries written by v0.5.9 are
forwards-compatible: missing `source` key returns all-None and the
regex fallback fires. 13 new tests cover the source-feature parser
(qualifier present / absent), cache-sharing with protein QC,
legacy-cache forward compatibility, all six branches of the CLI
helper's gating, the strain-as-isolate fallback, and the bool config
validation. 437 offline tests total.

`v0.5.9` adds a **`repseq doctor`** self-test subcommand for
bench-scientist debugging. Emits a grouped report
(Python packages / external tools / network / configuration) with
`[OK]`/`[WARN]`/`[FAIL]` tags and a one-line summary; exits non-zero
only on `[FAIL]`. Policy: required Python packages (biopython, click,
PyYAML, requests) are `FAIL` if missing; optional `[viz]` extras
(umap-learn, matplotlib) are `WARN` (only needed for `--plot`); every
external binary (mmseqs, cd-hit, cd-hit-est, mafft, FastTree) is
`WARN` if missing (none is strictly required — pick a backend you
have, or use a diversity-only mode); NCBI Entrez + UniProt REST are
pinged with a 5s timeout and `WARN` on unreachable (you can run with
`--no-resolve`); cache directory unwritable is `FAIL`; missing
`taxonomy.ncbi_email` is `WARN` (works but rate-limited). The actual
import is attempted (not just `find_spec`) so a broken install
— e.g. a numpy/scipy ABI mismatch — is reported as missing rather
than silently passing. `--no-network` skips the database pings.
`_package_version` reads version via `importlib.metadata` so click 9.x
losing `__version__` is a non-issue. 17 new tests cover the WARN-vs-
FAIL policy, the network unreachable/timeout branches, cache-dir
write failure, config validation surfacing, and the click-CLI exit
codes. 424 offline tests total.

`v0.5.8` adds an optional **phylogeny step** behind the new
`--phylo` flag (works on every mode command). Builds an MSA with MAFFT
(`--auto`) and an approximate-ML tree with FastTree on the final
representative sequences (which in segmented mode are the concatenated
per-isolate sequences, so the tree's leaves are isolates not segments).
Every rep gets a deterministic short id (`S0001`…) before the MSA step
because long names, whitespace, and pipes break many phylo tools; the
final phyloXML restores each terminal clade's name to `seq.id` via
`Bio.Phylo`, while the intermediate Newick keeps the short ids
(decodable from `{prefix}_tree_id_map.tsv`). FastTree's substitution
model is auto-picked from the rep alphabet (`-nt -gtr` for nucleotide,
default JTT for protein). The step is fail-soft: skipped with a stderr
`[phylo skipped]` when there are `<3` reps, when `mafft` or `FastTree`
are missing, or when either subprocess errors — the rest of the run's
outputs are always written. New tests cover short-id round-trip, name
restoration in phyloXML (including pipe / non-ASCII originals), the
`<3` skip rule, NT-vs-AA model selection, and the `[phylo skipped]`
stderr path. 407 offline tests total.

`v0.5.7` adds an optional **cd-hit clustering backend** alongside
the existing MMseqs2 one. `cfg["clustering"]["backend"]` selects
`"mmseqs2"` (default, unchanged behaviour) or `"cdhit"`. The cd-hit
wrapper (`clustering/cdhit.py`) auto-picks `cd-hit` for protein input,
`cd-hit-est` for nucleotide; auto-picks `-n` (word size) from the
threshold per cd-hit's required table; writes input FASTAs through
`_write_id_fasta` (so cd-hit's whitespace-truncation cannot corrupt the
round-trip) and parses `.clstr` output via the `>id...` token marked
with `*`. A shared dispatcher (`clustering/__init__.py:run_clustering`)
routes modes between backends — modes now import from `..clustering`
instead of `..clustering.mmseqs2`. Identity floors differ:
mmseqs2 = 0.0, cd-hit (protein) = 0.40, cd-hit-est (nucleotide) = 0.80;
`clustering.min_threshold(cfg, sequences)` returns the active floor and
`_binary_search_threshold` clamps `lo` to it so the search never asks
the backend for a value it would refuse. New regression tests cover the
`.clstr` parser, auto-binary selection, threshold floor, dispatch, the
binary-search floor-clamp, and config validation of the new
`clustering.cdhit` block. 394 offline tests total.

`v0.5.6`. All 8 modes structurally complete, optional protein-annotation
QC (batched GenBank fetch + per-segment count check), a protein-FASTA
output reconstructed from cached records, per-segment nucleotide length
bounds (`segment_lengths` in virus config, applied after completeness
filter), and an optional UMAP scatter of the clustering result (`--plot`,
behind the `[viz]` extras: `matplotlib` + `umap-learn`). In segmented mode
the whole-pool QC length filter is skipped automatically (mixed segment
lengths make a pooled median meaningless and would drop short segments,
leaving every isolate incomplete) — use `segment_lengths` instead.
Exact-duplicate removal is likewise skipped on the segment pool and
applied instead to the concatenated per-isolate sequences in
`_handle_segmented` — a segment conserved across two distinct isolates
must not delete either isolate by leaving it "incomplete". Every
run ends with a one-line CLI summary (representatives/clusters selected, QC
pass rate) or, when nothing is selected, a stderr warning naming the most
likely cause (`cli._final_summary`). Every mode also records per-group
before/after counts (`RunResult.group_stats`) written to
`{prefix}_group_counts.tsv`. Taxonomic lineage is resolved via NCBI
`efetch` XML (`ncbi._parse_taxonomy_xml`) — the taxonomy *esummary*
endpoint carries no lineage for viruses, which used to group every
viral sequence under "Unknown". Offline regression tests cover IO,
QC, selector, segmented logic (including per-segment length filtering),
cache TTL, config validation, diversity selection, resolver fallback, mode
dispatch (clustering mocked), protein parser + filter, FASTA writer, segment
aliases, the closing CLI summary, per-group count reporting, and the
clustering plot (auto-skipped if `umap-learn` is missing). NCBI Entrez paths have been
live-tested against the H1N1 RefSeq genome (8 segments, 11 proteins).

A full pipeline audit (see git history around this note) fixed, with
regression tests: clustering ID round-trip for UniProt/CONCAT inputs,
RefSeq-nucleotide accession misrouting, the `MAG:` keyword filter, an
**inverted binary-search direction** (the threshold↔cluster-count
relationship was backwards, so `n-per-group` searches walked away from
the target), NCBI host/country/date now harvested from esummary
`subtype`/`subname` with the DB taking precedence over header heuristics,
a length-robust containment distance for diversity selection, thread-safe
caches and rate limiters, and assorted QC/parsing corrections.

`v0.5.5` hardens every place where header or identifier text crosses a
format boundary. The trigger was a Peribunyaviridae run where a 73-seq
genus came out of clustering as 73 reps at `cutoff = 1.0000`: isolate
names like `yaba-7 virus strain yaba 7` contain whitespace, so the
normalised id (and the derived `CONCAT|<isolate>` `seq.id`) carried
spaces; MMseqs2 truncates a FASTA header token at the first whitespace,
the parser then could not match cluster-TSV rows back to inputs, and
the binary search misread the empty cluster list as "≤ target" — falling
through to the final `best_reps = list(sequences)` fallback at `hi = 1.0`.
Fixes, with regression tests:
`_normalise_isolate_id` replaces both whitespace **and** the pipe
separator with `_` so `seq.id` survives MMseqs2 truncation and
`split("|")[1]` CONCAT parsing; `run_clustering` raises
`MMseqs2Error("Cluster round-trip mismatch...")` if the parsed cluster
membership does not account for every input; both FASTA writers
(`io/fasta.write_fasta`, `clustering/mmseqs2._write_id_fasta`)
collapse embedded `\r`/`\n` in headers to a space so a malformed
header cannot inject a phantom record; every TSV writer routes free-
text fields through a new `output/report._tsv_safe` helper that
neutralises tab and line-break characters so a stray value cannot
shift columns or rows. A 200-case brittleness probe
(`tests/test_unusual_characters.py`) exercises every (surface,
character) pair across whitespace and the punctuation classes most
common in viral metadata. 362 offline regression tests total.

`v0.5.6` enriches `{prefix}_isolate_proteins.tsv`. The previous six
columns (`isolate_id, segment, accession, protein_id, product, length`)
are reordered and extended to sixteen: `protein_id, product, length,
isolate_id, segment, segment_length, accession, species, subgenus,
genus, subfamily, family, suborder, order, subclass, class`. The
protein columns lead so the table reads as "what gene, on which
segment of which isolate, in what organism." `segment_length` is the
nucleotide length of the parent segment (`seq.length`). Taxonomy is
read via `TaxonomyInfo.get_rank(rank)`, which checks the standard
fields first and falls back to the resolver's lineage map (NCBI
`LineageEx`); the four sub-ranks (`subgenus`, `subfamily`, `suborder`,
`subclass`) have no standard field and come exclusively from the
lineage map — they are commonly blank for viruses, where ICTV often
skips intermediate ranks. Two regression tests cover the lineage-
backed sub-ranks and the missing-taxonomy fallback (all 9 rank cells
blank). 364 offline regression tests total.
