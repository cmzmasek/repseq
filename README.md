# repseq

**Pick a small, clean, representative set of sequences out of a big messy FASTA file.**

You downloaded 80,000 influenza sequences from NCBI. You need 200 good ones that
still cover the real diversity — for a tree, a figure, a training set, or to seed
a curated reference database. Doing that by hand is miserable. `repseq` does it
for you: it cleans up the data, looks up what each sequence actually is, groups
similar sequences together, and keeps one good example from each group.

It works on protein **or** nucleotide FASTA files, from UniProt, NCBI, or NCBI
Virus, with a strong focus on viral sequences (including segmented viruses like
influenza).

---

## What "representative selection" means here

If you have 5,000 nearly-identical H3N2 sequences and 3 unusual ones, a random
sample of 200 will be 200 H3N2 sequences and you'll lose the unusual ones.
`repseq` instead:

1. **Cleans** — drops duplicates, truncated junk, sequences full of `N`s, and
   records labelled "hypothetical", "synthetic", "partial", etc.
2. **Identifies** — looks up each sequence's organism, host, country, and
   collection date from NCBI/UniProt (results are cached, so it's only slow once).
3. **Groups** — buckets the sequences (by similarity, or by genus, host, year,
   country… your choice).
4. **Keeps the best example from each group** — preferring curated records
   (RefSeq, reviewed UniProt) over random ones, and longer over shorter.

The result is a FASTA file small enough to work with but still spanning the
diversity that was in the original.

---

## Installation

You need **Python 3.10 or newer**. Then:

```bash
git clone https://github.com/cmzmasek/repseq.git
cd repseq
pip install -e .
```

That's it for most uses. Four optional pieces:

**A sequence-clustering program** — needed only for the similarity-based modes;
if you group by genus/host/year and your groups are small, or you just want N
diverse sequences, you can skip it. Pick **one** backend; `repseq` finds the
binary (`mmseqs`, or `cd-hit` / `cd-hit-est`) on your `PATH` automatically.

```bash
# MMseqs2 (default — fast, scales to very large datasets)
brew install mmseqs2                       # macOS
conda install -c bioconda mmseqs2          # Linux (conda)

# cd-hit (classic alternative — set clustering.backend: cdhit). Slower on big
# inputs, has identity floors (0.40 protein / 0.80 nucleotide), but produces
# tight all-vs-all clusters many groups prefer for reference work.
brew install cd-hit                        # macOS
conda install -c bioconda cd-hit           # Linux (conda)
```

**Plots & tree figures** — the `[viz]` extra (matplotlib only, installs
everywhere) powers both the optional clustering scatter plot (`--plot`) and the
graphical tree-figure PDFs/PNGs that `--phylo` renders by default. For the
optional UMAP embedding upgrade (falls back to classical MDS if absent) add
`[viz-umap]`:

```bash
pip install -e '.[viz]'
pip install -e '.[viz-umap]'   # optional UMAP upgrade
```

**MSA + tree (`--phylo`)** — install MAFFT plus a tree builder. `repseq`
auto-picks **IQ-TREE** for protein alignments (ModelFinder + UFBoot, best
topology quality) and **FastTree** for nucleotide (much faster on large NT
MSAs). Pin one with `phylo.tool: iqtree` (or `fasttree`). If MAFFT or the tree
builder is missing, the run still finishes — the tree step is skipped with a
clear stderr message.

```bash
brew install mafft iqtree fasttree              # macOS
conda install -c bioconda mafft iqtree fasttree # Linux (conda)
```

---

## Quickstart

**Step 1 — make a config file** (cleaning thresholds, your NCBI email, etc.). A
short wizard asks a few questions, then writes a **complete, fully-commented**
config — every setting present at its default with an inline explanation, your
answers pre-filled:

```bash
repseq init-config -o my_config.yaml
```

**Step 2 — run a mode.** Point it at your FASTA and pick how to select:

```bash
# Keep one sequence per cluster of 90%-identical sequences
repseq global -c my_config.yaml -i seqs.fasta -T 0.90

# Keep up to 5 sequences per genus
repseq taxonomic1 -c my_config.yaml -i seqs.fasta -r genus -n 5

# Keep up to 10 sequences per host species
repseq host -c my_config.yaml -i seqs.fasta -n 10
```

When it finishes you'll see a one-line summary — how many sequences passed
cleaning and how many representatives were selected. If *nothing* came out, it
tells you the most likely reason (see [Troubleshooting](#troubleshooting)).

> **Tip:** add `--no-resolve` while experimenting. It skips the NCBI/UniProt
> lookups, so runs are fast — at the cost of not knowing each sequence's
> organism/host/country. Drop it for the real run.

---

## Choosing a mode

A "mode" is just *how you want the sequences grouped before one is picked from
each group*. All modes take `-i` (input), `-c` (config), and most take `-n`
(how many to keep per group).

| Command | What it does | Use it when… |
| --- | --- | --- |
| `global` | One big pass over everything. Either cluster at a similarity threshold `-T` (e.g. `-T 0.95` = group sequences ≥95% identical), or just ask for `-n` maximally-different sequences. | You want a flat, even sampling of the whole dataset. |
| `taxonomic1` | Up to `-n` sequences per taxonomic rank — `--rank genus`, `family`, `species`, etc. | You want even coverage across the tree of life (or of viruses). |
| `taxonomic2` | Like `taxonomic1` but **nested**: e.g. 20 per family, then 5 per genus within each. | One rank isn't enough — you want a hierarchy. |
| `host` | Up to `-n` sequences per host organism. | You care which host the virus came from. |
| `time` | Up to `-n` per time window — `--window year`, `decade`, or a number like `5` for 5-year bins. | You want even coverage across collection dates. |
| `geographic` | Up to `-n` per country. | You want even geographic coverage. |
| `custom` | Up to `-n` per *any* field — a sequence attribute, a taxonomy rank, or a column in a metadata spreadsheet you provide. | Your grouping isn't one of the built-ins. |
| `hybrid` | Up to `-n` per *combination* of fields, e.g. `--fields genus,host,decade`. | You want a balanced grid across several variables at once. |

**Within each group**, if it already has `-n` or fewer sequences, all are kept.
If it's bigger, `repseq` clusters it down to about `-n` representatives. Add
`--overflow trim` if you need *exactly* `-n`.

The `-n` "maximally-different" selection (`global -n`, and the `--overflow trim`
step) is an alignment-free **MaxMin** sampler over k-mer distance, run on the
**same representation your clustering uses**: with
`--alphabet-for-clustering protein` (the default) it compares each sequence's
**marker protein**, with `nucleotide` it compares the genome. Nucleotide k-mer
diversity is inherently weak for highly diverged viral families (short k-mers
saturate on whole genomes, long ones share nothing across genera), so
**protein is the recommended alphabet** when you want `-n` to track real
divergence.

Every mode also accepts: `--input/-i`, `--output-dir/-o`, `--config/-c`,
`--threads`, `--seed`, `--segmented`, `--dry-run`, `--no-resolve`,
`--source {auto,uniprot,ncbi,ncbi_virus}`, `--overflow {keep,trim}`, `--plot`,
`--phylo`, `--per-protein-phylo`, `--per-segment-phylo` (segmented only),
`--pre-cluster-tree`, `--fast`, `--newick/--no-newick`, `--pdf/--no-pdf`,
`--verbose`, `--alphabet-for-clustering {protein,nucleotide}`,
`--concatenate-markers/--no-concatenate-markers`,
`--protect-ids FILE` (force-keep a list through QC), `--pin-ids FILE`
(force-select a list as representatives), `--exclude-ids FILE` (blocklist a list
out of the input) — see
[Sequences of special importance](#sequences-of-special-importance--overrides).

Two diagnostic/utility subcommands:

| Command | What it does |
| --- | --- |
| `repseq stats -i x.fasta` | Pre-flight inspection of an input FASTA — count, source breakdown, NT length distribution, per-field metadata coverage, RefSeq/reviewed counts, top-N taxonomy per rank. No clustering, no output files. Default `--no-resolve` (fast); pass `--resolve` to fetch taxonomy. |
| `repseq replay <lockfile.json> -o <out>` | Re-materialise representatives from a `{prefix}_lockfile.json`: re-fetches the listed accessions (via the cache) and re-emits the same representative FASTAs into a fresh directory. QC/clustering/selection are **not** re-run. Trees skipped unless `--rebuild-trees`. |
| `repseq doctor` | Self-test the install (dependencies, external tools, network, config). |

---

## What happens during a run

```
your FASTA file(s)
        │
        ▼
  read it, figure out the format (UniProt / NCBI / NCBI Virus)
        │
        ▼
  look up organism, host, country, date          (from NCBI & UniProt, cached locally)
        │
        ▼
  clean (QC):
    drop duplicates → length filter → too many ambiguous chars
    → bad-keyword annotations → (optional) wrong protein count
    → (optional) bad-residue protein quality
        │
        ▼
  (segmented mode only)
    populate /isolate, /segment from GenBank
    → flag/drop strain-as-isolate collisions
    → drop isolates whose segments disagree on species (taxonomy_consistency)
    → keep only isolates with ALL expected segments
    → per-segment length bounds (drop the whole isolate if any segment is out of range)
        │
        ▼
  HMM-based identity QC                                          (when any spec has hmms:)
    scan every CDS with hmmscan → drop isolates/sequences whose
    expected markers don't have a satisfying CDS (one bad segment
    fails the whole isolate)
        │
        ▼
  (segmented mode only)
    fetch the marker protein per segment, concatenate per isolate
        │
        ▼
  group by your chosen mode (rank / host / decade / similarity / …)
  cluster each group down to ~N reps (binary search over identity threshold)
  pick the best representative per cluster (RefSeq > reviewed UniProt > longest)
        │
        ▼
  write representatives (FASTA + per-segment + concatenated), all proteins,
  metadata spreadsheets, per-marker / extra / polyprotein FASTAs, taxonomic
  reports (+ tidy TSV companions), a Methods-section starter, a run log, a
  config snapshot, and (optionally) clustering plot + trees
```

The run log (`{prefix}_run.log`) records exactly what settings were used and
what got dropped at each step — keep it with your results so the selection is
reproducible.

---

## Output files

Everything is written to a single output directory (`./repseq_output/` by
default; set via `output.dir` or `-o/--output-dir`). `repseq` **refuses to write
into a non-empty directory** so a run can never silently overwrite or mix into a
previous one — empty it, delete it, or point somewhere new.

The exact set depends on whether you ran in **segmented** mode and which optional
flags you passed. At a glance:

| File family | Always written? | What it is |
| --- | --- | --- |
| `{prefix}_run.log` | yes | Settings used and per-step counts. |
| `{prefix}_qc_removed.tsv` | yes | Every dropped sequence and the reason. |
| `{prefix}_overrides.tsv` | when any sequence was force-kept | Force-keep audit: each sequence `overrides.protect_qc` rescued from a QC stage, with the reason it would otherwise have been dropped. |
| `{prefix}_force_selected.tsv` | when any sequence was force-selected | Force-select audit: each `overrides.force_select` / `--pin-ids` pin and the action taken. |
| `{prefix}_excluded.tsv` | when any sequence was blocklisted | Input-blocklist audit: each `overrides.exclude` / `--exclude-ids` record dropped on load, plus any exclude id that matched nothing. |
| `{prefix}_group_counts.tsv` | yes | One row per stratum: in / out / clustered / cutoff. |
| `{prefix}_clusters.tsv` | yes | Which sequences ended up grouped, who represents whom. |
| `{prefix}_representatives.fasta` | non-segmented | Selected representative sequences. |
| `{prefix}_representative_sequences.tsv` | non-segmented | Per-representative metadata spreadsheet. |
| `{prefix}_concatenated.fasta` | segmented | Per-isolate head-to-tail concat of all segments. |
| `{prefix}_segment_<NAME>.fasta` | segmented | One file per expected segment, just the representative isolates. |
| `{prefix}_representative_isolates.tsv` | segmented | Per-representative-isolate metadata spreadsheet. |
| `{prefix}_isolate_proteins.tsv` | segmented + GenBank | Every protein of every isolate that survived QC, with a `representative` TRUE/FALSE column. |
| `{prefix}_representative_isolate_proteins.tsv` | segmented + GenBank | Same schema, row-filtered to representatives. |
| `{prefix}_sequence_proteins.tsv` | non-segmented + GenBank | Per-CDS counterpart of `_isolate_proteins.tsv`. |
| `{prefix}_representative_sequence_proteins.tsv` | non-segmented + GenBank | Reps-only filtered view of the above. |
| `{prefix}_representative_isolate_proteins.fasta` | segmented + GenBank | AA FASTA of every protein of every representative isolate. |
| `{prefix}_representative_sequence_proteins.fasta` | non-segmented + GenBank | AA FASTA of every protein of every representative sequence. |
| `{prefix}_per_protein_fasta/{prefix}_<family>.fasta` | any run with `cluster_protein` / `segment_markers` declared | Unaligned per-marker protein FASTA, one CDS per rep carrying that marker. |
| `{prefix}_extra_protein_fasta/{prefix}_<name>.fasta` | any run with `extra_protein:` declared | Same shape, for accessory proteins not required everywhere. |
| `{prefix}_polyprotein/{prefix}_<spec>_<peptide>.fasta` + `{prefix}_<spec>_peptides.tsv` | `polyprotein:` declared + HMM tier active | Mature peptides sliced from each rep's polyprotein CDS, one FASTA per peptide, plus an audit TSV. |
| `{prefix}_polyprotein/{prefix}_<spec>_polyprotein.fasta` | `polyprotein:` declared + HMM tier active | Unaligned **whole-polyprotein** FASTA — the entire parent CDS per rep. |
| `{prefix}_polyprotein/{prefix}_<spec>_polyprotein_tree.{xml,pdf,png}` (+ `_msa.fasta`, `_tree_id_map.tsv`) | `--per-protein-phylo` + `phylo.per_protein.whole_polyprotein_tree` | Whole-polyprotein tree; its phyloXML `<domain_architecture>` shows every HMM domain end-to-end. Off by default. |
| `{prefix}_representatives_protein.fasta` | when `alphabet_for_clustering: protein` (default) | The AA strings actually fed into the clusterer. |
| `{prefix}_taxonomic_report.txt` | every run | Per-rank diversity table: distinct taxa before vs after clustering. |
| `{prefix}_protein_taxonomic_report.txt` | any run with `cluster_protein` / `segment_markers` / `extra_protein` declared | Per-rank protein coverage + AA length statistics. |
| `{prefix}_nucleotide_taxonomic_report.txt` | every run | Per-rank NT length statistics: per-segment + `total` (segmented) or single `genome` column. |
| `{prefix}_polyprotein_taxonomic_report.txt` | `polyprotein:` declared + HMM tier active | Per-rank peptide coverage + AA length statistics. |
| `{prefix}_subtype_report.txt` / `.tsv` | only when the representatives span >1 viral subtype | Per-serotype distribution (e.g. influenza H5N1) before vs after clustering. |
| `{prefix}_*_taxonomic_report.tsv` | parallel to the matching `.txt` reports | Machine-readable tidy long-format companions (8-column schema). |
| `{prefix}_clustering.png` | only with `--plot` | Diagnostic scatter of the clustering. |
| `{prefix}_msa.fasta`, `_tree.xml`, `_tree_id_map.tsv` | only with `--phylo` | Alignment + annotated tree (phyloXML) + short-id↔accession map. |
| `{prefix}_tree.pdf`, `_tree.png` | `--phylo` (default; `phylo.pdf: true`) | Graphical tree figure (ladderized phylogram, taxonomy-coloured). On by default; `--no-pdf` disables. Emitted for *every* tree (`*_tree.pdf` / `*_tree.png`). Needs matplotlib (`[viz]`). |
| `{prefix}_tree.nwk` | `--phylo` **and** `phylo.newick: true` (or `--newick`) | Plain-text Newick. **Off by default** — the phyloXML is a topological superset. Same for every tree. |
| `{prefix}_partition.nex`, `_msa_<family>.fasta` | `--phylo`, protein + IQ-TREE | NEXUS partition file + per-family alignments (partitioned supermatrix). |
| `{prefix}_msa_untrimmed.fasta` | `--phylo` + `phylo.trimal.enabled` | Raw MAFFT alignment retained when trimAl trimming ran. |
| `{prefix}_taxonomy_review.tsv` | `--phylo` + `phylo.taxonomy_review.enabled` | Phylogeny-based taxonomy verdicts: blank ranks imputed, disagreeing ranks flagged. |
| `{prefix}_representative_*_corrected.tsv`, `_representative_*_proteins_corrected.fasta` | above + `write_corrected` | Copies with high-confidence imputed blanks filled (originals kept). |
| `{prefix}_iqtree_summary.txt` | `--phylo` + IQ-TREE | IQ-TREE ModelFinder report. |
| `{prefix}_iqtree_model.txt` | `--phylo` + IQ-TREE | Grep-friendly `<label>: <model>` sidecar, one line per partition. |
| `{prefix}_per_protein/` | only with `--per-protein-phylo` | One tree (MSA + phyloXML + id map + PDF/PNG; Newick opt-in) per marker; plus `_incongruence.tsv` of pairwise Robinson-Foulds distances. |
| `{prefix}_pre_cluster_tree.xml`, `_id_map.tsv` (+ `.pdf`/`.png`) | `--pre-cluster-tree` | Rough overview tree of every post-QC sequence; `[repr] ` prefix + `repseq:is_representative` boolean `<property>` on the leaves. |
| `{prefix}_extra_protein/` | `--per-protein-phylo` + `extra_protein:` declared | Accessory-protein trees (kept out of the incongruence table by design). |
| `{prefix}_per_segment/` | `--per-segment-phylo` (segmented only) | One **nucleotide** tree per declared segment, from the reps' raw per-segment NT. |
| `{prefix}_msa_conservation.tsv` | any phylo step that wrote an MSA (default on) | One conservation number per alignment (JSD-to-background, Henikoff-weighted, gap-penalised). |
| `{prefix}_monophyly.tsv` | any phylo step that built a tree (always on) | Per-taxon monophyly status off every annotated `*_tree.xml`. |
| `{prefix}_segment_status_matrix.tsv` | whenever `_monophyly.tsv` is written | Cross-tree pivot of the monophyly report; `single_marker_break` localises which taxon breaks on which one marker tree (reassortment). |
| `{prefix}_flags.txt` | when any conflict table exists | Plain-English synthesis of the conflict tables into one skimmable list. |
| `{prefix}_report.html` | any run with figures or a conflict table | Single-file HTML report: flags + embedded tree-figure gallery + output-file index. |
| `{prefix}_hmm_diagnostic.tsv` | HMM tier ran AND a spec declares `hmms:` | One row per (rep, marker spec, HMM profile): hit / best-E / coverage / cutoff / passing. |
| `{prefix}_summary.md` | every run | Auto-generated Methods-section starter (prose + numbers + tool citations), led by an "Analysis at a glance" table. |
| `{prefix}_config_repseq<ver>.yaml` | every run | Sanitized effective-config snapshot (every setting resolved, comments stripped, credentials blanked). Re-runnable as-is. |
| `{prefix}_lockfile.json` | every run | Machine-readable reproducibility record (rep IDs, accessions, post-run config, tool versions, input + HMM-DB sha256). Replayable via `repseq replay`. |

"Segmented + GenBank" means: segmented mode is on **and** the GenBank source
features are reachable (cached or fetched) so the per-isolate CDS list could be
built — true by default; falls away under `--no-resolve` or when
`segmented.use_genbank_metadata: false` and no protein QC ran.

The detailed sections below tell you what's in each file, what columns you'll
see, and what each is good for downstream.

### Sequence files (FASTA)

#### `{prefix}_representatives.fasta` — non-segmented mode

The main result for a non-segmented run: the selected representatives in
whatever alphabet you fed in. Header format matches the auto-detected input
format (NCBI, UniProt, or NCBI Virus). One record per representative. Use for:
the curated set, BLAST/DIAMOND database, downstream alignment.

#### `{prefix}_concatenated.fasta` — segmented mode

One record per **representative isolate**, with each isolate's segments joined
head-to-tail in the declared canonical order (`segmented.viruses.<v>.segments`,
e.g. `[L, M, S]`). Header is `>CONCAT|<isolate_id>|<acc1>|<acc2>|<acc3>`. Use
for: whole-isolate phylogenetics, gene-content overviews — anywhere you want
each isolate as a single entity.

#### `{prefix}_segment_<NAME>.fasta` — segmented mode

One file per expected segment (`_segment_L.fasta`, …). Each holds **only the
representative isolates' copies of that one segment**, with original headers
preserved. The per-isolate accessions are also listed in
`{prefix}_representative_isolates.tsv` under `accessions`, in segment order. Use
for: per-segment BLAST databases, one tree per segment (reassortment), per-segment
HMM builds.

#### `{prefix}_representative_isolate_proteins.fasta` *(segmented)* / `{prefix}_representative_sequence_proteins.fasta` *(non-segmented)*

Amino-acid FASTA covering **every annotated protein of every representative**,
reconstructed from the cached GenBank records (no extra network calls). Headers
carry NCBI-style bracket tags so each protein stays self-describing:

```
>NP_123456.1 [organism=Hantaan orthohantavirus] [ncbi_taxon_id=11599] \
   [species=Orthohantavirus hantanense] [genus=Orthohantavirus] \
   [family=Hantaviridae] [order=Elliovirales] [class=Bunyaviricetes] \
   [isolate=76-118] [segment=L] [host=Apodemus agrarius] \
   [country=South Korea] [collection_date=1976] [length=2151] \
   [parent=NC_005222.1]
```

Empty fields are omitted; bracket characters in values are scrubbed so
organism names like `Foo virus [strain X]` can't break the tag syntax. A
`[subtype=H5N1]` tag is added after `[host=...]` when the GenBank `/serotype`
qualifier is present. **This is usually "the" output file for protein-centric
workflows** — pre-curated and pre-annotated, ready for BLAST/DIAMOND/HMMER/MMseqs2.

#### `{prefix}_per_protein_fasta/{prefix}_<family>.fasta` — one file per declared marker

Unaligned protein FASTA, one file per marker spec under
`clustering.cluster_protein` (non-segmented) or `segment_markers` /
`cluster_protein` (segmented). Each holds **one record per representative
carrying that marker** — the CDS that satisfied the marker's HMM (or the alias
chain if no `hmms:` is declared). `<family>` is the spec's `name:`
(e.g. `Spike`), else the single token (`Bunya_nucleocap`), else the first token
+ `_altN`; segment-prefixed in segmented mode (`M_Spike`).

Headers are **byte-identical** to the all-protein FASTA above, so a Spike record
in either file is the same string. Written on every run (independent of
`--per-protein-phylo`) so you always have the per-marker input FASTAs. Specs no
representative satisfies are skipped (no empty file); an HMM-off run against
HMM-only specs emits one stderr note.

#### `{prefix}_extra_protein_fasta/{prefix}_<name>.fasta` — when `extra_protein:` declares entries

The **accessory-protein** analogue of the per-protein FASTA: same selection
chain, same bracket-tag headers. Driven by `clustering.extra_protein:`
(non-segmented) or `virus.extra_protein:` (segmented). Filename uses the spec's
`name:` verbatim (no segment prefix; names must be unique). Sparse coverage is
expected, so unsatisfied specs are silently skipped. See
[Accessory proteins](#accessory-proteins-extra_protein).

#### `{prefix}_representatives_protein.fasta` — when `alphabet_for_clustering: protein` fired

The **AA strings fed into the clusterer** — the per-isolate marker concat
(segmented) or per-rep marker (non-segmented; the multi-marker concat when
`clustering.concatenate_markers: true`). A *diagnostic* for reproducing/auditing
the clustering input; not written for `alphabet_for_clustering: nucleotide`. For
a clean per-protein set use `_representative_*_proteins.fasta` or
`_per_protein_fasta/`.

### Spreadsheets (TSV)

All TSVs use a harmonised vocabulary: the canonical identifier column is always
`accession`; booleans are `TRUE`/`FALSE`; length columns carry their alphabet
(`length_nt`, `length_aa`, `segment_length_nt`, `total_length_nt`); taxonomic
ranks always appear in the same nine-rank ladder (`species`, `subgenus`, `genus`,
`subfamily`, `family`, `suborder`, `order`, `subclass`, `class`). Sub-ranks come
from the NCBI lineage and are commonly blank for viruses.

#### `{prefix}_representative_sequences.tsv` — non-segmented mode

One row per representative sequence. **Schema-identical to
`_representative_isolates.tsv`** so the same script reads both modes. Columns:

`isolate_id`, `isolate_id_source`, `organism`, `strain`, `host`, `subtype`,
`collection_date`, `country`, `n_segments`, `segments`, `accessions`,
`total_length_nt`, `is_refseq`, `is_reviewed`, `ncbi_taxon_id`, then the
nine-rank taxonomic ladder. (`subtype` is the viral serotype, e.g. `H5N1`, from
the GenBank `/serotype` qualifier.)

In non-segmented mode the isolate-only cells are blanked (`isolate_id`,
`isolate_id_source`, `n_segments`, `segments`) or remapped: `accessions` is the
single accession, `total_length_nt` is the sequence's NT length. Use
`_sequence_proteins.tsv` for per-CDS detail.

#### `{prefix}_representative_isolates.tsv` — segmented mode

One row per representative **isolate**. Columns:

`isolate_id`, `isolate_id_source`, `organism`, `strain`, `host`, `subtype`,
`collection_date`, `country`, `n_segments`, `segments` (comma-joined names in
concat order), `accessions` (comma-joined per-segment accessions in concat
order), `total_length_nt`, `is_refseq`, `is_reviewed`, `ncbi_taxon_id`, then the
nine-rank ladder.

`isolate_id_source` records where the `isolate_id` came from: `isolate` (GenBank
`/isolate` qualifier), `strain` (`/isolate` absent, `/strain` used as a fallback
grouping key — see [Strain-as-isolate fallback](#strain-as-isolate-fallback-and-collision-detection)),
or `regex` (header-regex fallback — typically `--no-resolve`, UniProt input, or
no accession). The per-sequence columns (`accession`, `segment`, `description`,
`molecule_type`, `length_nt`) are absent — they have no isolate-level meaning;
use `_isolate_proteins.tsv` for per-segment detail.

#### `{prefix}_isolate_proteins.tsv` — segmented mode (when proteins are reachable)

One row per **CDS** of every isolate that survived QC, whether or not it was
picked. Columns: `protein_id`, `product`, `length_aa`, `isolate_id`,
`isolate_id_source`, `segment`, `segment_length_nt`, `accession`,
`representative` (`TRUE` if the isolate made the final set), `hmmscan`, then the
nine-rank ladder. Use for: a full per-protein QC audit, and to join a downstream
hit back to its isolate.

#### `{prefix}_representative_isolate_proteins.tsv` — segmented mode (when proteins are reachable)

Same schema as `_isolate_proteins.tsv`, filtered to representative isolates only
(every row `representative=TRUE`).

#### `{prefix}_sequence_proteins.tsv` — non-segmented mode (when proteins are reachable)

The non-segmented counterpart of `_isolate_proteins.tsv`: one row per **CDS** of
every sequence that survived QC, **same column schema** so one script reads both
modes. The isolate-only columns are blanked; `segment_length_nt` holds the
parent sequence's NT length and `accession` is the parent accession.
`representative` is `TRUE` when the parent sequence ended up a representative.

#### `{prefix}_representative_sequence_proteins.tsv` — non-segmented mode (when proteins are reachable)

Row-filtered companion to `_sequence_proteins.tsv` (reps only). The
non-segmented analogue of `_representative_isolate_proteins.tsv`.

#### `{prefix}_clusters.tsv` — always

One row per cluster (the row shows the cluster's representative). Columns:

| Column | What it is |
| --- | --- |
| `cluster_id` | Stable cluster identifier. |
| `accession` | Accession of the representative. |
| `organism` | Representative's organism. |
| `cluster_size` | Total members (representative + others). |
| `is_refseq` / `is_reviewed` | Representative's quality flags. |
| `species_purity` / `species_distinct` | Fraction of populated species labels matching the most common one (1.000 = pure) and the count of distinct labels. |
| `genus_purity` / `genus_distinct` | Same at genus. **`genus_distinct > 1` is the actionable "threshold too lax for this stratum" signal.** |
| `host_purity` / `host_distinct` | Same against `/host` — flags cross-host contamination / mis-annotation. |
| `member_length_min` / `_max` / `_median` | NT-length distribution across members; a wide range can flag a length-outlier inclusion. |

Purity is computed over populated values only (a member with no resolved genus
doesn't count against genus purity). Blank purity (with `_distinct = 0`) means no
member carries a value at that rank. Use for: tracing which sequences clustered
together, or re-running selection with a different priority without re-clustering.

#### `{prefix}_group_counts.tsv` — always

One row per stratum. Columns: `stratified_by`, `stratum`,
`stratum_size_before`, `stratum_size_after`, `clustered` (did the
binary-search clusterer run?), `cutoff` (the identity threshold it settled on).

When `clustering.diversity_curve_cutoffs` is non-empty (default
`[0.99, 0.95, 0.9, 0.8, 0.7]`), every clustered stratum also carries trailing
`n_clusters_<c>` columns — the cluster count at that fixed threshold. This is a
**diagnostic only** (selection is unaffected): it lets you read off how
conserved a group is. Cells are `NA` for cutoffs below the backend's identity
floor and empty for unclustered strata. **The fastest way to see where the
reduction happened.**

#### `{prefix}_qc_removed.tsv` — always

Every sequence dropped during cleaning, with the reason. Columns: `accession`,
`reason`. Typical reasons:

- `duplicate` — exact duplicate of another record
- `length:<n><[<>]<m>` — failed the whole-pool length filter (non-segmented)
- `segment_length:<seg>:<n><[<>]<m>` — failed per-segment length bounds (whole isolate dropped, so all its segments appear)
- `ambiguous_fraction:<f>>0.05` — too many `N`/`X`/ambiguous characters
- `annotation_keyword:<word>` — description contains a blacklisted keyword
- `protein_count:<n><[<>]<m>` — wrong number of annotated CDS
- `segmented_filter:could_not_identify_isolate_or_segment` — GenBank lookup and `isolate_regex` both failed
- `incomplete_isolate:missing_segments:<list>` — missing one or more expected segments
- `taxonomy_mismatch:<rank>` — segments of one isolate disagreed at the configured rank
- `strain_collision:<segment>` — distinct accessions shared a strain-derived isolate id (when `strain_collision_action: drop`)
- `hmm_failed:<seg>:<token>` — no CDS satisfied the marker's tokens (HMM-QC)

**Check this file first when more was dropped than expected.** Sort by `reason`
to see which filter did the damage; the same info is summarised in `_run.log`.

#### `{prefix}_overrides.tsv` — when a sequence was force-kept

The flip side of `_qc_removed.tsv`: every sequence the **force-keep** list
(`overrides.protect_qc`) rescued from a QC removal stage. Columns: `id`, `stage`,
`would_be_reason`. Only written when at least one sequence was protected. See
[Sequences of special importance](#sequences-of-special-importance--overrides).

#### `{prefix}_force_selected.tsv` — when a sequence was force-selected

The selection-side audit for `overrides.force_select` / `--pin-ids`. Columns:
`id`, `action`, `detail`. `action` ∈ `elected_representative` (won its cluster),
`split_singleton` (split out because another pin won the same cluster),
`added_representative` (clustering had deselected it), `already_representative`,
or `unavailable` (no surviving sequence matched). Only written when at least one
sequence was pinned.

#### `{prefix}_excluded.tsv` — when a sequence was blocklisted

The input-blocklist audit for `overrides.exclude` / `--exclude-ids`. Columns:
`id`, `action`, `detail`. `action` ∈ `excluded` (a record dropped on load, with
`detail` naming the field it matched on) or `unavailable` (an exclude id that
matched no input sequence — a typo-catcher). Only written when at least one
record was excluded. See
[Excluding bad sequences](#excluding-bad-sequences).

#### `{prefix}_hmm_diagnostic.tsv` — when the HMM tier ran

One row per (representative, marker spec, individual HMM profile). Multidomain
tokens like `"CoV_S1--CoV_S2"` expand to one row per component HMM. Columns:

| Column | What it is |
| --- | --- |
| `isolate_id` / `accession` / `segment` | Identity of the rep (segment blank in non-segmented mode). |
| `spec_name` | The marker spec's `name:` (or its single token). |
| `profile` | One declared HMM profile from any token in the spec. |
| `hit` | `TRUE` if any CDS in scope hit this profile. |
| `best_dom_evalue` / `best_dom_score` / `best_coverage` | The lowest-E hit's stats (blank when `hit=FALSE`). |
| `ga_cutoff` | The profile's Pfam GA threshold (blank when none). |
| `cutoff_used` | `GA=<value>` when available, else `default_evalue=<value>`. |
| `passing` | `TRUE` if the best hit cleared BOTH the similarity and coverage gates. |

The **diagnostic that surfaces near-misses the binary gate hides** ("65 isolates
barely missed RdRP_4 at E=2e-5"). Sort by `passing`; group by `profile` to see
whether an HMM is too strict for your dataset.

### Run record

#### `{prefix}_run.log`

Plain-text record: version, date, exact command line, every active config
setting, per-step QC counts, the mode's results, and the output file list. Save
it with your data. The "QC SUMMARY" block also shows the per-segment
length-filter breakdown (`L too short : N`, …), so you can tell *which* segment
caused isolates to fall out.

#### `{prefix}_config_repseq<ver>.yaml` — sanitized effective-config snapshot

A clean, shareable copy of the configuration the run actually used — written on
every run. Because repseq deep-merges your YAML over its built-in `DEFAULTS`,
this shows **every setting at the value it ran with**. Sanitized three ways:
comments stripped; credentials (`taxonomy.ncbi_email` / `ncbi_api_key`) blanked
to `null` (keys kept, so it's still valid); runtime state (any `_`-prefixed key)
dropped. The filename embeds the version (`cov_config_repseq0_56_0.yaml`), so a
directory of runs self-documents which version produced each. Re-loadable
verbatim: `repseq <mode> -c <this-file> -i <input>`.

#### `{prefix}_lockfile.json` — structured reproducibility record

A machine-readable companion: where `_summary.md` records *what was done* in
prose, the lockfile records *what came out* in JSON — the elected representative
IDs, their parent accessions, the post-run config (credentials blanked), the
version of every external tool, an sha256 of every input FASTA, and (when the
HMM tier ran) an sha256 of the HMM database. Always written.

```jsonc
{
  "schema_version": 1,
  "repseq_version": "0.56.0",
  "created_utc": "2026-06-07T12:00:00Z",
  "mode": "taxonomic1",
  "command": "repseq taxonomic1 ...",
  "python_version": "3.11.6",
  "config": { ... },                       // post-run cfg
  "inputs": [{"path": "in.fasta", "sha256": "…", "size_bytes": 12345}],
  "tools": {"mafft": "v7.520", "iqtree2": "2.3.2", ...},
  "hmm_db": {"path": "...", "sha256": "...", "bundled": true, "n_profiles": 26},
  "representatives": [
    {"id": "AB1234", "kind": "sequence", "accession": "AB1234", "organism": "Virus X"},
    {"id": "CONCAT|iso1", "kind": "isolate", "isolate_id": "iso1",
     "organism": "Virus Y", "segment_accessions": {"L": "AB1", "M": "AB2", "S": "AB3"}}
  ]
}
```

Sequences are **not** embedded — they're re-fetched by `repseq replay`. Cluster
membership is **not** recorded (the authoritative source is
`{prefix}_clusters.tsv`). Written with sorted keys + 2-space indent, so two
lockfiles from identical pipelines diff cleanly on `created_utc` and nothing else.

#### `repseq replay <lockfile.json>` — re-materialise the representatives

Reads a lockfile and re-emits the representative FASTAs into a fresh directory,
without re-running QC, clustering, or selection:

```bash
repseq replay run1/run1_lockfile.json -o replay_of_run1/
repseq replay run1/run1_lockfile.json -o replay_of_run1/ --rebuild-trees   # slower
```

It validates `schema_version` / `repseq_version` (major mismatch is a hard
error, minor a stderr note), verifies the HMM-DB sha256 if the original used the
HMM tier, re-fetches each accession from NCBI via the cache (source
`ncbi_nuc_seq`; for segmented reps fetches all segments and rebuilds the CONCAT),
re-emits the representative FASTAs, optionally rebuilds trees, and writes a short
`{prefix}_replay.md`. When NCBI can't return an accession, it **warns loudly and
continues with survivors** — the dropped accessions land in
`{prefix}_replay_missing.tsv`.

### Taxonomic reports

Plain-text reports that turn "what kind of diversity did I select?" into a
one-glance answer. Each `_*_taxonomic_report.txt` begins with a `#`-commented
provenance line (version, mode, dataset type, clustering substrate + tool, rep
count), so an isolated report is self-describing. The TSV companions carry the
same facts as rows.

#### `{prefix}_taxonomic_report.txt` — diversity before vs after

Per-rank counts of **distinct taxa** in the pool fed to selection versus in the
final representatives. Counting unit is **isolates** (segmented) or **sequences**
otherwise. Section 1 is a 9-rank ladder with two columns (distinct before /
after). Section 2 is, per rank with ≥1 populated taxon, the per-taxon breakdown
(name, before-count, after-count, sorted by before-count desc; truncated to the
top `output.protein_report.max_breakdown`, default 20). Blank rank values are
excluded.

#### `{prefix}_subtype_report.txt` / `.tsv` — serotype distribution

The subtype analogue of one rank of the taxonomic report, over `seq.subtype`
(the viral **serotype**, e.g. influenza `H5N1`, from the GenBank `/serotype`
qualifier in segmented mode or the esummary path otherwise). **Written only
when the final representatives carry more than one distinct subtype** — a
single-serotype (or serotype-less) selection produces nothing. A header
coverage line says how many of the pool / reps carry a serotype at all (blank
/ unknown values are not a subtype and never appear as a row), then a
distinct-count line, then a per-subtype table (`subtype`, before-count,
after-count, sorted by before-count desc, truncated to the top
`output.protein_report.max_breakdown`). Counting unit is **isolates**
(segmented) or **sequences** otherwise. The `.tsv` is the tidy long-format
companion (same 8-column schema as the other report TSVs, `report=subtype`).
Use for: "how are the H/N subtypes distributed in my influenza set, and did
clustering preserve the rare ones?"

#### `{prefix}_protein_taxonomic_report.txt` — per-marker coverage + length

Per-rank tables for **each declared protein** (every `cluster_protein` /
`segment_markers` / `extra_protein` spec). For each rank from `subgenus` to
`class` (skipping `species`, dominated by annotation noise), four sub-tables:

1. **Coverage (post-QC pool)** — cells `<count> <pct>%`: how many in that taxon had a CDS satisfying the marker.
2. **Coverage (representatives)** — same, on the final set.
3. **Length statistics (post-QC pool)** — cells `min, max, median, Q3-Q1, n` (amino acids; Q3-Q1 = IQR; n = items contributing the length).
4. **Length statistics (representatives)** — same on the reps.

Marker columns are headed by the spec's `name:`; **cluster-driving markers carry
a trailing `*`** so you can tell load-bearing proteins from accessory ones. A
`== HMM marker architectures ==` block lists each HMM spec's token alternatives
(joined with `OR`) and the cutoff policy that gated it. Taxa are truncated by
member count (top `max_breakdown`) then sorted alphabetically. Use for: "is my
rep set still covering polymerase / glycoprotein / nucleocapsid across the
families I started with, and is the length distribution still sensible?"

#### `{prefix}_nucleotide_taxonomic_report.txt` — per-rank NT length

The nucleotide analogue: per-rank length-statistics tables on the input
nucleotide sequences. Two sub-tables per rank (post-QC pool, representatives) in
the same `min, max, median, Q3-Q1, n` format. Columns: **segmented** = one per
segment in configured order + a trailing **`total`** (sum per isolate);
**non-segmented** = a single **`genome`** column. Same rank scope and truncation
as the protein report. **No coverage tables** — every passing entity carries
every required nucleotide unit by construction (completeness QC guarantees it).
Use for: catching a rep set that's lost length distribution against the pool.

#### `{prefix}_polyprotein_taxonomic_report.txt` — per-peptide coverage + length

The **sliced-peptide analogue** of the protein report — emitted only when a
`polyprotein:` spec is declared and the HMM tier is active. One H2 section per
spec; columns are the spec's declared peptides in N→C order (a declared peptide
stays a column even at 0% coverage — itself an audit signal). Four sub-tables per
rank (same ladder, same truncation): coverage pool / coverage reps / peptide
length pool / peptide length reps. A `== Peptide architectures ==` block lists
each peptide's token alternatives and `cleavage_motif`. "Covered" counts peptides
with `status` `ok` or `overlap` (the same statuses that produce per-peptide FASTA
records). Use for: spotting genus-specific peptide loss (e.g. coronavirus NSP1
silently dropping out of one subgenus).

#### `{prefix}_*_taxonomic_report.tsv` — tidy long-format companions

Each of the four `.txt` reports has a **machine-readable tidy long-format TSV**
companion — one row per observation, no banners — designed for Excel pivot
tables, R (`readr` + `tidyr`), or pandas. All four share one 8-column schema, so
they concatenate and key on `report`:

| column | values |
| --- | --- |
| `report` | `diversity`, `protein`, `nucleotide`, `polyprotein` |
| `rank` | the 9 `_TAX_RANKS` |
| `pool` | `post_qc` or `reps` |
| `taxon` | taxon name at that rank; `*ALL*` for rank-level diversity rows |
| `taxon_count` | items in that (rank, pool, taxon) cell — the denominator of coverage_pct |
| `spec` | marker / extra-protein name (cluster-driving carry `*`); `<polyprotein>:<peptide>`; `genome` / `<segment>` / `total`; `_diversity` |
| `metric` | `distinct_taxa`, `member_count`, `coverage_count`, `coverage_pct`, `length_min`, `length_max`, `length_median`, `length_iqr`, `length_n` |
| `value` | numeric (integers as integers; percentages up to 4 decimals, trailing zeros trimmed) |

```r
library(tidyverse)
df <- read_tsv("coronavirus_protein_taxonomic_report.tsv")
df %>% filter(rank == "genus", pool == "reps", metric == "coverage_pct")
```
```python
import pandas as pd
df = pd.read_csv("coronavirus_polyprotein_taxonomic_report.tsv", sep="\t")
df.query("metric == 'length_median' and pool == 'reps'")
```

**Zero-coverage handling.** When a marker/peptide/segment is absent from a taxon
the row emits `coverage_count=0, coverage_pct=0, length_n=0` and **no other
length metrics** — empty distributions aren't fabricated. Filter on
`length_n > 0` before computing across-taxon length summaries. The architecture /
banner sections of the `.txt` reports are humans-only and have no TSV equivalent.

### Diagnostic plot

#### `{prefix}_clustering.png` — only if `--plot` is passed

Two-panel scatter of the clustered sequences (k-mer Jaccard distance):

- **Left** — every (sub)sampled sequence as a dot, coloured by genus.
- **Right** — the same dots, coloured by cluster, sized by cluster population,
  with faint lines from each member to its representative.

The embedding uses **UMAP** when `umap-learn` imports cleanly, else a numpy-only
**classical MDS** (PCoA) on the same distance matrix; the figure labels its axes
with whichever ran. Only `matplotlib` is required (`'.[viz]'`); the UMAP upgrade
is `'.[viz-umap]'`. Bounded by a default 2000-point subsample (representatives
always kept); skipped when the run produced no clusters (`global -n`).

### Phylogeny outputs

Written only with `--phylo`. The pipeline: short-id remap → MAFFT (`--auto`) →
*(optional)* trimAl → tree builder → root → label internal nodes by LCA →
phyloXML. The tree builder is auto-picked from the clustering alphabet:
**IQ-TREE for protein** (ModelFinder + UFBoot) and **FastTree for nucleotide**
(`-nt -gtr`). Pin one with `phylo.tool: fasttree` (or `iqtree`).

MAFFT / IQ-TREE / FastTree each print one start and one finish line so the
terminal stays readable. Pass **`--verbose`** to also stream live child stderr
(prefixed `[mafft]` / `[iqtree]` / `[fasttree]`) — useful for a step that looks
stuck. The stderr is buffered either way, so an on-failure error still contains
the full subprocess output.

##### Plain-text Newick is opt-in: `phylo.newick` / `--newick`

Each tree is written as annotated **phyloXML** (`*_tree.xml`) and, when asked,
as plain-text **Newick** (`*_tree.nwk`). The Newick is **off by default** to keep
directories tidy — the phyloXML carries the same topology plus taxonomy,
metadata, LCA labels, and colours. Set `phylo.newick: true` (or `--newick`) to
retain the `*_tree.nwk` for **every** tree of the run; `--no-newick` forces it
off. Turning it off costs nothing scientifically: the Newick is still produced
internally (the phyloXML is parsed from it, and the incongruence RF table reads
it) and only deleted at the end; `_tree_id_map.tsv` and `_msa.fasta` are kept
regardless.

##### Graphical tree figures: `phylo.pdf` / `--pdf` / `--no-pdf`

Every phyloXML tree is **also rendered as a figure** — a `*_tree.pdf` (vector)
and `*_tree.png` (150 dpi) sibling next to each `*_tree.xml`. The figure is a
ladderized rectangular phylogram (matplotlib + Bio.Phylo): taxonomy-coloured leaf
labels, a genus/subfamily colour legend, internal-node labels for the
genus→family ranks, and branch-support labels for nodes with support ≥ 50.
Everything is reconstructed **from the phyloXML itself**, so figure colours match
the tree annotation exactly, and one end-of-run sweep covers *every* tree type.
**On by default** (`phylo.pdf: true`); pass **`--no-pdf`** to suppress (the
phyloXML is retained either way and opens in e.g.
[Archaeopteryx](https://sites.google.com/view/archaeopteryx)). Needs matplotlib
(`[viz]`); soft-fails with a single note if missing.

##### Preliminary-run shortcut: `--fast`

When iterating on a config, the slow/strict tree pipeline is overkill. `--fast`
overrides the `phylo:` block of your YAML in memory (your file is untouched) and
forces: **FastTree** for every tree (no IQ-TREE/ModelFinder/UFBoot); **MAFFT
`--retree 1`** (single-pass FFT-NS-1, drops `--auto` and any L-INS-i); **no
trimAl**; **no partitioned supermatrix**; **midpoint rooting**. A banner is
printed at startup and `_summary.md` describes the fast pipeline. **Switch
`--fast` off for the final publication run** — trimAl, IQ-TREE/ModelFinder,
UFBoot, and the partitioned supermatrix all change topology meaningfully on real
data.

##### Alignment trimming (trimAl, optional, off by default)

Enable `phylo.trimal.enabled: true` to trim poorly-aligned / gap-rich columns
**before** tree inference with [trimAl](http://trimal.cgenomics.org/) (it sits
between MAFFT and the tree builder). `mode` is the trimAl method (default
`automated1`; also `gappyout`/`strict`/`strictplus`/`nogaps`/`noallgaps`), and
`extra_args` is a raw passthrough for threshold trimming (`["-gt", "0.8"]`). The
per-protein trees have an independent `phylo.per_protein.trimal`. When trimming
runs, `{prefix}_msa.fasta` becomes the trimmed (tree-input) alignment and the
raw MAFFT output is kept as `{prefix}_msa_untrimmed.fasta`; in partitioned mode
each per-family alignment is trimmed before concatenation. If trimal is missing
or strips the alignment to nothing, repseq warns loudly and builds on the
**untrimmed** alignment rather than failing. The trimal version + mode are
recorded in the phyloXML `<description>`.

##### Partitioned supermatrix (default for protein + IQ-TREE)

For protein + IQ-TREE runs, the whole-genome tree is built as a **partitioned
supermatrix** rather than by gluing markers into one string. Each declared marker
family (one per `hmms:` spec) is aligned **separately** with MAFFT, the per-family
alignments are concatenated **column-wise** into a supermatrix, and IQ-TREE fits a
substitution model **per partition** (ModelFinder). This is the statistically
correct multi-marker analysis: MAFFT never aligns an L-polymerase against an
M-glycoprotein, and proteins no longer share one model.

**On by default** (`phylo.partition.enabled: true`); applies only when the run
resolves to protein + IQ-TREE **and** the HMM tier resolved ≥2 marker families.
Otherwise it falls back transparently to the single-alignment, single-model path.
Knobs under `phylo.partition`:

- `linkage` — IQ-TREE branch-length linkage: `proportional` (`-p`, the default),
  `equal` (`-q`), or `unlinked` (`-Q`, independent per partition).
- `models` — optional per-family model pin keyed by family label
  (`{L_RdRP_4: "LG+G4", S_Bunya_nucleocap: "WAG+G4"}`). **All-or-nothing**: pins
  are honoured only when *every* family is pinned (the NEXUS then carries a
  `charpartition` and IQ-TREE runs without `-m`); any unpinned family falls back
  to per-partition ModelFinder with a warning. With no pins (the default), the
  NEXUS is charsets-only and IQ-TREE runs with `-m MFP`.

Extra outputs: per-family alignments `{prefix}_msa_<family>.fasta` and the NEXUS
partition file `{prefix}_partition.nex`. `{prefix}_msa.fasta` is then the
concatenated supermatrix.

#### `{prefix}_msa.fasta`

The MAFFT alignment (the concatenated supermatrix in partitioned mode). Headers
are `>S0001 <pretty-label>` — the short `SNNNN` id stays the first
whitespace-separated token (safe for any phylo tool), and the descriptive label
(from `phylo.labeling.format` / `segmented_format`) is appended so AliView /
Jalview / MEGA show recognisable names.

#### `{prefix}_tree.nwk` — off by default

Tree-builder Newick, short-id leaves (`S0001`, …). Use the id map to recover real
names. **Not retained by default** (`phylo.newick: false`); see "Plain-text
Newick is opt-in" above.

#### `{prefix}_tree.xml` — the tree you'll usually open

phyloXML with rich, browseable annotation:

- Every leaf gets a formatted `<name>`, a `<taxonomy>` block with NCBI taxon id,
  a `<sequence>` block with the GenBank accession + title, and one
  repseq-namespaced `<property>` per leaf attribute (host, subtype,
  collection_date, country, strain, isolate_id, year, plus the full 9-rank
  ladder — empties dropped).
- Leaf labels are **coloured by taxonomy** (default on, by genus): each carries a
  `<property ref="style:font_color">` hex colour. Configure under
  `phylo.coloring` — one rank gives each value a distinct hue; two ranks
  (`[genus, subgenus]`) shade each subgenus within its genus colour. Unresolved
  taxa show grey (`missing_color`). The same palette is shared across `--phylo`
  and every `--per-protein-phylo` tree, so a genus is the same colour
  everywhere — which is what makes cross-tree incongruence (reassortment) jump
  out visually. `saturation` and `value` (both 0–1, defaults `0.65`/`0.90`) are
  the HSV vividness and brightness held constant across the palette. Set
  `phylo.coloring.enabled: false` to turn it off.
- Every annotated internal clade gets a `<name>` and a `<taxonomy>` block with
  the LCA's scientific name + rank (`min_rank=genus` by default; each clade is
  labelled at its crown, and a non-monophyletic taxon split across disjoint
  clades is labelled on each crown).
- The `<phylogeny>` carries a `<name>` and a `<description>` that leads with a
  plain-English "what this tree is based on" sentence — e.g. *"This tree is based
  on a per-isolate concatenation of one marker protein per segment (L:RdRp,
  M:GPC, S:N); amino-acid; each leaf is one representative isolate."* — followed
  by MAFFT/IQ-TREE/FastTree versions, the selected model, bootstrap settings, and
  the rooting method that fired. The same substrate is also machine-readable as
  phylogeny-level `<property>` elements (`repseq:tree_basis`,
  `repseq:analysis_mode`, `repseq:substrate`, `repseq:alphabet`,
  `repseq:leaf_unit`). Every tree type carries its own basis.
- Confidence values are normalised to 0–100 integers (`sh_like` for FastTree,
  `ufboot` for IQ-TREE). The tree is ladderized.

Opens directly in [Archaeopteryx](https://sites.google.com/view/aptxjs)
([source](https://github.com/cmzmasek/forester)), which makes use of the rich
annotation; also readable in any phyloXML-aware tool (Dendroscope, etc.).

#### `{prefix}_tree.pdf`, `{prefix}_tree.png` — the figure you'll usually paste

A ready-to-use graphical rendering of the tree above, so you don't have to open a
viewer to glance at the topology. PDF is vector (clean at any zoom); PNG is a
150-dpi raster. Both are a ladderized rectangular phylogram (matplotlib +
Bio.Phylo) showing: leaf labels coloured by taxonomy (matching the phyloXML)
with a genus/subfamily legend; internal-node labels for genus→family ranks;
branch-support labels for nodes with support ≥ 50; a title + one-line provenance
caption. Produced for *every* tree type. **On by default** (`phylo.pdf`); see
"Graphical tree figures" above for `--no-pdf` and the matplotlib requirement.

#### `{prefix}_tree_id_map.tsv`

Two-column `short_id` ↔ `accession`. Use it to trace a leaf in the MSA (or the
opt-in `_tree.nwk`) back to a real sequence id.

#### `{prefix}_iqtree_summary.txt` — only when IQ-TREE ran

The IQ-TREE ModelFinder report — which model was selected and why, plus the
log-likelihood and bootstrap settings.

#### `{prefix}_iqtree_model.txt` — only when IQ-TREE ran

The same picks distilled to one `<label>: <model>` line per partition. A
non-partitioned run emits a single `GENOME: <model>` line; partitioned runs emit
one line per charset. The picks also land in the phyloXML `<description>`'s
`model=` field and the `_summary.md` Methods section.

> **Fail-soft:** if MAFFT, IQ-TREE, or FastTree are missing, or fewer than 3
> representatives survived, the phylogeny step is skipped with a stderr message
> and the rest of the run's outputs are still written. IQ-TREE refuses UFBoot
> with fewer than 4 sequences; the wrapper drops bootstrap automatically but
> still produces the tree.

### Per-protein trees — `--per-protein-phylo`

`--phylo` builds **one** tree from the whole representatives. `--per-protein-phylo`
instead builds **one tree per declared HMM marker** — one per `hmms:` spec. A
spec's `hmms:` list holds **alternative architectures** (OR; see
[Token notation](#token-notation)), so a marker like coronavirus Spike
(`["CoV_S1--CoV_S2", "bCoV_S1_N--bCoV_S1_RBD--CoV_S2"]`) builds one tree spanning
every rep whose Spike matches *either* architecture. repseq picks the satisfying
CDS on each carrying rep and runs the *same* MAFFT → IQ-TREE/FastTree → root →
LCA pipeline. Alignments use MAFFT `--auto` by default; for a high-accuracy
publication run set `phylo.per_protein.mafft.extra_args: ["--maxiterate", "1000",
"--localpair"]` (L-INS-i — slower, opt-in; a non-empty list is passed without
`--auto`).

Why you'd want it: comparing single-marker trees side by side reveals
**topological incongruence** — an L-segment polymerase tree disagreeing with an
M-segment glycoprotein tree is the classic signature of **reassortment**, which
the concatenated `--phylo` tree hides.

Outputs land in `{prefix}_per_protein/`, one set per built family. `<family>` is
the spec's `name:`, else its single token, else first token + `_altN`;
segment-prefixed in segmented mode:

```
{prefix}_per_protein/
  Spike_msa.fasta
  Spike_tree.xml
  Spike_tree.pdf
  Spike_tree.png
  Spike_tree_id_map.tsv
  S_Bunya_nucleocap_msa.fasta
  …
  {prefix}_incongruence.tsv
```

(`*_tree.nwk` only when `phylo.newick` / `--newick` is on; `*_tree.pdf` / `.png`
by default, `--no-pdf` to skip. The phyloXML, MSA, id-map, and incongruence table
are always written.)

**Extra-protein trees:** when you've also declared
[`extra_protein:`](#accessory-proteins-extra_protein), a **separate**
`{prefix}_extra_protein/` directory is emitted with the same per-tree shape,
built by the same engine. The split is intentional — accessory proteins are
sparse, so they're **excluded from `_incongruence.tsv`** to keep the reassortment
signal in the required-marker pairs from being drowned out by `NA` rows. For RF
distances involving an accessory tree, run the math off its `.nwk` yourself.

Each tree file has the same rich annotation as its `--phylo` counterpart, with
two deliberate differences. First, a leaf shows **only the CDS that fed that
tree** as `<sequence type="protein">` (the `CoV_nucleocap` tree shows just the
nucleocapsid). The `<sequence type="dna">` element and the `repseq:nuc_acc` /
`repseq:protein_acc` / `repseq:protein_names` summary properties still describe
the leaf's full gene content. Second, that protein carries its **HMM domain
architecture** as a phyloXML `<domain_architecture>` — one `<domain>` per hit,
with protein coordinates and the hit's E-value as `confidence`:

```xml
<sequence type="protein">
  <accession source="ncbi">QHD43416</accession>
  <name>spike glycoprotein</name>
  <domain_architecture length="1273">
    <domain from="13" to="305" confidence="2e-40">CoV_S1</domain>
    <domain from="334" to="1273" confidence="1e-120">CoV_S2</domain>
  </domain_architecture>
</sequence>
```

[Archaeopteryx](https://sites.google.com/view/aptxjs) draws these as domain boxes
and lets you filter them with its interactive E-value slider — so *all* hits are
emitted (not just those that cleared repseq's cutoffs), giving the slider its full
range. Turn the block off with `phylo.per_protein.domain_architecture: false`.
The flag runs alone or alongside `--phylo`.

#### `{prefix}_incongruence.tsv` — incongruence as a number

repseq scores the **pairwise unrooted Robinson-Foulds (RF) distance** between
every pair of marker trees (and, when `--phylo` also ran, the whole-genome tree
as a `GENOME` row):

| tree_a | tree_b | rf | norm_rf | n_common_taxa |
|--------|--------|----|---------|---------------|
| Spike  | N      | 4  | 0.3333  | 18            |
| Spike  | GENOME | 0  | 0.0000  | 21            |
| N      | GENOME | 4  | 0.2857  | 18            |

- **`rf`** — bipartitions that differ between the two trees, scored on their
  **shared** taxa only. `0` = identical unrooted topology. Rooting is ignored.
- **`norm_rf`** — `rf / 2·(n−3)`, comparable across pairs with different overlap.
  `NA` when fewer than 4 taxa are shared.
- **`n_common_taxa`** — how many representatives the pair has in common.

A `GENOME` column tells you which marker departs most from the consensus history.
Turn the table off with `phylo.per_protein.incongruence: false`.

> **Requirements / fail-soft:** needs the HMM tier (`hmm.enabled` + configured
> `hmms:`) to have run, and `mafft` + a tree builder on PATH. A family carried by
> fewer than `phylo.per_protein.min_taxa` (default 3) reps is skipped with a log
> note; if no family qualifies the whole step is skipped with one stderr line.

### Per-segment NT trees — `--per-segment-phylo` (segmented only)

A complement to `--per-protein-phylo`. Where that builds one **protein** tree per
marker (one CDS per segment), `--per-segment-phylo` builds one **nucleotide** tree
per declared segment — every representative isolate contributes its raw
per-segment NT (from `concat_segments`). Outputs land in `{prefix}_per_segment/`
with the same per-tree shape (`<segment>_msa.fasta`, `<segment>_tree.xml`,
`<segment>_tree_id_map.tsv` per segment, plus `<segment>_tree.pdf` / `.png` by
default; `.nwk` only when `phylo.newick`).

Why on top of the per-marker trees: reassortment can show up at the segment level
without appearing in any single marker tree (a marker is one CDS within a
segment). The per-segment NT trees catch reassortment the per-marker AA trees can
miss, especially when a marker is conserved across the parent species. Same
MAFFT / tree-builder / rooting / LCA / colouring / `min_taxa` / shared palette as
the other trees. Runs alone or alongside `--phylo` and `--per-protein-phylo`;
raises a stderr note on a non-segmented run.

#### Outgroup rooting

Set `phylo.rooting.method: outgroup` and name the outgroup by accession or taxon:

```yaml
phylo:
  rooting:
    method: outgroup
    outgroup: "AB123456"
    # …or a list (rooted at the MRCA of all matching leaves):
    # outgroup: ["AB123456", "CD789012"]
    # …or a whole clade by taxonomy rank:
    # outgroup_rank: {family: "Hantaviridae"}
```

Accession matching is case-insensitive against each rep's `accession`, `id`, or
`isolate_id` (so a CONCAT isolate id works too). A multi-leaf spec roots at the
MRCA of all matching leaves. `outgroup` and `outgroup_rank` are unioned. When the
spec matches no representative the rooter falls through to midpoint with a stderr
note, so a typo never voids the tree.

The default is `phylo.rooting.method: auto` (taxonomy → MAD → midpoint chain); the
other accepted values are `taxonomy`, `mad`, `midpoint`, and `none`.

### Pre-cluster overview tree — `--pre-cluster-tree`

A diagnostic tree over **every post-QC sequence**, built BEFORE clustering would
have collapsed redundancy — so you can see where the elected representatives land
in the broader diversity of the input pool. The phyloXML `<name>` of each
representative leaf is prefixed with `[repr] `, and **every** leaf additionally
carries a `repseq:is_representative` boolean `<property>` (`true`/`false`) — a
machine-readable flag you can filter, select, or colour by in a viewer such as
Archaeopteryx. This is the only repseq tree that emits it (every other tree
contains only representatives).

The leaf set is **capped at `phylo.pre_cluster_tree.max_leaves` (default 5000)** —
above the cap, all representatives are kept and the background is randomly
subsampled to fill it (`0` = no cap; see the scale note). The pipeline is
otherwise **hard-coded for speed**: **MAFFT `--retree 1`** (or **`--parttree`**
above `phylo.pre_cluster_tree.parttree_threshold` when uncapped), **FastTree**,
**midpoint rooting**, no LCA / trimAl / bootstrap. The alphabet mirrors
`clustering.alphabet_for_clustering` (AA when the run clustered on protein, else
NT); in segmented mode the leaves are CONCAT isolates.

Output files:

- `{prefix}_pre_cluster_tree.nwk` — short-id leaves. **Opt-in** (`phylo.newick`).
- `{prefix}_pre_cluster_tree.xml` — phyloXML with the full taxonomy `<property>`
  enrichment + per-leaf colour (same palette as the other trees), the `[repr] `
  prefix on representative leaves, and a `repseq:is_representative` boolean
  `<property>` on every leaf.
- `{prefix}_pre_cluster_tree.pdf` / `.png` — figure, rendered by default
  (`phylo.pdf`; `--no-pdf` to skip).
- `{prefix}_pre_cluster_tree_id_map.tsv` — three columns (`short_id`,
  `accession`, `is_rep`).

The MAFFT MSA is written to a temp dir and deleted after FastTree consumes it —
only the tree + id_map land in the output dir (run `--phylo` if you want the
MSA). Trigger via `--pre-cluster-tree` or `phylo.pre_cluster_tree.enabled: true`.
Soft-fails like the other phylo steps.

> **Scale note:** on a very large pool (tens of thousands of isolates — e.g. a
> global influenza run) FastTree memory scales ~linearly with leaf count (~1.7 GB
> at 1,000 protein leaves → ~180 GB at 100,000). The **`max_leaves` cap (default
> 5000)** keeps FastTree near ~7 GB and the figure legible — a few minutes to
> build. Only if you uncap (**`max_leaves: 0`**) does the MAFFT O(N²) matrix
> matter: at/above `phylo.pre_cluster_tree.parttree_threshold` (default 10000)
> MAFFT switches to **`--retree 2 --parttree`** (PartTree guide, scales to 10⁵⁺;
> `0` = always, very large = never).

### One conservation number per MSA — `{prefix}_msa_conservation.tsv`

repseq reduces **every alignment a run produces** to a single conservation score
and collects them into one table. Whenever any phylogeny step writes an MSA — the
whole-genome tree, the partitioned per-family alignments and supermatrix, the
per-protein marker trees, and any extra-protein, polyprotein-peptide, or
per-segment trees — a post-hoc sweep scores each `*_msa*.fasta` under the output
dir, one row per alignment. Runs by default; set `phylo.conservation.enabled:
false` to turn it off. Being a decoupled sweep, a scoring problem never affects a
tree, and if no tree was built the file isn't written.

**The metric.** Mean per-column **Jensen-Shannon divergence to a residue
background** (Capra & Singh 2007), with two standard corrections:

- **Henikoff & Henikoff (1994) position-based sequence weighting**, so a clade of
  near-identical sequences plus one outlier doesn't read as artificially
  conserved.
- a **gap penalty**: each column's JSD is multiplied by
  `(1 − weighted_gap_fraction)`, so a mostly-gap column contributes little.

Each column score is bounded `[0, 1]` (JSD, log base 2, symmetric mixture
`m = (p+q)/2`). Interpretation: `~0` looks like the background (unrelated
sequences); `~0.6–0.8` is typical for a real protein family; a perfectly
conserved column tops out at ~0.85–0.95 for protein and ~0.55 for nucleotide.
JSD-to-background **never reaches exactly 1** — it measures *divergence from a
background distribution*, not raw uniformity, so an invariant column's value
depends on how rare the conserved residue is (rarer ⇒ closer to 1). This is a
deliberate property of the Capra & Singh formulation and is why published JSD
conservation values don't saturate. Scores are therefore comparable only *within*
an alphabet — the `alphabet` column flags which is which.

**Columns:** `msa` (path relative to the output dir), `role` (`genome` /
`partition_family` / `marker` / `extra_protein` / `peptide` / `segment_nt`),
`label` (family / segment / peptide name), `alphabet` (`protein` / `nucleotide`,
auto-detected), `trimmed` (`TRUE` only when trimAl actually trimmed *this*
alignment), `n_seqs`, `n_sites`, `mean_conservation` (the headline all-column
mean), and `mean_conservation_core` (mean over well-occupied — ≥ 50% non-gap —
columns, ignoring the gap penalty).

**References.**

- **Capra JA & Singh M (2007).** Predicting functionally important residues from
  sequence conservation. *Bioinformatics* 23(15): 1875–1882.
  doi:[10.1093/bioinformatics/btm270](https://doi.org/10.1093/bioinformatics/btm270)
- **Henikoff S & Henikoff JG (1994).** Position-based sequence weights.
  *J Mol Biol* 243(4): 574–578.
  doi:[10.1016/0022-2836(94)90032-9](https://doi.org/10.1016/0022-2836(94)90032-9)
- **Valdar WSJ (2002).** Scoring residue conservation. *Proteins* 48(2): 227–241.
  doi:[10.1002/prot.10146](https://doi.org/10.1002/prot.10146) — survey of
  conservation schemes and the substitution-matrix-weighted alternative to JSD.

### Per-taxon monophyly — `{prefix}_monophyly.tsv`

For **every tree a run builds** (whole-genome, per-protein/extra/segment/
polyprotein, pre-cluster), repseq classifies how each taxon at each rank sits on
that rooted tree. Like the conservation and incongruence reports it's a decoupled,
always-on sweep over the annotated `*_tree.xml` files; it needs no flag and writes
nothing if no tree was built.

**Status (per taxon, per rank, per tree):**

- **monophyletic** — the taxon's members form exactly one clade (their MRCA
  subtends only members).
- **paraphyletic** — non-monophyletic, and the foreign leaves inside the members'
  MRCA span form a **single** excluded clade.
- **polyphyletic** — non-monophyletic, and the foreign leaves form **two or more**
  separate clades.

Each rank is judged **only among leaves classified at that rank** — a leaf with no
genus is neither member nor intruder when genus monophyly is assessed, so an
annotation gap never masquerades as non-monophyly. Taxa with fewer than two
classified leaves are omitted.

**Ranks assessed: subgenus → class by default; species is opt-in.** Viral
species labels are annotation-noisy (the same reason the taxonomic reports skip
species), so species rows are off by default. Set
`phylo.monophyly.include_species: true` to add them — species is the rank at
which **intra-genus reassortment** and most **misannotation** live: a species
monophyletic on the whole-genome tree but broken on a single marker tree is the
localised reassortant call that the coarser ranks cannot see.

**Why it's useful.** A taxon monophyletic on the whole-genome tree but **not** on
a marker tree is the per-marker **reassortment / recombination** signal — the
taxon-resolved companion to `_incongruence.tsv`. It also operationalises the
non-monophyly the tree's internal labels now show.

**The para/poly tag is a documented heuristic.** Monophyletic-vs-not is an exact
MRCA test, but the para/poly split is definition-dependent; repseq keys it on
`intruder_clusters` (1 → paraphyletic, ≥2 → polyphyletic), the standard
topology-only convention, and reports the unambiguous counts alongside the tag.

**Support-aware by default.** A single weakly-supported branch can make a taxon
look non-monophyletic when the tree doesn't actually resolve it that way. So
branches with support below `phylo.monophyly.min_support` (**default 70**) are
collapsed into polytomies first, and a taxon is called monophyletic when no
*well-supported* node contradicts it. Only **confident** non-monophyly is flagged.
Set `min_support: 0` to trust every branch (topology-only). The threshold used is
recorded in the report.

**Columns:** `tree` (path relative to the output dir), `rank`, `taxon`, `n_leaves`
(classified members), `status`, `n_clusters` (separate in-group blocks),
`n_intruders` (foreign classified leaves in the MRCA span), `intruder_clusters`
(separate foreign clades — the para/poly key), `intruder_taxa` (`;`-joined), and
`min_support` (the threshold applied).

### Segment-status matrix — `{prefix}_segment_status_matrix.tsv`

`_monophyly.tsv` has one row per (taxon, rank, **tree**); the biological
question is **cross-tree**. This file pivots it so each (rank, taxon) is one
row, and answers *which taxon breaks on which one tree*. It's a pure pivot of
the monophyly report (no trees are re-read), so it inherits the support-aware
collapse and — with `phylo.monophyly.include_species` on — the species rows
where intra-genus reassortment lives.

The payoff column is **`single_marker_break`**: the label of the *one* marker
tree a taxon breaks on, populated **only** when it is monophyletic on the
whole-genome tree and broken on exactly that one tree — the clean
single-segment-discordance (reassortment) signal. Filtering this column to
non-empty gives the candidate list directly; rows with a value sort to the top
of each rank. Other columns: `n_leaves`, `n_trees` (trees that assessed the
taxon), `n_nonmono`, `genome_status` (blank if no whole-genome tree was built),
and `nonmono_trees` (`;`-joined labels of every tree it breaks on). Where
`_incongruence.tsv` says *how much* two trees disagree as one global number,
this says *which taxon, on which tree*.

### Plain-English flags — `{prefix}_flags.txt`

The interesting anomalies a run can find are spread across several tables
(`_monophyly.tsv`, `_incongruence.tsv`, and — when the opt-in review ran —
`_taxonomy_review.tsv`). `{prefix}_flags.txt` distils them into one skimmable,
plain-English list so you don't have to join TSVs by hand:

- **Reassortment / recombination signals** — taxa monophyletic on the
  whole-genome tree but broken on a marker tree, and marker-tree pairs with a high
  normalised Robinson-Foulds distance (≥ 0.2).
- **Taxonomy ↔ tree conflicts** — taxa that are para- or polyphyletic on the
  whole-genome tree.
- **Per-isolate taxonomy conflicts** — isolates whose label the tree
  neighbourhood disagrees with (from the taxonomy review).

It's a pure synthesis (no new computation), written whenever any source table
exists; a clean run gets a `No flags raised ✓` file so you can tell "nothing
flagged" from "report not produced". The flags are heuristics — the source tables
remain authoritative.

### One-page HTML report — `{prefix}_report.html`

A single self-contained HTML file you can open in any browser or e-mail to a
collaborator. It bundles the run's **analysis flags**, a **gallery of every tree
figure** (the `*_tree.png` images embedded directly), and an **index of all
output files** (with sizes and relative links). Pure-Python, no extra dependency;
written last so its file index captures the whole run.

---

## Cleaning (QC) — what gets dropped and how to control it

All cleaning settings live under `qc:`. Defaults are sensible; loosen them if
you're losing sequences you want. Every drop is logged in
`{prefix}_qc_removed.tsv` with a precise reason.

```yaml
qc:
  remove_duplicates: true        # drop byte-identical sequences (keeps the curated copy)

  length_filter:
    mode: median_percent         # judge length relative to the dataset's median…
    min_percent: 50              #   …drop anything shorter than 50% of the median
    # max_percent: 200           #   …optionally also drop anything over 200%
    # ── or ──
    # mode: min_max              # judge length against fixed numbers instead
    # min_length: 1000
    # max_length: 20000

  ambiguous_threshold: 0.05      # drop sequences that are >5% N / X / other ambiguous letters

  annotation_filter:
    enabled: true
    keywords: ["MAG:", synthetic, partial, hypothetical, fragment, uncultured, ...]
    # any sequence whose description contains one of these words is dropped
```

A few things worth knowing:

- **`median_percent` compares every sequence to the median length of the whole
  file.** Perfect for a single gene, but **wrong for a mixed file** (several genes,
  or a whole genome plus its individual genes) — the median is meaningless and
  you'll drop things unfairly. For mixed files, use `min_max` with explicit
  numbers.
- **In segmented mode the whole-file length filter is skipped automatically** — a
  file of influenza segments mixes 2,300-nt and 890-nt sequences, so a single
  median can't work. Use per-segment bounds (`segmented.viruses.<v>.segment_lengths`).
  The QC summary breaks per-segment drops down (`L too short : N`, …), so you can
  see which cutoff was the bottleneck.
- **Segmented mode adds a taxonomy-consistency check.** Any isolate whose segments
  resolve to *different* taxa at the configured rank (default `species`) is
  dropped — usually a reassortant, a contaminated record, or two viruses sharing
  an `/isolate` qualifier. Missing labels are ignored (only *populated*
  disagreement counts). Configure under `segmented.taxonomy_consistency.{enabled,
  rank}`; `enabled: false` skips it.
- **Segmented mode also flags strain-as-isolate collisions.** When the GenBank
  pre-pass falls back to `/strain` for records without `/isolate`, two accessions
  can share a grouping key. The detector warns on colliding `(isolate_id, segment)`
  pairs; set `segmented.strain_collision_action: drop` to remove them. See
  [Strain-as-isolate fallback](#strain-as-isolate-fallback-and-collision-detection).

---

## Sequences of special importance — `overrides`

Sometimes a handful of records *must* be handled specially: a type strain, a
vaccine reference, a genome reviewers will look for by name. The `overrides` block
names them once (an id list and/or a file) and applies **two independent
capabilities** over that one list:

- **`protect_qc`** — **force-keep**: the named sequences bypass the QC removal
  stages they would otherwise fail (a whitelist), without loosening thresholds for
  everything else.
- **`force_select`** — **force-select**: the named sequences are guaranteed to
  appear among the representatives, regardless of how clustering would have
  collapsed them.

A third capability, **`exclude`**, is the mirror image — a blocklist that
*removes* named sequences from the input before anything runs (see
[Excluding bad sequences](#excluding-bad-sequences) below). It carries its own
separate id list, not the shared one above.

```yaml
overrides:
  ids: ["NC_045512.2", "MN908947"]   # accessions (or isolate_ids, segmented)
  ids_file: vip.txt                  # OR a file, one id per line; unioned with `ids`
  protect_qc: true                   # survive QC
  protect_stages: all                # or a subset (see below)
  force_select: true                 # appear in the representatives
```

Or pass lists on the command line — `--protect-ids FILE` enables force-keep,
`--pin-ids FILE` enables force-select; both union into the one id list:

```bash
repseq global -c my_config.yaml -i seqs.fasta -T 0.9 \
    --protect-ids vip.txt --pin-ids vip.txt
```

Each CLI flag reads one id per line (`#` comments and blanks ignored), unions with
any `overrides.ids` / `overrides.ids_file`, and turns its flag on for that run.

- **Matching** is by `accession` (non-segmented) or `isolate_id` (segmented),
  **case- and version-insensitive** — `NC_045512` matches `NC_045512.2`. In
  segmented mode you usually list `isolate_id`s; note the *early* QC stages
  (`ambiguous`, `annotation`) run before `isolate_id` is populated, so those only
  match on accession.
- **`protect_stages`** is `all` or a list naming any of: `duplicates`, `length`,
  `ambiguous`, `annotation`, `protein_count`, `taxonomy_consistency`,
  `protein_quality`, `hmm`. Per-stage on purpose: you can keep a slightly noisy
  reference through the `hmm` identity gate while still dropping it if its
  translation is genuine garbage (`protein_quality`).
- **Segmented completeness is deliberately *not* protectable.** Protection can't
  synthesize a missing segment, so a protected isolate missing a segment still
  drops (it can't be concatenated or treed) — with a distinct message rather than
  a silent fake-keep.

### How force-select decides (the hybrid policy)

When `force_select` is on, after clustering each pinned sequence is guaranteed into
the representatives by the rule that disturbs the rest of the selection least:

- **Wins its cluster** — a pinned member beats the configured
  refseq/reviewed/longest priority and becomes that cluster's representative (the
  displaced rep drops back into the members).
- **Splits on collision** — if two or more pinned members land in the *same*
  cluster, the best one wins the slot and the others split into their own
  singleton clusters, so all survive.
- **Added as a singleton** — a pin clustering deselected entirely (the
  diversity-only `global -n` path) is appended as its own singleton.
- **Unavailable** — a pinned id no surviving sequence matched (dropped in QC, or
  never in the input) is reported as `unavailable` with a warning.

> **Force-keep needs to pair with force-select for a hard guarantee.**
> `force_select` can only promote a sequence that *survived QC*. If a pinned record
> would be dropped by QC, enable **both** (`protect_qc` + `force_select`, or
> `--protect-ids` + `--pin-ids` on the same list).

- **Nothing is silent.** Force-kept records (and why each would have been removed)
  go to `{prefix}_overrides.tsv`; force-select actions go to
  `{prefix}_force_selected.tsv`; both are summarised in `{prefix}_summary.md`.
- **Count reconciliation.** Force-select runs *after* clustering, so the
  per-stratum counts in `{prefix}_group_counts.tsv` reflect clustering **before**
  any pins were added — they can sum to fewer than the final representative count.
  The pinned additions show up in the representative files,
  `{prefix}_force_selected.tsv`, and the final counts in `{prefix}_summary.md`.

### Excluding bad sequences

The inverse of force-keep: a **blocklist** for records you know are bad
(chimeric, mislabelled, contaminated). Listed sequences are dropped the moment
the FASTA is read — **before** metadata resolution, QC, and clustering — exactly
as if you had deleted them from the input file. They cost no network lookups and
appear in no QC counter.

```yaml
overrides:
  exclude:
    enabled: true
    ids: ["KJ642623.1"]    # accessions/ids; OR a file below
    ids_file: bad.txt      # one id per line; unioned with `ids`
```

Or on the command line — `--exclude-ids FILE` reads one id per line (`#` comments
and blanks ignored), unions with any `overrides.exclude.ids` /
`overrides.exclude.ids_file`, and turns the blocklist on for that run:

```bash
repseq global -c my_config.yaml -i seqs.fasta -T 0.9 --exclude-ids bad.txt
```

- **Matching** is by `accession` or `id` only, **case- and version-insensitive**
  (`NC_222` matches `NC_222.3`). Unlike force-keep/force-select it does *not* match
  `isolate_id` — exclusion runs before the GenBank lookup that populates it, so the
  rule is deliberately header-id-only.
- **Its own id list.** The `exclude` block has a dedicated `ids` / `ids_file`,
  separate from the shared keep/pin list. Listing the same id under both `exclude`
  and `protect_qc`/`force_select` is a **hard error** — an id cannot be both
  removed and guaranteed.
- **Audited, never silent.** Every dropped record goes to `{prefix}_excluded.tsv`
  (`id`, `action`, `detail`); an exclude id that matched nothing is reported there
  as `unavailable` plus a stderr warning (a typo-catcher). The count is summarised
  in `{prefix}_summary.md`.

---

## Segmented viruses (influenza, etc.)

Segmented viruses store their genome in several separate pieces. NCBI gives you
one FASTA record per segment, so a single isolate is spread across multiple
records. `repseq` stitches them back together:

1. **Group records by isolate** — by default from the `/isolate` (or `/strain`)
   qualifier of each record's GenBank source feature. The strain-name regex below
   is the fallback for sequences without an NCBI accession (e.g. UniProt) or whose
   record lacks the qualifier.
2. **Identify each record's segment** — by default from the GenBank `/segment`
   qualifier, or by name / number / a synonym you define (`hemagglutinin` → `HA`).
3. **Keep only complete isolates** — one missing any expected segment is dropped.
4. **(Optional) length-check each segment.**
5. **Concatenate** the segments of each complete isolate into one sequence, so the
   normal grouping/selection runs on whole isolates.

The GenBank lookup reuses the same cache as protein-annotation QC, so with both
on, repseq fetches each record once. To turn the lookup off entirely (relying on
the regex for *every* record), set `segmented.use_genbank_metadata: false`.

Configure under `segmented:` and turn it on with `--segmented` (or `enabled:
true`):

```yaml
segmented:
  enabled: false
  virus: influenza_a              # which entry below to use
  viruses:
    influenza_a:
      expected_segments: 8
      segments: [PB2, PB1, PA, HA, NP, NA, M, NS]   # canonical order
      isolate_regex: "(?P<isolate>[AB]/[^/(\\s]+/[^/(\\s]+/[^/(\\s]+/\\d{4})"
      segment_aliases:            # optional: words in headers that mean a segment
        HA: [hemagglutinin]
        NA: [neuraminidase]
        NP: [nucleoprotein, "nucleocapsid protein"]
      segment_lengths:            # optional: drop an isolate if a segment is out of range
        HA: {min: 1600, max: 1800}
        NS: {min: 800,  max: 1000}
```

`config/examples/influenza_a.yaml` is a complete, commented example you can copy.

> The `isolate_regex` is a *fallback* — by default
> (`use_genbank_metadata: true`) repseq prefers the GenBank `/isolate` qualifier.
> The regex still needs to match for sequences without an NCBI accession or whose
> record has no isolate qualifier, and it must capture the isolate either as a
> group named `isolate` or as the first parenthesised group.

### Strain-as-isolate fallback and collision detection

GenBank source features have two distinct qualifiers: `/isolate` (a specific
collection event, usually unique to one biological sample) and `/strain` (a named
variant, often shared across many samples). When a record has `/strain="L99"` but
no `/isolate=`, repseq uses the strain value as the isolate grouping key —
otherwise strain-only records (common in older bunyaviridae / hantaviridae
submissions) would drop out of segmented mode entirely.

The provenance of every `isolate_id` is written to the `isolate_id_source` column
in `_representative_isolates.tsv` and `_isolate_proteins.tsv`: `isolate` (from
`/isolate`, preferred), `strain` (from `/strain` as a fallback), or `regex` (the
header-regex fallback).

**The risk of `strain`-derived ids is over-merging.** Two submitters can deposit
independent samples under the same strain name; repseq would group them as one
isolate and dedup-drop the duplicate segment. After the GenBank pre-pass, repseq
runs a **strain-collision detector** that flags any `(isolate_id, segment)` pair
where two or more distinct accessions share a strain-derived id:

```
  Strain-collision check: found 2 (isolate, segment) pair(s) where /strain is shared across distinct accessions.
    isolate 'L99' segment 'S': 3 accessions (ACC1, ACC2, ACC3)
    isolate 'M77' segment 'L': 2 accessions (ACC4, ACC5)
    Action: warn (no records dropped). Set segmented.strain_collision_action: drop to remove them.
```

The default action is `warn` (the pipeline keeps the longest per-segment accession
and dedup-drops the rest downstream). Switch to `drop` to remove every accession
involved in any collision before the completeness filter runs:

```yaml
segmented:
  strain_collision_action: drop   # warn (default) | drop
```

Dropped records appear in `_qc_removed.tsv` with reason `strain_collision:<segment>`.

---

## Optional: protein-annotation QC

This step asks NCBI how many protein-coding genes (CDS features) each record has,
and drops records with too few — or, for segmented viruses, the wrong number per
segment. Off by default:

```yaml
qc:
  protein_annotation:
    enabled: true
    min_proteins: 1     # drop any record with fewer annotated proteins than this

segmented:
  viruses:
    influenza_a:
      expected_proteins_per_segment:
        HA: 1
        M:  2           # M1 + M2
        PB1: [1, 2]     # PB1 alone, or PB1 + PB1-F2 — a list means "any of these"
        NS: [1, 2]      # NS1 alone, or NS1 + NEP
```

Records are fetched in batches and cached locally, so a second run needs no
network. Skipped automatically under `--no-resolve`.

---

## Phylogeny-based taxonomy review (`phylo.taxonomy_review`)

Many viral records carry **missing** taxonomy (a blank genus, "unclassified") or a
**wrong/stale** label that pre-dates an ICTV reclassification. Once `--phylo` has
built a well-supported, LCA-labelled tree, the tree itself is evidence: a
representative nested deep inside a strongly-supported, taxonomically *pure* clade
is very likely the same taxon as its neighbours.

Enable it (off by default — it's an inference step):

```yaml
phylo:
  taxonomy_review:
    enabled: true
    ranks: [family, genus, subgenus]   # coarse→fine; species omitted on purpose
    min_support: 90        # enclosing-clade branch support (0-100)
    min_purity: 0.9        # fraction of labelled neighbours agreeing
    min_agreeing: 3        # min labelled neighbours backing the call
    require_refseq_anchor: true
    write_corrected: true
```

For each representative leaf and rank, it walks up to the **smallest enclosing
clade** that clears the gates (support, purity, neighbour count, RefSeq/reviewed
anchor) and then **imputes** a *blank* rank from the clade's majority value, or
**flags** a *populated* rank that disagrees (a suggestion only, **never**
auto-changed). Ranks are evaluated coarse→fine and imputations are kept
**hierarchy-consistent**. Confidence is tiered (high/medium) from support × purity
× neighbour count × anchor.

**Outputs.** `{prefix}_taxonomy_review.tsv` lists every verdict (impute + conflict)
with evidence columns — the authoritative per-cell imputation ledger. With
`write_corrected: true`, the **high-confidence imputations** are filled into
`{prefix}_representative_*_corrected.tsv` and the corresponding corrected protein
FASTA (clean values; originals kept untouched). Conflicts are *never* written into
the corrected copies.

**Scope & cautions.** Representatives only. The tree topology, colouring, and LCA
labels are **not** modified — no feedback loop. The support + purity + anchor
gates guard against rooting bias and recombination artifacts. Species is
intentionally excluded (viral species monophyly is too often violated). Soft-fails
(a stderr note, no files) if `--phylo` didn't build a tree.

---

## Optional: clustering plot

Pass `--plot` (and install with `'.[viz]'`) for a two-panel scatter plot,
`{prefix}_clustering.png`, to eyeball whether the clustering looks sensible: left
panel coloured by genus, right panel coloured by cluster (bigger dots for bigger
clusters, faint lines from each sequence to its representative). For big datasets
it's drawn from a subsample (representatives always included); skipped for
`global -n` runs, which produce no clusters.

---

## The config file

`repseq init-config` writes a **complete, fully-commented** config: every
option present at its default with an inline explanation, so the generated file
*is* the reference — read the comments in place and edit what you need. (The
same canonical content is bundled at `repseq/data/default_config.yaml`.) The
most-changed settings:

```yaml
qc:
  remove_duplicates: true
  ambiguous_threshold: 0.05
  # genome_length_filter is non-segmented-only and opt-in (default off):
  # genome_length_filter: {enabled: true, min: 1000, max: 50000}

taxonomy:
  ncbi_email: you@institute.org   # NCBI asks for this; without it you'll be rate-limited
  ncbi_api_key: null              # optional — get one from NCBI for faster lookups

clustering:
  backend: mmseqs2                # or "cdhit"
  alphabet_for_clustering: protein  # protein (default) or nucleotide
  mmseqs2_mode: easy-linclust     # fast; use easy-cluster for tighter, slower clustering
  coverage: 0.8
  # cd-hit options (only used when backend == cdhit) live under
  # `clustering.cdhit:` — see that block in your generated config (it's
  # there with comments) for the full set.

representative:
  # Ordered preference for picking each cluster's "best" sequence: the first
  # criterion wins, later ones break ties, sequence length is the final
  # tiebreaker. The order is honoured — [reviewed_uniprot, refseq, longest]
  # prefers a Swiss-Prot-reviewed entry over a RefSeq one. Drop a criterion
  # to deactivate that preference. (`longest` is required.)
  priority: [refseq, reviewed_uniprot, longest]
```

You can also set your NCBI email/key via the environment variables
`REPSEQ_NCBI_EMAIL` and `REPSEQ_NCBI_API_KEY` instead of putting them in the file.

### Unknown keys are rejected

repseq validates your config before doing any work and **aborts with a
`[config error]`** if it finds a key it doesn't recognise — a typo
(`ambiguous_threshhold`), a renamed/removed key, or a setting placed in the wrong
section (e.g. a `phylo.per_protein` knob written onto a `clustering.polyprotein`
spec). The message names the offending key, its full dotted path, and — when
there's a close match — the key you probably meant. This is deliberate: a
silently-ignored setting is worse than a loud stop.

To keep your own annotation/reminder key in the file, prefix it with an underscore
(`_my_note: ...`) and repseq leaves it alone. To check a config without launching
a run, use `repseq <mode> --dry-run -c your.yaml`.

### About `clustering.alphabet_for_clustering`

**This setting *only* changes what the clustering backend sees.** It does **not**
disable GenBank CDS download, protein-count QC, or the per-segment
`virus.expected_proteins_per_segment` check — those run on every isolate
regardless. Pick it purely on what kind of identity threshold makes sense for your
data. Only `protein` and `nucleotide` are accepted.

| Value | Clustering input | When to pick it |
|---|---|---|
| `protein` *(default)* | Non-segmented: the per-sequence marker protein (longest CDS by default, overridable via `cluster_protein` aliases). Segmented: per-isolate concatenation of each segment's marker. Triggers a one-shot GenBank CDS fetch if not cached. | Diverged viral families: synonymous substitutions inflate NT divergence by 30–40% with no biological signal, and protein homology stays reliable down to ~25–30% identity vs ~50–60% for nucleotide. |
| `nucleotide` | Non-segmented: the input FASTA as-is. Segmented: concatenation of all segments in `segments` order. | Tight species-level reference sets where genome-identity targets are what you want; non-coding input; or fully offline runs (`--no-resolve`), which require this value. |

#### Clustering on more than one marker (`clustering.concatenate_markers`)

By default a non-segmented protein-alphabet run clusters on a **single** marker —
the CDS from the first `cluster_protein` spec a record satisfies (e.g. Spike alone
for coronaviruses). Set:

```yaml
clustering:
  concatenate_markers: true
```

to instead select the marker CDS from **every** `cluster_protein` spec and cluster
on their amino-acid **concatenation** in declared spec order (e.g.
`Spike+Nucleocapsid`). This mirrors the segmented per-segment concat. A record
missing any required marker is dropped (when every spec declares `hmms:`, the
AND-gate already removes those upstream). `extra_protein` specs are never
concatenated.

This flag changes **only the clustering input** (and the
`{prefix}_representatives_protein.fasta` it writes). The whole-genome tree is
*not* affected. For a multi-gene whole-genome tree, leave the trees to the
**partitioned supermatrix** path (`phylo.tool: auto` or `iqtree`): each marker is
aligned separately and concatenated column-wise with a per-gene model. Gluing the
proteins into one string and aligning the chimera (what you'd get under FastTree)
is biologically dubious and not recommended.

Override at run time without editing the YAML:

```bash
repseq global -c my.yaml -i x.fasta -T 0.95 --alphabet-for-clustering nucleotide
repseq global -c my.yaml -i x.fasta -T 0.95 --concatenate-markers      # or --no-concatenate-markers
```

---

## HMM-based identity QC

Annotation of viral CDSes in GenBank is famously inconsistent
("RNA-dependent RNA polymerase" vs "polymerase" vs "L protein" vs "RNA polymerase
L"), and pseudogenes / misannotations can carry plausible-looking `/product`
strings. repseq's primary use of **HMMER hmmscan** is as a **QC step** — it
verifies that each segment/sequence actually carries the expected proteins by
structural identity rather than by name, and drops anything that fails the gate.
Marker selection for protein-alphabet clustering is a downstream consumer of the
same hit cache.

The HMM tier is **off-by-default at the marker level** — it activates only for
markers that carry an `hmms:` list — so you can phase HMMs in one marker at a
time. It **runs regardless of clustering alphabet**: if you cluster on
nucleotides, HMM-QC still fires and drops isolates whose segments don't carry the
expected proteins.

### Token notation

Each element of `hmms:` is a **token** string:

| Form | Meaning |
| --- | --- |
| `"Name"`    | Single HMM. A CDS satisfies it if it has a passing hit to `Name`. |
| `"A--B"`    | Multidomain. Satisfied only when the CDS has passing hits to **both** `A` and `B`, with `A` N-terminal to `B` (forward-progressing endpoints, ≤ `hmm.multidomain_overlap_tolerance` aa overlap at the seam — default 30 aa). |
| `"A--B--C"` | Same, three domains: `A` most N-terminal, `C` most C-terminal. |

**N-to-C order.** The first HMM in a multidomain token is the most N-terminal
domain — the direction molecular biology writes a protein — so a token mirrors the
architecture you'd draw (coronavirus Spike is `"CoV_S1--CoV_S2"`: S1 N-terminal,
S2 C-terminal). The named domains are compared against the hmmscan
`ali_from`/`ali_to` columns in that order.

**Why the overlap tolerance.** Pfam profile boundaries rarely align exactly to a
real cleavage site (the CoV Spike S1/S2 furin seam, bunyaviridae G1/G2). A strict
non-overlap rule would silently drop those representatives; the default 30-aa
tolerance accepts the fuzz (`hmm.multidomain_overlap_tolerance: 0` for strict
non-overlap). **Full containment is always rejected** regardless of tolerance —
both `ali_from` and `ali_to` must progress forward between consecutive named
domains. **Extra domains are fine** — a CDS annotated `A--B--HMMX` still satisfies
`"A--B"`.

**Alternative architectures (OR).** Multiple tokens in one marker's `hmms:` list
are *alternatives* — a CDS satisfying **any one** satisfies the marker. This lets
one marker span divergent forms of a protein, e.g. Spike across alpha- and
beta-coronaviruses:

```yaml
- {name: Spike, aliases: ["spike", "surface glycoprotein"],
   hmms: ["CoV_S1--CoV_S2", "bCoV_S1_N--bCoV_S1_RBD--CoV_S2"]}
```

Across *different* markers (separate dict entries) the rule is AND: each is an
independent marker that must be present. The OR applies only within one marker's
token list.

### How the gate works

For a marker spec that includes `hmms: [<token>, ...]`:

1. Every CDS is scanned against the configured HMM database (once per run, batched
   into a single hmmscan call, cached on AA sequence so re-runs are free).
2. A hit **passes** when both gates clear: **similarity** — the curated Pfam GA
   threshold (when present) OR `hmm.default_evalue` (default `1e-5`); **coverage**
   — the alignment span covers ≥ `hmm.relative_length_cutoff` (default `0.5`) of
   the HMM model length.
3. A **segment / sequence passes a marker** when at least one CDS satisfies **any
   one** of that marker's tokens (OR within a spec; AND across separate specs).
4. If a marker has no satisfying CDS, the segment fails. In segmented mode the
   whole isolate is dropped (one bad segment fails the isolate, counted under
   `removed_hmm_failed`, the reason naming the unmatched alternatives joined with
   `|`). In non-segmented mode the single sequence is dropped.
5. For protein-alphabet clustering, the marker CDS of each survivor is the
   **longest** CDS satisfying any token in the spec.

Markers WITHOUT an `hmms:` list keep the legacy alias → longest chain and are not
gated by HMM-QC.

### Bundled vs user-supplied database

`hmm.database: null` (default) uses the bundled set at
`repseq/data/hmms/repseq_viral_core.hmm` — 26 hand-picked Pfam-A profiles covering
the most common viral marker proteins (RdRp, nucleocapsid, glycoprotein, helicase,
protease, the coronavirus Spike S1/S2 architectures, the small accessory proteins)
across the main families. Pfam-A is CC0, so redistribution is unrestricted. The
set was assembled via `scripts/build_bundled_hmms.sh`; re-run it to refresh against
a newer Pfam release.

To use your own profiles, set `hmm.database` to one of:

- **a single `.hmm` file** — the classic user-supplied database;
- **a directory of `.hmm` files** — every `*.hmm` is concatenated into one
  combined database (cached, rebuilt only when a member changes);
- **a bare bundled-set name** — a subdirectory under `repseq/data/hmms/`, e.g.
  `hmm.database: Filoviridae` selects `repseq/data/hmms/Filoviridae/`. Drop `.hmm`
  files into such a directory to add a family-specific set.

Whatever resolves is auto-`hmmpress`-ed on first use if its `.h3*` indexes are
missing or stale. `repseq doctor` reports DB status (path, profile count, indexed
Y/N, GA-cutoff coverage).

> **HMM hits are cached**, keyed on the CDS amino-acid sequence **and** a
> fingerprint of the database file. Swapping in a different DB (or rebuilding the
> bundled one) invalidates old entries automatically. Editing cutoffs needs no
> clear — cutoffs are reapplied each run. The one case where a manual clear helps
> is editing a `.hmm` file in place without changing its size or mtime, or
> reclaiming disk after a DB swap: `repseq cache clear --source hmmscan` wipes
> only the HMM hits. See [The local cache](#the-local-cache).

### Configuration examples

Non-segmented marker spec (top-level `clustering.cluster_protein`). Each entry is
either an alias string (legacy, alias-only) or a dict with `name`, optional
`aliases`, optional `hmms`:

```yaml
clustering:
  cluster_protein:
    - "polymerase"                                                          # legacy: alias-only
    - {name: "RdRp",  aliases: ["polymerase"], hmms: ["RdRP_4"]}            # single-HMM token
    - {name: "Spike", aliases: ["spike"],      hmms: ["CoV_S1--CoV_S2"]}    # multidomain token
```

Segmented marker spec (per-virus `segment_markers` — the recommended HMM-aware
form). Two flavours of the M-segment spec illustrate multidomain vs separate-token
notation:

```yaml
segmented:
  enabled: true
  virus: peribunyaviridae
  viruses:
    peribunyaviridae:
      expected_segments: 3
      segments: [S, M, L]
      isolate_regex: "..."
      segment_markers:
        S: {hmms: ["Bunya_nucleocap"]}            # single HMM
        # Permissive: passes when the M segment has Bunya_G1 hit AND
        # Bunya_G2 hit, possibly on different CDSes (post-cleavage Gn/Gc case).
        M: {hmms: ["Bunya_G1", "Bunya_G2"]}
        # OR — strict polyprotein form: a single CDS must carry both domains
        # in N-to-C order. Use this to reject post-cleavage annotations.
        # M: {hmms: ["Bunya_G2--Bunya_G1"]}
        L: {hmms: ["RdRP_4"]}
```

When both `segment_markers` and the legacy `cluster_protein` define a marker for
the same segment, `segment_markers` wins.

### HMM-related config keys

```yaml
hmm:
  enabled: true                       # master switch
  database: null                      # null → bundled core; .hmm file, dir of .hmm files, or bundled-set name
  default_evalue: 1.0e-5              # used when a profile has no curated GA
  use_ga_when_available: true         # prefer curated GA over default_evalue
  relative_length_cutoff: 0.5         # ali_span / hmm_len ≥ this
  multidomain_overlap_tolerance: 30   # max aa overlap at a multidomain seam (0 = strict)
  threads: null                       # null → falls through to cfg.threads
```

### Soft-fail posture

The HMM tier never aborts a run on its own. If `hmmscan` is missing, the database
is unreadable, or hmmscan returns non-zero, repseq emits one stderr warning and
falls through to the alias / longest-CDS chain. Run `repseq doctor` to check the
wiring.

---

## Accessory proteins (`extra_protein`)

`cluster_protein` (non-segmented) and `segment_markers` (segmented) declare the
**required** marker proteins: every isolate must carry one, the whole-genome tree
is built from them, and missing one fails HMM-QC. That's right for an RdRp or a
structural glycoprotein.

It's the **wrong** model for sparse accessory proteins. The coronavirus ORF7a,
ORF8, ORF9b, or the betacoronavirus-only HE appear in some isolates and not
others; adding them to `cluster_protein` would force HMM-QC to drop every
representative that doesn't carry them.

`extra_protein:` is the home for those. The spec shape is identical to
`cluster_protein` — `{name, aliases?, hmms?}` — but extras are **not** used for
clustering, **not** required for HMM-QC, and **not** part of the whole-genome
tree. They are reported where it makes sense:

- Per-marker protein FASTA in `{prefix}_extra_protein_fasta/` (always-on when any
  extra is declared).
- A separate per-protein tree in `{prefix}_extra_protein/` with
  `--per-protein-phylo`.
- Coverage and length columns in `{prefix}_protein_taxonomic_report.txt`, without
  the trailing `*`.
- **Excluded** from `{prefix}_incongruence.tsv` by design — sparse coverage would
  produce a lot of `NA` rows that drown out the real reassortment signal. For RF
  distances against an extra, run the math off its `<name>_tree.nwk` yourself.

```yaml
segmented:
  viruses:
    sarbecovirus:
      expected_segments: 1
      segments: [genome]
      isolate_regex: "..."
      segment_markers:
        genome:
          - {name: Spike, hmms: ["CoV_S1--CoV_S2",
                                 "bCoV_S1_N--bCoV_S1_RBD--CoV_S2"]}
          - {name: N,     hmms: ["CoV_nucleocap"]}
      extra_protein:
        genome:
          - {name: ORF7a, hmms: ["Corona_7"]}
          - {name: ORF3a, hmms: ["Corona_NS3b"]}
          - {name: HE,    aliases: ["hemagglutinin-esterase"]}
```

(The non-segmented form is a top-level `clustering.extra_protein: [...]` list
using the same dict shape.) Names must be unique across all `extra_protein`
entries for a virus.

---

## Polyprotein cutting (`polyprotein`)

Some virus families express their CDS as one giant precursor that is cleaved into
mature peptides: picornavirus P1/P2/P3, coronavirus ORF1ab (NSP1..NSP16),
flavivirus polyprotein (C/prM/E/NS1/…). Representative-level analysis on the whole
polyprotein loses signal — diversity comparisons, alignments, and trees are far
more informative on the **mature peptides**.

`polyprotein:` declares one or more polyproteins, each with the list of mature
peptides it cleaves into. Each peptide carries an **HMM token** (single name like
`CoV_NSP8`, or a multidomain N-to-C architecture like `CoV_RPol_N--RdRP_1` — same
syntax as `cluster_protein.hmms[]`). After representatives are elected, repseq
finds each rep's polyprotein CDS (the CDS satisfying the most peptide tokens),
slices it into mature peptides, and emits them as accessory artifacts.
**Clustering and the whole-genome tree are untouched** — the polyprotein still
drives those; mature peptides are additive output.

### Config shape

Non-segmented (one polyprotein per virus, common for Picornaviridae and
Coronaviridae):

```yaml
clustering:
  polyprotein:
    - name: ORF1ab                       # required, unique; keys the output filenames
      peptides:                          # required, ordered N→C, ≥2 entries
        - {name: NSP8,  hmm: CoV_NSP8,                            cleavage_motif: "LQ"}
        - {name: NSP9,  hmm: CoV_NSP9,                            cleavage_motif: "LQ"}
        - {name: NSP10, hmm: CoV_NSP10,                           cleavage_motif: "LQ"}
        - {name: NSP12, hmm: "CoV_RPol_N--RdRP_1",                cleavage_motif: "LQ"}
        - {name: NSP13, hmm: "CoV_NSP13_ZBD--CoV_NSP13_stalk--CoV_NSP13_1B", cleavage_motif: "LQ"}
      cut_strategy: motif                # boundary | bisect | motif (see below)
      motif_window_aa: 10                # ±aa around bisect to search for the motif
      min_peptides_hit: 2                # parent-CDS identification threshold
```

`hmm:` accepts a single profile name or a multidomain N→C token (HMMs joined with
`--`, most N-terminal first). For a multidomain peptide, the token is satisfied
only when every named HMM hits the parent CDS in N→C order; the peptide's
footprint is the synthetic span from the first domain's start to the last domain's
end (Pfam-boundary fuzz tolerated up to `hmm.multidomain_overlap_tolerance` aa).

**Alternative architectures (OR) — `hmms:`.** When a peptide's architecture
differs across genera, declare a list of alternative tokens with `hmms:` instead
of the singular `hmm:`. The peptide is located when **any one** token is
satisfied; if several match, the slicer picks the one with the best (lowest)
worst-domain E-value. The chosen architecture is recorded in a
`matched_architecture` column of the audit TSV and a `[matched_arch=...]` tag in
the peptide FASTA header.

```yaml
clustering:
  polyprotein:
    - name: ORF1ab
      peptides:
        # NSP1 is non-homologous between alpha- and beta-CoV — list both.
        - {name: NSP1, hmms: [aCoV_NSP1, bCoV_NSP1], cleavage_motif: "LQ"}
        - {name: NSP2, hmm: CoV_NSP2, cleavage_motif: "LQ"}
        - {name: NSP12, hmm: "CoV_RPol_N--RdRP_1", cleavage_motif: "LQ"}
```

Each peptide must set **exactly one** of `hmm` or `hmms` — both is an error
(ambiguous), neither is an error (peptide has no locator).

Segmented mode uses a per-segment dict (parity with `extra_protein`):

```yaml
segmented:
  viruses:
    flavivirus:
      ...
      polyprotein:
        genome:
          - name: polyprotein
            peptides:
              - {name: C, hmm: Flavi_C}
              - {name: M, hmm: Flavi_M}
              - {name: E, hmm: Flavi_E}
              - {name: NS1, hmm: Flavi_NS1}
              ...
```

Spec names must be unique across all `polyprotein:` entries within a virus.

### Cut strategies (the load-bearing science choice)

1. **`boundary`** — each peptide spans its HMM hit's `ali_from..ali_to` verbatim.
   Deterministic, but loses inter-domain residues at the seams.
2. **`bisect`** — cuts at the midpoint between adjacent peptide HMM hits.
   Endpoints extend to the protein N-/C-term **only when the first/last declared
   peptide is itself located** — a missing leading/trailing peptide leaves the gap
   unassigned. Cut sites are geometric rather than biological.
3. **`motif`** (default when any peptide declares `cleavage_motif`) — start from
   `bisect`, then snap each inter-peptide cut to the last occurrence of the
   downstream peptide's `cleavage_motif` within ±`motif_window_aa` of the bisect
   point (falls back to `bisect` for that cut if no motif is in the window).
   `cleavage_motif:` is the residues just **N-terminal** of the cut that liberates
   that peptide — `LQ` for coronavirus 3CL, `Q` for picornavirus 3C, `RR`/`R[RK]`
   for flavivirus NS2B-3.

The default is `motif` if any peptide declared `cleavage_motif`, else `bisect`.
`boundary` is opt-in.

### Edge cases

- **Peptide HMM doesn't hit** the parent CDS → peptide skipped, audit row marks
  `missing`. Flanking peptides keep their HMM-hit boundary on the missing side
  (`cut_method_actual=hit-boundary`); the missing peptide's territory is left
  unassigned rather than absorbed by a neighbour.
- **Peptide HMMs hit out of N→C order** → spec fails for that rep, all its rows
  land as `out_of_order` (almost always an HMM-database problem).
- **Adjacent peptide HMM hits overlap** → cut placed at the midpoint of the hit
  centres; sequences still emitted but rows flagged `overlap`.
- **No CDS has hits from ≥ `min_peptides_hit` peptide HMMs** → no parent
  identified; no peptide FASTAs for that rep, and one `no_parent_cds` audit row
  per peptide so the user sees *why*.

### Outputs

`{prefix}_polyprotein/{prefix}_<spec>_<peptide>.fasta` — one FASTA per peptide of
each spec, with all clean (`ok` / `overlap`) slices across representatives.
Headers carry the standard protein-FASTA bracket tags plus polyprotein-specific
ones: `[polyprotein=<spec>]`, `[peptide=<peptide>]`, `[peptide_range_aa=<from>-<to>]`,
`[cut_method=<actual>]`, `[matched_arch=<token>]`.

`{prefix}_polyprotein/{prefix}_<spec>_polyprotein.fasta` — the **whole
polyprotein** FASTA: the entire parent CDS (amino acids) for every representative
carrying it, one record each, same bracket-tag headers. Always written (parallel
to the per-peptide FASTAs); independent of the `whole_polyprotein_tree` knob below.

`{prefix}_polyprotein/{prefix}_<spec>_peptides.tsv` — audit TSV, one row per
(representative × peptide) attempt: `isolate_id`, `peptide_name`,
`parent_accession`, `parent_protein_id`, `range_aa_from`, `range_aa_to`,
`length_aa`, `cut_method_actual`, `matched_architecture`, `status`, `note`. Use it
to spot heuristic vs motif-supported cuts and which OR-alternative fired on each
rep.

**Per-peptide phylogenetic trees** (opt-in via `--per-protein-phylo`). Each
peptide of each spec gets its own ML tree built the same way the marker trees are
— MAFFT, IQ-TREE or FastTree per `phylo.tool`, rooting + LCA + shared colour
palette, optional trimAl. Settings come from `phylo.per_protein` verbatim
(`min_taxa`, `mafft.extra_args`, `trimal.*`, `domain_architecture`). Files sit
alongside the peptide FASTAs:

* `{prefix}_<spec>_<peptide>_msa.fasta`
* `{prefix}_<spec>_<peptide>_tree.xml`
* `{prefix}_<spec>_<peptide>_tree.pdf` / `.png` (figure, default-on `phylo.pdf`)
* `{prefix}_<spec>_<peptide>_tree_id_map.tsv`
* `{prefix}_<spec>_<peptide>_tree.nwk` (opt-in, `phylo.newick`)

Sparse peptides (fewer than `phylo.per_protein.min_taxa` reps with an `ok` slice)
are skipped. Peptide trees are **not** included in `_incongruence.tsv` (that
compares whole-genome markers). The phyloXML `<domain_architecture>` on each leaf
is **peptide-local** (boxes on the peptide's 1..length coordinate space), so
multidomain peptides like NSP12 (`CoV_RPol_N--RdRP_1`) show both subdomain boxes
properly scoped to the cut peptide.

**Whole-polyprotein tree** (opt-in via `phylo.per_protein.whole_polyprotein_tree`,
default off). In addition to the per-peptide trees, build **one tree per spec on
the entire polyprotein CDS** — not the sliced peptides. Each rep's parent CDS feeds
the same `_build_tree` engine, so the files mirror the peptide-tree schema with a
`_polyprotein` basename:

* `{prefix}_<spec>_polyprotein_msa.fasta`
* `{prefix}_<spec>_polyprotein_tree.xml`
* `{prefix}_<spec>_polyprotein_tree.pdf` / `.png`
* `{prefix}_<spec>_polyprotein_tree_id_map.tsv`
* `{prefix}_<spec>_polyprotein_tree.nwk` (opt-in, `phylo.newick`)

The point is the phyloXML `<domain_architecture>`: where a peptide tree shows one
peptide-local box, the whole-polyprotein tree's leaf shows **every HMM domain
across the entire polyprotein, end-to-end**, so Archaeopteryx draws the complete
domain diagram on each leaf. Off by default (one extra large alignment + tree per
spec); only fires under `--per-protein-phylo` with the HMM tier active. Like the
peptide trees, **not** added to `_incongruence.tsv`.

### Soft-fail posture

When the HMM tier didn't run (HMMER not installed, no HMMs declared), polyprotein
cutting prints a single stderr line and emits no outputs. Specs whose every rep
yields `no_parent_cds` still write the audit TSV (so the user sees why no FASTAs
appeared) but no FASTAs.

---

## The local cache

Every NCBI/UniProt lookup is saved to a small SQLite database (`~/.repseq/cache/`
by default) so you only pay the network cost once. Several "kinds" of lookups
share this cache, each under its own **source** name:

| Source | What's cached |
|---|---|
| `ncbi_taxonomy` | Taxonomy lineages for accessions / taxon IDs. |
| `ncbi_proteins` | GenBank CDS records (`/isolate`, `/segment`, protein translations). |
| `ncbi_nuc_seq` | Nucleotide bodies fetched by `repseq replay`. |
| `uniprot` | UniProt entries pulled by accession. |
| `hmmscan` | Per-CDS HMM hit lists from `hmmscan` (see below). |

```bash
repseq cache stats -c your_config.yaml                          # how big, what's in it
repseq cache purge-expired -c your_config.yaml                  # remove stale entries (TTL)
repseq cache clear -c your_config.yaml                          # wipe everything (asks first)
repseq cache clear -c your_config.yaml --source ncbi_taxonomy   # wipe just one source
```

`--source` accepts any of the source names above.

### When (and when not) to clear the `hmmscan` cache

The `hmmscan` cache stores, for every CDS amino-acid sequence the pipeline has
seen, the raw HMM hit list against your database. Because viral CDSes are very
repetitive across related isolates, this pays for itself within a single run.

1. **Cutoff changes do NOT require a clear.** The cache stores raw hits; the
   pass/fail flags (E-value, GA, coverage, multidomain order, overlap tolerance)
   are recomputed every run from your live config.
2. **DB changes usually do NOT require a clear either.** Each entry is keyed on
   `(AA sequence, database fingerprint)`, the fingerprint derived from the `.hmm`
   file's size and mtime — so pointing at a different file, or rebuilding the
   bundled DB, automatically invalidates old entries. They fall out via TTL
   (`taxonomy.cache_ttl_days`, 90 days by default).

You only need to manually clear it when: reclaiming disk after a DB swap; you
edited the `.hmm` file in place without changing its size or mtime; or you suspect
cache corruption. In those cases:

```bash
repseq cache clear -c your_config.yaml --source hmmscan
```

This wipes only the `hmmscan` rows; NCBI taxonomy, GenBank CDS, and UniProt caches
stay intact. Add `--yes` to skip the confirmation.

---

## Troubleshooting

**"WARNING: no representative sequences were selected."** The run finished but
nothing came out. `repseq` prints the most likely cause; the usual ones:

- *No sequences were loaded* — wrong input path, empty file, or unrecognised
  header format. Try `--source ncbi_virus` (or `ncbi` / `uniprot`) to force it.
- *QC removed everything* — thresholds too strict. Look at `_qc_removed.tsv` to
  see which step did it, then loosen that setting. A common one: `median_percent`
  length filtering on a mixed-gene file — switch to `min_max`.
- *The segmented step dropped everything* — no isolate had all its segments. Most
  often the `isolate_regex` doesn't match your headers; also check segment
  names/aliases and any `segment_lengths` bounds.

**`MMseqs2Error` / "mmseqs not found"** (or **`CDHitError` / "cd-hit not found"**
with `clustering.backend: cdhit`) — the clustering program isn't installed or
isn't on `PATH`. Install it (see [Installation](#installation)), or use a mode
that doesn't need it (`global -n`, or a stratified mode where every group is
already small).

**`cd-hit identity threshold X is below the supported floor`** — cd-hit refuses
identities below 0.40 (protein) or 0.80 (nucleotide). Raise your threshold, or
switch the backend to `mmseqs2` (no floor).

**Everything is grouped under "Unknown"** — the metadata lookups didn't run or
found nothing. Don't use `--no-resolve` for the real run, and set `ncbi_email`.

**Lookups are slow the first time** — expected; they're cached, so the *second*
run is fast. An NCBI API key speeds up the first run.

**Not sure if everything's installed?** Run `repseq doctor`. It checks every
required and optional dependency, the external tools (`mmseqs`, `cd-hit`,
`cd-hit-est`, `mafft`, `FastTree`, `iqtree2`, `trimal`, `hmmscan`), reaches NCBI
and UniProt to confirm the network, and verifies your config — then tells you in
plain English what needs fixing. Add `--no-network` to skip the database pings.

**`[phylo skipped]` / `[phylo failed]`** — the `--phylo` step is fail-soft: if
fewer than 3 representatives survived, or `mafft` / `iqtree2` / `FastTree` are
missing or errored, the message is printed to stderr and the rest of the run's
outputs are still written. IQ-TREE refuses UFBoot with fewer than 4 sequences; the
wrapper falls back to no-bootstrap automatically.

**`CDHitError: Cluster round-trip mismatch: N sequences in, M accounted for…`** —
the cd-hit wrapper does a strict round-trip check on its `.clstr` output and
refuses to silently undercount. The error names the seq IDs that fell out (both
directions). The usual cause is a seq id with characters cd-hit transforms —
internal `...`, very long IDs, or whitespace.

**`taxonomy_mismatch:<rank>` rows in `_qc_removed.tsv`** — isolates dropped because
their segments resolved to different taxa at the configured rank (reassortants or
`/isolate` collisions). To keep them, set
`segmented.taxonomy_consistency.enabled: false` (or relax the rank to e.g.
`genus`).

**Per-segment length-filter drops surprised you** — open `_run.log` and scroll to
the "QC SUMMARY" block, which prints a per-segment, isolate-level breakdown
(`L too short : 257` etc.). Widen the `segment_lengths` bounds for the bottleneck
segment.

**`Strain-collision check: found N (isolate, segment) pair(s)…`** — the detector
found two or more distinct accessions sharing a strain-derived `isolate_id` *and*
a segment, the over-merge signature of the `/strain → isolate_id` fallback. By
default the run continues (action: `warn`); to remove them set
`segmented.strain_collision_action: drop` (they then appear in `_qc_removed.tsv`
with reason `strain_collision:<segment>`). The `isolate_id_source` column shows
which records relied on the fallback.

---

## Testing

```bash
pip install pytest
pytest tests/
```

The default suite runs fully offline — all network calls are simulated — and
finishes in a couple of seconds. An opt-in integration test runs the real
mafft / FastTree / mmseqs end-to-end; run it with `pytest -m integration`.

---

## License

`repseq` is free software, released under the **GNU General Public License,
version 3** (GPLv3). You may redistribute it and/or modify it under the terms of
that license; it comes with **no warranty**. The full license text is in the
[`LICENSE`](LICENSE) file, and at <https://www.gnu.org/licenses/gpl-3.0.html>.

---

## Status

Current: **`v0.58.0`**. All 8 selection modes; protein-alphabet clustering by
default (`alphabet_for_clustering: protein`); MMseqs2 and cd-hit backends;
optional protein-annotation and protein-quality QC; per-isolate
taxonomy-consistency QC and strain-as-isolate provenance + collision detection for
segmented viruses; **HMM-based identity QC** with multidomain-architecture token
notation (HMMER hmmscan + a bundled 26-profile viral Pfam-A subset; user-supplied
DB supported); **rich phyloXML output** with taxonomy-driven leaf colouring,
partitioned-supermatrix trees for protein + IQ-TREE (default), optional trimAl
trimming, graphical PDF/PNG tree figures (default on), per-marker
domain-architecture trees with a pairwise Robinson-Foulds incongruence table,
per-segment NT trees, a pre-cluster overview tree, an opt-in phylogeny-based
taxonomy review, per-MSA conservation scoring, and an always-on per-taxon
monophyly report; **`extra_protein`** accessory-protein and **`polyprotein`**
mature-peptide channels (FASTAs, trees, reports); a force-keep / force-select /
input-blocklist `overrides` system; an optional UMAP/MDS clustering plot;
plain-English flags and a single-file HTML report; an auto-generated
Methods-section starter (`_summary.md`), a reproducibility lockfile
(`repseq replay`), and a sanitized config snapshot on every run. **1,419 offline
regression tests pass**; the
NCBI-backed paths have been validated end-to-end against live influenza-A,
peribunyaviridae, hantaviridae, and coronaviridae datasets. See `git log` for the
full release history.
