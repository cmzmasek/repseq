# repseq

**Pick a small, clean, representative set of sequences out of a big messy FASTA file.**

You downloaded 80,000 influenza sequences from NCBI. You need 200 good ones that
still cover the real diversity — for a tree, for a figure, for a training set, or
to seed a curated reference database. Doing that by hand is miserable. `repseq`
does it for you: it cleans up the data, looks up what each sequence actually is,
groups similar sequences together, and keeps one good example from each group.

It works on protein **or** nucleotide FASTA files, from UniProt, NCBI, or NCBI
Virus, and it has a strong focus on viral sequences (including segmented viruses
like influenza).

---

## What "representative selection" means here

If you have 5,000 nearly-identical H3N2 sequences and 3 unusual ones, a random
sample of 200 will be 200 H3N2 sequences and you'll lose the unusual ones.
`repseq` instead:

1. **Cleans** — drops duplicates, truncated junk, sequences full of `N`s, and
   records labelled "hypothetical", "synthetic", "partial", etc.
2. **Identifies** — looks up each sequence's organism, host, country, and
   collection date from NCBI/UniProt (results are cached, so it's only slow once).
3. **Groups** — buckets the sequences (by how similar they are, or by genus, host,
   year, country… your choice).
4. **Keeps the best example from each group** — preferring curated records
   (RefSeq, reviewed UniProt) over random ones, and longer over shorter.

The result is a FASTA file that's small enough to work with but still spans the
diversity that was in the original.

---

## Installation

You need **Python 3.10 or newer**. Then:

```bash
git clone https://github.com/cmzmasek/repseq.git
cd repseq
pip install -e .
```

That's it for most uses. Two optional pieces:

**A sequence-clustering program** — `repseq` calls one to group sequences
*by similarity*. You only need it for the similarity-based modes (see the
table below); if you group by genus/host/year and your groups are small, or
you just want N diverse sequences, you can skip it. Two backends are
supported; pick **one**.

*MMseqs2* (default — fast, scales to very large datasets):

```bash
brew install mmseqs2                       # macOS
conda install -c bioconda mmseqs2          # Linux (conda)
```

*cd-hit* (classic alternative — set `clustering.backend: cdhit` in your
config). Slower than MMseqs2 on big inputs and has minimum-identity floors
(0.40 for protein, 0.80 for nucleotide), but produces tight, all-vs-all
clusters that many groups prefer for reference-set work:

```bash
brew install cd-hit                        # macOS
conda install -c bioconda cd-hit           # Linux (conda)
```

`repseq` finds whichever binary it needs (`mmseqs`, or `cd-hit` /
`cd-hit-est`) on your `PATH` automatically.

**Plots & tree figures** — the `[viz]` extra (matplotlib only — installs
cleanly everywhere) powers both the optional diagnostic scatter plot of the
clustering (`--plot`) **and** the graphical tree-figure PDFs/PNGs that
`--phylo` renders by default (`phylo.pdf`):

```bash
pip install -e '.[viz]'
```

For the optional UMAP embedding upgrade (falls back to classical MDS if
not installed):

```bash
pip install -e '.[viz-umap]'
```

**MSA + tree (optional, `--phylo`)** — to build a multiple-sequence
alignment and a maximum-likelihood phylogeny over the final
representatives, install MAFFT plus a tree builder. `repseq` auto-picks
**IQ-TREE** for protein alignments (ModelFinder + UFBoot — best topology
quality) and **FastTree** for nucleotide alignments (much faster on
large NT MSAs). Install whichever pair you'll use:

```bash
brew install mafft iqtree fasttree              # macOS
conda install -c bioconda mafft iqtree fasttree # Linux (conda)
```

If MAFFT or the relevant tree builder is missing the rest of the run
still finishes; `repseq` just skips the tree step with a clear
stderr message. Pin one tool with `phylo.tool: iqtree` (or `fasttree`).

---

## Quickstart

**Step 1 — make a config file.** This is a small text file with your settings
(cleaning thresholds, your NCBI email, etc.). A wizard asks you the questions:

```bash
repseq init-config -o my_config.yaml
```

**Step 2 — run a mode.** Point it at your FASTA file and pick how to select:

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

> **Tip:** add `--no-resolve` while you're experimenting. It skips the NCBI/UniProt
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
| `host` | Up to `-n` sequences per host organism. | You care about which host the virus came from. |
| `time` | Up to `-n` per time window — `--window year`, `decade`, or a number like `5` for 5-year bins. | You want even coverage across collection dates. |
| `geographic` | Up to `-n` per country. | You want even geographic coverage. |
| `custom` | Up to `-n` per *any* field — a sequence attribute, a taxonomy rank, or a column in a metadata spreadsheet you provide. | Your grouping isn't one of the built-ins. |
| `hybrid` | Up to `-n` per *combination* of fields, e.g. `--fields genus,host,decade`. | You want a balanced grid across several variables at once. |

**Within each group**, if the group already has `-n` or fewer sequences, all of
them are kept. If it's bigger, `repseq` clusters it down to about `-n`
representatives. Add `--overflow trim` if you need *exactly* `-n` and not "about
`-n`".

The `-n` "maximally-different" selection (`global -n`, and the `--overflow trim`
step in every mode) is an alignment-free **MaxMin** sampler over k-mer distance.
It runs on the **same representation your clustering uses**: with
`--alphabet-for-clustering protein` (the default) it compares each sequence's
**marker protein**, and with `nucleotide` it compares the genome. Note that
nucleotide k-mer diversity is inherently weak for highly diverged viral
families (short k-mers saturate on whole genomes, long ones share nothing
across genera), so **protein is the recommended alphabet** when you want `-n`
to track real divergence.

Every mode also accepts: `--input/-i`, `--output-dir/-o`, `--config/-c`,
`--threads`, `--seed`, `--segmented`, `--dry-run`, `--no-resolve`,
`--source {auto,uniprot,ncbi,ncbi_virus}`, `--overflow {keep,trim}`, `--plot`,
`--phylo`, `--per-protein-phylo`, `--per-segment-phylo` (segmented only),
`--pre-cluster-tree`, `--fast`, `--verbose`,
`--alphabet-for-clustering {protein,nucleotide}`,
`--protect-ids FILE` (force-keep a list of special-importance sequences
through QC), `--pin-ids FILE` (force-select a list as representatives) —
see [Sequences of special importance](#sequences-of-special-importance--overrides).

In addition to the selection modes, two diagnostic/utility subcommands:

| Command | What it does |
| --- | --- |
| `repseq stats -i x.fasta` | Pre-flight inspection of an input FASTA — count, source breakdown, NT length distribution, per-field metadata coverage (organism, host, subtype, country, segment, taxonomy, …), RefSeq/reviewed flag counts, and a top-N taxonomy breakdown per rank. No clustering, no output files; use this before committing to a long run. Default `--no-resolve` (fast, network-free); pass `--resolve` to fetch taxonomy. |
| `repseq replay <lockfile.json> -o <out>` | Re-materialise representatives from a `{prefix}_lockfile.json` written by an earlier run: re-fetches the listed accessions from NCBI (via the cache) and re-emits the same representative FASTAs into a fresh output directory. QC, clustering, and selection are **not** re-run — the lockfile is authoritative. Trees are skipped unless `--rebuild-trees`. v0.30.0. |
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
  write:
    • selected representatives (FASTA, per-segment FASTAs, concatenated FASTA)
    • all proteins of each representative (one big AA FASTA + per-CDS TSV)
    • one AA FASTA per declared marker        →  {prefix}_per_protein_fasta/
    • one AA FASTA per declared extra_protein →  {prefix}_extra_protein_fasta/
    • mature peptides for declared polyproteins →  {prefix}_polyprotein/        (peptide + whole-polyprotein FASTAs + audit TSV)
    • per-rep metadata spreadsheet
    • per-stratum + per-cluster + per-drop TSVs
    • taxonomic diversity report              →  {prefix}_taxonomic_report.txt
    • per-marker coverage + AA length stats   →  {prefix}_protein_taxonomic_report.txt
    • per-rank NT length stats                →  {prefix}_nucleotide_taxonomic_report.txt
    • per-peptide coverage + AA length stats  →  {prefix}_polyprotein_taxonomic_report.txt
                                                  (when polyprotein: declared)
    • tidy TSV companions to all 4 reports    →  {prefix}_*_taxonomic_report.tsv
                                                  (8-col tidy schema; Excel-pivot / R-pandas friendly)
    • Methods-section starter                 →  {prefix}_summary.md
    • sanitized effective-config snapshot      →  {prefix}_config_repseq<ver>.yaml
    • plain-text run log                      →  {prefix}_run.log
    • (optional) UMAP/MDS clustering plot                        (--plot)
    • (optional) MAFFT MSA + IQ-TREE / FastTree + phyloXML tree  (--phylo)
    •            (protein + IQ-TREE: partitioned-supermatrix tree, per-marker model)
    •            (optional trimAl column trimming; off by default)
    • (optional) one tree per HMM marker      →  {prefix}_per_protein/      (--per-protein-phylo)
    •            one tree per extra_protein   →  {prefix}_extra_protein/    (--per-protein-phylo)
    • (optional) one NT tree per segment      →  {prefix}_per_segment/      (--per-segment-phylo)
    • per-(rep, HMM profile) diagnostic       →  {prefix}_hmm_diagnostic.tsv (when HMM tier ran)
```

The run log (`{prefix}_run.log`) records exactly what settings were used
and what got dropped at each step — keep it with your results so the
selection is reproducible.

---

## Output files

Everything is written to a single output directory (`./repseq_output/` by
default, set via `output.dir` in the config or `-o/--output-dir` on the
command line). `repseq` **refuses to write into a non-empty directory** so a
run can never silently overwrite or mix into a previous one — empty it,
delete it, or point somewhere new.

The exact set of files depends on whether you ran in **segmented** mode and
which optional flags you passed (`--plot`, `--phylo`). At a glance:

| File family | Always written? | What it is |
| --- | --- | --- |
| `{prefix}_run.log` | yes | Settings used and per-step counts. Keep with your results. |
| `{prefix}_qc_removed.tsv` | yes | Every dropped sequence and the reason. |
| `{prefix}_overrides.tsv` | when any sequence was force-kept | Force-keep audit: each sequence the `overrides.protect_qc` whitelist rescued from a QC stage, with the reason it would otherwise have been dropped. |
| `{prefix}_force_selected.tsv` | when any sequence was force-selected | Force-select audit: each `overrides.force_select` / `--pin-ids` pin and the action taken (elected / split / added / already-a-rep / unavailable). |
| `{prefix}_group_counts.tsv` | yes | One row per stratum: in / out / clustered / cutoff. |
| `{prefix}_clusters.tsv` | yes | Which sequences ended up grouped, who represents whom. |
| `{prefix}_representatives.fasta` | non-segmented | Selected representative sequences. |
| `{prefix}_representative_sequences.tsv` | non-segmented | Per-representative metadata spreadsheet (since v0.21.0: same schema as `_representative_isolates.tsv`). |
| `{prefix}_concatenated.fasta` | segmented | Per-isolate head-to-tail concat of all segments. |
| `{prefix}_segment_<NAME>.fasta` | segmented | One file per expected segment, just the representative isolates. |
| `{prefix}_representative_isolates.tsv` | segmented | Per-representative-isolate metadata spreadsheet. |
| `{prefix}_isolate_proteins.tsv` | segmented + GenBank | Every protein of every isolate that survived QC, with a `representative` TRUE/FALSE column. |
| `{prefix}_representative_isolate_proteins.tsv` | segmented + GenBank | Same schema as above, row-filtered to representatives only. |
| `{prefix}_sequence_proteins.tsv` | non-segmented + GenBank | Per-CDS counterpart of `_isolate_proteins.tsv` (since v0.21.0). |
| `{prefix}_representative_sequence_proteins.tsv` | non-segmented + GenBank | Reps-only filtered view of the above. |
| `{prefix}_representative_isolate_proteins.fasta` | segmented + GenBank | AA FASTA of every protein of every representative isolate. |
| `{prefix}_representative_sequence_proteins.fasta` | non-segmented + GenBank | AA FASTA of every protein of every representative sequence. |
| `{prefix}_per_protein_fasta/{prefix}_<family>.fasta` | any run with `cluster_protein` / `segment_markers` declared | Unaligned per-marker protein FASTA, one CDS per rep carrying that marker (always-on since v0.22.0). |
| `{prefix}_extra_protein_fasta/{prefix}_<name>.fasta` | any run with `extra_protein:` declared | Same shape, for accessory proteins that aren't required everywhere (v0.22.0). |
| `{prefix}_polyprotein/{prefix}_<spec>_<peptide>.fasta` + `{prefix}_<spec>_peptides.tsv` | any run with `polyprotein:` declared and HMM tier active | Mature peptides sliced from each representative's polyprotein CDS, one FASTA per peptide of each spec, plus an audit TSV recording every (rep × peptide) attempt with status (`ok` / `missing` / `out_of_order` / `overlap` / `no_parent_cds`). v0.33.0. |
| `{prefix}_polyprotein/{prefix}_<spec>_polyprotein.fasta` | any run with `polyprotein:` declared and HMM tier active | Unaligned **whole-polyprotein** FASTA — the entire parent CDS per representative, parallel to the per-peptide FASTAs. v0.50.0. |
| `{prefix}_polyprotein/{prefix}_<spec>_polyprotein_tree.{xml,pdf,png}` (+ `_msa.fasta`, `_tree_id_map.tsv`) | `--per-protein-phylo` + `phylo.per_protein.whole_polyprotein_tree` | **Whole-polyprotein tree** — one tree per spec on the entire polyprotein CDS; its phyloXML `<domain_architecture>` shows every HMM domain across the whole polyprotein end-to-end. Off by default. v0.50.0. |
| `{prefix}_representatives_protein.fasta` | when `alphabet_for_clustering: protein` (default) | The AA strings actually fed into the clusterer. |
| `{prefix}_taxonomic_report.txt` | every run | Per-rank diversity table: distinct taxa before vs after clustering. |
| `{prefix}_protein_taxonomic_report.txt` | any run with `cluster_protein` / `segment_markers` / `extra_protein` declared | Per-rank protein coverage + AA length statistics (v0.22.0). |
| `{prefix}_nucleotide_taxonomic_report.txt` | every run | Per-rank NT length statistics: per-segment + `total` (segmented) or single `genome` column (non-segmented). v0.23.0. |
| `{prefix}_polyprotein_taxonomic_report.txt` | any run with `polyprotein:` declared and HMM tier active | Per-rank peptide coverage + AA length statistics (one column per declared peptide of each spec) — sliced-peptide analogue of `_protein_taxonomic_report.txt`. v0.35.0. |
| `{prefix}_taxonomic_report.tsv`, `_protein_taxonomic_report.tsv`, `_nucleotide_taxonomic_report.tsv`, `_polyprotein_taxonomic_report.tsv` | parallel to the matching `.txt` reports | **Machine-readable tidy long-format TSV** companions (8-column schema: `report, rank, pool, taxon, taxon_count, spec, metric, value`). Excel-pivot / R-pandas friendly. v0.35.0. |
| `{prefix}_clustering.png` | only with `--plot` | Diagnostic scatter of the clustering. |
| `{prefix}_msa.fasta`, `_tree.xml`, `_tree_id_map.tsv` | only with `--phylo` | Alignment + annotated tree (phyloXML) + short-id↔accession map. |
| `{prefix}_tree.pdf`, `_tree.png` | `--phylo` (default; `phylo.pdf: true`) | **Graphical tree figure** — a ladderized rectangular phylogram (matplotlib + Bio.Phylo) with taxonomy-coloured leaf labels, a genus/subfamily legend, and branch-support labels, rendered from the phyloXML. **On by default**; `--no-pdf` (or `phylo.pdf: false`) disables it. Emitted for *every* tree below (`*_tree.pdf` / `*_tree.png`). Needs matplotlib (`[viz]`); soft-fails if missing. v0.49.0. |
| `{prefix}_tree.nwk` | `--phylo` **and** `phylo.newick: true` (or `--newick`) | Plain-text Newick. **Off by default** (`phylo.newick: false`) to cut clutter — the phyloXML is a topological superset. Same for every tree below (`*_tree.nwk`). |
| `{prefix}_partition.nex`, `_msa_<family>.fasta` | `--phylo`, protein + IQ-TREE | NEXUS partition file + per-family alignments (partitioned-supermatrix tree). |
| `{prefix}_msa_untrimmed.fasta` | `--phylo` + `phylo.trimal.enabled` | Raw MAFFT alignment retained when trimAl trimming ran (`_msa.fasta` is then the trimmed tree input). |
| `{prefix}_taxonomy_review.tsv` | `--phylo` + `phylo.taxonomy_review.enabled` | Phylogeny-based taxonomy verdicts: blank ranks imputed from the enclosing clade (`impute_missing`) and populated-but-disagreeing ranks flagged (`conflict_flag`), with evidence + confidence. v0.39.0. |
| `{prefix}_representative_*_corrected.tsv`, `_representative_*_proteins_corrected.fasta` | above + `write_corrected` | Copies of the rep TSV + protein FASTA with **high-confidence imputed blanks filled** (clean values; originals kept). v0.39.0. |
| `{prefix}_iqtree_summary.txt` | only with `--phylo` + IQ-TREE | IQ-TREE ModelFinder report. |
| `{prefix}_iqtree_model.txt` | only with `--phylo` + IQ-TREE | Grep-friendly sidecar with the ModelFinder pick(s): one `<label>: <model>` line per partition (or a single `GENOME: <model>` line for the non-partitioned path). v0.28.0. |
| `{prefix}_per_protein/` | only with `--per-protein-phylo` | One tree (MSA + phyloXML + id map + PDF/PNG figure by default; Newick opt-in via `phylo.newick`) per marker; plus `_incongruence.tsv` of pairwise Robinson-Foulds distances. |
| `{prefix}_pre_cluster_tree.xml`, `_id_map.tsv` (+ `.pdf` / `.png`) | only with `--pre-cluster-tree` (or `phylo.pre_cluster_tree.enabled`) | Rough overview tree of every post-QC sequence (one leaf per CONCAT isolate in segmented mode), with `[repr] ` prefix on representative leaves in the phyloXML `<name>`. A PDF/PNG figure is rendered by default (`phylo.pdf`); the `.nwk` is opt-in (`phylo.newick`). v0.32.0. |
| `{prefix}_extra_protein/` | `--per-protein-phylo` + `extra_protein:` declared | Same shape, accessory-protein trees (kept out of the incongruence table by design). |
| `{prefix}_per_segment/` | only with `--per-segment-phylo` (segmented only) | One **nucleotide** tree per declared segment (each with a PDF/PNG figure by default), from the representative isolates' raw per-segment NT. Complement to the per-marker AA trees: reassortment may show up here even when no single marker tree captures it. v0.26.0. |
| `{prefix}_msa_conservation.tsv` | any phylo step that wrote an MSA (default on; `phylo.conservation.enabled`) | One conservation number per alignment written this run (genome, partition families + supermatrix, per-protein, extra-protein, polyprotein peptides, per-segment NT). Mean per-column **Jensen-Shannon divergence to a residue background** (Capra & Singh 2007), **Henikoff-weighted** + **gap-penalised**. Columns: `msa, role, label, alphabet, trimmed, n_seqs, n_sites, mean_conservation, mean_conservation_core`. v0.43.0. |
| `{prefix}_monophyly.tsv` | any phylo step that built a tree (always on) | Per-taxon monophyly status, one row per (tree, rank, taxon): `monophyletic` / `paraphyletic` / `polyphyletic` plus the counts behind the call. Swept off every annotated `*_tree.xml` (genome + per-protein/extra/segment/polyprotein/pre-cluster). A taxon monophyletic on the genome tree but not on a marker tree is the per-marker reassortment/recombination signal — the taxon-resolved companion to `_incongruence.tsv`. Support-aware by default (`phylo.monophyly.min_support`, default 70 — collapses weak branches so only confident non-monophyly is flagged). Columns: `tree, rank, taxon, n_leaves, status, n_clusters, n_intruders, intruder_clusters, intruder_taxa, min_support`. v0.54.0 (support-aware v0.56.0). |
| `{prefix}_flags.txt` | when any conflict table exists (a tree was built) | **Plain-English synthesis** of the conflict tables (`_monophyly.tsv`, `_incongruence.tsv`, `_taxonomy_review.tsv`) into one skimmable list of reassortment/recombination signals, non-monophyletic taxa, and per-isolate taxonomy conflicts. A clean run gets a "No flags raised ✓" file. v0.55.0. |
| `{prefix}_report.html` | any run with figures or a conflict table | **Single-file HTML report** bundling the flags, a gallery of every tree figure (embedded), and an index of all output files — one page to open in a browser or send to a collaborator. v0.55.0. |
| `{prefix}_hmm_diagnostic.tsv` | when the HMM tier ran AND a spec declares `hmms:` | One row per (representative, marker spec, declared HMM profile): `hit` T/F, `best_dom_evalue`, `best_dom_score`, `best_coverage`, `ga_cutoff`, `cutoff_used`, `passing`. Surfaces near-miss patterns the binary gate hides ("65 isolates barely missed the cutoff at E=2e-5"). v0.26.0. |
| `{prefix}_summary.md` | every run | Auto-generated Methods-section starter (prose + numbers + tool citations). Leads with an **"Analysis at a glance"** table (v0.41.0) stating the run's actual decisions — dataset type, input→representative counts, clustering substrate + alphabet + tool, selection mode, whole-genome tree substrate, per-marker/per-segment tree presence, taxonomy review. The **Clustering substrate** and **Whole-genome tree** cells name the actual marker(s) — the canonical `name:`, not the aliases/synonyms — and, when the HMM tier ran, each marker's domain architecture (`name (HMM: tok OR tok)`); single-marker mode lists them in priority order, concat joins with `+`, segmented prefixes each with its segment. v0.42.0. |
| `{prefix}_config_repseq<ver>.yaml` | every run | **Sanitized effective-config snapshot** (v0.42.0): the full run-time config — every setting resolved, defaults filled in — with all comments stripped and the NCBI credentials (`ncbi_email` / `ncbi_api_key`) blanked to `null`. Re-runnable as-is via `repseq <mode> -c <this-file>`. Filename embeds the repseq version (periods→underscores), e.g. `cov_config_repseq0_42_0.yaml`. |
| `{prefix}_lockfile.json` | every run | Machine-readable reproducibility record: the elected representative IDs, their parent accessions, the post-mutation config (NCBI credentials blanked, v0.42.0), tool versions, input-FASTA + HMM-DB sha256. Replayable via `repseq replay`. v0.30.0. |

"Segmented + GenBank" means: segmented mode is on **and** the GenBank source
features are reachable (either cached or fetched on demand) so the per-isolate
CDS list could be built — true by default; falls away under `--no-resolve` or
when `segmented.use_genbank_metadata: false` and no protein QC ran.

The detailed sections below tell you what's actually in each file, what
columns you'll see, and what each one is good for downstream.

### Sequence files (FASTA)

#### `{prefix}_representatives.fasta` — non-segmented mode

The main result for a non-segmented run: the selected representative
sequences in whatever alphabet you fed in (nucleotide or amino acid).
Header format matches the input format that was auto-detected (NCBI,
UniProt, or NCBI Virus). One record per representative.

Use it for: the new curated set, BLAST/DIAMOND database build, downstream
alignment, anything that needs the original sequences.

#### `{prefix}_concatenated.fasta` — segmented mode

One record per **representative isolate**, with each isolate's segments
joined head-to-tail in the canonical order you declared
(`segmented.viruses.<v>.segments`, e.g. `[L, M, S]`). The header is
`>CONCAT|<isolate_id>|<acc1>|<acc2>|<acc3>` and the body is the literal
concatenation of all segment sequences.

Use it for: whole-isolate phylogenetics on synonymous sites, gene-content
overviews, dotplots — anywhere you want each isolate as a single entity.

#### `{prefix}_segment_<NAME>.fasta` — segmented mode

One file per expected segment (e.g. `_segment_L.fasta`, `_segment_M.fasta`,
`_segment_S.fasta`). Each file holds **only the representative isolates'
copies of that one segment**, with the original GenBank/UniProt headers
preserved. The per-isolate accessions for the representatives are also
listed in `{prefix}_representative_isolates.tsv` under the `accessions`
column, in segment order.

Use it for: per-segment BLAST databases, building one tree per segment
(useful for reassortment analysis), feeding a per-segment HMM build.

#### `{prefix}_representative_isolate_proteins.fasta` *(segmented)* / `{prefix}_representative_sequence_proteins.fasta` *(non-segmented)*

Amino-acid FASTA covering **every annotated protein of every representative**.
Reconstructed from the same cached GenBank records the pipeline already
fetched, so no extra network calls are made.

Headers carry NCBI-style bracketed tags so each protein remains
self-describing in any downstream pipeline:

```
>NP_123456.1 [organism=Hantaan orthohantavirus] [ncbi_taxon_id=11599] \
   [species=Orthohantavirus hantanense] [genus=Orthohantavirus] \
   [family=Hantaviridae] [order=Elliovirales] [class=Bunyaviricetes] \
   [isolate=76-118] [segment=L] [host=Apodemus agrarius] \
   [country=South Korea] [collection_date=1976] [length=2151] \
   [parent=NC_005222.1]
```

Empty fields are omitted; bracket characters in values are scrubbed so
organism names like `Foo virus [strain X]` can't break the tag syntax.
When the GenBank `/serotype` qualifier is present (e.g. influenza A), a
`[subtype=H5N1]` tag is added after `[host=...]` — omitted otherwise.

Use it for: **the protein set you'll almost certainly hand to BLAST,
DIAMOND, HMMER, MMseqs2 search, or any sequence-search tool** — it's
both pre-curated and pre-annotated. This is usually "the" output file
for protein-centric workflows.

#### `{prefix}_per_protein_fasta/{prefix}_<family>.fasta` *(always-on since v0.22.0, one file per declared marker)*

Unaligned protein FASTA, one file per marker spec you declared under
`clustering.cluster_protein` (non-segmented) or `segment_markers` /
`cluster_protein` (segmented). Each file holds **one record per
representative that carries that marker** — exactly the CDS that
satisfied the marker's HMM (or the legacy alias chain if no `hmms:` is
declared on that spec).

`<family>` is the spec's `name:` when given (e.g. `Spike`), or the
single token (`Bunya_nucleocap`), or the first token plus `_altN` for
an unnamed multi-architecture spec. In segmented mode the segment is
prefixed (`M_Spike`, `S_Bunya_nucleocap`).

Headers are **byte-identical** to the all-protein FASTA above
(`_representative_isolate_proteins.fasta` / `_representative_sequence_proteins.fasta`),
so a Spike record in `_per_protein_fasta/<prefix>_M_Spike.fasta` and
in the all-protein file are the same string. That means you can pull
one marker, hand it to BLAST/HMMER/MAFFT, and the metadata travels
with it.

Use it for: per-marker downstream analysis (alignment, profile build,
ML tree off-pipeline) without having to re-extract the CDS from the
all-protein FASTA. This is independent of `--per-protein-phylo` —
it's written on every run so you always have the input FASTAs even if
you skip the tree step.

Specs that no representative satisfies are silently skipped (empty
file isn't written). A run where the HMM tier didn't fire and every
spec is HMM-only emits one stderr note explaining why nothing came
out.

#### `{prefix}_extra_protein_fasta/{prefix}_<name>.fasta` *(when `extra_protein:` declares entries, v0.22.0+)*

The **accessory-protein** analogue of the per-protein FASTA above.
Same selection chain, same bracket-tag header format. Driven by
`clustering.extra_protein:` (non-segmented) or `virus.extra_protein:`
(segmented, per-segment dict). See the
[Accessory proteins](#accessory-proteins-extra_protein) section for
what these are for and when you'd want one.

The filename uses the spec's `name:` verbatim (no segment prefix even
in segmented mode — `extra_protein` names must be unique across all
segments). Sparse coverage is the expected case here, so specs that no
representative satisfies are silently skipped without comment.

#### `{prefix}_representatives_protein.fasta` *(when `clustering.alphabet_for_clustering: protein` actually fired)*

The **AA strings that were fed into the clusterer** — the per-isolate
marker-protein concat in segmented mode, or the per-rep marker in
non-segmented mode (the multi-marker concat, e.g. Spike+Nucleocapsid, when
`clustering.concatenate_markers: true`). Written only if every
representative ended up with a populated `protein_sequence` (i.e. the
protein-alphabet path completed; not written for
`alphabet_for_clustering: nucleotide`).

This is a *diagnostic*, not a primary output: useful if you want to
reproduce or audit the clustering input. For a clean per-protein set, use
`_representative_*_proteins.fasta` (all proteins) or
`_per_protein_fasta/` (split by marker).

### Spreadsheets (TSV)

All TSVs use a harmonised vocabulary (since v0.8.0): the canonical
identifier column is always `accession`; booleans are always `TRUE` /
`FALSE`; length columns carry their alphabet (`length_nt`, `length_aa`,
`segment_length_nt`, `total_length_nt`); taxonomic ranks always appear in
the same nine-rank ladder (`species`, `subgenus`, `genus`, `subfamily`,
`family`, `suborder`, `order`, `subclass`, `class`). Sub-ranks come from
the NCBI lineage and are commonly blank for viruses.

#### `{prefix}_representative_sequences.tsv` — non-segmented mode

One row per representative sequence. **Schema-identical to
`_representative_isolates.tsv` (below)** since v0.21.0, so the same
analysis script reads both modes. Columns:

`isolate_id`, `isolate_id_source`, `organism`, `strain`, `host`,
`subtype`, `collection_date`, `country`, `n_segments`, `segments`,
`accessions`, `total_length_nt`, `is_refseq`, `is_reviewed`,
`ncbi_taxon_id`, then the nine-rank taxonomic ladder. (`subtype` is the
viral serotype, e.g. influenza A `H5N1`, from the GenBank `/serotype`
qualifier — blank when absent.)

In non-segmented mode the isolate-only cells are blanked
(`isolate_id`, `isolate_id_source`, `n_segments`, `segments`) or
remapped to their per-sequence meaning: `accessions` is the single
accession, `total_length_nt` is the sequence's NT length. The
per-sequence-only columns the old schema carried (`description`,
`segment`, `molecule_type`, `length_nt` under that name) are absent
— they have no slot in the shared schema. Use
`_sequence_proteins.tsv` for per-CDS detail.

Open in Excel / Numbers / your scripting language — this is the
spreadsheet you'll usually hand to a collaborator.

#### `{prefix}_representative_isolates.tsv` — segmented mode

One row per representative **isolate** (not per sequence — segmented
representatives are whole isolates, not single segments). Columns:

`isolate_id`, `isolate_id_source`, `organism`, `strain`, `host`,
`subtype`, `collection_date`, `country`, `n_segments`, `segments` (comma-joined
segment names in concat order), `accessions` (comma-joined per-segment
GenBank accessions in concat order), `total_length_nt`, `is_refseq`,
`is_reviewed`, `ncbi_taxon_id`, then the nine-rank taxonomic ladder.

`isolate_id_source` records where the `isolate_id` value came from:
`isolate` (GenBank `/isolate` qualifier — submitter-asserted unique per
biological sample), `strain` (`/isolate` was absent and `/strain` was
used as a fallback grouping key — see [Strain-as-isolate fallback](#strain-as-isolate-fallback-and-collision-detection)),
or `regex` (the header-regex fallback fired — typically `--no-resolve`,
UniProt input, or no accession). Blank for synthetic records without an
upstream provenance.

The per-sequence columns (`accession`, `segment`, `description`,
`molecule_type`, `length_nt`) are deliberately absent — they have no
isolate-level meaning. Use `_isolate_proteins.tsv` (or the per-segment
FASTA files) for per-segment / per-CDS detail.

#### `{prefix}_isolate_proteins.tsv` — segmented mode (when proteins are reachable)

One row per **CDS** of every isolate that survived QC, whether or not it
was picked as a representative. Columns:

`protein_id`, `product`, `length_aa`, `isolate_id`, `isolate_id_source`,
`segment`, `segment_length_nt`, `accession`, `representative` (`TRUE` if
the isolate made it into the final set, `FALSE` otherwise), then the
nine-rank taxonomic ladder.

`isolate_id_source` has the same meaning as in
`_representative_isolates.tsv` — see that section for the legend.

Use it for: a full per-protein audit of what made it through QC, and to
join back from a hit in a downstream analysis to the isolate it came from.

#### `{prefix}_representative_isolate_proteins.tsv` — segmented mode (when proteins are reachable)

Same exact schema as `_isolate_proteins.tsv`, but **filtered to the
representative isolates only** (i.e. every row has `representative=TRUE`).
Easier on the eye when you just want "the proteins in my reduced set".

#### `{prefix}_sequence_proteins.tsv` — non-segmented mode (when proteins are reachable, v0.21.0+)

The non-segmented counterpart of `_isolate_proteins.tsv`. One row per
**CDS** of every sequence that survived QC, with the exact same column
schema as the segmented version so the same analysis script reads both
modes:

`protein_id`, `product`, `length_aa`, `isolate_id`, `isolate_id_source`,
`segment`, `segment_length_nt`, `accession`, `representative` (`TRUE` if
the parent sequence ended up as a representative, `FALSE` otherwise),
`hmmscan`, then the nine-rank taxonomic ladder.

The isolate-only columns are blanked in non-segmented mode
(`isolate_id`, `isolate_id_source`, `segment`); `segment_length_nt`
holds the parent sequence's NT length and `accession` is the parent
accession.

Use it the same way as `_isolate_proteins.tsv`: a per-protein audit of
what passed QC, and the table to join against for downstream hits.

#### `{prefix}_representative_sequence_proteins.tsv` — non-segmented mode (when proteins are reachable, v0.21.0+)

Row-filtered companion to `_sequence_proteins.tsv` (reps only — every
row has `representative=TRUE`). Same column schema; the non-segmented
analogue of `_representative_isolate_proteins.tsv`.

#### `{prefix}_clusters.tsv` — always

One row per cluster (NOT per member — the row shows the cluster's
representative). Columns:

| Column | What it is |
| --- | --- |
| `cluster_id` | Stable cluster identifier. |
| `accession` | Accession of the representative. |
| `organism` | Representative's organism. |
| `cluster_size` | Total members in the cluster (representative + others). |
| `is_refseq` / `is_reviewed` | Representative's quality flags. |
| `species_purity` / `species_distinct` | **v0.26.0 diagnostic.** Fraction of populated species labels in the cluster that match the most common one (1.000 = pure) and the count of distinct species labels (1 = pure, >1 = mixed). |
| `genus_purity` / `genus_distinct` | Same at the genus rank. **`genus_distinct > 1` is the actionable "threshold too lax for this stratum" signal** — a single cluster crossing genus boundaries is almost always worth investigating. |
| `host_purity` / `host_distinct` | Same against the `/host` qualifier — flags cross-host contamination / mis-annotation. |
| `member_length_min` / `_max` / `_median` | NT-length distribution across the cluster's members. A wide range can flag a length-outlier inclusion (partial sequence dragged into a full-genome cluster). |

Purity is computed across populated values only — a member with no
resolved genus does **not** count against the cluster's genus purity.
Blank purity (with `_distinct = 0`) means no member carries a value at
that rank.

Use it for: tracing which sequences ended up in the same cluster, or for
re-running selection with a different priority (e.g. "I want the reviewed
UniProt entry not the longest one") without re-clustering. The purity
columns are your fastest tell that a stratum was clustered too aggressively.

#### `{prefix}_group_counts.tsv` — always

One row per stratum (whatever your mode grouped on: genus, host, decade,
country, custom field, …). Columns: `stratified_by`, `stratum`,
`stratum_size_before`, `stratum_size_after`, `clustered` (TRUE/FALSE —
did the binary-search clusterer run, or was the group already small
enough to keep whole?), `cutoff` (the identity threshold the clusterer
settled on, if it ran).

When `clustering.diversity_curve_cutoffs` is non-empty (default
`[0.99, 0.95, 0.9, 0.8, 0.7]`), every clustered stratum also carries
trailing `n_clusters_<c>` columns — the cluster count obtained by
re-running the backend at that fixed identity threshold. This is a
**diagnostic only — representative selection is not affected**. It
lets you read off how conserved a group is at a glance: a stratum that
collapses to 5 clusters at 0.99 vs one that needs 0.70 to collapse to
the same 5 are very different beasts. Cells are `NA` for cutoffs below
the backend's identity floor (cd-hit-est < 0.80, cd-hit protein < 0.40),
and empty for unclustered strata.

**The fastest way to see where the reduction happened.** If a single
genus ate most of your "selected" budget, you'll see it here.

#### `{prefix}_qc_removed.tsv` — always

Every sequence dropped during cleaning, with the reason. Two columns:
`accession` and `reason`. Typical reasons:

- `duplicate` — exact duplicate of another record
- `length:<n><[<>]<m>` — failed the whole-pool length filter (non-segmented mode)
- `segment_length:<seg>:<n><[<>]<m>` — failed the per-segment length bounds (segmented mode; whole isolate is dropped, so all of its segments appear with this reason)
- `ambiguous_fraction:<f>>0.05` — too many `N`/`X`/other ambiguous characters
- `annotation_keyword:<word>` — description contains a blacklisted keyword
- `protein_count:<n><[<>]<m>` — wrong number of annotated CDS (when protein QC is on)
- `segmented_filter:could_not_identify_isolate_or_segment` — the GenBank lookup and `isolate_regex` fallback both failed
- `incomplete_isolate:missing_segments:<list>` — an isolate that was missing one or more expected segments
- `taxonomy_mismatch:<rank>` — segments of one isolate disagreed on the configured taxonomy rank (default `species`), so the whole isolate was dropped (new in v0.9.0)

**Check this file first when more was dropped than you expected.** Sort
by `reason` to see which filter did the damage; the same information
shows up summarised in `_run.log`.

#### `{prefix}_overrides.tsv` — when a sequence was force-kept (v0.45.0+)

The flip side of `_qc_removed.tsv`: every sequence that the **force-keep
override list** (`overrides.protect_qc`) rescued from a QC removal stage.
Three columns — `id`, `stage` (which QC stage it bypassed), and
`would_be_reason` (exactly the reason it would have appeared with in
`_qc_removed.tsv` had it not been protected). Only written when at least
one sequence was actually protected, so a normal run leaves no empty
file behind.

See **[Sequences of special importance](#sequences-of-special-importance--overrides)**
below for how to configure the list.

#### `{prefix}_force_selected.tsv` — when a sequence was force-selected (v0.46.0+)

The selection-side audit for `overrides.force_select` / `--pin-ids`. One
row per pinned sequence — `id`, `action`, and `detail`. `action` is one of:
`elected_representative` (won its cluster), `split_singleton` (split out
because another pin won the same cluster; `detail` names the cluster),
`added_representative` (clustering had deselected it; added as a singleton),
`already_representative` (it was already a rep — nothing to do), or
`unavailable` (no surviving sequence matched the id — it didn't survive QC,
or was never in the input). Only written when at least one sequence was
pinned.

#### `{prefix}_hmm_diagnostic.tsv` — when the HMM tier ran (v0.26.0+)

One row per (representative, declared marker spec, individual HMM
profile). Emitted only when the HMM tier was active this session AND at
least one marker spec declared `hmms:`. Multidomain tokens like
`"CoV_S1--CoV_S2"` expand to **one row per component HMM**
(`CoV_S1`, `CoV_S2`) so you see per-profile granularity.

Columns:

| Column | What it is |
| --- | --- |
| `isolate_id` / `accession` / `segment` | Identity of the rep (segment is blank in non-segmented mode). |
| `spec_name` | The marker spec's `name:` (or its single token). |
| `profile` | One declared HMM profile from any token in the spec. |
| `hit` | `TRUE` if any CDS in scope produced a hit on this profile, else `FALSE`. |
| `best_dom_evalue` / `best_dom_score` / `best_coverage` | The lowest-E hit's stats (blank when `hit=FALSE`). |
| `ga_cutoff` | The profile's Pfam GA threshold from the .hmm headers (blank when the profile has no curated GA). |
| `cutoff_used` | `GA=<value>` when `use_ga_when_available` AND a curated GA exists, else `default_evalue=<value>`. |
| `passing` | `TRUE` if the best hit cleared BOTH the similarity gate AND the coverage gate; blank when no hit. |

This is the **diagnostic that surfaces silent failures the binary
HMM-QC gate hides**. Sort by `passing` to see borderline misses; group
by `profile` to see whether a particular HMM is too strict for your
dataset ("65 isolates barely missed RdRP_4 at E=2e-5, consider relaxing
`hmm.default_evalue`"). Rows with `hit=FALSE` tell you when a profile
is genuinely absent rather than just below the cutoff.

### Run record

#### `{prefix}_run.log`

Plain-text record of the run: version, date, exact command-line, every
config setting that was active, the per-step QC counts, the mode's
results, and the list of output files. Save this with your data — it
makes the selection fully reproducible.

The "QC SUMMARY" block also shows the new **per-segment length-filter
breakdown** (since v0.9.1) so you can tell *which* segment caused
isolates to fall out, not just that "some did".

#### `{prefix}_config_repseq<ver>.yaml` — sanitized effective-config snapshot (v0.42.0+)

A clean, shareable copy of the configuration the run actually used —
written on every run, no flag. Because repseq deep-merges your YAML over
its built-in `DEFAULTS`, this file shows **every setting at the value it
ran with**: anything you set, plus every default you didn't. It is
*sanitized* in three ways:

- **Comments stripped** — the dump is pure settings (no schema
  commentary), so it's easy to diff and re-use.
- **Credentials blanked** — `taxonomy.ncbi_email` and
  `taxonomy.ncbi_api_key` are set to `null`. The keys stay (so the
  structure is complete and the file is a valid config) but your secret
  is never written to disk. Supply your own when re-running, or via the
  `REPSEQ_NCBI_EMAIL` / `REPSEQ_NCBI_API_KEY` environment variables.
- **Runtime state dropped** — pipeline-internal keys (anything starting
  with `_`, e.g. computed HMM cutoffs) are removed; they're recomputed
  on every run and aren't configuration.

The filename embeds the repseq version with periods replaced by
underscores (`cov_config_repseq0_42_0.yaml`), so a directory of runs
self-documents which version produced each. The file is re-loadable
verbatim: `repseq <mode> -c cov_config_repseq0_42_0.yaml -i <input>`
reproduces the same settings. (The machine-readable
`{prefix}_lockfile.json` below also embeds the config — likewise with
credentials blanked since v0.42.0 — but as JSON, with the elected
representatives and tool fingerprints, for `repseq replay`.)

#### `{prefix}_lockfile.json` — structured reproducibility record (v0.30.0+)

A machine-readable companion to `_summary.md`: where the summary
records *what was done* in prose, the lockfile records *what came out*
in JSON — the elected representative IDs, their parent accessions,
the post-mutation config the pipeline actually ran under, the version
of every external tool that ran, an sha256 fingerprint of every input
FASTA, and (when the HMM tier ran) an sha256 of the HMM database.
Always written, no flag.

Top-level keys (`schema_version: 1`):

```jsonc
{
  "schema_version": 1,
  "repseq_version": "0.30.0",
  "created_utc": "2026-05-25T12:00:00Z",
  "mode": "taxonomic1",
  "command": "repseq taxonomic1 ...",
  "python_version": "3.11.6",
  "config": { ... },                       // post-mutation cfg
  "inputs": [{"path": "in.fasta", "sha256": "…", "size_bytes": 12345}],
  "tools": {"mafft": "v7.520", "iqtree2": "2.3.2", ...},
  "hmm_db": {"path": "...", "sha256": "...", "bundled": true, "n_profiles": 47},
  "representatives": [
    {"id": "AB1234", "kind": "sequence", "accession": "AB1234", "organism": "Virus X"},
    {"id": "CONCAT|iso1", "kind": "isolate",
     "isolate_id": "iso1", "organism": "Virus Y",
     "segment_accessions": {"L": "AB1", "M": "AB2", "S": "AB3"}}
  ]
}
```

Sequences are **not** embedded — they're re-fetched on demand by
`repseq replay`. Cluster membership is **not** recorded; the
authoritative source for that is `{prefix}_clusters.tsv` (recording it
twice would let the two go out of sync). The JSON is written with
sorted keys + 2-space indent, so two lockfiles from identical pipelines
diff cleanly on `created_utc` and nothing else.

#### `repseq replay <lockfile.json>` — re-materialise the representatives

Reads a lockfile and re-emits the representative FASTAs into a fresh
output directory, without re-running QC, clustering, or selection:

```bash
repseq replay run1/run1_lockfile.json -o replay_of_run1/
# Optional: re-build trees too (slower; bootstrap variance even at same seed)
repseq replay run1/run1_lockfile.json -o replay_of_run1/ --rebuild-trees
```

What it does:

1. Validate the lockfile's `schema_version` and `repseq_version` (major
   mismatch is a hard error; minor difference is a stderr note).
2. Verify the HMM database content sha256 if the original run used the
   HMM tier — a mismatch is a loud stderr warning so the user knows
   HMM-derived state may differ.
3. Re-fetch each accession from NCBI via the cache (new cache source
   `ncbi_nuc_seq`); for segmented reps, fetch all three segments and
   rebuild the CONCAT.
4. Re-emit `{prefix}_representatives.fasta` (non-segmented) or
   `{prefix}_concatenated.fasta` + `{prefix}_segment_<NAME>.fasta`
   (segmented).
5. If `--rebuild-trees`: run the phylogeny step against the
   re-fetched reps using the lockfile's `phylo:` config.
6. Write a short `{prefix}_replay.md` pointing at the original
   lockfile and listing what was emitted.

When NCBI can't return an accession (retired, network down), the replay
**warns loudly + continues with survivors** — better to deliver
partial output than crash. The dropped accessions land in
`{prefix}_replay_missing.tsv` (representative_id, accession) so the
discrepancy is auditable.

### Taxonomic reports

Three plain-text reports that turn "what kind of diversity did I just
select?" into a one-glance answer. All three are written on every run
that has the relevant inputs available; all are safe to open in any
text editor.

> **Provenance header (v0.41.0).** Each of the four
> `_*_taxonomic_report.txt` files below begins with a single
> `#`-commented line recording the run's repseq version, selection mode,
> dataset type (segmented/non-segmented), clustering substrate + tool,
> and representative count — so an isolated report file is
> self-describing without needing `_summary.md` alongside it. The TSV
> companions don't carry it (a comment line would break a tidy-data
> parser); the same facts live in their rows.

#### `{prefix}_taxonomic_report.txt` — diversity before vs after

Per-rank counts of **distinct taxa** in the pool fed to selection
versus in the final representatives. Counting unit is **isolates** in
segmented mode and **sequences** otherwise. Two sections:

- **Section 1** — a 9-rank ladder (species → class) with two numeric
  columns: distinct taxa before, distinct taxa after. A glance tells
  you which ranks the clustering compressed and which it preserved.
- **Section 2** — for each rank that has at least one populated taxon,
  the per-taxon breakdown: taxon name, before-count, after-count,
  sorted by before-count desc. Ranks with more than
  `output.protein_report.max_breakdown` (default 20) populated taxa
  show only the top 20, with a note in the rank label.

Blank rank values are excluded from every count, so subfamily /
suborder / subclass cells (commonly empty for viruses) don't inflate
diversity.

#### `{prefix}_protein_taxonomic_report.txt` — per-marker coverage + length (v0.22.0+)

Per-rank tables for **each declared protein** (every
`cluster_protein` / `segment_markers` / `extra_protein` spec). For
each rank from `subgenus` to `class` (skipping `species` because the
annotation noise there dominates), four sub-tables are written:

1. **Coverage (post-QC pool)** — cells are `<count> <pct>%`: how many
   isolates/sequences in that taxon had a CDS satisfying the marker.
2. **Coverage (representatives)** — same, restricted to the final
   selected set.
3. **Length statistics (post-QC pool)** — cells are
   `min, max, median, Q3-Q1, n` in amino acids (Q3-Q1 = IQR; n = the
   number of items contributing the length, so you can judge whether
   the quartiles are trustworthy).
4. **Length statistics (representatives)** — same on the reps.

Marker columns are headed by the spec's `name:`. **Cluster-driving
markers (the `cluster_protein` / `segment_markers` ones — i.e. the
proteins that drove both clustering and the whole-genome tree) carry a
trailing `*`**, so you can tell at a glance which are the load-bearing
proteins vs accessory `extra_protein` ones.

A `== HMM marker architectures ==` block at the bottom lists each
HMM-bearing spec's token alternatives (joined with `OR`) and the
cutoff policy that gated it (GA when curated, else
`default_evalue`, plus `relative_length_cutoff`) — so the Methods
section of a paper can quote the gate verbatim.

Taxa within each sub-table are **truncated first by member count**
(top `max_breakdown` by population, default 20) **then sorted
alphabetically** for display — so a long-tail rank never silently
drops a high-count taxon for being late in the alphabet.

Use it for: a one-shot answer to "is my representative set still
covering polymerase / glycoprotein / nucleocapsid across the
families I started with, and is the protein length distribution still
sensible?" — the bench-scientist QC that catches a marker quietly
dropping out of one genus.

#### `{prefix}_nucleotide_taxonomic_report.txt` — per-rank NT length (v0.23.0+)

The nucleotide analogue of `_protein_taxonomic_report.txt`: per-rank
**length-statistics** tables on the input nucleotide sequences
themselves, not on the encoded proteins. Two sub-tables per rank —
post-QC pool and representatives — with cells in the same
`min, max, median, Q3-Q1, n` format (Q3-Q1 = IQR; n = number of items
contributing the length).

Columns depend on the mode:

- **Segmented:** one column per **segment** in the configured order
  (S/M/L for hanta; HA/NA/NS/… for influenza), plus a trailing
  **`total`** column (sum of segment lengths per isolate, i.e. the
  concatenated genome length). One row per isolate. Lets you see e.g.
  whether the L segment shrank in your reps compared to the pool
  (a sign your clustering disproportionately kept short / fragmented
  records), or whether one genus has a systematically larger total
  genome than another.
- **Non-segmented:** a single **`genome`** column with NT lengths.
  Useful for sanity-checking that the reps don't skew toward the
  short or long end of the distribution within each taxon.

Same rank scope as the protein report (`subgenus` to `class`,
skipping `species`), the same top-20-by-member-count truncation knob
(`output.protein_report.max_breakdown`), and the same "blank ranks
are excluded" rule.

There are **no coverage tables** here — every passing entity carries
every required nucleotide unit by construction (segmented
completeness QC guarantees all `expected_segments` are present), so a
coverage column would always be 100% and add no signal. The
protein-coverage tables exist precisely because proteins are
variably present; segments and whole genomes are not.

Use it for: catching a representative set that's lost length
distribution against the pool — a frequent and easy-to-miss artifact
when clustering at an aggressive identity threshold or when one
sub-clade is over-represented by short submissions.

#### `{prefix}_polyprotein_taxonomic_report.txt` — per-peptide coverage + length (v0.35.0+)

The **sliced-peptide analogue** of
`_protein_taxonomic_report.txt` — only emitted when at least one
`clustering.polyprotein:` / `virus.polyprotein:` spec is declared and
the HMM tier is active (the same gate as the per-peptide FASTAs under
`{prefix}_polyprotein/`).

One H2 section per declared polyprotein spec
(`========== ORF1ab ==========`). Within each section, columns are the
spec's declared peptides in N→C order — and **every declared peptide
stays as a column even if no item carries it**: a 0% coverage row for
a peptide you declared is itself a useful audit signal (it tells you
that on this dataset the slicer never located that peptide on any
parent CDS).

For each rank from `subgenus` to `class` (same ladder as the protein
report, same `output.protein_report.max_breakdown` top-N truncation),
four sub-tables are written:

1. **Coverage (post-QC pool)** — cells are `<count> <pct>%`: how many
   isolates/sequences in that taxon had a peptide successfully sliced
   from the polyprotein parent CDS.
2. **Coverage (representatives)** — same, restricted to the final
   selected set.
3. **Peptide length statistics (post-QC pool)** — cells are
   `min, max, median, Q3-Q1, n` in amino acids on the sliced peptide
   sequence.
4. **Peptide length statistics (representatives)** — same on the reps.

A trailing `== Peptide architectures ==` block per spec lists each
peptide's HMM token alternatives (joined with `OR` when more than one)
and its `cleavage_motif:` when set, so the Methods section of a paper
can quote the locator for each peptide verbatim.

**What counts as "covered":** peptides with `status` of `ok` or
`overlap` in the per-spec audit TSV — the **same statuses that produce
records in the per-peptide FASTAs**. Peptides with status `missing`,
`out_of_order`, or `no_parent_cds` are not counted as covered; check
`{prefix}_polyprotein/{prefix}_<spec>_peptides.tsv` to see why an
individual representative fell into each bucket.

Use it for: spotting genus-specific peptide loss in a polyprotein
(e.g. coronavirus NSP1 silently dropping out of one
subgenus — non-homologous between alpha- and beta-CoV, easy to miss
in the per-isolate audit, immediately visible here as 0% in one row
and 100% in another).

#### `{prefix}_*_taxonomic_report.tsv` — tidy long-format companions (v0.35.0+)

Each of the four `_*_taxonomic_report.txt` files above (diversity,
protein, nucleotide, polyprotein) is accompanied by a **machine-readable
tidy long-format TSV** with the same basename and `.tsv` extension.
These are designed to drop straight into Excel pivot tables, R
(`readr::read_tsv` + `tidyr`), pandas, or any other analysis tool —
one row per observation, no multi-table layout, no header banners.

All four TSVs share one canonical 8-column schema, so they can be
concatenated and keyed on the `report` column:

| column | values |
| --- | --- |
| `report` | `diversity`, `protein`, `nucleotide`, `polyprotein` |
| `rank` | `species`, `subgenus`, `genus`, …, `class` (the 9 `_TAX_RANKS`) |
| `pool` | `post_qc` (the input pool fed to clustering) or `reps` (the elected representatives) |
| `taxon` | taxon name at that rank (e.g. `Alphacoronavirus`); `*ALL*` for rank-level diversity rows |
| `taxon_count` | number of items in that (rank, pool, taxon) cell — the denominator of coverage_pct |
| `spec` | marker / extra protein name (cluster-driving markers carry a trailing `*`, matching the `.txt`); `<polyprotein_name>:<peptide_name>` for polyprotein rows (e.g. `ORF1ab:NSP1`); `genome` / `<segment_name>` / `total` for nucleotide rows; `_diversity` for diversity rows |
| `metric` | one of: `distinct_taxa`, `member_count`, `coverage_count`, `coverage_pct`, `length_min`, `length_max`, `length_median`, `length_iqr`, `length_n` |
| `value` | numeric; integers as integers, percentages with up to 4 decimals (trailing zeros trimmed, so `100.0%` is written `100`) |

**Excel:** Data → From Text/CSV → tab delimiter → Insert →
PivotTable. Then drag `taxon` to rows, `metric` to columns, `value` to
the cell aggregator — you get any wide view the `.txt` showed and
more.

**R / pandas:**
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

**Zero-coverage handling.** When a marker / peptide / segment is
absent from a taxon, the row emits `coverage_count=0,
coverage_pct=0, length_n=0` and **no other length metrics** — we
won't fabricate `length_min`/`length_max`/`length_median` for an empty
distribution, since downstream aggregators would treat them as real
zeros and skew averages. Filter on `length_n > 0` before computing
across-taxon length summaries.

The architecture / banner sections of the `.txt` reports (HMM marker
architectures, peptide architectures, header preamble) are humans-only
and have no TSV equivalent — consult the matching `.txt` when you need
to know which HMM tokens defined a `spec`.

### Diagnostic plot

#### `{prefix}_clustering.png` — only if `--plot` is passed

Two-panel scatter of the clustered sequences (k-mer Jaccard distance):

- **Left** — every (sub)sampled sequence as a dot, coloured by genus.
- **Right** — the same dots, coloured by cluster, sized by cluster
  population, with faint lines from each member to its representative.

The embedding uses **UMAP** when `umap-learn` is installed and imports
cleanly, and otherwise falls back to a **numpy-only classical MDS**
(PCoA) on the same distance matrix; the figure labels its axes and title
with whichever method ran. Only `matplotlib` is required —
`pip install -e '.[viz]'` — so `--plot` works out of the box. For the
UMAP upgrade (nicer separation on large diverse sets) add it with
`pip install -e '.[viz-umap]'`; if that package fails to build/import in
your environment, `--plot` still produces the MDS plot.

Bounded by a default 2000-point subsample for very large runs
(representatives are always kept). Skipped when the run produced no
clusters (`global -n` mode).

### Phylogeny outputs

Written only when you pass `--phylo`. The full pipeline:
short-id remap → MAFFT (`--auto`) → *(optional)* trimAl trimming → tree
builder → root → label internal
nodes by LCA → phyloXML. Tree builder is auto-picked from your
clustering alphabet: **IQ-TREE for protein** (ModelFinder + UFBoot
bootstrap by default) and **FastTree for nucleotide** (with `-nt -gtr`).
Set `phylo.tool: fasttree` (or `iqtree`) to pin one.

Since **v0.27.0**, MAFFT / IQ-TREE / FastTree print one start line and
one finish line each (`[phylo] starting MAFFT (...)` / `[phylo] MAFFT
finished (Ns)`), so the terminal stays readable. Pass **`--verbose`**
to also stream the live child stderr (prefixed with `[mafft]` /
`[iqtree]` / `[fasttree]`) — useful for diagnosing a step that looks
stuck or for picking up an IQ-TREE ModelFinder banner mid-run. The
stderr is buffered either way, so the on-failure error message still
contains the full subprocess output regardless of `--verbose`.
(v0.26.0 introduced the streaming; v0.27.0 made it opt-in because the
firehose was too noisy for a successful run.)

##### Plain-text Newick is opt-in: `phylo.newick` / `--newick` (v0.48.0+)

Each tree is written as annotated **phyloXML** (`*_tree.xml`) and, when
you ask for it, as plain-text **Newick** (`*_tree.nwk`). The Newick is
**off by default** (`phylo.newick: false`) to keep output directories
tidy — the phyloXML carries the same topology plus taxonomy, metadata,
LCA labels and colours, so the `.nwk` is redundant for most workflows.

Set `phylo.newick: true` (or pass `--newick`) to retain the `*_tree.nwk`
for **every** tree of the run — whole-genome, per-protein, extra-protein,
per-segment, pre-cluster, and partitioned. `--no-newick` forces it off
even if the YAML turned it on. Turning it off costs nothing scientifically:
the Newick is still produced internally (the phyloXML is parsed from it,
and the per-protein incongruence RF table reads it) and only deleted at
the very end of the run; the `_tree_id_map.tsv` (short-id↔accession) and
`_msa.fasta` are kept regardless, so the retained alignment stays
decodable.

##### Graphical tree figures: `phylo.pdf` / `--pdf` / `--no-pdf` (v0.49.0+)

Every phyloXML tree is **also rendered as a graphical figure** — a
`*_tree.pdf` (vector) and a `*_tree.png` (150 dpi raster) sibling next to
each `*_tree.xml`. The figure is a ladderized rectangular phylogram drawn
with **matplotlib + Bio.Phylo**: taxonomy-coloured leaf labels, a
genus/subfamily colour legend, internal-node labels for the genus→family
ranks, and branch-support labels for nodes with support ≥ 50. Everything in
the figure (labels, colours, the legend, ranks, support) is reconstructed
**from the phyloXML itself**, so the figure colours match the tree
annotation exactly — and a single end-of-run sweep covers *every* tree type
(whole-genome, per-protein, extra-protein, per-segment, pre-cluster,
partition) with no per-tree configuration.

Rendering is **on by default** (`phylo.pdf: true`). Pass **`--no-pdf`** (or
set `phylo.pdf: false`) to suppress it — the annotated phyloXML is retained
either way and can be opened later in a dedicated viewer (e.g.
[Archaeopteryx](https://sites.google.com/view/archaeopteryx)). `--pdf`
forces it on even if the YAML turned it off. Rendering needs **matplotlib**
(the `[viz]` extra: `pip install 'repseq[viz]'`); if it is missing, the run
prints a single note and skips the figures — the same soft-fail posture as
`--plot`. (The approach is ported from the sister project
[`vfam_trees`](https://github.com/cmzmasek), which renders viral-family
reference trees the same way.)

##### Preliminary-run shortcut: `--fast` (v0.24.0+)

When you are iterating on a config (tuning QC thresholds, HMM
architectures, cluster cut-offs, etc.), the slow/strict tree pipeline
is overkill — every run rebuilds full IQ-TREE bootstraps over a fresh
MAFFT `--auto` alignment. Pass `--fast` and repseq overrides the
`phylo:` block of your YAML for that one run and forces:

- **FastTree** for every tree (whole-genome 2E **and** per-protein 2F).
  No IQ-TREE, no ModelFinder, no UFBoot.
- **MAFFT `--retree 1`** — single-pass FFT-NS-1 progressive alignment.
  Drops `--auto`, ignores `phylo.mafft.extra_args` /
  `phylo.per_protein.mafft.extra_args`. (Skips any L-INS-i you may
  have set for publication runs.)
- **No trimAl** (`phylo.trimal.enabled` and
  `phylo.per_protein.trimal.enabled` are forced off).
- **No partitioned supermatrix** (`phylo.partition.enabled` is forced
  off — FastTree can't fit per-partition models anyway).
- **Midpoint rooting** for every tree (`phylo.rooting.method` is forced
  to `midpoint`, overriding any configured method/outgroup). Skips the
  slower taxonomy-guided / MAD rooting chain, which would otherwise
  re-score internal-clade specificity against NCBI lineages on every
  tree.

The override is in-memory only — your YAML file is untouched. A
banner is printed to stderr at startup so you always know a fast-mode
run is happening, and `{prefix}_summary.md` describes the fast
pipeline (not whatever your YAML said). **Switch `--fast` off for
the final publication run.** trimAl trimming, IQ-TREE / ModelFinder,
UFBoot bootstrap, and the partitioned supermatrix all change the
topology meaningfully on real data; the fast pipeline trades that
accuracy for speed.

##### Alignment trimming (trimAl, optional, off by default)

Enable `phylo.trimal.enabled: true` to trim poorly-aligned / gap-rich
columns from the alignment **before** tree inference with
[trimAl](http://trimal.cgenomics.org/) — it sits between MAFFT and the
tree-builder. `mode` is the trimAl method (default `automated1`, its
heuristic best for ML trees; also `gappyout`/`strict`/`strictplus`/
`nogaps`/`noallgaps`), and `extra_args` is a raw passthrough for threshold
trimming (`["-gt", "0.8"]`). The per-protein trees have an independent
`phylo.per_protein.trimal` with the same shape. When trimming runs,
`{prefix}_msa.fasta` becomes the trimmed (tree-input) alignment and the
raw MAFFT output is kept as `{prefix}_msa_untrimmed.fasta`; in partitioned
mode each per-family alignment is trimmed before concatenation (so the
NEXUS charset ranges stay correct). If trimal is missing or strips the
alignment to nothing, repseq emits a loud warning and builds the tree on
the **untrimmed** alignment rather than failing — so you never lose the
tree over a trimming tool. `repseq doctor` reports trimal as an optional
binary. The trimal version + mode are recorded in the phyloXML
`<description>`.

##### Partitioned supermatrix (default for protein + IQ-TREE)

For protein runs that use IQ-TREE, the whole-genome tree is built as a
**partitioned supermatrix** rather than by gluing every marker into one
string. Each declared marker family (one per `hmms:` spec you already use
for QC — alternative architectures within a spec collapse to one family)
is aligned **separately** with MAFFT, the
per-family alignments are concatenated **column-wise** into a supermatrix,
and IQ-TREE fits a substitution model **per partition** (ModelFinder).
This is the statistically correct multi-marker analysis: MAFFT is never
asked to align an L-polymerase against an M-glycoprotein at a segment
seam, and a polymerase and a glycoprotein no longer share one model.

This is **on by default** (`phylo.partition.enabled: true`); it
applies only when the run resolves to protein + IQ-TREE **and** the HMM
tier resolved at least two marker families. Otherwise — FastTree, no HMM
tier, or a single family — it transparently falls back to the legacy
single-alignment, single-model path. Knobs under `phylo.partition`:

- `linkage` — IQ-TREE branch-length linkage across partitions:
  `proportional` (`-p`, edge-linked + per-partition rate; the default),
  `equal` (`-q`, shared branch lengths), or `unlinked` (`-Q`, independent
  branch lengths per partition — most flexible, useful when segments may
  have different histories).
- `models` — optional per-family model pin keyed by family label
  (`{L_RdRP_4: "LG+G4", S_Bunya_nucleocap: "WAG+G4"}`). This is **all-or-
  nothing**: pins are honoured only when *every* family is pinned (the
  NEXUS then carries a `charpartition` with those concrete models and
  IQ-TREE runs without `-m`). If any family is left unpinned, all
  partitions fall back to per-partition ModelFinder (a `charpartition`
  cannot mix concrete models with ModelFinder), with a warning. With no
  pins (the default), the NEXUS is charsets-only and IQ-TREE is run with
  `-m MFP` so ModelFinder picks a model per partition. (We deliberately do
  *not* write `MFP` into the `charpartition` — IQ-TREE 2.x treats it as a
  model-file path and aborts with `File not found MFP`.)

Extra outputs in this mode: the per-family alignments
`{prefix}_msa_<family>.fasta` and the NEXUS partition file
`{prefix}_partition.nex` (the column ranges, plus per-partition models
when pinned). `{prefix}_msa.fasta` is then the **concatenated
supermatrix** the tree was actually inferred on.

#### `{prefix}_msa.fasta`

The MAFFT alignment (the concatenated supermatrix in partitioned mode).
Headers are `>S0001 <pretty-label>` — the short
`SNNNN` id stays the first whitespace-separated token (safe for any
phylo tool), and the descriptive label (built from `phylo.labeling.format`
/ `segmented_format`) is appended as the FASTA description so AliView /
Jalview / MEGA show recognisable names without confusing the underlying
tools.

#### `{prefix}_tree.nwk` — **off by default**

Tree-builder Newick output. Leaf names are the **short ids** (`S0001`,
`S0002`, …) so the file works with any downstream tool that expects
short, safe leaf labels. Use the id map below to recover real names.

**Not retained by default** (`phylo.newick: false`) — the phyloXML below
is a topological superset, so the plain-text Newick is redundant for most
workflows and was removed to cut output clutter. Set `phylo.newick: true`
(or pass `--newick`) to keep the `*_tree.nwk` for *every* tree of the run
(whole-genome, per-protein, extra-protein, per-segment, pre-cluster,
partition); `--no-newick` forces it off regardless of config. The Newick
is still generated internally during the run — the phyloXML is parsed from
it and the per-protein incongruence table reads it — so turning it off
costs nothing but the files; the `_tree_id_map.tsv` and `_msa.fasta` are
kept either way.

#### `{prefix}_tree.xml` — **the tree you'll usually open**

phyloXML with rich, browseable annotation:

- Every leaf gets a formatted `<name>`, a `<taxonomy>` block with NCBI
  taxon id, a `<sequence>` block with the GenBank accession + title, and
  one repseq-namespaced `<property>` per per-leaf attribute (host,
  subtype, collection_date, country, strain, isolate_id, year, plus the
  full 9-rank taxonomy ladder: species, subgenus, genus, subfamily, family,
  suborder, order, subclass, class — empties dropped).
- Leaf labels are **coloured by taxonomy** (default on, by genus): each
  leaf carries a `<property ref="style:font_color" applies_to="node">`
  hex colour that Archaeopteryx renders on the label. Configure under
  `phylo.coloring` — one rank (`ranks: [genus]`) gives each value a
  distinct hue; two ranks (`ranks: [genus, subgenus]`) shade each
  subgenus within its genus's colour. Unresolved taxa show grey. The
  same palette is shared across `--phylo` and every `--per-protein-phylo`
  tree, so a genus is the same colour everywhere — which is what makes
  cross-tree incongruence (reassortment) jump out visually. Set
  `phylo.coloring.enabled: false` to turn it off.
  - `saturation` and `value` (both 0–1) are the **HSV** vividness and
    brightness held constant across the palette — only the hue changes
    per taxon. `saturation` is how vivid vs. washed-out (`0` = grey,
    `1` = fully vivid); `value` is how bright vs. dark (`0` = black,
    `1` = brightest). The defaults (`0.65` / `0.90`) give strong but
    not garish colours that stay legible on a white background. Lower
    `saturation` for pastels; lower `value` for darker, muted tones.
  - `missing_color` is the `#RRGGBB` grey used for unresolved taxa.
- Every annotated internal clade gets a `<name>` and a `<taxonomy>`
  block holding the LCA's scientific name + rank (`min_rank=genus` by
  default; the labeller keeps each monophyletic clade labelled at its
  crown and suppresses obvious duplications like a 2-leaf same-species
  pair).
- The `<phylogeny>` element carries a `<name>` and `<description>`. The
  `<description>` now **leads with a plain-English "what this tree is
  based on" sentence** (v0.41.0) — e.g. *"This tree is based on a
  per-isolate concatenation of one marker protein per segment (L:RdRp,
  M:GPC, S:N); amino-acid; each leaf is one representative isolate."* —
  followed by the MAFFT/IQ-TREE/FastTree versions, the selected
  substitution model, the bootstrap settings, and the rooting method that
  actually fired. The same substrate is also emitted as machine-readable
  phylogeny-level `<property>` elements (`repseq:tree_basis`,
  `repseq:analysis_mode` = segmented/non_segmented, `repseq:substrate`,
  `repseq:alphabet`, `repseq:leaf_unit` = isolate/sequence), so a script
  can read what a tree is built on without parsing prose. Every tree
  type carries its own basis: the whole-genome tree, each per-marker
  tree (*"the Spike marker protein only … HMM architecture …"*), the
  per-segment NT trees, the polyprotein-peptide trees, the
  accessory-protein trees, and the pre-cluster overview tree.
- Confidence values are normalised to 0–100 integers (`sh_like` for
  FastTree, `ufboot` for IQ-TREE). The tree is ladderized.

Opens directly in [Archaeopteryx](https://sites.google.com/view/aptxjs)
([web version](https://sites.google.com/view/aptxjs))
([source](https://github.com/cmzmasek/forester)) which makes use of the
the rich annotation.
The annotation is also visible in any other phyloXML-aware tool (Dendroscope, etc).

#### `{prefix}_tree.pdf`, `{prefix}_tree.png` — **the figure you'll usually paste** (v0.49.0+)

A ready-to-use graphical rendering of the tree above, so you don't have to
open a viewer just to glance at the topology. The PDF is vector (clean at any
zoom, drops straight into a figure panel); the PNG is a 150-dpi raster for
quick previews and slide decks. Both are a **ladderized rectangular
phylogram** drawn with matplotlib + Bio.Phylo, showing:

- **Leaf labels coloured by taxonomy** (the same `style:font_color` the
  phyloXML carries — so the figure matches what Archaeopteryx would show),
  with a **genus/subfamily colour legend** when colours resolve. Unresolved
  taxa (e.g. an offline `--no-resolve` run) fall back to grey, no legend.
- **Internal-node labels** for the genus→family ranks (species-level and
  unranked internals are suppressed to cut clutter).
- **Branch-support labels** for nodes with support ≥ 50.
- A title + a one-line provenance caption pulled from the phyloXML.

Rendered by a single end-of-run sweep over every `*_tree.xml`, so the same
figure is produced for *every* tree type (whole-genome, per-protein,
extra-protein, per-segment, pre-cluster, partition). **On by default**
(`phylo.pdf`); see **"Graphical tree figures"** under
[Phylogeny outputs](#phylogeny-outputs) for the `--no-pdf` switch and the
matplotlib requirement.

#### `{prefix}_tree_id_map.tsv`

Two-column mapping `short_id` ↔ `accession`. Use it when you need to
trace a leaf in `_tree.nwk` or the MSA back to a real sequence id.

#### `{prefix}_iqtree_summary.txt` — only when IQ-TREE ran

The IQ-TREE ModelFinder report — which substitution model was selected
and why, plus the log-likelihood and bootstrap settings. Worth a glance
before quoting a model in a methods section.

#### `{prefix}_iqtree_model.txt` — only when IQ-TREE ran (v0.28.0)

The same ModelFinder pick(s), distilled to one `<label>: <model>` line
per partition so you don't have to grep through the long `.iqtree`
report. Non-partitioned runs emit a single `GENOME: <model>` line (e.g.
`GENOME: LG+I+G4`); partitioned-supermatrix runs emit one line per
charset, in partition order (e.g. `CoV_S1: LG+I+G4`, `CoV_M: JTT+G4`).
The same picks also land in the phyloXML `<description>`'s `model=`
field (so a phyloXML viewer can show them) and in the `_summary.md`
Methods section (so the auto-generated prose says which model fit
rather than the generic "ModelFinder for substitution-model selection").

> **Fail-soft:** if MAFFT, IQ-TREE, or FastTree are missing, or fewer
> than 3 representatives survived, the phylogeny step is skipped with a
> stderr message and the rest of the run's outputs are still written.
> IQ-TREE additionally refuses UFBoot with fewer than 4 sequences; the
> wrapper drops bootstrap automatically in that case but still produces
> the tree.

### Per-protein trees — `--per-protein-phylo`

`--phylo` builds **one** tree from the whole representatives (the marker
concat in segmented mode). `--per-protein-phylo` instead builds **one
tree per declared HMM marker** — one tree per `hmms:` spec you set for QC
under `segment_markers` / `cluster_protein`. A spec's `hmms:` list holds
**alternative domain architectures** (OR; see [Token
notation](#token-notation)), so a marker such as
coronavirus Spike — `["CoV_S1--CoV_S2", "bCoV_S1_N--bCoV_S1_RBD--CoV_S2"]`
— builds a single tree spanning every representative whose Spike matches
*either* architecture (alpha- and beta-CoV together). repseq picks the CDS
that satisfies any of the spec's tokens on each carrying representative and
runs the *same* MAFFT → IQ-TREE/FastTree → root → LCA pipeline on those
protein translations. The alignments use MAFFT
`--auto` by default (fast). For a high-accuracy publication run, set
`phylo.per_protein.mafft.extra_args: ["--maxiterate", "1000",
"--localpair"]` (L-INS-i) — affordable on these small single-gene
alignments but noticeably slower, so it's opt-in; a non-empty list is
passed to MAFFT without `--auto` so the strategy takes effect.

Why you'd want it: comparing the single-marker trees side by side reveals
**topological incongruence** — an L-segment polymerase tree disagreeing
with an M-segment glycoprotein tree is the classic signature of
**reassortment**, which the concatenated `--phylo` tree hides.

Outputs land in a `{prefix}_per_protein/` subdirectory, one set per built
family. `<family>` is the spec's `name:` when given (e.g. `Spike`), else
its single token (e.g. `Bunya_nucleocap`), else the first token + `_altN`
for an unnamed multi-architecture spec — prefixed with the segment in
segmented mode (`M_Spike`, `S_Bunya_nucleocap`):

```
{prefix}_per_protein/
  Spike_msa.fasta
  Spike_tree.nwk
  Spike_tree.xml
  Spike_tree.pdf
  Spike_tree.png
  Spike_tree_id_map.tsv
  S_Bunya_nucleocap_msa.fasta
  …
  {prefix}_incongruence.tsv
```

(The `*_tree.nwk` files shown above are written only when `phylo.newick`
/ `--newick` is on — off by default; see **"Plain-text Newick is opt-in"**
under [Phylogeny outputs](#phylogeny-outputs). The `*_tree.pdf` / `*_tree.png`
figures are written by default — see **"Graphical tree figures"** in the same
section; pass `--no-pdf` to skip them. The phyloXML, MSA, id-map, and
incongruence table are always written.)

**Extra-protein trees (v0.22.0+):** when you've also declared
[`extra_protein:`](#accessory-proteins-extra_protein) accessory
proteins, a **separate** `{prefix}_extra_protein/` directory is emitted
with the same per-tree shape (one `<name>_msa.fasta` / `<name>_tree.nwk`
/ `<name>_tree.xml` / `<name>_tree_id_map.tsv` per spec, built by the
same engine). The split is intentional: accessory proteins are sparse
by design, so they are **excluded from the `_incongruence.tsv` table**
to keep the reassortment signal in the required-marker pairs from being
drowned out by `NA`/low-`n_common_taxa` rows. If you want RF distances
involving an accessory tree, run the math off its `.nwk` file
yourself.

Each tree file has the same format and rich annotation as its `--phylo`
counterpart, with two deliberate differences. First, a leaf shows **only
the CDS that fed that tree** as `<sequence type="protein">` (the
`CoV_nucleocap` tree shows just the nucleocapsid protein, not every CDS
of the genome). The `<sequence type="dna">` element encoding it, and the
`repseq:nuc_acc` / `repseq:protein_acc` / `repseq:protein_names` summary
properties, still describe the leaf's full gene content. Second, that
protein carries its **HMM domain architecture** as a phyloXML
`<domain_architecture>` — one `<domain>` per hit, with the protein
coordinates and the hit's E-value as `confidence`:

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

[Archaeopteryx](https://sites.google.com/view/aptxjs) draws these as
domain boxes on the leaf and lets you filter them with its interactive
E-value slider — so *all* hits are emitted (not just the ones that
cleared repseq's cutoffs), giving the slider its full range. Turn the
block off with `phylo.per_protein.domain_architecture: false`. The flag
runs alone or alongside `--phylo`.

#### `{prefix}_incongruence.tsv` — incongruence as a number

So you don't have to eyeball it, repseq scores the **pairwise unrooted
Robinson-Foulds (RF) distance** between every pair of marker trees (and,
when `--phylo` also ran, the whole-genome tree as a `GENOME` row):

| tree_a | tree_b | rf | norm_rf | n_common_taxa |
|--------|--------|----|---------|---------------|
| Spike  | N      | 4  | 0.3333  | 18            |
| Spike  | GENOME | 0  | 0.0000  | 21            |
| N      | GENOME | 4  | 0.2857  | 18            |

- **`rf`** — number of bipartitions that differ between the two trees,
  scored on their **shared** taxa only (marker trees have different leaf
  sets). `0` = identical unrooted topology; higher = more disagreement.
  Rooting is ignored, so a marker that merely *roots* differently isn't
  counted as incongruent.
- **`norm_rf`** — `rf` divided by the maximum possible RF for that many
  shared taxa (`2·(n−3)`), so values are comparable across pairs with
  different overlap. `NA` when fewer than 4 taxa are shared.
- **`n_common_taxa`** — how many representatives the pair has in common
  (the basis for that row's score).

A `GENOME` column tells you which marker departs most from the consensus
history. Turn the table off with `phylo.per_protein.incongruence: false`.

> **Requirements / fail-soft:** needs the HMM tier (`hmm.enabled` plus
> configured `hmms:`) to have run, and `mafft` + a tree builder on PATH.
> A family carried by fewer than `phylo.per_protein.min_taxa` (default 3)
> representatives is skipped with a log note; if no family qualifies — or
> the HMM tier didn't run — the whole step is skipped with one stderr
> line, leaving the rest of the run intact.

### Per-segment NT trees — `--per-segment-phylo` (v0.26.0+, segmented only)

A complement to `--per-protein-phylo`. Where 2F builds one **protein**
tree per declared marker (one CDS per segment), `--per-segment-phylo`
builds one **nucleotide** tree per declared segment — every
representative isolate contributes its raw per-segment NT (from
`concat_segments`) to one tree per segment. Outputs land in a
`{prefix}_per_segment/` subdirectory with the same per-tree shape
(`<segment>_msa.fasta`, `<segment>_tree.xml`,
`<segment>_tree_id_map.tsv` per segment, plus a `<segment>_tree.pdf` /
`<segment>_tree.png` figure by default; the `<segment>_tree.nwk` only
when `phylo.newick` / `--newick` is on).

Why a per-**segment** tree on top of the per-**marker** tree set:
reassortment can show up at the segment level without appearing in any
single marker tree (a marker is one CDS within a segment, so per-marker
incongruence is a lower-resolution view). The per-segment NT trees
catch reassortment events the per-marker AA trees can miss, especially
when a marker is conserved across the parent species. Same MAFFT /
tree-builder / rooting / LCA / colouring as 2E and 2F; same `min_taxa`
floor; same shared colour palette so a taxon stays the same colour
across every tree of the run. Runs alone or alongside `--phylo` and
`--per-protein-phylo`. Segmented-mode only; raises a stderr note on a
non-segmented run.

#### Outgroup rooting (v0.26.0+)

Set `phylo.rooting.method: outgroup` and name the outgroup by accession
or by taxon — the most common request from bench scientists who already
know which lineage is basal:

```yaml
phylo:
  rooting:
    method: outgroup
    # Single accession:
    outgroup: "AB123456"
    # …or a list (the tree is rooted at the MRCA of all matching leaves):
    # outgroup: ["AB123456", "CD789012"]
    # …or a whole clade by taxonomy rank:
    # outgroup_rank: {family: "Hantaviridae"}
```

Accession matching is case-insensitive against each rep's `accession`,
`id`, or `isolate_id` (so a CONCAT isolate id works too). A multi-leaf
spec roots at the MRCA of all matching leaves. `outgroup` and
`outgroup_rank` are unioned — you can mix "this specific reference
accession plus every member of family X". When the spec matches no
representative the rooter falls through to midpoint with a stderr
note, so a typo never voids the tree. The summary's
`{prefix}_summary.md` names the outgroup target honestly (it no longer
falsely claims "MAD/midpoint fallbacks" for non-`auto` methods).

The default remains `phylo.rooting.method: auto` (taxonomy → MAD →
midpoint chain); the other accepted values are `taxonomy`, `mad`,
`midpoint`, and `none` (no rooting — for an input tree already rooted
upstream).

### Pre-cluster overview tree — `--pre-cluster-tree` (v0.32.0+)

A diagnostic / sanity-check tree over **every post-QC sequence**,
built BEFORE clustering would have collapsed redundancy. The point
is to see at a glance where the elected representatives land in the
broader diversity of the input pool: the phyloXML `<name>` of each
representative leaf is prefixed with `[repr] ` so it stands out in
Archaeopteryx, FigTree, or any phyloXML viewer.

Pipeline is intentionally **hard-coded for speed**, regardless of
what the rest of `phylo:` looks like:

- **MAFFT `--retree 1`** (FFT-NS-1, no `--auto`, no L-INS-i)
- **FastTree** (no IQ-TREE / ModelFinder / UFBoot) — regardless of
  alphabet; the model selection IQ-TREE does isn't worth the budget
  on a many-thousand-leaf rough tree
- **Midpoint rooting** only (no taxonomy-guided / MAD)
- **No LCA annotation, no trimAl, no bootstrap**

The alphabet mirrors `clustering.alphabet_for_clustering`: if the
run clustered on protein and every sequence carries a populated
`protein_sequence`, the pre-cluster tree is built on AA; otherwise
it falls back to NT. This keeps the overview tree's distance space
consistent with the selection it depicts. In segmented mode the
leaves are CONCAT isolates (same unit as the post-cluster tree).

Output files:

- `{prefix}_pre_cluster_tree.nwk` — Newick, short-id leaves
  (`S0001`, `S0002`, …). **Opt-in** — written only when `phylo.newick`
  / `--newick` is on (off by default, like every other `*_tree.nwk`)
- `{prefix}_pre_cluster_tree.xml` — phyloXML with the full
  taxonomy `<property>` enrichment + per-leaf taxonomy colour
  (same palette as `--phylo` / `--per-protein-phylo`, so taxa stay
  the same colour across every tree of a run), and the `[repr] `
  prefix on representative leaves
- `{prefix}_pre_cluster_tree.pdf` / `.png` — graphical figure of the
  overview tree, rendered by default (`phylo.pdf`; `--no-pdf` to skip)
- `{prefix}_pre_cluster_tree_id_map.tsv` — three columns
  (`short_id`, `accession`, `is_rep`), so a user grepping for
  reps without opening the XML can find them

The MSA produced by MAFFT is written to the run's temp directory
and deleted after FastTree consumes it — only the tree + id_map
land in the output dir. This is by design (the lock-in answer to
"output file set?"); if you want to keep the MSA, run `--phylo`
instead.

Trigger either by passing `--pre-cluster-tree` on the command line
or by setting `phylo.pre_cluster_tree.enabled: true` in the YAML
config. Soft-fails like the other phylo steps: missing
MAFFT/FastTree, fewer than 3 post-QC sequences, or any subprocess
error surfaces as a single stderr line and the rest of the run
continues.

> **Scale note:** runtime grows roughly linearly with the post-QC
> pool size. MAFFT `--retree 1` + FastTree on 5000 leaves typically
> finishes in a few minutes; tens of thousands of leaves are
> achievable but the resulting tree may be visually unreadable
> without subsampling.

### One conservation number per MSA — `{prefix}_msa_conservation.tsv` (v0.43.0+)

repseq reduces **every alignment a run produces** to a single
conservation score and collects them all into one table. Whenever any
phylogeny step writes an MSA — the
whole-genome tree, the partitioned per-family alignments and supermatrix,
the per-protein marker trees, and any extra-protein, polyprotein-peptide,
or per-segment trees — a post-hoc sweep scores each `*_msa*.fasta` under
the output directory and writes `{prefix}_msa_conservation.tsv`, one row
per alignment. It runs by default (no flag); set
`phylo.conservation.enabled: false` to turn it off. Being a decoupled
sweep, a scoring problem never affects a tree, and if no tree was built
the file simply isn't written.

**The metric.** Mean per-column **Jensen-Shannon divergence to a residue
background** (Capra & Singh, *Bioinformatics* 2007), with two standard
corrections:

- **Henikoff & Henikoff (1994) position-based sequence weighting**, so a
  clade of near-identical sequences plus one outlier doesn't read as
  artificially conserved (the redundant sequences share the weight one
  real sequence would carry).
- a **gap penalty**: each column's JSD is multiplied by `(1 −
  weighted_gap_fraction)`, so a mostly-gap column contributes little.

Each column score is bounded `[0, 1]` (JSD with log base 2 and the
symmetric mixture `m = (p+q)/2`). Interpretation:

- `~0` — columns look like the background, i.e. unrelated sequences;
- `~0.6–0.8` — typical for a real protein family of close homologs;
- a **perfectly conserved** column tops out at ~0.85–0.95 for protein
  (the exact ceiling depends on the conserved residue's background
  frequency — a rarer residue lands closer to 1) and ~0.55 for
  nucleotide (only four symbols, flat background). JSD-to-background
  **never reaches exactly 1**, which is why published JSD conservation
  values don't saturate. Scores are therefore comparable only *within*
  an alphabet — don't compare a protein MSA's number against a
  nucleotide MSA's; the `alphabet` column flags which is which.

**Why a perfectly conserved column isn't 1.** The score measures
*divergence from a background distribution*, not raw column uniformity. An
invariant column is a point mass on one residue; its Jensen-Shannon
divergence from the background is large but finite, and its exact value is
set by how rare that residue is in the background (rarer ⇒ more surprising
⇒ closer to 1). With only four symbols and a flat background, nucleotide
columns have a lower ceiling still (~0.55). This is a deliberate property
of the Capra & Singh formulation — it rewards columns that are *both*
invariant *and* enriched for an otherwise-rare residue — and it's why
published JSD conservation values never saturate at 1. If you want a
strict 0→1 "fraction identical" instead, that's a different (and less
informative) statistic; this column reports conservation in the
divergence-from-background sense.

**Columns:** `msa` (path relative to the output dir), `role` (`genome` /
`partition_family` / `marker` / `extra_protein` / `peptide` /
`segment_nt`), `label` (family / segment / peptide name), `alphabet`
(`protein` / `nucleotide`, auto-detected), `trimmed` (`TRUE` only when
trimAl actually trimmed *this* alignment — i.e. an `_untrimmed`
companion exists beside it; `FALSE` for an untrimmed tree input, which
is the default since trimAl is off by default, and `FALSE` for the
retained `_untrimmed` companion row itself), `n_seqs`, `n_sites`,
`mean_conservation` (the headline all-column mean), and
`mean_conservation_core` (mean over well-occupied — ≥ 50% non-gap —
columns, ignoring the gap penalty; a view of the well-aligned core that
ignores ragged ends).

**References.**

- **Capra JA & Singh M (2007).** Predicting functionally important
  residues from sequence conservation. *Bioinformatics* 23(15):
  1875–1882. doi:[10.1093/bioinformatics/btm270](https://doi.org/10.1093/bioinformatics/btm270)
  — the Jensen-Shannon-divergence conservation score this metric
  implements (including the gap-penalty term and the
  divergence-from-background ceiling behaviour above).
- **Henikoff S & Henikoff JG (1994).** Position-based sequence weights.
  *J Mol Biol* 243(4): 574–578.
  doi:[10.1016/0022-2836(94)90032-9](https://doi.org/10.1016/0022-2836(94)90032-9)
  — the position-based sequence weighting used to down-weight redundant
  near-identical sequences.
- **Valdar WSJ (2002).** Scoring residue conservation. *Proteins*
  48(2): 227–241. doi:[10.1002/prot.10146](https://doi.org/10.1002/prot.10146)
  — a thorough survey of conservation-scoring schemes and the
  chemistry-aware (substitution-matrix-weighted) alternative to JSD.
  repseq reports the JSD score rather than Valdar's, but this is the
  reference for the broader method space if you need a substitution-
  similarity-aware measure.

### Per-taxon monophyly — `{prefix}_monophyly.tsv` (v0.54.0+)

For **every tree a run builds** (whole-genome, per-protein/extra/segment/
polyprotein, pre-cluster), repseq classifies how each taxon at each rank sits
on that rooted tree and collects the results into one table. Like the
conservation and incongruence reports it's a decoupled, always-on sweep — it
reads the annotated `*_tree.xml` files, so it needs no flag and simply writes
nothing if no tree was built.

**Status (per taxon, per rank, per tree):**

- **monophyletic** — the taxon's members form exactly one clade (their MRCA
  subtends only members).
- **paraphyletic** — non-monophyletic, and the foreign leaves inside the
  members' MRCA span form a **single** excluded clade (the taxon is an
  ancestral grade with one derived group carved out).
- **polyphyletic** — non-monophyletic, and the foreign leaves form **two or
  more** separate clades (the taxon is interrupted in several places).

Each rank is judged **only among leaves classified at that rank** — a leaf
with no genus is neither a member nor an intruder when genus monophyly is
assessed, so an annotation gap never masquerades as non-monophyly. Taxa with
fewer than two classified leaves are omitted (a singleton has no meaningful
status).

**Why it's useful.** A taxon that's monophyletic on the whole-genome tree but
**not** on a marker tree is the per-marker **reassortment / recombination**
signal — this is the taxon-resolved companion to `_incongruence.tsv` (which
gives whole-tree-vs-tree RF distances). It also operationalises the
non-monophyly that the tree's internal labels now show (a split taxon is
labelled on each of its clades).

**The para/poly tag is a documented heuristic.** Monophyletic-vs-not is an
exact MRCA test, but the para/poly split is definition-dependent; repseq keys
it on `intruder_clusters` (1 → paraphyletic, ≥2 → polyphyletic), the standard
topology-only convention. A taxon with several nested carve-outs can be
paraphyletic under the strict cladistic definition yet flagged polyphyletic
here, so the unambiguous counts are reported alongside the tag for you to
judge.

**Support-aware by default.** A single weakly-supported branch can make a
taxon look non-monophyletic when the tree doesn't actually resolve it that way.
So the assessment is **support-aware**: branches with support below
`phylo.monophyly.min_support` (**default 70**) are collapsed into polytomies
first, and a taxon is called monophyletic when no *well-supported* node
contradicts it (it could be a clade once the weak branches are resolved). Only
**confident** non-monophyly is flagged — a taxon broken solely by weak branches
reads as monophyletic. Set `min_support: 0` to trust every branch
(topology-only). On a real coronavirus tree this correctly turns an
apparently-paraphyletic *Alphacoronavirus* (split by a weak branch) back into a
clean monophyletic genus. The threshold used is recorded in the report.

**Columns:** `tree` (path relative to the output dir), `rank`, `taxon`,
`n_leaves` (classified members), `status`, `n_clusters` (separate in-group
blocks), `n_intruders` (foreign classified leaves in the MRCA span),
`intruder_clusters` (separate foreign clades — the para/poly key),
`intruder_taxa` (the intruding taxon names, `;`-joined), and `min_support`
(the support threshold applied).

### Plain-English flags — `{prefix}_flags.txt` (v0.55.0+)

The interesting anomalies a run can find are spread across several tables
(`_monophyly.tsv`, `_incongruence.tsv`, and — when the opt-in review ran —
`_taxonomy_review.tsv`). `{prefix}_flags.txt` distils them into one skimmable,
plain-English list so you don't have to join four TSVs by hand:

- **Reassortment / recombination signals** — taxa monophyletic on the
  whole-genome tree but broken on a marker tree, and marker-tree pairs with a
  high normalised Robinson-Foulds distance (≥ 0.2).
- **Taxonomy ↔ tree conflicts** — taxa that are para- or polyphyletic on the
  whole-genome tree.
- **Per-isolate taxonomy conflicts** — isolates whose label the tree
  neighbourhood disagrees with (from the taxonomy review).

It's a pure synthesis (no new computation), written whenever any source table
exists; a clean run gets a `No flags raised ✓` file so you can tell "nothing
flagged" from "report not produced". The flags are heuristics — the source
tables remain authoritative.

### One-page HTML report — `{prefix}_report.html` (v0.55.0+)

A single self-contained HTML file you can open in any browser or e-mail to a
collaborator. It bundles the run's **analysis flags**, a **gallery of every
tree figure** (the `*_tree.png` images embedded directly so the visuals travel
with the file), and an **index of all output files** (with sizes and relative
links to the full machine-readable detail). Pure-Python, no extra dependency;
written last so its file index captures the whole run.

---

## Cleaning (QC) — what gets dropped and how to control it

All cleaning settings live under `qc:` in your config file. Defaults are
sensible; loosen them if you're losing sequences you want to keep. Every
drop is logged in `{prefix}_qc_removed.tsv` with a precise reason — that
file is your first stop when QC ate more than expected.

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

- **`median_percent` compares every sequence to the median length of the
  whole file.** That's perfect for a single gene, but **wrong for a mixed
  file** (e.g. several different genes, or a whole genome plus its
  individual genes) — the median is meaningless and you'll drop things
  unfairly. For mixed files, use `min_max` with explicit numbers instead.
- **In segmented-virus mode, the whole-file length filter is skipped
  automatically** — a file of influenza segments mixes 2,300-nt and 890-nt
  sequences, so a single median can't work. Use per-segment length bounds
  instead (`segmented.viruses.<v>.segment_lengths`, see below). Since
  v0.9.1 the QC summary breaks per-segment drops down as
  `L too short : N`, `M too long : N`, etc., so you can see exactly which
  segment cutoff was the bottleneck.
- **Segmented mode adds a taxonomy-consistency check (since v0.9.0).** Any
  isolate whose segments resolve to *different* taxa at the configured
  rank (default `species`) is dropped — usually a reassortant, a
  contaminated record, or two unrelated viruses sharing an `/isolate`
  qualifier. Missing labels are ignored (only *populated disagreement*
  counts), so isolates with sparse taxonomy survive. Configure under
  `segmented.taxonomy_consistency.{enabled, rank}` and set
  `enabled: false` to skip it entirely.
- **Segmented mode now also flags strain-as-isolate collisions
  (since v0.9.3).** When the GenBank pre-pass falls back to `/strain` for
  records without `/isolate`, two distinct accessions can end up sharing
  the same grouping key. The new detector warns on `(isolate_id, segment)`
  pairs that collide; set `segmented.strain_collision_action: drop` to
  remove the colliders instead of just warning. See
  [Strain-as-isolate fallback](#strain-as-isolate-fallback-and-collision-detection).

---

## Sequences of special importance — `overrides`

Sometimes a handful of records *must* be handled specially: a type strain,
a vaccine reference, a genome your reviewers will look for by name. The
`overrides` block names them once (an id list and/or a file) and applies
**two independent capabilities** over that one list:

- **`protect_qc`** — **force-keep**: the named sequences bypass the QC
  removal stages they would otherwise fail (a whitelist), without
  loosening thresholds for everything else.
- **`force_select`** — **force-select**: the named sequences are
  guaranteed to appear among the representatives, regardless of how
  clustering would have collapsed them.

```yaml
overrides:
  ids: ["NC_045512.2", "MN908947"]   # accessions (or isolate_ids, segmented)
  ids_file: vip.txt                  # OR a file, one id per line; unioned with `ids`
  protect_qc: true                   # survive QC
  protect_stages: all                # or a subset (see below)
  force_select: true                 # appear in the representatives
```

Or skip the config and pass lists on the command line — handy for a one-off
run. `--protect-ids FILE` enables force-keep, `--pin-ids FILE` enables
force-select; both union into the one id list:

```bash
# keep these through QC AND guarantee they're in the output
repseq global -c my_config.yaml -i seqs.fasta -T 0.9 \
    --protect-ids vip.txt --pin-ids vip.txt
```

Each CLI flag reads one accession/isolate-id per line (`#` comments and
blank lines ignored), **unions** with any `overrides.ids` /
`overrides.ids_file` in the config, and turns its respective flag **on**
for that run.

- **Matching** is by `accession` (non-segmented) or `isolate_id`
  (segmented), **case- and version-insensitive** — `NC_045512` matches
  `NC_045512.2` and vice versa. In segmented mode you usually list
  `isolate_id`s; note the *early* QC stages (`ambiguous`, `annotation`)
  run before `isolate_id` is populated, so those only match on accession.
- **`protect_stages`** is either `all` or a list naming any of:
  `duplicates`, `length`, `ambiguous`, `annotation`, `protein_count`,
  `taxonomy_consistency`, `protein_quality`, `hmm`. This is per-stage on
  purpose: you can keep a slightly noisy reference through the `hmm`
  identity gate while *still* dropping it if its translation is genuine
  garbage (`protein_quality`).
- **Segmented completeness is deliberately *not* protectable.** Protection
  can't synthesize a missing segment, so a protected isolate that is
  missing a segment still drops (it can't be concatenated or treed) — with
  a distinct message rather than a silent fake-keep.

### How force-select decides (the hybrid policy)

When `force_select` is on, after clustering each pinned sequence is
guaranteed into the representatives by the rule that disturbs the rest of
the selection least:

- **Wins its cluster** — a pinned member beats the configured
  refseq/reviewed/longest priority and becomes that cluster's
  representative (the displaced representative drops back into the
  cluster's members).
- **Splits on collision** — if two or more pinned members land in the
  *same* cluster, the best one (same priority) wins the slot and the
  others are split out into their own singleton clusters, so all survive.
- **Added as a singleton** — a pin that clustering deselected entirely
  (the diversity-only `global -n` path, where non-selected sequences
  belong to no cluster) is appended as its own singleton representative.
- **Unavailable** — a pinned id that no surviving sequence matched (it
  was dropped in QC, or was never in the input) can't be selected; it's
  reported as `unavailable` with a warning.

> **Force-keep needs to pair with force-select for a hard guarantee.**
> `force_select` can only promote a sequence that *survived QC*. If a
> pinned record would be dropped by QC, enable **both** (`protect_qc` +
> `force_select`, or `--protect-ids` + `--pin-ids` on the same list) so it
> survives QC *and* lands in the output.

- **Nothing is silent.** Force-kept records (and the reason each would
  have been removed) go to `{prefix}_overrides.tsv`; force-select actions
  (elected / split / added / already-a-rep / unavailable) go to
  `{prefix}_force_selected.tsv`; both are summarised in
  `{prefix}_summary.md`. A reviewer can see exactly which records bypassed
  which filter and which were pinned into the output.
- **Count reconciliation.** Force-select runs *after* clustering, so the
  per-stratum counts in `{prefix}_group_counts.tsv` reflect clustering
  **before** any pins were added — they can sum to fewer than the final
  representative count. The pinned additions show up in the representative
  files, `{prefix}_force_selected.tsv`, and the final counts in
  `{prefix}_summary.md`.

---

## Segmented viruses (influenza, etc.)

Segmented viruses store their genome in several separate pieces. NCBI gives you
one FASTA record per segment, so a single isolate is spread across multiple
records. `repseq` can stitch them back together:

1. **Group records by isolate** — by default, `repseq` reads the `/isolate`
   (or `/strain`) qualifier from each record's GenBank source feature.
   The strain-name regex below is the fallback for sequences that don't have
   an NCBI accession (e.g. UniProt input) or whose record lacks the qualifier.
2. **Identify each record's segment** — by default from the GenBank `/segment`
   qualifier, or by its name / number / a synonym you define
   (e.g. `hemagglutinin` → `HA`).
3. **Keep only complete isolates** — an isolate missing any expected segment is
   dropped.
4. **(Optional) length-check each segment** — drop an isolate if, say, its HA
   segment is suspiciously short.
5. **Concatenate** the segments of each complete isolate into one sequence, so the
   normal grouping/selection can run on whole isolates.

The GenBank lookup reuses the same cache as the protein-annotation QC, so if
you have both turned on, repseq only fetches each record once. To turn the
lookup off entirely (and rely on the regex below for *every* record), set
`segmented.use_genbank_metadata: false`.

Configure it under `segmented:` and turn it on with `--segmented` on the command
line (or `enabled: true` in the config):

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

> The `isolate_regex` is a *fallback* — by default (`use_genbank_metadata: true`),
> repseq prefers the GenBank source-feature `/isolate` qualifier. The regex still
> needs to match for sequences without an NCBI accession or whose record has no
> isolate qualifier. If isolate identification fails for those records, it has to
> match the strain identifier as it appears in *your* headers, and it must capture
> it either as a group named `isolate` or as the first parenthesised group.

### Strain-as-isolate fallback and collision detection

GenBank source features have two distinct qualifiers: `/isolate` (a
specific collection event, usually unique to a single biological sample)
and `/strain` (a named variant, often shared across many samples). When
a record has `/strain="L99"` but no `/isolate=`, repseq uses the strain
value as the isolate grouping key — otherwise strain-only records (very
common in older bunyaviridae / hantaviridae submissions) would drop out
of segmented mode entirely.

The provenance of every `isolate_id` is written to the new
`isolate_id_source` column in `_representative_isolates.tsv` and
`_isolate_proteins.tsv`:

- `isolate` — value came from the `/isolate` qualifier (preferred).
- `strain` — value came from `/strain` as a fallback (`/isolate` was absent).
- `regex` — value came from the header-regex fallback (UniProt input,
  `--no-resolve`, or no accession).

**The risk of `strain`-derived ids is over-merging.** Two submitters
can deposit independent samples under the same strain name; repseq
would group them as one isolate and dedup-drop the duplicate segment.
After the GenBank pre-pass, repseq runs a **strain-collision detector**
that flags any `(isolate_id, segment)` pair where two or more distinct
accessions share a strain-derived id. You'll see a stderr block like:

```
  Strain-collision check: found 2 (isolate, segment) pair(s) where /strain is shared across distinct accessions.
    isolate 'L99' segment 'S': 3 accessions (ACC1, ACC2, ACC3)
    isolate 'M77' segment 'L': 2 accessions (ACC4, ACC5)
    Action: warn (no records dropped). Set segmented.strain_collision_action: drop to remove them.
```

The default action is `warn` (informational only — the pipeline keeps
the longest per-segment accession and dedup-drops the rest downstream).
Switch to `drop` to remove every accession involved in any collision
before the completeness filter runs:

```yaml
segmented:
  strain_collision_action: drop   # warn (default) | drop
```

Dropped records appear in `_qc_removed.tsv` with reason
`strain_collision:<segment>` and are counted in `QCReport.removed_strain_collisions`.

---

## Optional: protein-annotation QC

This step asks NCBI how many protein-coding genes (CDS features) each record has,
and drops records with too few — or, for segmented viruses, the wrong number per
segment. It's off by default; turn it on in the config:

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

Records are fetched from NCBI in batches and cached locally, so a second run on
the same data needs no network. Skipped automatically under `--no-resolve`.

---

## Phylogeny-based taxonomy review (`phylo.taxonomy_review`, v0.39.0)

Many viral records carry **missing** taxonomy (a blank genus, "unclassified")
or a **wrong/stale** label that pre-dates an ICTV reclassification. Once
`--phylo` has built a well-supported, LCA-labelled tree, the tree itself is
evidence: a representative nested deep inside a strongly-supported,
taxonomically *pure* clade is very likely the same taxon as its neighbours.

Enable it (off by default — it's a new inference step):

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
anchor) and then:

- **imputes** a *blank* rank from the clade's majority value, or
- **flags** a *populated* rank that disagrees — a suggestion only, **never**
  auto-changed.

Ranks are evaluated coarse→fine and imputations are kept **hierarchy-consistent**
(a subgenus is only imputed from neighbours that agree on genus; a blank genus
is back-filled when the agreeing neighbours share one). Confidence is tiered
(high/medium) from support × purity × neighbour count × anchor.

**Outputs.** `{prefix}_taxonomy_review.tsv` lists every verdict (impute +
conflict) with evidence columns — this is the authoritative **per-cell
imputation ledger**. With `write_corrected: true`, the **high-confidence
imputations** are additionally filled into
`{prefix}_representative_*_corrected.tsv` and the corresponding corrected
protein FASTA (clean values; the originals are kept untouched, and the review
TSV records exactly which cells were filled). Conflicts are *never* written
into the corrected copies.

**Scope & cautions (v1).** Representatives only (the tree's leaves; a
non-representative inherits nothing here). The tree topology, colouring, and LCA
labels are **not** modified — no feedback loop. The walk uses the same rooted
tree shown in the phyloXML; the support + purity + anchor gates are what guard
against rooting bias and recombination artifacts. Species is intentionally
excluded (viral species monophyly is too often violated). Soft-fails (a stderr
note, no files) if `--phylo` didn't build a tree.

---

## Optional: clustering plot

Pass `--plot` (and install with `'.[viz]'`) to get a two-panel scatter plot,
`{prefix}_clustering.png`, that lets you eyeball whether the clustering looks
sensible:

- **Left** — every sequence as a dot, positioned so similar sequences sit close
  together, coloured by genus.
- **Right** — the same dots, coloured by cluster, with bigger dots for bigger
  clusters and faint lines from each sequence to its chosen representative.

For big datasets the plot is drawn from a subsample (the representatives are
always included). It's skipped for `global -n` runs, which produce no clusters.

---

## The config file

`config/default_config.yaml` is fully commented and documents every option — read
that file as the reference. The `repseq init-config` wizard writes a starter
config for you. The most-changed settings:

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
  # `clustering.cdhit:` — see default_config.yaml for the full block.

representative:
  # Ordered preference for picking each cluster's "best" sequence: the first
  # criterion wins, later ones break ties, sequence length is the final
  # tiebreaker. The order is honoured — [reviewed_uniprot, refseq, longest]
  # prefers a Swiss-Prot-reviewed entry over a RefSeq one. Drop a criterion
  # to deactivate that preference. (`longest` is required.)
  priority: [refseq, reviewed_uniprot, longest]
```

### Unknown keys are rejected

repseq validates your config before doing any work and **aborts with a
`[config error]`** if it finds a key it doesn't recognise — a typo
(`ambiguous_threshhold`), a renamed/removed key, or a setting placed in the
wrong section (for example a `phylo.per_protein` knob written onto a
`clustering.polyprotein` spec). The message names the offending key, its full
dotted path, and — when there's a close match — the key you probably meant.
This is deliberate: a silently-ignored setting (the run quietly using the
default for something you thought you'd changed) is worse than a loud stop.

If you want to keep your own annotation/reminder key in the file, prefix it
with an underscore (`_my_note: ...`) and repseq will leave it alone. To check a
config without launching a run, use `repseq <mode> --dry-run -c your.yaml`.

### About `clustering.alphabet_for_clustering`

**This setting *only* changes what the clustering backend sees.** It does
**not** disable GenBank CDS download, protein-count QC
(`qc.protein_annotation.enabled`), or the per-segment
`virus.expected_proteins_per_segment` check — those run on every isolate
regardless of which alphabet you cluster on. Pick this purely on what kind
of identity threshold makes sense for your data.

| Value | Clustering input | When to pick it |
|---|---|---|
| `protein` *(default)* | Non-segmented: the per-sequence marker protein (longest CDS by default, overridable via `cluster_protein` aliases). Segmented: per-isolate concatenation of each segment's marker. Triggers a one-shot GenBank CDS fetch if proteins aren't already cached. | Diverged viral families: synonymous substitutions inflate NT divergence by 30–40% with no biological signal, and protein homology stays reliable down to ~25–30% identity vs ~50–60% for nucleotide. |
| `nucleotide` | Non-segmented: the input FASTA sequence as-is. Segmented: concatenation of all segments in `segments` order. | Tight species-level reference sets where genome-identity targets are what you want; non-coding input; or fully offline runs (`--no-resolve`), which require this value. |

`auto` was removed in v0.10.0 — pick `protein` or `nucleotide` explicitly.

#### Clustering on more than one marker (`clustering.concatenate_markers`)

By default a non-segmented protein-alphabet run clusters on a **single**
marker — the CDS from the first `cluster_protein` spec a record satisfies
(e.g. Spike alone for coronaviruses, since Spike is listed first). Set

```yaml
clustering:
  concatenate_markers: true
```

to instead select the marker CDS from **every** `cluster_protein` spec and
cluster on their amino-acid **concatenation** in declared spec order
(e.g. `Spike+Nucleocapsid`). This mirrors the segmented per-segment concat.
A record missing any required marker is dropped (when every spec declares
`hmms:`, the AND-gate already removes those upstream). `extra_protein`
specs are never concatenated — they are accessory/sparse by design.

This flag changes **only the clustering input** (and the
`{prefix}_representatives_protein.fasta` it writes). The whole-genome tree
is *not* affected. For a multi-gene whole-genome tree, leave the trees to
the **partitioned supermatrix** path by setting `phylo.tool: auto` (or
`iqtree`): each marker is aligned separately and concatenated column-wise
with a per-gene model — the correct way to build a Spike+N tree. Gluing the
proteins into one string and aligning the chimera (which is what you'd get
under FastTree) is biologically dubious and not recommended.

Override at run time without editing the YAML:

```bash
repseq global -c my.yaml -i x.fasta -T 0.95 --alphabet-for-clustering nucleotide
repseq global -c my.yaml -i x.fasta -T 0.95 --concatenate-markers      # or --no-concatenate-markers
```

`--concatenate-markers` / `--no-concatenate-markers` overrides
`clustering.concatenate_markers` for a single run (see above).

You can also set your NCBI email/key via the environment variables
`REPSEQ_NCBI_EMAIL` and `REPSEQ_NCBI_API_KEY` instead of putting them in the file.

---

## HMM-based identity QC

Annotation of viral CDSes in GenBank is famously inconsistent
("RNA-dependent RNA polymerase" vs "polymerase" vs "L protein" vs
"RNA polymerase L") and pseudogenes / misannotations can carry plausible-
looking `/product` strings. Since v0.14.0, repseq's primary use of
**HMMER hmmscan** is as a **QC step** — it verifies that each
segment/sequence actually carries the expected proteins by structural
identity rather than by name, and drops anything that fails the gate.
Marker selection for protein-alphabet clustering is a downstream
consumer of the same hit cache.

The HMM tier is **off-by-default at the marker level** — it activates
only for markers that carry an `hmms:` list — so you can phase HMMs in
one marker at a time without touching the rest of your pipeline.

**Runs regardless of clustering alphabet.** If you cluster on
nucleotides, HMM-QC still fires and drops isolates whose segments
don't carry the expected proteins (this is a change from v0.13.0,
which only ran HMMs on alphabet=protein runs).

### Token notation

Each element of `hmms:` is a **token** string. A token is either:

| Form | Meaning |
| --- | --- |
| `"Name"`         | Single HMM. A CDS satisfies the token if it has a passing hit to `Name`. |
| `"A--B"`         | Multidomain. A CDS satisfies the token only when it has passing hits to **both** `A` and `B`, with `A` lying N-terminal to `B` (forward-progressing endpoints, with at most `hmm.multidomain_overlap_tolerance` aa of overlap at the seam — default 30 aa). |
| `"A--B--C"`      | Same idea, three domains: `A` most N-terminal, `C` most C-terminal. |

**N-to-C order**. The first HMM in a multidomain token is the most
**N-terminal** domain on the protein — the same direction molecular
biology writes a protein sequence, so a token mirrors the domain
architecture as you'd draw it (e.g. the coronavirus Spike is
`"CoV_S1--CoV_S2"`: S1 N-terminal, S2 C-terminal). The named
domains are compared against the hmmscan `ali_from`/`ali_to` columns
in that order.

**Why the overlap tolerance.** Pfam profile boundaries rarely align
exactly to a real cleavage site. The classic case is the coronavirus
Spike S1/S2 furin seam, where the two Pfam profiles can overlap by
20-25 aa even though biology cuts a single peptide bond; the
bunyaviridae G1/G2 glycoproteins behave the same way. A strict
non-overlap rule would silently drop those representatives. The
default 30-aa tolerance accepts that fuzz; set
`hmm.multidomain_overlap_tolerance: 0` for strict non-overlap.
**Full containment is always rejected** regardless of the tolerance —
both `ali_from` and `ali_to` must progress forward between consecutive
named domains, so a single profile hit fully inside another can never
satisfy the multidomain token.

**Extra domains are fine.** A CDS annotated as `A--B--HMMX` still
satisfies the token `"A--B"` because `A` remains N-terminal to `B`.

**Alternative architectures (OR).** Multiple tokens in one marker's
`hmms:` list are *alternatives* — a CDS that satisfies **any one** of
them satisfies the marker. This lets one marker span divergent forms of
the same protein. For example, to verify Spike across alpha- and
beta-coronaviruses in one run:

```yaml
- {name: Spike, aliases: ["spike", "surface glycoprotein"],
   hmms: ["CoV_S1--CoV_S2", "bCoV_S1_N--bCoV_S1_RBD--CoV_S2"]}
```

A sequence passes the Spike marker when its Spike CDS matches *either*
the `CoV_S1--CoV_S2` architecture *or* the
`bCoV_S1_N--bCoV_S1_RBD--CoV_S2` architecture. (Across *different*
markers — separate dict entries — the rule is AND: each is an
independent marker that must be present. The OR applies only within one
marker's token list.)

### How the gate works

For a marker spec that includes `hmms: [<token>, <token>, ...]`:

1. Every CDS in every input sequence is scanned against the configured
   HMM database (once per sequence, batched into a single hmmscan call
   per run, cached on AA sequence so re-runs are free).
2. A hit **passes** when both gates clear:
   - **Similarity:** the curated Pfam GA threshold (when the profile
     carries one) OR `hmm.default_evalue` (default `1e-5`).
   - **Coverage:** the alignment span covers at least
     `hmm.relative_length_cutoff` (default `0.5`) of the HMM model
     length.
3. A **segment / sequence passes a marker** when at least one CDS
   satisfies **any one** of that marker's tokens. (Tokens in the same
   `hmms:` list are OR — alternative architectures; satisfying one is
   enough. Across *separate* marker specs the rule is AND — each is an
   independent marker that must be present.)
4. If a marker has no satisfying CDS (none of its alternative tokens
   matched), the segment fails. In segmented mode the entire isolate is
   then dropped (one bad segment fails the whole isolate, counted under
   `removed_hmm_failed` with a per-marker breakdown in
   `removed_hmm_by_marker`; the reason names the unmatched alternatives
   joined with `|`). In non-segmented mode the single sequence is
   dropped.
5. For protein-alphabet clustering, the marker CDS of each surviving
   segment / sequence is then the **longest** CDS that satisfies any
   token in the spec.

Markers WITHOUT an `hmms:` list keep the legacy alias → longest chain
unchanged and are not gated by HMM-QC.

### Bundled vs user-supplied database

`hmm.database: null` (default) uses the bundled set at
`repseq/data/hmms/repseq_viral_core.hmm` — 26 hand-picked Pfam-A
profiles covering the most common viral marker proteins (RdRp,
nucleocapsid, glycoprotein, helicase, protease, plus the coronavirus
Spike S1/S2 architectures and the small accessory proteins M/E/3a/7a/
viroporin) across the main families. Pfam-A is licensed CC0 so redistribution is unrestricted.
The bundled set was assembled from Pfam-A via
`scripts/build_bundled_hmms.sh`; you can re-run that script to refresh
it against a newer Pfam release.

To use your own profiles (e.g. VOGdb, a custom curated set, or the
full Pfam-A.hmm), set `hmm.database` to one of:

- **a single `.hmm` file** — the classic user-supplied database;
- **a directory of `.hmm` files** — every `*.hmm` in it is concatenated
  into one combined database (cached under `cache_dir`, rebuilt only when
  a member file is added/removed/changed). Handy for keeping
  family-specific profiles as separate files instead of `cat`-ing them
  yourself;
- **a bare bundled-set name** — a subdirectory under
  `repseq/data/hmms/`, e.g. `hmm.database: Filoviridae` selects the
  bundled `repseq/data/hmms/Filoviridae/` set (combined as above). Drop
  `*.hmm` files into such a directory to add a family-specific bundled
  set; see that directory's `README.md` for the convention.

Whatever resolves is auto-`hmmpress`-ed on first use if its `.h3*`
indexes are missing or stale.

`repseq doctor` reports the DB status (path, profile count, indexed
Y/N, GA-cutoff coverage).

> **HMM hits are cached.** Each CDS is scanned against the HMM
> database once and the result is kept in the local cache, keyed on
> the CDS amino-acid sequence **and** a fingerprint of the database
> file. Swapping in a different DB (or rebuilding the bundled one)
> invalidates old entries automatically — no manual step needed.
> Editing cutoffs (`default_evalue`, `relative_length_cutoff`,
> `multidomain_overlap_tolerance`, …) likewise needs no clear; cutoffs
> are reapplied each run. The one case where a manual clear helps is
> editing a `.hmm` file in place without changing its size or mtime,
> or wanting to reclaim disk after a DB swap:
> `repseq cache clear -c your_config.yaml --source hmmscan` wipes only
> the HMM hits and leaves the NCBI / UniProt caches intact. See
> [The local cache](#the-local-cache) for the full breakdown.

### Configuration examples

Non-segmented marker spec (top-level `clustering.cluster_protein`).
Each entry is either an alias string (legacy, alias-only) or a dict
with `name`, optional `aliases`, optional `hmms`:

```yaml
clustering:
  cluster_protein:
    - "polymerase"                                                          # legacy: alias-only
    - {name: "RdRp",  aliases: ["polymerase"], hmms: ["RdRP_4"]}            # single-HMM token
    - {name: "Spike", aliases: ["spike"],      hmms: ["CoV_S1--CoV_S2"]}    # multidomain token
```

Segmented marker spec (per-virus `segment_markers` — the recommended
HMM-aware form). Two flavours of the M-segment spec, illustrating the
difference between multidomain and separate-token notation:

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
        # Bunya_G2 hit, possibly on different CDSes (handles the
        # post-cleavage Gn/Gc annotation case).
        M: {hmms: ["Bunya_G1", "Bunya_G2"]}
        # OR — strict polyprotein form: a single CDS must carry both
        # domains in N-to-C order (Bunya_G2/Gn N-terminal, Bunya_G1/Gc
        # C-terminal). Use this when you want to reject post-cleavage
        # annotations as inconsistent.
        # M: {hmms: ["Bunya_G2--Bunya_G1"]}
        L: {hmms: ["RdRP_4"]}
```

When both `segment_markers` and the legacy `cluster_protein` define a
marker for the same segment, `segment_markers` wins.

### HMM-related config keys

```yaml
hmm:
  enabled: true                       # master switch
  database: null                      # null → bundled core; .hmm file, dir of .hmm files, or bundled-set name (e.g. Filoviridae)
  default_evalue: 1.0e-5              # used when a profile has no curated GA
  use_ga_when_available: true         # prefer curated GA over default_evalue
  relative_length_cutoff: 0.5         # ali_span / hmm_len ≥ this
  multidomain_overlap_tolerance: 30   # max aa overlap at the seam of a multidomain
                                      # token (0 = strict non-overlap; default 30
                                      # accommodates Pfam-boundary fuzz at e.g.
                                      # the CoV S1/S2 furin seam)
  threads: null                       # null → falls through to cfg.threads
```

### Soft-fail posture

The HMM tier never aborts a run on its own. If `hmmscan` is missing,
the database is unreadable, or hmmscan returns non-zero, repseq emits
one stderr warning and falls through to the alias / longest-CDS
chain. To check that everything is wired up, run `repseq doctor`.

---

## Accessory proteins (`extra_protein`)

`cluster_protein` (non-segmented) and `segment_markers` (segmented)
declare the **required** marker proteins: every isolate has to carry
one, the whole-genome tree is built from them, and missing one fails
HMM-QC. That's the right model for an RNA-dependent RNA polymerase or
a structural glycoprotein — proteins present in every member of the
family.

It's the **wrong** model for sparse accessory proteins. The
coronavirus ORF7a, ORF8, ORF9b, or the betacoronavirus-only HE
(haemagglutinin-esterase) appear in some isolates and not others.
Adding them to `cluster_protein` would force HMM-QC to drop every
representative that doesn't carry them — exactly the wrong outcome.

`extra_protein:` (since v0.22.0) is the home for those proteins. The
spec shape is identical to `cluster_protein` — `{name, aliases?,
hmms?}` — but extras are **not** used for clustering, are **not**
required for HMM-QC to pass, and are **not** part of the whole-genome
tree. They are reported wherever it makes sense to:

- Per-marker protein FASTA in `{prefix}_extra_protein_fasta/`
  (always-on when any extra is declared).
- A separate per-protein tree in `{prefix}_extra_protein/` when you
  also pass `--per-protein-phylo` (built with the same MAFFT +
  IQ-TREE / FastTree engine as the required-marker trees).
- Coverage and length columns in `{prefix}_protein_taxonomic_report.txt`,
  without the trailing `*` that marks the cluster-driving markers.
- **Excluded** from `{prefix}_per_protein/{prefix}_incongruence.tsv` by
  design — sparse coverage would produce a lot of `NA`/low-`n_common_taxa`
  rows that would drown out the real reassortment signal across required
  markers. If you want RF distances against an extra, run the math
  off the `<name>_tree.nwk` file yourself.

Configuration shape (segmented, one entry per segment as a dict; the
non-segmented form is a top-level `clustering.extra_protein: [...]`
list using the same dict shape):

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

Names must be unique across all `extra_protein` entries for a virus
(the name is used as the filename and also has to identify the protein
unambiguously in the protein-taxonomic report).

---

## Polyprotein cutting (`polyprotein`) — v0.33.0+

Some virus families express their CDS as one giant precursor protein
that is then post-translationally cleaved into mature peptides:
picornavirus P1/P2/P3 (VP4/VP2/VP3/VP1, 2A/2B/2C, 3A/3B/3C/3D),
coronavirus ORF1ab (NSP1..NSP16), flavivirus polyprotein (C/prM/E/NS1/
NS2A/NS2B/NS3/…). Representative-level analysis on the polyprotein as a
whole loses signal — diversity comparisons, alignments, and trees are
far more biologically informative on the **mature peptides**.

`polyprotein:` declares one or more polyproteins, each with the list of
mature peptides it cleaves into. Each peptide carries an **HMM token**
(single profile name like `CoV_NSP8`, or a multidomain N-to-C
architecture like `CoV_RPol_N--RdRP_1` — same token syntax as
`cluster_protein.hmms[]`). After representatives are elected, repseq
finds each rep's polyprotein CDS (the CDS that satisfies the most
peptide tokens), slices it into mature peptides, and emits them as
accessory artifacts. **Clustering and the whole-genome tree are
untouched** — the polyprotein still drives those; mature peptides are
additive output.

### Config shape

Non-segmented (one polyprotein per virus, common case for
Picornaviridae and Coronaviridae):

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

`hmm:` accepts either a single profile name or a multidomain N→C
token (HMMs joined with `--`, most N-terminal first). For a
multidomain peptide, the token is satisfied only when every named HMM
hits the parent CDS in N→C order; the peptide's footprint is the
synthetic span from the first domain's start to the last domain's
end. Pfam-boundary fuzz at the seam is tolerated up to
`hmm.multidomain_overlap_tolerance` aa (default 30), the same setting
the `cluster_protein` HMM gate uses.

**Alternative architectures (OR) — `hmms:`** (v0.34.0+). When a mature
peptide's architecture differs across virus genera, declare a list of
alternative tokens with `hmms:` instead of the singular `hmm:`. The
peptide is located when **any one** of the listed tokens is satisfied;
if multiple alternatives are satisfied on a CDS, the slicer picks the
one with the best (lowest) worst-domain E-value (mirrors
`cluster_protein.hmms[]` ranking). The chosen architecture is recorded
in a `matched_architecture` column of the audit TSV and a
`[matched_arch=...]` tag in the peptide FASTA header so the bench
scientist can see which fold fired on each rep.

```yaml
clustering:
  polyprotein:
    - name: ORF1ab
      peptides:
        # Coronavirus NSP1 is non-homologous between alpha- and
        # beta-CoV — list both architectures and either one will fire.
        - {name: NSP1, hmms: [aCoV_NSP1, bCoV_NSP1], cleavage_motif: "LQ"}
        # NSP2 has one architecture; legacy singular `hmm:` form is fine.
        - {name: NSP2, hmm: CoV_NSP2, cleavage_motif: "LQ"}
        # NSP12 is multidomain; each entry of an `hmms:` list can be
        # multidomain too — alternative whole architectures.
        - {name: NSP12, hmm: "CoV_RPol_N--RdRP_1", cleavage_motif: "LQ"}
```

Each peptide must set **exactly one** of `hmm` or `hmms` — both is
an error (ambiguous), neither is an error (peptide has no locator).

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

Spec names must be unique across all `polyprotein:` entries within a
virus (they key the output filenames).

### Cut strategies (the load-bearing science choice)

1. **`boundary`** — each peptide spans its HMM hit's `ali_from..ali_to`
   verbatim. Deterministic, but loses the inter-domain residues at the
   seams.
2. **`bisect`** — cuts placed at the midpoint between adjacent peptide
   HMM hits. Endpoints extend to the protein N-/C-term **only when the
   first/last declared peptide is itself located** — a missing leading
   or trailing peptide leaves the gap unassigned (the v0.37.0 missing-
   neighbour-aware fix). Cut sites are geometric rather than biological.
3. **`motif`** (default when any peptide declares `cleavage_motif`) —
   start from `bisect`, then snap each inter-peptide cut to the last
   occurrence of the downstream peptide's `cleavage_motif` within
   ±`motif_window_aa` of the bisect point. Falls back to `bisect` for
   that cut if no motif lies in the window. `cleavage_motif:` is the
   residues found just **N-terminal** of the cut that liberates that
   peptide — `LQ` for coronavirus 3CL (Q is the recognised residue),
   `Q` for picornavirus 3C, `RR`/`R[RK]` for flavivirus NS2B-3.

The default cut strategy is `motif` if any peptide declared
`cleavage_motif`, else `bisect`. `boundary` is opt-in.

### Edge cases

- **Peptide HMM doesn't hit** the parent CDS → peptide skipped, audit
  row marks `missing`. Flanking peptides keep their HMM-hit boundary on
  the missing side (`cut_method_actual=hit-boundary`); the missing
  peptide's territory is left unassigned rather than absorbed by a
  neighbour. The v0.37.0 fix; pre-v0.37.0 the bisect blindly split the
  missing peptide's gap evenly between its neighbours, which on
  SARS-CoV-2 (NSP2 untraceable by the bundled HMMs) inflated NSP1 from
  ~180 aa to 505 aa and bled 313 aa of NSP2 territory into NSP3.
- **Peptide HMMs hit out of N→C order** on the parent → spec fails for
  that rep, all rows for it land as `out_of_order`. (Real polyproteins
  are linear; this is almost always an HMM-database problem.)
- **Adjacent peptide HMM hits overlap** on the parent → cut placed at
  the midpoint of the hit centres; sequences still emitted but rows
  flagged `overlap` so the user can audit.
- **No CDS has hits from ≥ `min_peptides_hit` peptide HMMs** → no
  parent identified; the spec emits no peptide FASTAs for that rep, and
  the audit TSV records one `no_parent_cds` row per peptide so the user
  sees *why*.

### Outputs

`{prefix}_polyprotein/{prefix}_<spec>_<peptide>.fasta` — one FASTA per
peptide of each spec, with all clean (`ok` / `overlap`) slices across
representatives. Headers carry the standard protein-FASTA bracket-tag
set (organism, taxonomy ladder, isolate, segment, host, subtype, country,
date, length, parent) plus polyprotein-specific tags:
`[polyprotein=<spec>]`, `[peptide=<peptide>]`,
`[peptide_range_aa=<from>-<to>]`, `[cut_method=<actual>]`,
`[matched_arch=<token>]`.

`{prefix}_polyprotein/{prefix}_<spec>_polyprotein.fasta` — the **whole
polyprotein** FASTA: the entire parent CDS (amino acids) for every
representative carrying it, one record each, with the same bracket-tag
headers as the per-peptide FASTAs. Always written (parallel to the
per-peptide FASTAs) so you can re-align / re-tree the whole-polyprotein
set yourself; independent of the `whole_polyprotein_tree` knob below.

`{prefix}_polyprotein/{prefix}_<spec>_peptides.tsv` — audit TSV, one row
per (representative × peptide) attempt: `isolate_id`, `peptide_name`,
`parent_accession`, `parent_protein_id`, `range_aa_from`, `range_aa_to`,
`length_aa`, `cut_method_actual`, `matched_architecture`, `status`,
`note`. Use this to spot peptides whose cuts are purely heuristic vs.
motif-supported, and which OR-alternative architecture fired on each rep.

**Per-peptide phylogenetic trees** (v0.34.0+, opt-in via
`--per-protein-phylo`). Each peptide of each spec gets its own ML tree
built the same way 2F builds the marker trees — MAFFT alignment,
IQ-TREE or FastTree per `phylo.tool`, rooting + LCA + colour palette
shared with every other tree of the run, optional trimAl trimming.
Settings come from `phylo.per_protein` verbatim — `min_taxa`,
`mafft.extra_args`, `trimal.*`, `domain_architecture`. Outputs sit
alongside the peptide FASTAs in `{prefix}_polyprotein/` with the
matching basename schema:

* `{prefix}_<spec>_<peptide>_msa.fasta`
* `{prefix}_<spec>_<peptide>_tree.nwk` (opt-in, `phylo.newick`)
* `{prefix}_<spec>_<peptide>_tree.xml`
* `{prefix}_<spec>_<peptide>_tree.pdf` / `.png` (figure, default-on `phylo.pdf`)
* `{prefix}_<spec>_<peptide>_tree_id_map.tsv`

Sparse peptides (fewer than `phylo.per_protein.min_taxa` reps with an
`ok`-status slice) are skipped with a log note. Peptide trees are
intentionally **not** included in the
`{prefix}_per_protein/{prefix}_incongruence.tsv` table — that one
compares whole-genome markers, and peptide-vs-peptide RF within a
polyprotein is a different scientific question. The phyloXML
`<domain_architecture>` on each leaf is **peptide-local** (boxes drawn
on the peptide's 1..length coordinate space, not on the parent
polyprotein), so multidomain peptides like NSP12 (`CoV_RPol_N--RdRP_1`)
show both subdomain boxes properly scoped to the cut peptide.

**Whole-polyprotein tree** (v0.50.0+, opt-in via
`phylo.per_protein.whole_polyprotein_tree`, default off). In addition to
the per-peptide trees, build **one tree per spec on the entire
polyprotein CDS** — not the sliced peptides. Each rep's parent CDS feeds
the same `_build_tree` engine (MAFFT / IQ-TREE-or-FastTree / rooting /
LCA / colour palette / `min_taxa` / the v0.49.0 PDF+PNG figure), so the
files mirror the peptide-tree schema with a `_polyprotein` basename:

* `{prefix}_<spec>_polyprotein_msa.fasta`
* `{prefix}_<spec>_polyprotein_tree.xml`
* `{prefix}_<spec>_polyprotein_tree.pdf` / `.png`
* `{prefix}_<spec>_polyprotein_tree_id_map.tsv`
* `{prefix}_<spec>_polyprotein_tree.nwk` (opt-in, `phylo.newick`)

The point is the phyloXML `<domain_architecture>`: where a peptide tree
shows one peptide-local box, the whole-polyprotein tree's leaf shows
**every HMM domain across the entire polyprotein, end-to-end** (the full
nsp / peptide layout), so Archaeopteryx draws the complete domain diagram
on each leaf. It's the natural companion to the peptide trees — the
peptide trees zoom in on one cleavage product's phylogeny; the
whole-polyprotein tree gives the full-length view with the whole domain
architecture. Off by default because it's one extra (large) alignment +
tree per spec; only fires under `--per-protein-phylo` with the HMM tier
active. Like the peptide trees, it is **not** added to the
`_incongruence.tsv` table.

### Soft-fail posture

When the HMM tier didn't run this session (HMMER not installed, no
HMMs declared, etc.), polyprotein cutting prints a single stderr line
and emits no outputs. Specs whose every rep yields `no_parent_cds`
still write the audit TSV (so the user sees why no FASTAs appeared)
but no FASTAs.

---

## The local cache

Every NCBI/UniProt lookup is saved to a small SQLite database
(`~/.repseq/cache/` by default) so you only pay the network cost once.
A few different "kinds" of lookups share this cache, each under its own
**source** name:

| Source              | What's cached                                                          |
|---------------------|------------------------------------------------------------------------|
| `ncbi_taxonomy`     | Taxonomy lineages for accessions / taxon IDs.                          |
| `ncbi_proteins`     | GenBank CDS records (the `/isolate`, `/segment`, protein translations). |
| `ncbi_nuc_seq`      | Nucleotide bodies fetched by `repseq replay` (separate from `ncbi_proteins` because old cached entries don't carry the body). v0.30.0. |
| `uniprot`           | UniProt entries pulled by accession.                                   |
| `hmmscan`           | Per-CDS HMM hit lists from `hmmscan` (see below).                       |

Manage it with:

```bash
repseq cache stats -c your_config.yaml                # how big is it, what's in it
repseq cache purge-expired -c your_config.yaml        # remove stale entries (TTL)
repseq cache clear -c your_config.yaml                # wipe everything (asks first)
repseq cache clear -c your_config.yaml --source ncbi_taxonomy   # wipe just one source
```

`--source` accepts any of the source names from the table above.

### When (and when not) to clear the `hmmscan` cache

The `hmmscan` cache stores, for every CDS amino-acid sequence the
pipeline has ever seen, the **raw HMM hit list** returned by `hmmscan`
against your HMM database. Because viral CDSes are very repetitive
across closely related isolates, this cache pays for itself within a
single run — and on re-runs the HMM step is essentially free.

Two things to know:

1. **Cutoff changes do NOT require a clear.** The cache stores the raw
   hits; the pass/fail flags (E-value, GA threshold, coverage,
   multidomain order, overlap tolerance) are recomputed every run from
   your live config. Editing `hmm.default_evalue`,
   `hmm.relative_length_cutoff`, `hmm.use_ga_when_available`, or
   `hmm.multidomain_overlap_tolerance` is fine — the next run will
   apply the new cutoffs to the cached hits.

2. **DB changes usually do NOT require a clear either.** Each cache
   entry is keyed on `(AA sequence, database fingerprint)`. The
   fingerprint is derived from the `.hmm` file's size and modification
   time, so pointing `hmm.database` at a different file — or
   rebuilding the bundled DB with `scripts/build_bundled_hmms.sh` —
   automatically invalidates the old entries: they simply become
   unreachable and the next run rescans against the new DB. They'll
   eventually fall out of the cache via TTL (`taxonomy.cache_ttl_days`,
   30 days by default).

You only need to manually clear the `hmmscan` cache when:

- **You want to reclaim disk space** after a DB swap — the orphaned
  old-fingerprint rows are dead weight until the TTL evicts them.
- **You edited the `.hmm` file in place without changing its size or
  mtime** (uncommon, but possible with tools that preserve timestamps,
  or with content-only edits at the exact same length). The DB
  fingerprint will look unchanged, and stale hits will be reused.
- **You suspect cache corruption** — a partial / interrupted earlier
  run, or a manual edit to the SQLite file.

In those cases:

```bash
repseq cache clear -c your_config.yaml --source hmmscan
```

This wipes only the `hmmscan` rows. Your NCBI taxonomy, GenBank CDS,
and UniProt caches stay intact — so the next run still skips the
network roundtrips and only the HMM scan re-fires. Add `--yes` to skip
the confirmation prompt.

To see how much you're holding before / after, use `repseq cache stats
-c your_config.yaml` — it lists entry counts per source.

---

## Troubleshooting

**"WARNING: no representative sequences were selected."**
The run finished but nothing came out. `repseq` prints the most likely cause; the
usual ones are:

- *No sequences were loaded* — the input path is wrong, the file is empty, or its
  header format wasn't recognised. Try `--source ncbi_virus` (or `ncbi` /
  `uniprot`) to force it.
- *QC removed everything* — your cleaning thresholds are too strict for this data.
  Look at `{prefix}_qc_removed.tsv` to see which step did it, then loosen that
  setting. A common one: `median_percent` length filtering on a mixed-gene file —
  switch to `min_max`.
- *The segmented step dropped everything* — no isolate had all its segments. Most
  often the `isolate_regex` doesn't match your headers; also check the segment
  names/aliases and any `segment_lengths` bounds.

**`MMseqs2Error` / "mmseqs not found"** (or **`CDHitError` / "cd-hit not
found"** if you've set `clustering.backend: cdhit`) — the similarity-clustering
program isn't installed or isn't on your `PATH`. Install it (see
[Installation](#installation)), or use a mode that doesn't need it (`global -n`,
or a stratified mode where every group is already small).

**`cd-hit identity threshold X is below the supported floor`** — cd-hit refuses
identities below 0.40 (protein) or 0.80 (nucleotide). Either raise your
threshold to the floor, or switch the backend to `mmseqs2`, which has no
identity floor.

**Everything is grouped under "Unknown"** in a taxonomic/host/geographic run — the
metadata lookups didn't run or didn't find anything. Don't use `--no-resolve` for
the real run, and make sure your `ncbi_email` is set in the config.

**Lookups are slow the first time** — that's expected; they're cached, so the
*second* run on the same data is fast. An NCBI API key speeds up the first run.

**Not sure if everything's installed?** Run `repseq doctor`. It checks
every required and optional dependency, the external tools (`mmseqs`,
`cd-hit`, `cd-hit-est`, `mafft`, `FastTree`, `iqtree2`), reaches NCBI
and UniProt to confirm the network is working, and verifies your config
is valid — then tells you in plain English what (if anything) needs
fixing. Add `--no-network` to skip the database pings if you're offline.

**`[phylo skipped]` / `[phylo failed]`** — the `--phylo` step is
fail-soft: if fewer than 3 representatives survived, or
`mafft` / `iqtree2` / `FastTree` are missing or errored, the message is
printed to stderr and the rest of the run's outputs are still written.
To enable it, install MAFFT plus the tree builder you want
(see [Installation](#installation)). IQ-TREE additionally refuses
UFBoot bootstrap with fewer than 4 sequences; the wrapper falls back to
no-bootstrap automatically in that case.

**`CDHitError: Cluster round-trip mismatch: N sequences in, M accounted for…`** —
the cd-hit wrapper does a strict round-trip check on its `.clstr` output
and refuses to silently undercount. The error names the specific seq IDs
that fell out (both directions: input IDs absent from the `.clstr`, and
`.clstr` IDs not in the input). The usual cause is a seq id with
characters cd-hit transforms in its output — internal `...`, very long
IDs, or whitespace. If you hit a new one, the error message points
straight at the offending id.

**`taxonomy_mismatch:<rank>` rows in `_qc_removed.tsv`** — these are
isolates dropped because their segments resolved to different taxa at
the configured rank (usually reassortants or `/isolate` collisions). To
keep them, set `segmented.taxonomy_consistency.enabled: false` in your
config (or relax the rank to a higher level, e.g. `genus`).

**Per-segment length-filter drops surprised you** — open `_run.log` and
scroll to the "QC SUMMARY" block. Since v0.9.1, segmented runs print a
per-segment, isolate-level breakdown (`L too short : 257` etc.). If a
single segment ate most of your input, that's where to widen the
`segment_lengths` bounds.

**`Strain-collision check: found N (isolate, segment) pair(s)…`** — the
detector (since v0.9.3) found two or more distinct accessions sharing a
strain-derived `isolate_id` *and* a segment, which is the over-merge
signature of the `/strain → isolate_id` fallback. By default the run
continues (action: `warn`); the colliding records get the longest
per-segment kept and the rest dedup-dropped. To remove them outright,
set `segmented.strain_collision_action: drop` — they'll then appear in
`_qc_removed.tsv` with reason `strain_collision:<segment>`. The
`isolate_id_source` column in `_representative_isolates.tsv` /
`_isolate_proteins.tsv` shows which records relied on the fallback in
the first place.

---

## Testing

```bash
pip install pytest
pytest tests/
```

The tests run fully offline — all network calls are simulated — so they're safe
to run anywhere and finish in a couple of seconds.

---

## License

`repseq` is free software, released under the **GNU General Public License,
version 3** (GPLv3). You may redistribute it and/or modify it under the terms
of that license; it comes with **no warranty**. The full license text is in the
[`LICENSE`](LICENSE) file at the root of this repository, and is also available
at <https://www.gnu.org/licenses/gpl-3.0.html>.

---

## Status

Current: **`v0.22.0`**. All 8 selection modes, protein-alphabet
clustering by default (`alphabet_for_clustering: protein`), MMseqs2 and
cd-hit backends, optional protein-annotation and protein-quality QC,
per-isolate taxonomy-consistency QC for segmented viruses, **HMM-based
identity QC with multidomain-architecture token notation** (HMMER
hmmscan + bundled 26-profile viral Pfam-A subset; user-supplied DB
supported), strain-as-isolate provenance + collision detection,
segment-name synonyms, **rich phyloXML output with taxonomy-driven leaf
coloring, partitioned-supermatrix tree for protein + IQ-TREE (default),
optional trimAl trimming, per-marker domain-architecture trees with
pairwise Robinson-Foulds incongruence table**, an optional UMAP/MDS
plot of the clustering, an auto-generated Methods-section starter
(`_summary.md`) on every run, **per-marker protein FASTAs always-on, a
separate `extra_protein:` channel for sparse accessory proteins, and
two plain-text diversity reports (`_taxonomic_report.txt`,
`_protein_taxonomic_report.txt`)**. **983 offline regression tests
pass**; the NCBI-backed paths have been validated end-to-end against
live influenza-A, peribunyaviridae, hantaviridae, and coronaviridae
datasets.

Highlights of recent releases (newest first; `git log` for the
complete history):

- **`v0.22.0`** — protein-centric outputs. Always-on per-marker
  protein FASTAs in `{prefix}_per_protein_fasta/`. New top-level
  `extra_protein:` for sparse accessory proteins (CoV ORF7a / HE /
  …) — emits per-protein FASTAs in `{prefix}_extra_protein_fasta/`
  and, with `--per-protein-phylo`, separate trees in
  `{prefix}_extra_protein/` (deliberately kept out of the
  incongruence table). New per-rank coverage + AA-length report
  `{prefix}_protein_taxonomic_report.txt`, cluster-driving markers
  tagged `*` so they stand out from accessories. **Relaxed
  multidomain token rule** — new `hmm.multidomain_overlap_tolerance`
  (default 30 aa) lets Pfam-boundary fuzz at e.g. the CoV S1/S2
  furin seam pass; full containment is still rejected.

- **`v0.21.0`** — non-segmented output schema parity + 7 new bundled
  HMM profiles. `_representative_sequences.tsv` now uses the
  **identical column schema** as `_representative_isolates.tsv` so
  one analysis script reads both modes; two new TSVs
  `_sequence_proteins.tsv` / `_representative_sequence_proteins.tsv`
  give non-segmented mode the same per-CDS coverage the segmented
  side already had. Bundled HMM DB grew from 19 → 26 with the
  coronavirus accessory + alternative-Spike profiles. A stderr
  warning fires when a configured HMM token names a profile absent
  from the database.

- **`v0.20.x`** — phylogenetic stack hardening. **v0.20.0** added
  partitioned-supermatrix trees for protein + IQ-TREE runs (default
  on; per-family MSA → column-wise supermatrix → ModelFinder per
  partition; falls back transparently to single-alignment when only
  one marker family is present), trimAl alignment trimming (off by
  default; soft-fails to the untrimmed alignment), and two new
  bundled HMMs. **v0.20.1** introduced the N-to-C
  domain-architecture token notation (`"A--B"` means A is N-terminal
  to B, mirroring how biology draws a protein). **v0.20.2** added
  alternative architectures within one marker (`hmms: [token1,
  token2]` is OR), so e.g. one Spike spec spans both alpha-CoV
  `CoV_S1--CoV_S2` and beta-CoV `bCoV_S1_N--bCoV_S1_RBD--CoV_S2`
  architectures.

- **`v0.18.0` / v0.17.0 / v0.16.0 / v0.15.0`** — per-protein
  phylogenetics. **v0.15.0** introduced `--per-protein-phylo` (one
  tree per HMM marker spec, useful for spotting reassortment as
  topological incongruence between segment-marker trees). **v0.16.0**
  added taxonomy-driven leaf colouring shared across the
  whole-genome tree and every per-marker tree (so a genus is the
  same colour everywhere — what makes incongruence visually obvious
  in Archaeopteryx). **v0.17.0** added the
  `{prefix}_incongruence.tsv` pairwise unrooted-Robinson-Foulds table
  so you don't have to eyeball it. **v0.18.0** scoped each leaf in a
  per-marker tree to its satisfying CDS and made MAFFT strategy
  configurable (default `--auto`; opt-in L-INS-i for publication
  runs).

- **`v0.14.x`** — HMM-tier matured into universal QC. **v0.14.1**
  reframed the HMM scan from "alphabet=protein only" (v0.13.0) to
  "runs on every run that declares any `hmms:`" — segments / sequences
  failing the HMM gate are dropped regardless of clustering alphabet.
  **v0.14.2** added a protein-quality QC step (ambiguous-residue
  fraction on CDS translations — closes the gap where a garbled
  translation could pass a presence-only protein-count check).
  **v0.14.3** was a scientific-accuracy audit pass (HMM coverage,
  length filter, cluster-size accounting). **v0.14.4** hardened the
  `--plot` step so matplotlib is the only hard dep and `umap-learn`
  is best-effort with a classical-MDS fallback.

For releases before v0.14, see `git log` — the v0.6 to v0.13 line
brought protein-alphabet clustering, the harmonised TSV/FASTA
schemas, the rich phyloXML output, `_summary.md`, per-isolate
taxonomy-consistency QC, strain-as-isolate provenance + collision
detection, and the initial HMM-based marker selection (`v0.13.0`).
