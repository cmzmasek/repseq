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

**Plots** — if you want the optional diagnostic scatter plot of the clustering
result (matplotlib only — installs cleanly everywhere):

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

Every mode also accepts: `--input/-i`, `--output-dir/-o`, `--config/-c`,
`--threads`, `--seed`, `--segmented`, `--dry-run`, `--no-resolve`,
`--source {auto,uniprot,ncbi,ncbi_virus}`, `--overflow {keep,trim}`, `--plot`.

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
        │
        ▼
  (segmented mode only)
    populate /isolate, /segment from GenBank
    → drop isolates whose segments disagree on species (taxonomy_consistency)
    → keep only isolates with ALL expected segments
    → per-segment length bounds (drop the whole isolate if any segment is out of range)
    → fetch the marker protein per segment, concatenate per isolate
        │
        ▼
  group by your chosen mode (rank / host / decade / similarity / …)
  cluster each group down to ~N reps (binary search over identity threshold)
  pick the best representative per cluster (RefSeq > reviewed UniProt > longest)
        │
        ▼
  write:
    • selected representatives (FASTA, per-segment FASTAs, concatenated FASTA)
    • their proteins (per-CDS FASTA + per-CDS TSV)
    • per-rep metadata spreadsheet
    • per-stratum + per-cluster + per-drop TSVs
    • plain-text run log
    • (optional) UMAP/MDS clustering plot                        (--plot)
    • (optional) MAFFT MSA + IQ-TREE / FastTree + phyloXML tree  (--phylo)
    • (optional) one tree per HMM domain-architecture marker     (--per-protein-phylo)
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
| `{prefix}_group_counts.tsv` | yes | One row per stratum: in / out / clustered / cutoff. |
| `{prefix}_clusters.tsv` | yes | Which sequences ended up grouped, who represents whom. |
| `{prefix}_representatives.fasta` | non-segmented | Selected representative sequences. |
| `{prefix}_representative_sequences.tsv` | non-segmented | Per-representative metadata spreadsheet. |
| `{prefix}_concatenated.fasta` | segmented | Per-isolate head-to-tail concat of all segments. |
| `{prefix}_segment_<NAME>.fasta` | segmented | One file per expected segment, just the representative isolates. |
| `{prefix}_representative_isolates.tsv` | segmented | Per-representative-isolate metadata spreadsheet. |
| `{prefix}_isolate_proteins.tsv` | segmented + GenBank | Every protein of every isolate that survived QC, with a `representative` TRUE/FALSE column. |
| `{prefix}_representative_isolate_proteins.tsv` | segmented + GenBank | Same schema as above, row-filtered to representatives only. |
| `{prefix}_representative_isolate_proteins.fasta` | segmented + GenBank | AA FASTA of every protein of every representative isolate. |
| `{prefix}_representative_sequence_proteins.fasta` | non-segmented + GenBank | AA FASTA of every protein of every representative sequence. |
| `{prefix}_representatives_protein.fasta` | when `alphabet_for_clustering: protein` (default) | The AA strings actually fed into the clusterer. |
| `{prefix}_clustering.png` | only with `--plot` | Diagnostic scatter of the clustering. |
| `{prefix}_msa.fasta`, `_tree.nwk`, `_tree.xml`, `_tree_id_map.tsv` | only with `--phylo` | Alignment + tree + name mapping. |
| `{prefix}_iqtree_summary.txt` | only with `--phylo` + IQ-TREE | IQ-TREE ModelFinder report. |
| `{prefix}_summary.md` | every run | Auto-generated Methods-section starter (prose + numbers + tool citations). |

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

Use it for: **the protein set you'll almost certainly hand to BLAST,
DIAMOND, HMMER, MMseqs2 search, or any sequence-search tool** — it's
both pre-curated and pre-annotated. This is usually "the" output file
for protein-centric workflows.

#### `{prefix}_representatives_protein.fasta` *(when `clustering.alphabet_for_clustering: protein` actually fired)*

The **AA strings that were fed into the clusterer** — the per-isolate
marker-protein concat in segmented mode, or the per-rep marker in
non-segmented mode. Written only if every representative ended up with a
populated `protein_sequence` (i.e. the protein-alphabet path completed;
not written for `alphabet_for_clustering: nucleotide`).

This is a *diagnostic*, not a primary output: useful if you want to
reproduce or audit the clustering input. For a clean per-protein set, use
`_representative_*_proteins.fasta` instead.

### Spreadsheets (TSV)

All TSVs use a harmonised vocabulary (since v0.8.0): the canonical
identifier column is always `accession`; booleans are always `TRUE` /
`FALSE`; length columns carry their alphabet (`length_nt`, `length_aa`,
`segment_length_nt`, `total_length_nt`); taxonomic ranks always appear in
the same nine-rank ladder (`species`, `subgenus`, `genus`, `subfamily`,
`family`, `suborder`, `order`, `subclass`, `class`). Sub-ranks come from
the NCBI lineage and are commonly blank for viruses.

#### `{prefix}_representative_sequences.tsv` — non-segmented mode

One row per representative sequence. Columns:

`accession`, `organism`, `description`, `strain`, `host`,
`collection_date`, `country`, `segment`, `isolate_id`, `molecule_type`,
`length_nt`, `is_refseq`, `is_reviewed`, `ncbi_taxon_id`, then the
nine-rank taxonomic ladder.

Open in Excel / Numbers / your scripting language — this is the
spreadsheet you'll usually hand to a collaborator.

#### `{prefix}_representative_isolates.tsv` — segmented mode

One row per representative **isolate** (not per sequence — segmented
representatives are whole isolates, not single segments). Columns:

`isolate_id`, `isolate_id_source`, `organism`, `strain`, `host`,
`collection_date`, `country`, `n_segments`, `segments` (comma-joined
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

#### `{prefix}_clusters.tsv` — always

One row per cluster member. Columns: `cluster_id`, `accession` (of the
*representative*), `organism`, `cluster_size`, `is_refseq`, `is_reviewed`.

Use it for: tracing which sequences ended up in the same cluster, or for
re-running selection with a different priority (e.g. "I want the reviewed
UniProt entry not the longest one") without re-clustering.

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

### Run record

#### `{prefix}_run.log`

Plain-text record of the run: version, date, exact command-line, every
config setting that was active, the per-step QC counts, the mode's
results, and the list of output files. Save this with your data — it
makes the selection fully reproducible.

The "QC SUMMARY" block also shows the new **per-segment length-filter
breakdown** (since v0.9.1) so you can tell *which* segment caused
isolates to fall out, not just that "some did".

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
short-id remap → MAFFT (`--auto`) → tree builder → root → label internal
nodes by LCA → phyloXML. Tree builder is auto-picked from your
clustering alphabet: **IQ-TREE for protein** (ModelFinder + UFBoot
bootstrap by default) and **FastTree for nucleotide** (with `-nt -gtr`).
Set `phylo.tool: fasttree` (or `iqtree`) to pin one.

#### `{prefix}_msa.fasta`

The MAFFT alignment. Headers are `>S0001 <pretty-label>` — the short
`SNNNN` id stays the first whitespace-separated token (safe for any
phylo tool), and the descriptive label (built from `phylo.labeling.format`
/ `segmented_format`) is appended as the FASTA description so AliView /
Jalview / MEGA show recognisable names without confusing the underlying
tools.

#### `{prefix}_tree.nwk`

Tree-builder Newick output. Leaf names are the **short ids** (`S0001`,
`S0002`, …) so the file works with any downstream tool that expects
short, safe leaf labels. Use the id map below to recover real names.

#### `{prefix}_tree.xml` — **the tree you'll usually open**

phyloXML with rich, browseable annotation:

- Every leaf gets a formatted `<name>`, a `<taxonomy>` block with NCBI
  taxon id, a `<sequence>` block with the GenBank accession + title, and
  one repseq-namespaced `<property>` per per-leaf attribute (host,
  collection_date, country, strain, isolate_id, year, plus the full
  9-rank taxonomy ladder: species, subgenus, genus, subfamily, family,
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
- Every annotated internal clade gets a `<name>` and a `<taxonomy>`
  block holding the LCA's scientific name + rank (`min_rank=genus` by
  default; the labeller keeps each monophyletic clade labelled at its
  crown and suppresses obvious duplications like a 2-leaf same-species
  pair).
- The `<phylogeny>` element carries a `<name>` and `<description>` with
  MAFFT/IQ-TREE/FastTree versions, the selected substitution model, the
  bootstrap settings, and the rooting method that actually fired.
- Confidence values are normalised to 0–100 integers (`sh_like` for
  FastTree, `ufboot` for IQ-TREE). The tree is ladderized.

Opens directly in [Archaeopteryx](https://sites.google.com/view/aptxjs)
([web version](https://sites.google.com/view/aptxjs))
([source](https://github.com/cmzmasek/forester)) which makes use of the
the rich annotation.
The annotation is also visible in any other phyloXML-aware tool (Dendroscope, etc).

#### `{prefix}_tree_id_map.tsv`

Two-column mapping `short_id` ↔ `accession`. Use it when you need to
trace a leaf in `_tree.nwk` or the MSA back to a real sequence id.

#### `{prefix}_iqtree_summary.txt` — only when IQ-TREE ran

The IQ-TREE ModelFinder report — which substitution model was selected
and why, plus the log-likelihood and bootstrap settings. Worth a glance
before quoting a model in a methods section.

> **Fail-soft:** if MAFFT, IQ-TREE, or FastTree are missing, or fewer
> than 3 representatives survived, the phylogeny step is skipped with a
> stderr message and the rest of the run's outputs are still written.
> IQ-TREE additionally refuses UFBoot with fewer than 4 sequences; the
> wrapper drops bootstrap automatically in that case but still produces
> the tree.

### Per-protein trees — `--per-protein-phylo`

`--phylo` builds **one** tree from the whole representatives (the marker
concat in segmented mode). `--per-protein-phylo` instead builds **one
tree per declared HMM domain-architecture** — the same `hmms:` tokens you
set for QC under `segment_markers` / `cluster_protein`. For each token
(e.g. `Bunya_nucleocap`, or the multidomain `Bunya_G1--Bunya_G2`), repseq
picks the CDS that satisfies it on every representative carrying that
architecture and runs the *same* MAFFT → IQ-TREE/FastTree → root → LCA
pipeline on those protein translations.

Why you'd want it: comparing the single-marker trees side by side reveals
**topological incongruence** — an L-segment polymerase tree disagreeing
with an M-segment glycoprotein tree is the classic signature of
**reassortment**, which the concatenated `--phylo` tree hides.

Outputs land in a `{prefix}_per_protein/` subdirectory, one set per built
family (`<family>` = the sanitised token, prefixed with the segment in
segmented mode, e.g. `M_Bunya_G1--Bunya_G2`):

```
{prefix}_per_protein/
  M_Bunya_G1--Bunya_G2_msa.fasta
  M_Bunya_G1--Bunya_G2_tree.nwk
  M_Bunya_G1--Bunya_G2_tree.xml
  M_Bunya_G1--Bunya_G2_tree_id_map.tsv
  S_Bunya_nucleocap_msa.fasta
  …
```

Each file has the same format and rich annotation as its `--phylo`
counterpart. The flag runs alone or alongside `--phylo`.

> **Requirements / fail-soft:** needs the HMM tier (`hmm.enabled` plus
> configured `hmms:`) to have run, and `mafft` + a tree builder on PATH.
> A family carried by fewer than `phylo.per_protein.min_taxa` (default 3)
> representatives is skipped with a log note; if no family qualifies — or
> the HMM tier didn't run — the whole step is skipped with one stderr
> line, leaving the rest of the run intact.

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
  length_filter:
    mode: median_percent
    min_percent: 50
  ambiguous_threshold: 0.05

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
  priority: [refseq, reviewed_uniprot, longest]   # tie-break order for picking the "best"
```

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

Override at run time without editing the YAML:

```bash
repseq global -c my.yaml -i x.fasta -T 0.95 --alphabet-for-clustering nucleotide
```

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

### Token notation (v0.14.0)

Each element of `hmms:` is a **token** string. A token is either:

| Form | Meaning |
| --- | --- |
| `"Name"`         | Single HMM. A CDS satisfies the token if it has a passing hit to `Name`. |
| `"A--B"`         | Multidomain. A CDS satisfies the token only when it has passing hits to **both** `A` and `B`, with `A` lying C-terminal to `B` (strict non-overlap). |
| `"A--B--C"`      | Same idea, three domains: `A` most C-terminal, `C` most N-terminal. |

**C-to-N order**. The first HMM in a multidomain token is the most
**C-terminal** domain on the protein. Most molecular-biology notation
reads N-to-C left-to-right; this is the opposite. Check your token
order against the actual domain architecture before running. The
convention exists so multidomain tokens can be compared directly
against the hmmscan `ali_from`/`ali_to` columns.

**Extra domains are fine.** A CDS annotated as `HMMX--A--B` still
satisfies the token `"A--B"` because `A` remains C-terminal to `B`.

**Hard cutover from v0.13.0.** In v0.13.0, `hmms: [A, B]` meant "both
A and B must hit the same CDS." In v0.14.0 that same YAML now means
"two separate single-HMM tokens" (loose). To get the v0.13.0 strict
multidomain semantic, write `hmms: ["A--B"]`.

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
3. A **segment / sequence passes QC** when, for each token in the
   spec, at least one CDS in the segment satisfies that token. (Tokens
   in the same `hmms:` list are AND — every token must be satisfied.
   Tokens are independent of each other: the same CDS can satisfy
   multiple tokens, or different CDSes can each satisfy one.)
4. If any token has no satisfying CDS, the segment fails. In segmented
   mode the entire isolate is then dropped (one bad segment fails the
   whole isolate, counted under `removed_hmm_failed` with a per-marker
   breakdown in `removed_hmm_by_marker`). In non-segmented mode the
   single sequence is dropped.
5. For protein-alphabet clustering, the marker CDS of each surviving
   segment / sequence is then the **longest** CDS that satisfies any
   token in the spec.

Markers WITHOUT an `hmms:` list keep the legacy alias → longest chain
unchanged and are not gated by HMM-QC.

### Bundled vs user-supplied database

`hmm.database: null` (default) uses the bundled set at
`repseq/data/hmms/repseq_viral_core.hmm` — 17 hand-picked Pfam-A
profiles covering the most common viral marker proteins (RdRp,
nucleocapsid, glycoprotein, helicase, protease) across the main
families. Pfam-A is licensed CC0 so redistribution is unrestricted.
The bundled set was assembled from Pfam-A via
`scripts/build_bundled_hmms.sh`; you can re-run that script to refresh
it against a newer Pfam release.

To use your own profiles (e.g. VOGdb, a custom curated set, or the
full Pfam-A.hmm), set `hmm.database` to an absolute path. The first
run auto-`hmmpress`-es the file if the `.h3*` index files are missing.

`repseq doctor` reports the DB status (path, profile count, indexed
Y/N, GA-cutoff coverage).

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
        # domains in C-to-N order (Bunya_G1 C-terminal, Bunya_G2
        # N-terminal). Use this when you want to reject post-cleavage
        # annotations as inconsistent.
        # M: {hmms: ["Bunya_G1--Bunya_G2"]}
        L: {hmms: ["RdRP_4"]}
```

When both `segment_markers` and the legacy `cluster_protein` define a
marker for the same segment, `segment_markers` wins.

### HMM-related config keys

```yaml
hmm:
  enabled: true                  # master switch
  database: null                 # null → bundled; abs path → user-supplied
  default_evalue: 1.0e-5         # used when a profile has no curated GA
  use_ga_when_available: true    # prefer curated GA over default_evalue
  relative_length_cutoff: 0.5    # ali_span / hmm_len ≥ this
  threads: null                  # null → falls through to cfg.threads
```

### Soft-fail posture

The HMM tier never aborts a run on its own. If `hmmscan` is missing,
the database is unreadable, or hmmscan returns non-zero, repseq emits
one stderr warning and falls through to the alias / longest-CDS
chain. To check that everything is wired up, run `repseq doctor`.

---

## The local cache

Every NCBI/UniProt lookup is saved to a small database (`~/.repseq/cache/` by
default) so you only pay the network cost once. Manage it with:

```bash
repseq cache stats                           # how big is it, what's in it
repseq cache purge-expired                   # remove stale entries
repseq cache clear                           # wipe everything
repseq cache clear --source ncbi_taxonomy    # wipe just one kind of lookup
```

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

## Status

Current: **`v0.13.0`**. All 8 selection modes, protein-alphabet clustering
by default (`alphabet_for_clustering: protein`), MMseqs2 and cd-hit
backends, optional protein-annotation QC, per-isolate
taxonomy-consistency QC for segmented viruses, **HMM-based marker
selection (HMMER hmmscan + bundled Pfam-A subset)**, strain-as-isolate
provenance + collision detection, segment-name synonyms, rich phyloXML
output, an optional UMAP/MDS plot of the clustering, an auto-generated
Methods-section starter (`_summary.md`) on every run, and a per-stratum
diversity curve in `_group_counts.tsv`. **741 offline regression tests
pass**; the NCBI-backed paths have been validated end-to-end against
live influenza-A, peribunyaviridae, and hantaviridae datasets.

Highlights of recent releases (newest first):

- **`v0.13.0`** — HMM-based marker selection. When `hmm.enabled` is
  true (default) and a marker definition carries `hmms: [...]` (in
  `clustering.cluster_protein` or per-segment `segment_markers`),
  HMMER hmmscan is used as the AUTHORITATIVE gate — a CDS whose
  `/product` matches an alias but FAILS the HMM check (E-value or
  coverage) is rejected, and the segment / sequence is dropped.
  Multi-HMM lists are AND (every HMM must hit the same CDS). Bundled
  set of 17 viral Pfam-A profiles in `repseq/data/hmms/`; auto-
  `hmmpress`-ed on first run. User-supplied DB via `hmm.database`.
  Soft-fails (warns + falls back to alias/longest) when hmmscan is
  missing or the database is unavailable. Two new QC counters
  (`removed_hmm_failed`, `removed_hmm_by_marker`); new `hmmscan`
  column right of `representative` in `_isolate_proteins.tsv` and
  `_representative_isolate_proteins.tsv` showing passing hits
  (`Bunya_Gn(E=1e-30,cov=0.85);Bunya_Gc(E=3e-25,cov=0.78)` format,
  best-E first); same string emitted as `[hmmscan=...]` tag on
  proteins-FASTA headers. `repseq doctor` reports DB status (path,
  profile count, indexed Y/N, GA-cutoff coverage). The Methods-
  section summary names HMMER + the cutoffs when the tier fired. 63
  new tests bring the suite from 678 → 741. See the new
  "HMM-based marker selection" section above for configuration and
  examples.

- **`v0.12.0`** — `_group_counts.tsv` gains a per-stratum diversity
  curve. For each stratum where the binary-search clustering ran, the
  clustering backend is also invoked at each of the configured
  "standard" identity thresholds (default `[0.99, 0.95, 0.9, 0.8, 0.7]`)
  and the resulting cluster counts are reported as trailing
  `n_clusters_<c>` columns. Purely diagnostic — representative
  selection is unchanged. Cutoffs below the backend's identity floor
  (cd-hit-est < 0.80, cd-hit protein < 0.40) are reported as `NA`.
  Configurable via `clustering.diversity_curve_cutoffs`; set to `[]`
  to disable on very large runs where the extra clustering work would
  dominate runtime. Mentioned by the `_summary.md` writer so the
  Methods section is aware of it.
- **`v0.11.0`** — every successful run now writes `{prefix}_summary.md`,
  a Markdown Methods-section starter that reads as scientific prose
  with the live numbers from the run filled in (input count, per-stage
  QC counts, segmented-completeness counts, clustering input
  description, mode rationale, phylogenetic inference description if
  `--phylo` ran). Includes a Software & references table with versions
  auto-detected from `--version` probes on every external tool repseq
  could have used (cd-hit, cd-hit-est, MMseqs2, MAFFT, FastTree,
  IQ-TREE) plus full journal citations. Soft-fails on render error so a
  bug in the writer never voids a successful selection.
- **`v0.10.2`** — actual fix for cd-hit-est round-trip mismatches.
  v0.10.1 sanitized IUPAC ambiguity codes but a hantaviridae run still
  hit `91 in, 34 accounted for`. Root cause: the `.clstr` parser regex
  only matched cd-hit (protein) member lines (`at 99.93%`) and missed
  cd-hit-est (nucleotide) strand-prefixed ones (`at +/99.93%`,
  `at -/99.93%`). Every member line silently failed to parse, making
  every cluster look like a singleton. Regex now accepts the optional
  `+/`/`-/` prefix; protein output unaffected.
- **`v0.10.1`** — fixed a silent data-loss bug in cd-hit-est when
  clustering on nucleotides. cd-hit-est rejects sequences containing
  IUPAC ambiguity codes (W, R, Y, K, M, S, B, D, H, V) but writes the
  warning to stderr only — our wrapper swallowed it on a 0 exit, and a
  hantaviridae run hit 91 concat sequences in / 34 out / 57 silently
  dropped. Fixed by sanitizing the cd-hit FASTA input (non-ACGTN → N)
  before clustering; original `seq.sequence` and all downstream output
  FASTAs are untouched. MMseqs2 tolerates IUPAC codes natively so its
  path is unchanged. A one-time stderr notice prints when
  sanitization fires.
- **`v0.10.0`** — the clustering alphabet is now a user-facing choice
  end-to-end. The config knob `clustering.alphabet` was renamed to
  `clustering.alphabet_for_clustering` to make the scope unambiguous
  (it only chooses what the clustering backend sees — GenBank CDS
  download, protein-count QC, and the per-segment
  `expected_proteins_per_segment` check still run regardless). A new
  CLI flag `--alphabet-for-clustering [protein|nucleotide]` overrides
  the YAML at run time. The little-used `auto` value was dropped:
  pick `protein` or `nucleotide` explicitly (a `--no-resolve` run must
  pick `nucleotide`). **Breaking change** for any YAML that set
  `clustering.alphabet`.
- **`v0.9.4`** — the QC summary's `Passed QC` line was misleading: it
  only counted what survived the basic-QC stages (duplicates / length /
  ambiguous / annotation), not protein-QC, segmented completeness,
  taxonomy-consistency, strain-collision, or per-segment length. It's
  now relabeled `Passed basic QC`, and a new `Final survivors : N
  isolates (after every QC stage)` line is appended at the bottom of
  the QC Summary block (and the success summary in stderr) so the
  number the user actually cares about — what reached selection — is
  visible at a glance.
- **`v0.9.3`** — segmented mode now records the **provenance of every
  `isolate_id`** (new `isolate_id_source` column: `isolate` | `strain` |
  `regex`) and runs a **strain-collision detector** on the strain-derived
  ids. Default action is `warn`; set
  `segmented.strain_collision_action: drop` to remove colliders into
  `_qc_removed.tsv` with reason `strain_collision:<segment>`. See
  [Strain-as-isolate fallback](#strain-as-isolate-fallback-and-collision-detection).
- **`v0.9.2`** — hardened the cd-hit `.clstr` round-trip: the
  member-line regex used to silently drop seq ids containing internal
  `...`; now it's greedy + end-anchored, and the round-trip error
  message names the specific unmatched ids in both directions.
- **`v0.9.1`** — segmented runs now report a **per-segment, isolate-level
  length-filter breakdown** during the run and in the QC summary
  (`L too short : 257`, `M too long : 6`, …) instead of the misleading
  "skipped (segmented mode)" line.
- **`v0.9.0`** — new segmented-only QC step
  `segmented.taxonomy_consistency` (on by default). Drops any isolate
  whose segments disagree on the configured taxonomy rank (default
  `species`) — typically reassortants or `/isolate` collisions. Missing
  labels are ignored; dropped segments appear in `_qc_removed.tsv` with
  reason `taxonomy_mismatch:<rank>`.
- **`v0.8.0` / `v0.8.1`** — TSV and FASTA output schemas overhauled and
  harmonised. Canonical identifier column is always `accession`;
  booleans are always `TRUE` / `FALSE`; length columns carry their
  alphabet (`length_nt`, `length_aa`, `segment_length_nt`,
  `total_length_nt`); taxonomic ladder is the same nine ranks
  everywhere. Representative FASTA now splits by mode
  (`_representative_isolate_proteins.fasta` segmented vs
  `_representative_sequence_proteins.fasta` non-segmented) with NCBI
  bracket-tag headers carrying organism, full taxonomy, host, country,
  collection date, and protein length. **Breaking change** for any
  downstream script that parses output by column or filename.
- **`v0.7.0` / `v0.7.1`** — rich phyloXML output: every leaf carries
  taxonomy, accession, host/country/date as `<property>` elements; every
  annotated internal clade carries an LCA scientific name + rank;
  taxonomy-guided rooting (with MAD / midpoint fallbacks); deterministic
  short-id remap so phylo binaries can't lose track of identity.
- **`v0.6.0` / `v0.6.1`** — **protein-alphabet clustering is now the
  default.** Segmented runs cluster on the per-isolate concatenated
  marker proteins (longest CDS per segment, overridable). IQ-TREE
  becomes the protein-tree default (ModelFinder + UFBoot); FastTree
  stays the nucleotide default.
