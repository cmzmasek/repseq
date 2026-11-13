# Orthoherpesviridae HMM set (bundled, family-specific)

Drop one or more HMMER3 profile files (`*.hmm`) into this directory — one
protein per file is fine, or a few concatenated profiles per file. repseq
**combines every `*.hmm` in this directory into one pressed database** at
run time (cached, rebuilt only when a file changes), so you don't have to
`cat` and `hmmpress` them yourself.

Select this set from a config with the bare family name:

```yaml
hmm:
  enabled: true
  database: Orthoherpesviridae   # → repseq/data/hmms/Orthoherpesviridae/  (this directory)
```

(`database` also accepts an absolute path to a single `.hmm` file or to any
directory of `.hmm` files; `null` uses the bundled `repseq_viral_core.hmm`.)

Curated profiles that carry a `GA` gathering-threshold line are gated by it
automatically (`hmm.use_ga_when_available`); profiles without one fall back to
`hmm.default_evalue` + `hmm.relative_length_cutoff`. This README is ignored by
the loader (only `*.hmm` files are read).

**Name:** the ICTV renamed the family *Herpesviridae* → *Orthoherpesviridae*
and NCBI Taxonomy follows (taxid 3044472), so this set carries the current
name. It covers the three subfamilies *Alpha-*, *Beta-* and
*Gammaherpesvirinae* only — not *Alloherpesviridae* (fish, amphibians) or
*Malacoherpesviridae* (molluscs), whose proteins these profiles do not detect.

## Bundled profiles

All from Pfam-A 38.0 (CC0); **all 32 carry curated GA cutoffs**. They cover the
**30 proteins common to the whole family** (gB contributes 3 profiles, the
polymerase 2, the spliced terminase subunit 1 two), chosen from data — see
"How the set was chosen" below. Gene names are
the HSV-1 ones — orthologues carry a different name in every virus (the DNA
polymerase is UL30 in HSV, UL54 in HCMV, BALF5 in EBV, ORF9 in KSHV), which is
why herpesvirus markers must be identified by HMM rather than by `/product`
aliases.

| Pfam | Profile (`NAME`) | GA | Protein (HSV-1 gene) |
|------|------------------|----|----------------------|
| PF03104 | `DNA_pol_B_exo1`  | 23   | DNA polymerase catalytic subunit (UL30) — 3'-5' exonuclease domain |
| PF00136 | `DNA_pol_B`       | 23   | DNA polymerase catalytic subunit (UL30) — polymerase domain |
| PF02689 | `Herpes_Helicase` | 27   | helicase-primase helicase subunit (UL5) |
| PF03121 | `Herpes_UL52`     | 21   | helicase-primase primase subunit (UL52) |
| PF03122 | `Herpes_MCP`      | 25   | major capsid protein (UL19) |
| PF01802 | `Herpes_V23`      | 25   | capsid triplex subunit 2 (UL18, VP23) |
| PF01366 | `PRTP`            | 24.6 | DNA packaging terminase subunit 2 (UL28) |
| PF00747 | `Viral_DNA_bp`    | 25   | single-stranded DNA-binding protein (UL29, ICP8) |
| PF00606 | `Glycoprotein_B`  | 43   | envelope glycoprotein B (UL27) |
| PF01771 | `Viral_alk_exo`   | 20   | alkaline exonuclease (UL12) |
| PF03167 | `UDG`             | 24.3 | uracil-DNA glycosylase (UL2) |
| PF00716 | `Peptidase_S21`   | 27   | capsid maturation protease (UL26) |
| PF01528 | `Herpes_glycop`   | 25   | envelope glycoprotein M (UL10) |
| PF01499 | `Herpes_UL25`     | 25   | DNA packaging tegument protein UL25 |
| PF01646 | `Herpes_UL24`     | 27   | nuclear protein UL24 |
| PF03327 | `Herpes_VP19C`    | 25   | capsid triplex subunit 1 (UL38, VP19C) |
| PF01763 | `Herpes_UL6`      | 27   | capsid portal protein (UL6) |
| PF02718 | `Herpes_UL31`     | 25   | nuclear egress lamina protein (UL31) |
| PF04541 | `Herpes_U34`      | 26   | nuclear egress membrane protein (UL34) |
| PF04559 | `Herpes_UL17`     | 25   | DNA packaging tegument protein UL17 |
| PF17488 | `Herpes_glycoH_C` | 25   | envelope glycoprotein H (UL22) |
| PF01673 | `Herpes_env`      | 25   | DNA packaging protein UL32 |
| PF03044 | `Herpes_UL16`     | 25   | tegument protein UL16 |
| PF05459 | `Herpes_UL69`     | 25   | multifunctional expression regulator (UL54, ICP27) |
| PF01677 | `Herpes_UL7`      | 25   | tegument protein UL7 |
| PF04843 | `Herpes_teg_N`    | 25   | large tegument protein (UL36) |
| PF02867 | `Ribonuc_red_lgC` | 23.9 | ribonucleotide reductase subunit 1 (UL39) |
| PF03324 | `Herpes_HEPA`     | 25   | helicase-primase subunit (UL8) |
| PF03581 | `Herpes_UL33`     | 25   | DNA packaging protein UL33 |
| PF00692 | `dUTPase`         | 22   | deoxyuridine triphosphatase (UL50) |
| PF02500 | `DNA_pack_N`      | 25   | DNA packaging terminase subunit 1 (UL15) — N-terminal domain |
| PF02499 | `DNA_pack_C`      | 22.8 | DNA packaging terminase subunit 1 (UL15) — C-terminal domain |

Reference a profile in config by its `NAME`, e.g. `hmms: ["Herpes_MCP"]`, or
the polymerase as the multidomain token `"DNA_pol_B_exo1--DNA_pol_B"` (the
exonuclease domain lies N-terminal to the polymerase domain in every
herpesvirus polymerase checked).

## How the set was chosen (2026-11-12)

From data, not from profile names: the CDS translations of the **114 RefSeq
genome-length Orthoherpesviridae records** (55 alpha, 27 beta, 32 gamma;
10,739 proteins) were scanned against all of Pfam-A 38.0 at GA. 298 families
hit; about 35 hit at least half the genomes of every subfamily (the
herpesvirus core genes). 35 families hit at least half the genomes of every subfamily and collapse to
**exactly 30 proteins** — the herpesvirus core. The seven CLUSTERING markers are
the vfam_trees family reference tree's set with one substitution (UL28 for
UL15); the other 23 profiles cover the rest of the 30 for per-protein output
and conservation analysis.

Genomes with a CDS satisfying the token at GA:

| Protein | Token | RefSeq genomes |
|---|---|---|
| DNA polymerase | `DNA_pol_B_exo1--DNA_pol_B` | 111/114 |
| helicase | `Herpes_Helicase` | 114/114 |
| primase | `Herpes_UL52` | 114/114 |
| major capsid protein | `Herpes_MCP` | 112/114 |
| triplex subunit 2 | `Herpes_V23` | 114/114 |
| terminase subunit 2 (UL28) | `PRTP` | 114/114 |
| ssDNA-binding protein | `Viral_DNA_bp` | 110/114 |
| glycoprotein B | `Glycoprotein_B` | 114/114 |

Every miss above is a record with **no CDS annotated** for that gene (e.g.
the 2002 murine cytomegalovirus RefSeq NC_004065.1 annotates 161 CDS, 79 of
them "hypothetical protein", and has no polymerase, capsid or DNA-binding
protein CDS) — an annotation gap, not missing biology.

## Why the terminase ATPase subunit (UL15) is bundled but NOT a gated marker

The profiles are here (`DNA_pack_N--DNA_pack_C`) because UL15 is one of the 30
family-wide proteins and worth reporting, but **do not make it a
`cluster_protein`**: as a gated marker it fails on annotation, and gating on it
would drop a fifth of all genomes. the gene is **spliced** (two exons) across all three
subfamilies, and in 14 of 20 failing RefSeq genomes only the second exon was
annotated as the CDS (295–452 aa instead of ~700 aa). The C-terminal Pfam
domain (`DNA_pack_C`, PF02499) then hits strongly but the N-terminal one
(`DNA_pack_N`, PF02500) does not, so the two-domain token passes only 94 of
114 genomes. Terminase subunit 2 (UL28, `PRTP`) is unspliced and present in
all 114, and replaces it. If you add UL15 yourself, expect the same loss.

## Things to know when configuring

- **Lengths differ by subfamily** (medians, aa): polymerase 1,220 alpha /
  1,037 beta / 1,012 gamma; primase 1,071 / 917 / 842; UL28 783 / 768 / 684.
  A family-wide `length_aa` range therefore catches only gross truncations;
  `qc.marker_completeness.require_met_start` covers 5' truncation (none of
  the ~820 RefSeq marker proteins starts without Met).
- **Profile coverage.** A few strong hits are split in two by an insertion,
  so their best single domain covers only ~0.50 of the model (e.g. the
  chimpanzee CMV polymerase, 1,357 aa: 0.52; a bat gammaherpesvirus
  ssDNA-binding protein: 0.50). They pass the default
  `hmm.relative_length_cutoff: 0.5` only just.
- Adding profiles to this directory rescales every hmmscan E-value (E is
  proportional to database size) but leaves bit scores, and therefore every
  GA-gated pass/fail decision, unchanged.
