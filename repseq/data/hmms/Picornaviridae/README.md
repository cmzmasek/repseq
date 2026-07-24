# Picornaviridae HMM set (bundled, family-specific)

Drop one or more HMMER3 profile files (`*.hmm`) into this directory — one
protein per file is fine, or a few concatenated profiles per file. repseq
**combines every `*.hmm` in this directory into one pressed database** at
run time (cached, rebuilt only when a file changes), so you don't have to
`cat` and `hmmpress` them yourself.

Select this set from a config with the bare family name:

```yaml
hmm:
  enabled: true
  database: Picornaviridae    # → repseq/data/hmms/Picornaviridae/  (this directory)
```

(`database` also accepts an absolute path to a single `.hmm` file or to any
directory of `.hmm` files; `null` uses the bundled `repseq_viral_core.hmm`.)

Curated profiles that carry a `GA` gathering-threshold line are gated by it
automatically (`hmm.use_ga_when_available`); profiles without one fall back to
`hmm.default_evalue` + `hmm.relative_length_cutoff`. This README is ignored by
the loader (only `*.hmm` files are read).

## Bundled profiles

All from Pfam-A (CC0); **all 17 carry curated GA cutoffs**. Membership was
taken from InterPro's Pfam-entries-by-taxonomy query for *Picornaviridae*
(taxid 12058). That query returns 32 entries; the 15 excluded are host or
recombinant-construct artifacts (green fluorescent protein, KDPG aldolase,
AIG1, immunoglobulin C1-/V-set domains, SANTA, the telomere-complex families
STN1/TPP1/CST, LRAT, two DUFs) plus three families belonging to *other*
picorna-like families (calicivirus coat protein, waikavirus capsid, tungro
spherical virus peptidase).

| Pfam | Profile (`NAME`) | Polyprotein region |
|------|------------------|--------------------|
| PF05408 | `Peptidase_C28` | L (foot-and-mouth-disease leader protease) |
| PF11475 | `VP_N-CPKC` | Virion protein N-terminal domain |
| PF02226 | `Pico_P1A` | VP4 |
| PF08935 | `VP4_2` | VP4 (alternative family) |
| PF00073 | `Rhv` | Capsid jelly-roll — hits VP2, VP3 **and** VP1 |
| PF22663 | `Rhv_5` | Capsid, VP1-specific |
| PF12944 | `HAV_VP` | VP1 (*Hepatovirus*) |
| PF08762 | `CRPV_capsid` | VP1 (cripavirus-like) |
| PF00947 | `Pico_P2A` | 2A |
| PF01552 | `Pico_P2B` | 2B |
| PF20758 | `2B_soluble` | 2B soluble domain |
| PF00910 | `RNA_helicase` | 2C |
| PF08727 | `P3A` | 3A (poliovirus-like) |
| PF06363 | `Picorna_P3A` | 3A (alternative family) |
| PF06344 | `Parecho_VpG` | 3B / VPg (*Parechovirus*) |
| PF00548 | `Peptidase_C3` | 3C protease (picornain 3C) |
| PF00680 | `RdRP_1` | 3D polymerase |

## The capsid problem — why VP2 and VP3 cannot be sliced apart

`Rhv` (PF00073) is the picornavirus jelly-roll fold, and it matches **VP2,
VP3 and VP1** on the same polyprotein — averaging 2.2 hits per record.
Poliovirus (NC_002058, 2209 aa) scans as:

```
 [2-69 Pico_P1A] [93-308 Rhv] [370-529 Rhv] [626-852 Rhv_5] [637-799 Rhv]
 [899-1025 Pico_P2A] [1030-1128 Pico_P2B] [1243-1350 RNA_helicase]
 [1457-1515 P3A] [1566-1731 Peptidase_C3] [1775-2185 RdRP_1]
```

repseq resolves a peptide token to the **union span from its first matched
domain**, so on those hits:

| token | resolved span |
|---|---|
| `Rhv` | (93, 308) — VP2 |
| `Rhv--Rhv` | (93, 529) — VP2+VP3 |
| `Rhv--Rhv--Rhv` | (93, 799) |
| `Rhv_5` | (626, 852) — VP1 |

There is no token that yields *only* the second `Rhv` hit (370–529). Declaring
`VP2: {hmm: Rhv}` and `VP3: {hmm: "Rhv--Rhv"}` gives both peptides the same
span start (93), which `slicer.compute_cuts` rejects as out-of-N→C-order —
and that failure discards **the entire spec** for that representative, not
just the two peptides.

The bundled config therefore declares a single merged **`VP2-VP3`** peptide
using `Rhv--Rhv`, and takes VP1 from the distinct `Rhv_5` / `HAV_VP` /
`CRPV_capsid` families. Verified across 196 RefSeq polyproteins: **zero
N→C ordering violations**.

## Peptide coverage across the family

Computed with repseq's own `_satisfying_span_for_token`, on the longest CDS
of each of the 196 distinct *Picornaviridae* RefSeq organisms:

| Peptide | Token(s) | Organisms | % |
|---|---|---:|---:|
| L | `Peptidase_C28` | 5 | 3 % |
| VP4 | `Pico_P1A`, `VP4_2` | 43 | 22 % |
| VP2-VP3 | `Rhv--Rhv` | 176 | 90 % |
| VP1 | `Rhv_5`, `HAV_VP`, `CRPV_capsid` | 150 | 77 % |
| 2A | `Pico_P2A` | 33 | 17 % |
| 2B | `Pico_P2B`, `2B_soluble` | 39 | 20 % |
| 2C | `RNA_helicase` | 195 | 99 % |
| 3A | `P3A`, `Picorna_P3A` | 46 | 23 % |
| 3B | `Parecho_VpG` | 2 | 1 % |
| 3C | `Peptidase_C3` | 151 | 77 % |
| 3D | `RdRP_1` | 196 | 100 % |

The backbone (**VP2-VP3, VP1, 2C, 3C, 3D**) is family-wide; L, VP4, 2A, 2B,
3A and 3B are largely *Enterovirus*-specific because their Pfam families were
built on enteroviruses. The median record slices **5** peptides; 96 % slice at
least 3.

## Clustering marker

`RdRP_1` alone hits 196/196 organisms but is a broad family that also covers
caliciviruses, dicistroviruses and other picorna-like lineages. The
multidomain token **`RNA_helicase--RdRP_1`** (the canonical picornaviral
2C-helicase-upstream-of-3D-polymerase signature) covers 195/196 — 99.5 % —
while being genuinely picornavirus-specific, so the bundled config uses it.
Drop to plain `RdRP_1` if you see unexpected eliminations.

## Model-coverage cutoff

| Profile | median | 5th pct | min |
|---|---|---|---|
| `RdRP_1` | 0.97 | 0.81 | 0.76 |
| `RNA_helicase` | 1.00 | 0.98 | 0.98 |
| `Rhv` | 0.80 | 0.49 | 0.42 |
| `Rhv_5` | 1.00 | 0.68 | 0.38 |
| `Peptidase_C3` | 0.95 | 0.67 | 0.42 |

The gate profiles (`RdRP_1`, `RNA_helicase`) are tight, but the capsid and
protease families have a long tail. **`hmm.relative_length_cutoff: 0.35`** is
recommended so the `Rhv` (min 0.42) and `Rhv_5` (min 0.38) tails survive —
losing `Rhv` hits is expensive because `VP2-VP3` needs *two* of them. The
0.5 default would clip ~5 % of `Rhv` hits and silently shrink capsid coverage.
