# Phenuiviridae HMM set (bundled, family-specific)

Drop one or more HMMER3 profile files (`*.hmm`) into this directory — one
protein per file is fine, or a few concatenated profiles per file. repseq
**combines every `*.hmm` in this directory into one pressed database** at
run time (cached, rebuilt only when a file changes), so you don't have to
`cat` and `hmmpress` them yourself.

Select this set from a config with the bare family name:

```yaml
hmm:
  enabled: true
  database: Phenuiviridae     # → repseq/data/hmms/Phenuiviridae/  (this directory)
```

(`database` also accepts an absolute path to a single `.hmm` file or to any
directory of `.hmm` files; `null` uses the bundled `repseq_viral_core.hmm`.)

Curated profiles that carry a `GA` gathering-threshold line are gated by it
automatically (`hmm.use_ga_when_available`); profiles without one fall back to
`hmm.default_evalue` + `hmm.relative_length_cutoff`. This README is ignored by
the loader (only `*.hmm` files are read).

## Why a separate set — the bundled core does NOT cover this family

The core set ships four `Bunya_*` profiles, but three of them were built on
**Orthobunyavirus** (family *Peribunyaviridae*) and do **not** generalise to
*Phenuiviridae*. Scanning the Rift Valley fever / SFTS / Uukuniemi reference
proteomes against `repseq_viral_core.hmm` at curated GA gives:

| Core profile | Hits Phenuiviridae? |
|---|---|
| `Bunya_RdRp` (PF04196) | **yes** — all three L proteins, score 335–844, model coverage 0.92–0.94 |
| `Bunya_G1` (PF03557) | no |
| `Bunya_G2` (PF03563) | no |
| `Bunya_nucleocap` (PF00952) | no |

Gating the M segment on `Bunya_G1`/`Bunya_G2` would therefore drop **every**
phenuivirus isolate at the marker gate — the same silent-elimination failure
mode as the pre-v0.66.0 Amarillovirales/`Flavi_NS5` episode. The phenuivirus
glycoprotein and nucleocapsid have their own Pfam families, bundled here.

## Bundled profiles

All from Pfam-A (CC0); **all 15 carry curated GA cutoffs**. Family membership
was taken from InterPro's Pfam-entries-by-taxonomy query for *Phenuiviridae*
(taxid 1980418), not from name matching.

| Pfam | Profile (`NAME`) | Segment | Phenuivirus protein |
|------|------------------|---------|---------------------|
| PF04196 | `Bunya_RdRp` | L | RNA-dependent RNA polymerase (core domain) |
| PF15518 | `L_protein_N` | L | L protein N-terminus (endonuclease region) |
| PF12603 | `L_PA-C-like` | L | PA-C-like domain of L |
| PF07243 | `Phlebovirus_G1` | M | Glycoprotein Gn (older "G1" nomenclature) |
| PF07245 | `Phlebovirus_G2` | M | Glycoprotein Gc, fusion domain |
| PF19019 | `Phlebo_G2_C` | M | Glycoprotein Gc, C-terminal domain |
| PF07246 | `Phlebovirus_NSM` | M | Non-structural protein NSm |
| PF05733 | `Tenui_N` | S | Nucleocapsid N (Tenuivirus/Phlebovirus) |
| PF11073 | `NSs` | S | Non-structural protein NSs (RVFV-like) |
| PF03300 | `Tenui_NS4` | — | Tenuivirus NS4 movement protein |
| PF04876 | `Tenui_NCP` | — | Tenuivirus major non-capsid protein |
| PF05310 | `Tenui_NS3` | — | Tenuivirus movement protein |
| PF06656 | `Tenui_PVC2` | — | Tenuivirus PVC2 protein |
| PF25219 | `Znf_Tenuivirus` | — | Tenuivirus zinc finger |
| PF01107 | `MP` | — | Viral movement protein (plant-infecting genera) |

`Bunya_RdRp.hmm` is byte-identical to the copy in `repseq_viral_core.hmm` —
bundled sets are self-contained (`hmm.database` selects one set, it does not
add to the core), so the shared profile is duplicated here deliberately.

## The M-segment glycoprotein is one CDS, not two

Unlike *Peribunyaviridae*, the phenuivirus M segment encodes a **single
glycoprotein precursor polyprotein (GPC)** that is cleaved co-translationally
into Gn and Gc. All three domains therefore sit on **one** CDS, in strict
N→C order, with no overlap:

```
   RVFV  M (1197 aa):  [NSm 1-135] [G1 156-689] [G2 691-1012] [G2_C 1027-1196]
   SFTSV M (1073 aa):            [G1 198-555] [G2 563-872]  [G2_C  907-1071]
   Uuk   M (1008 aa):            [G1  80-494] [G2 514-834]  [G2_C  850-998 ]
```

That makes the multidomain token `Phlebovirus_G1--Phlebovirus_G2--Phlebo_G2_C`
a clean completeness check for a full-length GPC — the domain seams are
16–35 aa apart, well inside the default 30 aa
`hmm.multidomain_overlap_tolerance` (they do not overlap at all, so a
tolerance of `0` also works).

## Family coverage — read this before scoping a run

Scanned against all 473 *Phenuiviridae* RefSeq records (608 CDS, 154 distinct
organisms), the marker profiles cover the family **unevenly, and the gaps are
biological, not accidental**:

| Marker | Organisms covered |
|---|---|
| `Bunya_RdRp` (L / polymerase) | ~all |
| `Tenui_N` (S / nucleocapsid) | ~all |
| `Phlebovirus_G1--G2--G2_C` (M / glycoprotein) | **94 / 154 (61 %)** |

Of the 60 organisms the glycoprotein token misses:

* **34 have no glycoprotein CDS at all.** *Coguvirus*, *Rubodvirus*,
  *Laulavirus*, *Entovirus*, *Lentinuvirus*, *Bocivirus*, *Mechlorovirus* are
  bi-segmented plant viruses (L + S only); *Tenuivirus* has 4–8 segments and is
  non-enveloped. These are **correct** drops — a glycoprotein-based analysis
  cannot include a virus with no glycoprotein.
* **7 have an annotated glycoprotein that scores below GA** (*Mobuvirus*,
  *Horwuvirus*, some *Tenuivirus*) — genuinely too divergent for
  Phlebovirus-built profiles.
* **19 hit only part of the architecture** — mostly the insect-specific
  genera (*Phasivirus*, *Goukovirus*) and divergent *Uukuvirus*. A looser
  2-domain token `Phlebovirus_G1--Phlebovirus_G2` recovers 5 of these
  (94 → 99); a bare `Phlebovirus_G2` recovers more, at the cost of accepting
  partial GPCs.

Per-genus, the medically important tri-segmented genera are well covered:
*Phlebovirus* 67/68, *Bandavirus* 9/10, *Uukuvirus* 12/18.

**Practical consequence:** a family-wide run gated on the glycoprotein is
really a run over the *enveloped tri-segmented* phenuiviruses. That is a
defensible scope, but state it — and check the
`{prefix}_flags.txt` "Taxa eliminated entirely by QC" section, which will
name the dropped genera explicitly.

## Model-coverage cutoff

Model coverage of GA-passing hits is tight at the median but has a real tail
(fragmentary records and divergent genera):

| Profile | median | 5th pct | min |
|---|---|---|---|
| `Bunya_RdRp` | 1.00 | 0.54 | 0.05 |
| `Tenui_N` | 0.99 | 0.75 | 0.50 |
| `Phlebovirus_G1` | 0.99 | 0.43 | 0.20 |
| `Phlebovirus_G2` | 0.99 | 0.95 | 0.23 |
| `Phlebo_G2_C` | 0.94 | 0.74 | 0.50 |

`Phlebovirus_G1` is a long model (526 aa) and the divergent genera align to
less of it, so **`hmm.relative_length_cutoff: 0.4`** is recommended over the
0.5 default — 0.5 clips part of the legitimate `G1` tail. Raise it toward
0.8 if you want near-full-length markers only.
