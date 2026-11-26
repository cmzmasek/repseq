# Coronaviridae HMM set (bundled, family-specific)

Drop one or more HMMER3 profile files (`*.hmm`) into this directory — one
protein per file is fine, or a few concatenated profiles per file. repseq
**combines every `*.hmm` in this directory into one pressed database** at
run time (cached, rebuilt only when a file changes), so you don't have to
`cat` and `hmmpress` them yourself.

Select this set from a config with the bare family name:

```yaml
hmm:
  enabled: true
  database: Coronaviridae    # → repseq/data/hmms/Coronaviridae/  (this directory)
```

(`database` also accepts an absolute path to a single `.hmm` file or to any
directory of `.hmm` files; `null` uses the bundled `repseq_viral_core.hmm`.)

Curated profiles that carry a `GA` gathering-threshold line are gated by it
automatically (`hmm.use_ga_when_available`); profiles without one fall back to
`hmm.default_evalue` + `hmm.relative_length_cutoff`. This README is ignored by
the loader (only `*.hmm` files are read).

## Why this set exists

The generic `repseq_viral_core.hmm` carries only **four** of Pfam's **six**
coronavirus spike domains. It is missing `CoV_S1_C` (PF19209) and `CoV_S2_C`
(PF19214), and it has no hemagglutinin-esterase profile at all. The
consequence for spike work is that the C-terminal end of S1 and the
intravirion tail of S2 carry no box in any genus.

> **Corrected 2026-11-03.** An earlier version of this README called
> `CoV_S1_C` "the S1 CTD — the RBD in alpha/gamma/delta". That is wrong. Pfam
> describes PF19209 as **subdomain 2 (SD2)**, a small all-β domain named from
> the PEDV spike structure, and its hits confirm it: a ~60-residue box
> immediately N-terminal of the S1/S2 boundary in every genus (beta ~623–681,
> alpha ~647–703, gamma ~432–496, delta ~430–486 on full-length spikes). It
> does **not** split the alpha/delta `CoV_S1` box into NTD and RBD; that box
> still spans both. Gammacoronavirus (IBV) S1 has almost no GA-passing box
> other than `CoV_S1_C` (`CoV_S1` passes on 1 of 1,071 gamma spikes).

`CoV_S1_C` is not a marginal addition: in a full-Pfam domain-architecture
survey of Nidovirales it appears in **30 of 57** observed spike architectures,
against only 7 for the generic `CoV_S1`.

This set therefore = every coronavirus-specific profile from the core set, plus
the generic profiles coronavirus configs actually reference (`AAA_12`, `Macro`,
`Peptidase_C16`, `RdRP_1`), plus the three additions above.

**All 47 profiles carry GA cutoffs.** Because GA decisions are bit-score based
and bit scores do not depend on database size, switching a config from
`repseq_viral_core.hmm` (59 profiles) to this set (47) **cannot change any
gating decision** — only the reported E-values move.

## Bundled profiles (47)

All from Pfam-A (CC0). `GA` is the curated gathering threshold in bits;
`LENG` is the model length.

### Spike — the S1/S2 surface glycoprotein

| Pfam | Profile (`NAME`) | LENG | GA | Coronavirus protein |
|------|------------------|------|------|---------------------|
| PF01600 | `CoV_S1`        | 413 | 25.0 | Spike S1 (generic, receptor-binding subunit) |
| PF16451 | `bCoV_S1_N`     | 299 | 26.0 | Spike S1 N-terminal domain (NTD), betacoronavirus-like |
| PF09408 | `bCoV_S1_RBD`   | 163 | 25.0 | Spike S1 receptor-binding domain, betacoronavirus |
| **PF19209** | **`CoV_S1_C`**  | 59  | 27.0 | **Spike S1 C-terminal subdomain 2 (SD2), just before the S1/S2 boundary — not an RBD** |
| PF01601 | `CoV_S2`        | 525 | 40.0 | Spike S2 (fusion subunit) |
| **PF19214** | **`CoV_S2_C`**  | 42  | 27.0 | **Spike S2 intravirion / cytoplasmic tail** |

Architecture tokens these support (N-to-C; `--` separates domains, extra
unlisted domains are tolerated):

```
CoV_S1--CoV_S2                            generic spike
CoV_S1--CoV_S1_C--CoV_S2                  generic spike, SD2 resolved
bCoV_S1_N--CoV_S2                         NTD + S2
bCoV_S1_RBD--CoV_S2                       RBD + S2
bCoV_S1_N--bCoV_S1_RBD--CoV_S2            betacoronavirus spike, NTD + RBD
bCoV_S1_N--CoV_S1_C--CoV_S2               NTD + SD2
bCoV_S1_N--CoV_S1_C--CoV_S2--CoV_S2_C     ... with the intravirion tail
CoV_S1_C--CoV_S2                          SD2 + S2 (IBV form; tested as a gate
                                          token 2026-11-03 and rejected -- it
                                          drops 114 IBV records, see the config)
CoV_S2                                    S2 only (truncated / diverged S1)
```

### Other structural and accessory proteins

| Pfam | Profile (`NAME`) | LENG | GA | Coronavirus protein |
|------|------------------|------|------|---------------------|
| PF00937 | `CoV_nucleocap`   | 344 | 30.0 | Nucleocapsid (N) |
| PF01635 | `CoV_M`           | 206 | 27.0 | Membrane / matrix glycoprotein (M) |
| PF02723 | `CoV_E`           | 74  | 27.0 | Small envelope protein (E) |
| **PF03996** | **`Hema_esterase`** | 238 | 25.0 | **Hemagglutinin-esterase (HE) — embecoviruses, toroviruses** |
| PF11289 | `bCoV_viroporin`  | 273 | 27.0 | Betacoronavirus viroporin |
| PF03053 | `Corona_NS3b`     | 226 | 25.0 | ORF3b accessory protein |
| PF12133 | `bCoV_NS6`        | 62  | 27.0 | Betacoronavirus NS6 |
| PF02398 | `Corona_7`        | 78  | 27.0 | Protein 7 / ORF7a |

### ORF1ab replicase — mature peptides (NSP1–NSP16)

| Pfam | Profile (`NAME`) | LENG | GA | Peptide / domain |
|------|------------------|------|------|------------------|
| PF16688 | `CNV-Replicase_N`   | 108 | 27.1 | NSP1, N-terminal |
| PF11501 | `bCoV_NSP1`         | 138 | 25.0 | NSP1, betacoronavirus |
| PF11963 | `B-CoV_A_NSP1`      | 355 | 25.0 | NSP1, betacoronavirus lineage A |
| PF19211 | `CoV_NSP2_N`        | 201 | 29.7 | NSP2, N-terminal |
| PF19212 | `CoV_NSP2_C`        | 120 | 25.0 | NSP2, C-terminal |
| PF17896 | `NSP2_gammaCoV`     | 358 | 27.0 | NSP2, gammacoronavirus |
| PF12379 | `bCoV_NSP3_N`       | 149 | 26.7 | NSP3, N-terminal |
| PF01661 | `Macro`             | 118 | 22.0 | NSP3 macrodomain (ADRP) |
| PF11633 | `bCoV_SUD_M`        | 143 | 24.0 | NSP3 SUD middle domain |
| PF12124 | `bCoV_SUD_C`        | 66  | 25.0 | NSP3 SUD C-terminal domain |
| PF08715 | `CoV_peptidase`     | 320 | 26.0 | NSP3 papain-like protease (PL-pro) |
| PF01831 | `Peptidase_C16`     | 249 | 27.7 | PL-pro, peptidase C16 family |
| PF16251 | `bCoV_NAB`          | 128 | 28.7 | NSP3 nucleic-acid-binding domain |
| PF19218 | `CoV_NSP3_C`        | 457 | 25.0 | NSP3, C-terminal |
| PF19217 | `CoV_NSP4_N`        | 352 | 25.0 | NSP4, N-terminal |
| PF16348 | `CoV_NSP4_C`        | 96  | 27.0 | NSP4, C-terminal |
| PF05409 | `Peptidase_C30`     | 284 | 27.0 | NSP5 main protease (3CL-pro) |
| PF19213 | `CoV_NSP6`          | 263 | 25.0 | NSP6 |
| PF08716 | `CoV_NSP7`          | 83  | 31.8 | NSP7 |
| PF08717 | `CoV_NSP8`          | 199 | 27.0 | NSP8 |
| PF08710 | `CoV_NSP9`          | 110 | 27.0 | NSP9 |
| PF09401 | `CoV_NSP10`         | 122 | 28.8 | NSP10 |
| PF06478 | `CoV_RPol_N`        | 350 | 27.0 | NSP12 RdRp, N-terminal (NiRAN) |
| PF00680 | `RdRP_1`            | 450 | 23.9 | NSP12 RdRp core |
| PF20632 | `CoV_NSP13_ZBD`     | 93  | 27.7 | NSP13 helicase, zinc-binding domain |
| PF20633 | `CoV_NSP13_stalk`   | 48  | 27.7 | NSP13 helicase, stalk |
| PF20631 | `CoV_NSP13_1B`      | 82  | 28.4 | NSP13 helicase, 1B domain |
| PF13087 | `AAA_12`            | 195 | 27.0 | NSP13 helicase, AAA+ ATPase core |
| PF06471 | `CoV_ExoN`          | 515 | 27.0 | NSP14 proofreading exoribonuclease |
| PF19219 | `CoV_NSP15_N`       | 61  | 27.0 | NSP15, N-terminal oligomerisation |
| PF19216 | `CoV_NSP15_M`       | 118 | 25.0 | NSP15, middle domain |
| PF19215 | `CoV_NSP15_C`       | 154 | 23.9 | NSP15 uridylate-specific endoribonuclease |
| PF06460 | `CoV_Methyltr_2`    | 299 | 27.0 | NSP16 2'-O-methyltransferase |

## Notes and cautions

- **`CoV_S2` is the strongest single spike signal** (GA 40.0 over 525
  positions) and passed on 100 % of representatives in a 484-rep family-wide
  run. A bare `CoV_S2` token is therefore a recall-maximal spike gate — useful
  when you want every spike, but it is *not* an architecture check and it will
  admit S2-only fragments.
- **`CoV_S1_C` and `CoV_S2_C` are short models** (59 and 42 positions). With
  `hmm.relative_length_cutoff: 0.5` they need only ~30 and ~21 aligned
  residues, so they hit readily; rely on their GA cutoffs (27.0 each) rather
  than loosening the coverage floor.
- **`Hema_esterase` (HE) is lineage-restricted** — embecoviruses (OC43, HKU1,
  BCoV, MHV) and toroviruses carry it, most other coronaviruses do not.
  Declare it under `extra_protein:`, never as a required `cluster_protein`
  marker, or it will drop most of the family.
- **`Macro`, `AAA_12`, `Peptidase_C16` and `RdRP_1` are not
  coronavirus-specific.** They are included because coronavirus ORF1ab
  architectures reference them. Do not gate a whole run on one of them alone.
- **No torovirus spike profile is included.** The older `Spike_torovirin`
  model is not in current Pfam-A under that name; torovirus spikes will
  resolve only through `CoV_S2` (if at all). Add a profile here if torovirus
  coverage matters.
- Rename map for older Pfam names, if you are reusing pre-2022 domain
  architecture tables: `Corona_S1`→`CoV_S1`, `Corona_S2`→`CoV_S2`,
  `Spike_NTD`→`bCoV_S1_N`, `Spike_rec_bind`→`bCoV_S1_RBD`.

## Provenance

Profiles extracted from `repseq/data/hmms/repseq_viral_core.hmm` with
`hmmfetch` (44 of them), plus PF19209, PF19214 and PF03996 retrieved from the
InterPro/Pfam API on 2026-11-03. Pfam-A is CC0, so redistribution is
unrestricted.

Run `repseq doctor` to report the resolved database path, profile count,
indexed state and GA-cutoff coverage — use that rather than trusting the
count quoted here.
