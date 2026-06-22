# Amarillovirales HMM set (bundled)

This set covers the **order Amarillovirales** — the clade historically treated
as a single family *Flaviviridae*, which NCBI now splits into three families:

| Family | Representative genera |
|--------|-----------------------|
| Flaviviridae   | Orthoflavivirus (dengue, Zika, WNV, YFV, JEV, TBEV, …) |
| Hepaciviridae  | Hepacivirus (HCV), Pegivirus |
| Pestiviridae   | Pestivirus (BVDV, CSFV) |

The bundled profiles are Orthoflavivirus-tuned `Flavi_*` markers (which match
*only* family Flaviviridae) **plus** the pan-genus `RdRP_3` RdRp (which reaches
Hepaci-/Pegi-/Pestivirus). OR-ing `Flavi_NS5` with `RdRP_3` therefore gives a
marker that spans the whole order; using `Flavi_NS5` alone scopes a run to
family Flaviviridae (Orthoflavivirus).

Drop one or more HMMER3 profile files (`*.hmm`) into this directory — one
protein per file is fine, or a few concatenated profiles per file. repseq
**combines every `*.hmm` in this directory into one pressed database** at
run time (cached, rebuilt only when a file changes), so you don't have to
`cat` and `hmmpress` them yourself.

Select this set from a config with the bare set name:

```yaml
hmm:
  enabled: true
  database: Amarillovirales    # → repseq/data/hmms/Amarillovirales/  (this directory)
```

(`database` also accepts an absolute path to a single `.hmm` file or to any
directory of `.hmm` files; `null` uses the bundled `repseq_viral_core.hmm`.)

Curated profiles that carry a `GA` gathering-threshold line are gated by it
automatically (`hmm.use_ga_when_available`); profiles without one fall back to
`hmm.default_evalue` + `hmm.relative_length_cutoff`. This README is ignored by
the loader (only `*.hmm` files are read).

## Bundled profiles

All from Pfam-A (CC0), carrying curated GA cutoffs. The first block is the
**Orthoflavivirus** polyprotein, listed in N→C order (structural C–prM/M–E
first, then the non-structural NS1–NS5 cassette); these `Flavi_*` /
`Peptidase_S7` / `FtsJ` profiles are *Orthoflavivirus-specific* (family
Flaviviridae) and do **not** match the structurally distinct proteins of the
other two families (Hepaciviridae, Pestiviridae):

| Pfam | Profile (`NAME`) | Flaviviral protein / domain |
|------|------------------|------------------------------|
| PF01003 | `Flavi_capsid`      | Capsid protein C |
| PF01570 | `Flavi_propep`      | Polyprotein propeptide (pr of prM) |
| PF01004 | `Flavi_M`           | Membrane glycoprotein M |
| PF00869 | `Flavi_glycoprot`   | Envelope E — central + dimerisation domains |
| PF02832 | `Flavi_glycop_C`    | Envelope E — immunoglobulin-like domain |
| PF21659 | `Flavi_E_stem`      | Envelope E — stem/anchor domain |
| PF00948 | `Flavi_NS1`         | NS1 (non-structural) |
| PF01005 | `Flavi_NS2A`        | NS2A |
| PF01002 | `Flavi_NS2B`        | NS2B (protease cofactor) |
| PF00949 | `Peptidase_S7`      | NS3 — serine protease (peptidase S7) |
| PF07652 | `Flavi_DEAD`        | NS3 — DEAD-box helicase domain |
| PF20907 | `Flav_NS3-hel_C`    | NS3 — helicase C-terminal helical domain |
| PF01350 | `Flavi_NS4A`        | NS4A |
| PF01349 | `Flavi_NS4B`        | NS4B |
| PF01728 | `FtsJ`              | NS5 — FtsJ-like methyltransferase (N-terminal) |
| PF00972 | `Flavi_NS5`         | NS5 — RdRp fingers + palm domains |
| PF20483 | `Flavi_NS5_thumb`   | NS5 — RdRp thumb domain |

**Pan-genus RdRp** (added so the non-Orthoflavivirus families can be identified
at all — their RdRp is a different Pfam family, not `Flavi_NS5`):

| Pfam | Profile (`NAME`) | Flaviviral protein / domain |
|------|------------------|------------------------------|
| PF00998 | `RdRP_3`            | NS5B RdRp of **Hepacivirus / Pegivirus / Pestivirus** |

`RdRP_3` is the HCV NS5B family (representative structure 1gx6). Hepacivirus and
Pegivirus cover it near-fully; the more divergent Pestivirus NS5B clears the GA
bit-score by a wide margin but only matches its conserved core (~43 % of the
486-aa model), so admitting Pestivirus needs `hmm.relative_length_cutoff` ≤ 0.4.

## Scoping a run: order vs. family

- **Whole order Amarillovirales** (all genera) — OR the two RdRp profiles as the
  polyprotein marker: `hmms: ["Flavi_NS5", "RdRP_3"]` (Orthoflavivirus →
  `Flavi_NS5`, the other three genera → `RdRP_3`). See the worked example
  `repseq_config_amarillovirales.yaml`.
- **Family Flaviviridae only** (Orthoflavivirus) — use `hmms: ["Flavi_NS5"]`
  alone; the `RdRP_3`-only genera fail the gate and drop. See
  `repseq_config_flaviviridae.yaml`.

Reference a profile in config by its `NAME`, e.g. `hmms: ["Flavi_NS1"]` or a
multidomain token like
`"Peptidase_S7--Flavi_DEAD--Flav_NS3-hel_C"` (NS3 protease→helicase) or
`"FtsJ--Flavi_NS5--Flavi_NS5_thumb"` (NS5 MTase→RdRp).
