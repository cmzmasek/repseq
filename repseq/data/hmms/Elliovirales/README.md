# Elliovirales HMM set (bundled, order-specific)

Drop one or more HMMER3 profile files (`*.hmm`) into this directory — one
protein per file is fine, or a few concatenated profiles per file. repseq
**combines every `*.hmm` in this directory into one pressed database** at
run time (cached, rebuilt only when a file changes), so you don't have to
`cat` and `hmmpress` them yourself.

Select this set from a config with the bare order name:

```yaml
hmm:
  enabled: true
  database: Elliovirales     # → repseq/data/hmms/Elliovirales/  (this directory)
```

(`database` also accepts an absolute path to a single `.hmm` file or to any
directory of `.hmm` files; `null` uses the bundled `repseq_viral_core.hmm`.)

Curated profiles that carry a `GA` gathering-threshold line are gated by it
automatically (`hmm.use_ga_when_available`); profiles without one fall back to
`hmm.default_evalue` + `hmm.relative_length_cutoff`. This README is ignored by
the loader (only `*.hmm` files are read).

## What Elliovirales is

Elliovirales is one of the orders the old *Bunyavirales* was split into when
that order was elevated to class **Bunyaviricetes** (the other order is
Hareavirales). Its member families are the segmented, negative-sense
ssRNA viruses below. **Do not confuse it with Hareavirales** — Phenuiviridae,
Nairoviridae, Arenaviridae and Hantaviridae's sister families are in *that*
order, not this one. (Note: NCBI places **Hantaviridae in Elliovirales**;
different authorities disagree, and repseq follows NCBI.)

| Family | Segments | Hosts | RefSeq records |
|---|---|---|---|
| Peribunyaviridae | 3 (L,M,S) | vertebrates / arthropods | ~455 |
| Hantaviridae | 3 (L,M,S) | rodents / shrews | ~150 |
| Fimoviridae (emaraviruses) | **4–10** | plants | ~159 |
| Phasmaviridae | 3 (L,M,S) | insects | ~100 |
| Tospoviridae | 3 (L,M,S) | plants (thrips-borne) | ~82 |
| Cruliviridae | 3 (L,M,S) | crustaceans | ~9 |
| Tulasviridae | — | fungi | ~2 |

## The one universal marker: the L-segment RdRp

Scanning the entire order's RefSeq proteome (957 records, 1,120 CDS) shows a
single protein family that spans **every** member family: the L-segment
RNA-dependent RNA polymerase, Pfam **`Bunya_RdRp` (PF04196)**.

| Family | Organisms with a `Bunya_RdRp` hit at GA |
|---|---|
| Peribunyaviridae | 148 / 148 |
| Tospoviridae | 26 / 28 |
| Fimoviridae | 27 / 27 |
| Hantaviridae | 45 / 54 |
| Phasmaviridae | 10 / 17 |
| Cruliviridae | 3 / 3 |

The glycoprotein (M) and nucleocapsid (S), by contrast, are **strongly
family-partitioned** — each family has its own non-homologous Pfam family, and
several have no profile at all (Tospovirus glycoprotein hits only 2/28;
Phasmavirus and Crulivirus nucleocapsids hit 0). So `Bunya_RdRp` is the only
protein that can gate the whole order, exactly as `Flavi_NS5 OR RdRP_3` is for
Amarillovirales.

**The bundled `repseq_config_elliovirales.yaml` therefore gates only the L
segment on `Bunya_RdRp` and leaves M and S ungated** (longest-CDS fallback,
which reliably picks the glycoprotein and nucleocapsid respectively). See that
config for how to switch on family-specific M/S gating for a single-family run.

### RdRp model coverage — why `relative_length_cutoff` must be low

`Bunya_RdRp` is a 741-aa model built on the conserved polymerase core; the
L proteins are 2,000–2,400 aa. In the divergent families only the core motifs
clear GA, so **model coverage runs low**:

| Family | median model coverage | min |
|---|---|---|
| Tospoviridae | 1.00 | 1.00 |
| Peribunyaviridae | 1.00 | 0.35 |
| Hantaviridae | 1.00 | 0.15 |
| Fimoviridae | 0.94 | 0.20 |
| **Phasmaviridae** | **0.36** | 0.26 |

The GA **bit-score** is the real identity gate here; the length filter is
secondary. `repseq_config_elliovirales.yaml` sets
`hmm.relative_length_cutoff: 0.2` so Phasmaviridae (median 0.36) and the
Hantaviridae low-coverage tail survive. A cutoff of 0.5 would silently
eliminate the entire Phasmaviridae family.

## Bundled profiles

All from Pfam-A (CC0); **all 17 carry curated GA cutoffs**. Membership was
taken from InterPro's Pfam-entries-by-taxonomy query for *Elliovirales*
(taxid 3151837). The four `Bunya_*` profiles are byte-identical to the copies
in `repseq_viral_core.hmm` (bundled sets are self-contained, so shared profiles
are duplicated deliberately); the rest are the current Pfam release.

| Pfam | Profile (`NAME`) | Segment | Role |
|------|------------------|---------|------|
| PF04196 | `Bunya_RdRp` | L | **RdRp core — the order-wide gate** |
| PF15518 | `L_protein_N` | L | L protein N-terminus (endonuclease region) |
| PF21561 | `L_thumb_ring_vir` | L | L polymerase thumb/ring domain |
| PF21991 | `capSnatchArena` | L | Cap-snatching endonuclease |
| PF12426 | `DUF3674` | L | RdRp-associated (hantaviral L) |
| PF03557 | `Bunya_G1` | M | *Peribunyaviridae* glycoprotein Gn |
| PF03563 | `Bunya_G2` | M | *Peribunyaviridae* glycoprotein Gc |
| PF01567 | `Hanta_Gn-H` | M | *Hantaviridae* Gn head |
| PF20679 | `Hanta_Gn-B` | M | *Hantaviridae* Gn base |
| PF01561 | `Hanta_Gc_N` | M | *Hantaviridae* Gc N-terminal |
| PF20682 | `Hanta_Gc_C` | M | *Hantaviridae* Gc C-terminal |
| PF00952 | `Bunya_nucleocap` | S | *Peribunyaviridae* nucleocapsid N |
| PF00846 | `Hanta_nucleocap` | S | *Hantaviridae* nucleocapsid N |
| PF01533 | `Tospo_nucleocap` | S | *Tospoviridae* nucleocapsid N |
| PF25629 | `Fimo_NCAP` | S | *Fimoviridae* nucleocapsid N |
| PF01104 | `Bunya_NS-S` | S | *Peribunyaviridae* NSs (non-structural) |
| PF03231 | `Tospov_NS-S_N` | S | *Tospoviridae* NSs N-terminal |

Even in the RdRp-only order run, these M/S profiles earn their place: the HMM
scan runs the whole database against every CDS, so `{prefix}_hmm_diagnostic.tsv`
annotates each M/S CDS with its family-specific identity (Gn/Gc/N) — you get
the identification without the gating drop.

## Scope: Fimoviridae is excluded by the 3-segment config

The bundled config declares `expected_segments: 3` (L, M, S), which covers the
six tri-segmented families. **Fimoviridae (emaraviruses) carry 4–10 segments**
(median 5) and cannot be assembled by a 3-segment definition — they need their
own config with the appropriate `expected_segments` and per-RNA markers.
`Fimo_NCAP` is bundled here anyway so a Fimoviridae-specific config can reuse
this set.

## Reassortment analysis lives in the family configs

Because M and S are ungated in the order-wide run, it produces one whole-genome
(L+M+S concat) tree, not per-segment marker trees — so it does not compute the
cross-segment reassortment signal. That analysis is inherently within-family
(a hantavirus does not reassort with a tospovirus), and belongs in the
family-level configs (`repseq_config_peribunyaviridae.yaml`,
`repseq_config_hantaviridae.yaml`, `repseq_config_phasmaviridae.yaml`), which
gate all three segments and enable `--per-protein-phylo`.
