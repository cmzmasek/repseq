#!/bin/bash
# Build repseq/data/hmms/repseq_viral_core.hmm from individual Pfam-A profiles.
#
# Reproducible recipe: each profile is fetched from InterPro by accession,
# decompressed, validated (the NAME line must match the expected name),
# concatenated, and finally hmmpress-ed.
#
# Run from the repo root:
#     bash scripts/build_bundled_hmms.sh
#
# Pfam-A is licensed CC0 (public domain); redistribution is unrestricted.
# Pfam release pinned by the URL pattern (latest by default — change to
# /api/entry/pfam/<acc>/<release>/?annotation=hmm to pin a release).

set -euo pipefail

OUTDIR="repseq/data/hmms"
OUTFILE="$OUTDIR/repseq_viral_core.hmm"
mkdir -p "$OUTDIR"

# (Pfam accession, expected NAME, role)
#
# Substitutions made after audit against the InterPro REST API:
#   - PF02078 (originally guessed for "Filo_NP") was Synapsin (non-viral);
#     substituted PF05505 = Ebola_NP, which is the closest Pfam entry for
#     filovirus nucleoprotein (Pfam has no generic Filo_NP family).
#   - PF03537 (originally guessed for Bunya_G2) was Glyco_hydro_114
#     (non-viral); substituted the correct PF03563 = Bunya_G2.
#   - PF02477 (originally guessed for Bunya_nucleocap) was Nairo_nucleo
#     (nairovirus-specific); substituted PF00952, whose Pfam NAME is
#     literally "Bunya_nucleocap" (broader family coverage).
#   - PF00863 (originally guessed for Peptidase_C18) was Peptidase_C4
#     (picornaviral 2A protease). MEROPS family C18 has no dedicated
#     Pfam profile, so that exact slot stays dropped. The viral-protease
#     slot it was meant to fill is now occupied by PF05409
#     (Peptidase_C30, the coronavirus 3C-like main protease / Mpro) —
#     a real Pfam family with a curated GA cutoff that lets coronavirus
#     ORF1ab/replicase be HMM-gated on its protease. Other viral-protease
#     candidates if more are wanted: PF00548 (Picornain 3C),
#     PF00770 (Adenovirus endoprotease), PF00851 (Helper component
#     proteinase), PF00863 (Peptidase_C4), PF02902 (Ulp1 protease),
#     PF03290 (Vaccinia I7).
#   - PF04196 (Bunya_RdRp) added so bunyavirus L segments can be
#     HMM-gated on their RNA-dependent RNA polymerase (the bundle's
#     RdRP_4/RdRP_2 don't cover the bunyavirus L protein).
#
# Bundle ends up with 19 profiles.
PROFILES=(
    "PF02123 RdRP_4 RdRp_polymerase"
    "PF00978 RdRP_2 RdRp_polymerase"
    "PF00972 Flavi_NS5 Flaviviral_polymerase"
    "PF00946 Mononeg_RNA_pol Mononegavirales_polymerase"
    "PF04196 Bunya_RdRp Bunyavirus_polymerase"
    "PF01443 Viral_helicase1 Helicase"
    "PF00949 Peptidase_S7 Flavivirus_NS3_protease"
    "PF05409 Peptidase_C30 Coronavirus_3CL_main_protease"
    "PF00952 Bunya_nucleocap Bunyaviridae_nucleocapsid"
    "PF05505 Ebola_NP Filovirus_nucleoprotein"
    "PF00506 Flu_NP Orthomyxoviridae_nucleoprotein"
    "PF00937 CoV_nucleocap Coronavirus_nucleocapsid"
    "PF03557 Bunya_G1 Bunyavirus_glycoprotein_G1"
    "PF03563 Bunya_G2 Bunyavirus_glycoprotein_G2"
    "PF01600 CoV_S1 Coronavirus_spike_S1"
    "PF01601 CoV_S2 Coronavirus_spike_S2"
    "PF00509 Hemagglutinin Influenza_HA"
    "PF01107 MP Plant_viral_movement_protein"
    "PF05065 Phage_capsid Phage_capsid"
)

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

: > "$OUTFILE"
for entry in "${PROFILES[@]}"; do
    acc="$(echo "$entry" | awk '{print $1}')"
    expected_name="$(echo "$entry" | awk '{print $2}')"
    role="$(echo "$entry" | awk '{print $3}')"
    url="https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/${acc}?annotation=hmm"
    echo "[$acc] fetching ($expected_name — $role) ..."
    if ! curl -sf "$url" -o "$TMP/$acc.gz"; then
        echo "  FAILED — skipping $acc"
        continue
    fi
    gunzip -c "$TMP/$acc.gz" > "$TMP/$acc.hmm"
    # Validate the NAME line matches what we expected, just to catch
    # accession typos before they propagate into config snippets.
    actual_name="$(grep -m1 '^NAME' "$TMP/$acc.hmm" | awk '{print $2}' || true)"
    if [[ -z "$actual_name" ]]; then
        echo "  WARN — no NAME line in $acc.hmm; concatenating anyway"
    elif [[ "$actual_name" != "$expected_name" ]]; then
        echo "  NOTE — NAME is '$actual_name' (expected '$expected_name'); using actual"
    fi
    cat "$TMP/$acc.hmm" >> "$OUTFILE"
done

echo "Bundled: $OUTFILE"
n_profiles=$(grep -c '^NAME' "$OUTFILE")
echo "  $n_profiles profile(s)"

echo "Indexing with hmmpress ..."
hmmpress -f "$OUTFILE"
echo "Done. Index files (.h3f .h3i .h3m .h3p) are next to $OUTFILE."
