"""Plain-English "what is this tree based on" descriptions.

Every tree repseq writes (the whole-genome tree, each per-marker tree,
the per-segment nucleotide trees, the polyprotein-peptide trees, and the
pre-cluster overview tree) is inferred from a *different* slice of the
data. The phyloXML provenance fields record which tools and versions ran,
but on their own they never tell a reader **what the leaves are** (an
isolate? a sequence?) or **what sequence data the alignment was built
from** (a whole genome? one protein? a Spike+Nucleocapsid concatenation?
one segment's nucleotides?).

:func:`describe_tree_basis` closes that gap. It returns a one-sentence,
repseq-agnostic description of a tree's biological substrate plus a small
dict of machine-readable properties. The sentence leads the phyloXML
``<phylogeny><description>``; the properties are emitted as
phylogeny-level ``repseq:`` ``<property>`` elements so scripts can read
the basis without regexing prose.

The module is deliberately dependency-free (pure string assembly) so it
can be unit-tested in isolation and imported from both the phyloXML
writer and the summary renderer without cycles.

Substrate codes (stable, machine-readable — the ``substrate`` property):

* ``single_marker``        — one marker protein, amino-acid (2E, non-seg)
* ``marker_protein_concat``— concat of >1 marker protein, AA (2E, non-seg,
  ``concatenate_markers``)
* ``segment_marker_concat``— per-isolate concat of one marker per segment,
  AA (2E, segmented protein)
* ``genome_nt``            — whole-genome nucleotide (2E, non-seg NT)
* ``segment_nt_concat``    — per-isolate concat of segment NT (2E, seg NT)
* ``supermatrix``          — partitioned supermatrix of per-family MSAs (2E)
* ``marker``               — a single marker CDS (2F per-protein tree)
* ``segment_nt``           — a single segment's nucleotides (2H)
* ``accessory_protein``    — an ``extra_protein`` CDS, AA
* ``peptide``              — a sliced polyprotein peptide, AA
* ``overview``             — every post-QC sequence (pre-cluster overview)
"""
from __future__ import annotations

from typing import Optional

# Role tokens accepted by describe_tree_basis. Kept as plain strings (not
# an Enum) so callers stay terse and tests read naturally.
ROLES = (
    "genome",
    "genome_partitioned",
    "marker",
    "segment_nt",
    "extra_protein",
    "peptide",
    "pre_cluster",
)


def describe_tree_basis(
    role: str,
    *,
    alphabet: str,
    segmented: bool,
    markers: Optional[str] = None,
    families: Optional[list[str]] = None,
    family: Optional[str] = None,
    segment: Optional[str] = None,
    architecture: Optional[str] = None,
    parent: Optional[str] = None,
    concat_markers: bool = False,
) -> tuple[str, dict[str, str]]:
    """Describe the biological substrate of one tree.

    Args:
        role: one of :data:`ROLES`.
        alphabet: ``"protein"`` or ``"nucleotide"`` (the alphabet the tree
            was actually inferred on).
        segmented: whether the run is segmented (controls the leaf unit —
            isolate vs sequence).
        markers: a human marker summary (``_markers_summary`` output), e.g.
            ``"L:RdRp, M:GPC, S:N"`` (segmented) or ``"Spike, Nucleocapsid"``
            (non-segmented). Used by the ``genome``/``genome_partitioned``
            roles.
        families: per-family labels for the partitioned supermatrix.
        family: the marker / accessory-protein / peptide label (per-protein,
            extra-protein, peptide roles).
        segment: the segment name (segment_nt role).
        architecture: the HMM domain-architecture token alternatives
            (``"A--B OR C--D"``) for a marker/accessory tree, when known.
        parent: the parent polyprotein spec name (peptide role).
        concat_markers: whether ``clustering.concatenate_markers`` is on
            (distinguishes ``single_marker`` from ``marker_protein_concat``
            for the non-segmented genome tree).

    Returns:
        ``(sentence, properties)`` — ``sentence`` is the plain-English basis
        line; ``properties`` is a dict with keys ``tree_basis``,
        ``analysis_mode``, ``substrate``, ``alphabet``, ``leaf_unit``.
    """
    aa = alphabet == "protein"
    alpha_word = "amino-acid" if aa else "nucleotide"
    alpha_code = "amino_acid" if aa else "nucleotide"
    leaf_unit = "isolate" if segmented else "sequence"
    leaf_phrase = f"each leaf is one representative {leaf_unit}"

    def _with_markers(text: str) -> str:
        return f"{text} ({markers})" if markers else text

    if role == "genome":
        if segmented:
            if aa:
                substrate = "segment_marker_concat"
                what = _with_markers(
                    "a per-isolate concatenation of one marker protein "
                    "per segment"
                ) + "; amino-acid"
            else:
                substrate = "segment_nt_concat"
                what = (
                    "a per-isolate concatenation of the segment nucleotide "
                    "sequences"
                )
        else:
            if aa and concat_markers:
                substrate = "marker_protein_concat"
                what = _with_markers(
                    "a concatenation of marker proteins"
                ) + "; amino-acid"
            elif aa:
                substrate = "single_marker"
                what = f"the {markers or 'selected'} marker protein (amino-acid)"
            else:
                substrate = "genome_nt"
                what = "the whole-genome nucleotide sequence"
        sentence = f"This tree is based on {what}; {leaf_phrase}."

    elif role == "genome_partitioned":
        substrate = "supermatrix"
        fam_str = ", ".join(families) if families else (markers or "")
        fam_clause = f" ({fam_str})" if fam_str else ""
        what = (
            f"a partitioned supermatrix of separately-aligned marker "
            f"families{fam_clause}; amino-acid, one substitution model "
            f"per partition"
        )
        sentence = f"This tree is based on {what}; {leaf_phrase}."

    elif role == "marker":
        substrate = "marker"
        fam = family or "the marker"
        what = f"the {fam} marker protein only ({alpha_word})"
        if architecture:
            what += f"; HMM architecture {architecture}"
        sentence = f"This tree is based on {what}; {leaf_phrase}."

    elif role == "segment_nt":
        substrate = "segment_nt"
        seg = segment or "the"
        what = f"the {seg} segment nucleotide sequence only"
        sentence = f"This tree is based on {what}; {leaf_phrase}."

    elif role == "extra_protein":
        substrate = "accessory_protein"
        fam = family or "an accessory protein"
        what = f"the accessory protein {fam} ({alpha_word})"
        if architecture:
            what += f"; HMM architecture {architecture}"
        sentence = (
            f"This tree is based on {what}; {leaf_phrase}. Accessory "
            f"proteins do not drive clustering or the whole-genome tree."
        )

    elif role == "peptide":
        substrate = "peptide"
        fam = family or "a peptide"
        origin = f"the {parent} polyprotein" if parent else "a polyprotein"
        what = f"the {fam} mature peptide ({alpha_word}), sliced from {origin}"
        sentence = f"This tree is based on {what}; {leaf_phrase}."

    elif role == "pre_cluster":
        substrate = "overview"
        sentence = (
            f"This tree is an overview of every post-QC {leaf_unit} "
            f"({alpha_word}) — the full input diversity before "
            f"representative selection; the elected representatives are "
            f"marked with a [repr] prefix."
        )

    else:
        substrate = "unknown"
        sentence = f"This tree's basis is unspecified (role={role})."

    properties = {
        "tree_basis": sentence,
        "analysis_mode": "segmented" if segmented else "non_segmented",
        "substrate": substrate,
        "alphabet": alpha_code,
        "leaf_unit": leaf_unit,
    }
    return sentence, properties
