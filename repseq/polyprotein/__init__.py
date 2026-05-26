"""Polyprotein cutting: slice a polyprotein CDS into its mature peptides.

For viruses that express their CDS as one giant precursor protein that is
post-translationally cleaved into mature peptides — picornavirus
P1/P2/P3, coronavirus ORF1ab NSP1..16, flavivirus polyprotein —
representative-level analysis is far more biologically informative on the
mature peptides than on the polyprotein as a whole.

This module is **augment-only**: it does not affect clustering,
representative selection, or the whole-genome tree. After representatives
are elected, each declared :class:`PolyproteinSpec` is applied to every
representative; the spec's parent CDS is identified by counting peptide-
HMM hits, and the protein is sliced into mature peptides using one of
three cut strategies (``boundary``, ``bisect``, ``motif``). The peptides
are emitted as accessory artifacts (per-peptide FASTAs + an audit TSV)
alongside the existing per-protein / extra-protein outputs.

See :mod:`repseq.polyprotein.specs` for the config-to-dataclass parser
and :mod:`repseq.polyprotein.slicer` for the cut math.
"""

from .specs import (
    PeptideSpec,
    PolyproteinSpec,
    collect_polyprotein_specs,
)
from .slicer import (
    SlicedPeptide,
    slice_polyprotein,
)

__all__ = [
    "PeptideSpec",
    "PolyproteinSpec",
    "SlicedPeptide",
    "collect_polyprotein_specs",
    "slice_polyprotein",
]
