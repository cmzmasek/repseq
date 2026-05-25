"""Optional phylogeny step: MAFFT alignment + FastTree → phyloXML.

The orchestrator (``pipeline.run_phylogeny``) is the public surface; the
backend wrappers (``mafft``, ``fasttree``) are kept as separate modules
so they can be mocked in tests and replaced individually if we ever add
alternatives.
"""

from .pipeline import PhyloError, run_phylogeny
from .per_protein import run_per_protein_phylogeny, run_per_segment_phylogeny

__all__ = [
    "PhyloError",
    "run_phylogeny",
    "run_per_protein_phylogeny",
    "run_per_segment_phylogeny",
]

