"""End-to-end integration test through the REAL external tools.

Everything else in the suite mocks mafft / FastTree / mmseqs. This test runs
the actual binaries on a tiny bundled FASTA, because real-tool runs are what
caught the bugs the mocked unit tests passed (MAD rooting, the monophyly
miscounts). It's **opt-in** — registered under the ``integration`` marker and
skipped by default (`pytest -m "not integration"` via pyproject), and skipped
outright when the tools aren't on PATH. Run it with ``pytest -m integration``.
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

pytestmark = pytest.mark.integration

_TOOLS = ("mafft", "FastTree", "mmseqs")
_DATA = Path(__file__).parent / "data" / "integration_nt.fasta"


@pytest.mark.skipif(
    not all(shutil.which(t) for t in _TOOLS),
    reason="requires mafft, FastTree and mmseqs on PATH",
)
def test_global_nucleotide_phylo_end_to_end(tmp_path):
    cfg = tmp_path / "cfg.yaml"
    cfg.write_text(
        "clustering:\n"
        "  backend: mmseqs2\n"
        "  alphabet_for_clustering: nucleotide\n"
        "  mmseqs2_mode: easy-cluster\n"
        "hmm:\n"
        "  enabled: false\n"
        "phylo:\n"
        "  tool: fasttree\n"
        "segmented:\n"
        "  enabled: false\n"
    )
    out = tmp_path / "out"
    proc = subprocess.run(
        [
            "repseq", "global",
            "-c", str(cfg), "-i", str(_DATA), "-o", str(out),
            "-T", "0.90", "--no-resolve", "--phylo",
        ],
        capture_output=True, text=True, timeout=600,
    )
    assert proc.returncode == 0, f"repseq failed:\n{proc.stderr}\n{proc.stdout}"

    names = {p.name for p in out.iterdir() if p.is_file()}
    # Clustering produced representatives, and the real MAFFT→FastTree→phyloXML
    # path produced a tree + alignment.
    assert any("representative" in n and n.endswith(".fasta") for n in names), names
    assert any(n.endswith("_tree.xml") for n in names), names
    assert any(n.endswith("_msa.fasta") for n in names), names
    # The new reporting sweeps fired too.
    assert any(n.endswith("_report.html") for n in names), names
    assert any(n.endswith("_summary.md") for n in names), names

    # The phyloXML is well-formed and has the expected leaf count (3 clusters).
    import xml.etree.ElementTree as ET
    tree_xml = next(out / n for n in names if n.endswith("_tree.xml"))
    root = ET.parse(tree_xml).getroot()
    ns = "{http://www.phyloxml.org}"
    leaves = [
        c for c in root.iter(f"{ns}clade")
        if not c.findall(f"{ns}clade") and c.find(f"{ns}name") is not None
    ]
    assert len(leaves) == 3, f"expected 3 representative leaves, got {len(leaves)}"
