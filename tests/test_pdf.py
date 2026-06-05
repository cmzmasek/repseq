"""Graphical tree-figure rendering (2H): phylo/pdf.py.

The renderer reconstructs everything (leaf labels, colours, the genus legend,
internal ranks, support values) from the phyloXML alone, so these tests build a
real repseq phyloXML via ``write_phyloxml`` (with a colour scheme so the leaves
carry ``style:font_color`` + ``repseq:genus`` properties) and then exercise the
parser/draw path. matplotlib is in the test environment, so the happy path
actually rasterises; the missing-matplotlib branch is exercised by monkeypatch.
"""
from __future__ import annotations

from pathlib import Path

from Bio import Phylo

from repseq.models import Sequence, SequenceType, TaxonomyInfo
from repseq.phylo.coloring import build_color_scheme
from repseq.phylo.phyloxml_writer import write_phyloxml
from repseq.phylo import pdf as pdfmod


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _make_seq(seq_id, genus="Orthohantavirus", subfamily="Mammantavirinae") -> Sequence:
    return Sequence(
        id=seq_id,
        header=f"{seq_id} {genus} fake header",
        sequence="ACGT" * 10,
        seq_type=SequenceType.NUCLEOTIDE,
        accession=seq_id,
        organism=f"{genus} sp.",
        taxonomy=TaxonomyInfo(
            taxid=1980456,
            species=f"{genus} species",
            genus=genus,
            family="Hantaviridae",
            lineage={"subfamily": subfamily},
        ),
    )


def _write_phyloxml(tmp_path: Path, reps, *, name="x_tree.xml", color=True) -> Path:
    """Produce a real repseq phyloXML at ``tmp_path/name`` from a 3-leaf tree."""
    newick = tmp_path / "in.nwk"
    # 0.95 internal support exercises the confidence (branch-label) path.
    newick.write_text("((S0001:0.1,S0002:0.1)0.95:0.2,S0003:0.3);")
    id_map = {f"S{i + 1:04d}": r.id for i, r in enumerate(reps)}
    out = tmp_path / name
    # write_phyloxml builds its own colour scheme from cfg when none is
    # passed, so disabling colour means telling it via cfg, not just None.
    cfg: dict = {} if color else {"phylo": {"coloring": {"enabled": False}}}
    color_scheme = build_color_scheme(reps, cfg) if color else None
    write_phyloxml(
        newick, out, reps, id_map,
        cfg=cfg, prefix="x", alphabet="nucleotide",
        msa_tool="MAFFT", msa_version="v7.520",
        tree_tool="FastTree", tree_version="2.1.11",
        model="GTR", ufboot=None,
        color_scheme=color_scheme,
    )
    return out


def _reps():
    # Two genera so the reconstructed legend has >1 entry.
    return [
        _make_seq("A", genus="Orthohantavirus"),
        _make_seq("B", genus="Orthohantavirus"),
        _make_seq("C", genus="Mammarenavirus", subfamily="Arenavirinae"),
    ]


# ---------------------------------------------------------------------------
# Dependency check
# ---------------------------------------------------------------------------

def test_matplotlib_available_in_test_env():
    # The test environment ships matplotlib, so the reason must be None.
    assert pdfmod.matplotlib_unavailable_reason() is None


# ---------------------------------------------------------------------------
# phyloXML clade readers
# ---------------------------------------------------------------------------

def test_clade_readers_from_phyloxml(tmp_path):
    out = _write_phyloxml(tmp_path, _reps())
    tree = Phylo.read(str(out), "phyloxml")

    leaves = tree.get_terminals()
    assert leaves, "expected leaves in the parsed tree"
    # Every leaf carries a font colour and a genus property.
    for leaf in leaves:
        assert pdfmod._clade_property(leaf, "style:font_color")
        assert pdfmod._clade_property(leaf, "repseq:genus")
        # _internal_label returns the full leaf name on terminals.
        assert pdfmod._internal_label(leaf) == (leaf.name or "")

    # The internal node from the 0.95 Newick label has a numeric confidence.
    internals = [c for c in tree.find_clades() if not c.is_terminal()]
    confs = [pdfmod._confidence_value(c) for c in internals]
    assert any(c is not None and c >= 50 for c in confs)


def test_build_color_maps_reconstructs_from_leaves(tmp_path):
    out = _write_phyloxml(tmp_path, _reps())
    tree = Phylo.read(str(out), "phyloxml")
    display_to_color, genus_to_color, subfamily_to_genera = pdfmod._build_color_maps(tree)

    # Every leaf label got a colour, and both genera appear in the legend map.
    assert display_to_color
    assert all(v.startswith("#") for v in display_to_color.values())
    assert "Orthohantavirus" in genus_to_color
    assert "Mammarenavirus" in genus_to_color
    # Subfamily → genera grouping is reconstructed from the leaf properties.
    assert "Mammantavirinae" in subfamily_to_genera
    assert "Orthohantavirus" in subfamily_to_genera["Mammantavirinae"]


def test_build_color_maps_empty_without_colour_scheme(tmp_path):
    out = _write_phyloxml(tmp_path, _reps(), color=False)
    tree = Phylo.read(str(out), "phyloxml")
    display_to_color, genus_to_color, _ = pdfmod._build_color_maps(tree)
    # No style:font_color properties → no colours reconstructed.
    assert display_to_color == {}
    assert genus_to_color == {}


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------

def test_render_one_writes_pdf_and_png(tmp_path):
    out = _write_phyloxml(tmp_path, _reps())
    created = pdfmod.render_one(out, want_png=True)
    pdf = out.with_suffix(".pdf")
    png = out.with_suffix(".png")
    assert pdf in created and png in created
    assert pdf.exists() and pdf.stat().st_size > 0
    assert png.exists() and png.stat().st_size > 0


def test_render_one_pdf_only(tmp_path):
    out = _write_phyloxml(tmp_path, _reps())
    created = pdfmod.render_one(out, want_png=False)
    assert created == [out.with_suffix(".pdf")]
    assert out.with_suffix(".pdf").exists()
    assert not out.with_suffix(".png").exists()


def test_render_one_skips_tree_with_one_leaf(tmp_path):
    # A degenerate <2-leaf tree returns no figure (and no files).
    out = tmp_path / "tiny_tree.xml"
    out.write_text(
        '<phyloxml xmlns="http://www.phyloxml.org"><phylogeny rooted="true">'
        "<clade><name>only</name></clade></phylogeny></phyloxml>"
    )
    created = pdfmod.render_one(out)
    assert created == []
    assert not out.with_suffix(".pdf").exists()


def test_render_tree_pdfs_batch(tmp_path):
    a = _write_phyloxml(tmp_path, _reps(), name="a_tree.xml")
    b = _write_phyloxml(tmp_path, _reps(), name="b_tree.xml")
    created, skipped, failures = pdfmod.render_tree_pdfs([a, b])
    assert skipped is None
    assert failures == []
    # Two trees × (pdf + png) = four files.
    assert len(created) == 4
    assert a.with_suffix(".pdf").exists()
    assert b.with_suffix(".png").exists()


def test_render_tree_pdfs_empty_input():
    assert pdfmod.render_tree_pdfs([]) == ([], None, [])


def test_render_tree_pdfs_soft_skips_without_matplotlib(tmp_path, monkeypatch):
    out = _write_phyloxml(tmp_path, _reps())
    monkeypatch.setattr(
        pdfmod, "matplotlib_unavailable_reason",
        lambda: "matplotlib missing (test)",
    )
    created, skipped, failures = pdfmod.render_tree_pdfs([out])
    assert created == []
    assert skipped == "matplotlib missing (test)"
    assert failures == []
    # Nothing rendered.
    assert not out.with_suffix(".pdf").exists()


def test_render_tree_pdfs_collects_per_tree_failures(tmp_path):
    good = _write_phyloxml(tmp_path, _reps(), name="good_tree.xml")
    bad = tmp_path / "bad_tree.xml"
    bad.write_text("not valid phyloxml <<<")
    created, skipped, failures = pdfmod.render_tree_pdfs([good, bad])
    assert skipped is None
    # The good tree still rendered; the bad one is recorded, not raised.
    assert good.with_suffix(".pdf").exists()
    assert len(failures) == 1
    assert failures[0][0] == bad
