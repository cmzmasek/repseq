"""PhyloXML writer: rich per-leaf annotation, schema element order,
confidence normalisation, phylogeny header, ladderize, multi-sequence
emission (markers first), summary property lists.

The writer is XML, so most tests just parse the output with
``ElementTree`` and inspect the structure. The MAFFT/IQ-TREE/FastTree
binary calls are not exercised here — those are mocked in tests/test_phylo.py.
"""
from __future__ import annotations

import xml.etree.ElementTree as ET
from pathlib import Path

from repseq.models import Sequence, SequenceType, TaxonomyInfo
from repseq.phylo.phyloxml_writer import (
    _confidence_type_for,
    _normalize_confidence,
    write_phyloxml,
)


NS = "http://www.phyloxml.org"


def _ns(tag: str) -> str:
    return f"{{{NS}}}{tag}"


def _make_seq(seq_id, **overrides) -> Sequence:
    base = dict(
        id=seq_id,
        header=f"{seq_id} fake header",
        sequence="ACGT" * 10,
        seq_type=SequenceType.NUCLEOTIDE,
        accession=seq_id,
        organism="Hantaan virus",
        host="Apodemus agrarius",
        strain="76-118",
        collection_date="1976-04",
        country="South Korea",
        isolate_id="76-118",
        taxonomy=TaxonomyInfo(
            taxid=1980456,
            species="Hantaan orthohantavirus",
            genus="Orthohantavirus",
            family="Hantaviridae",
            lineage={"subfamily": "Mammantavirinae", "subgenus": "Hantaanvirus"},
        ),
    )
    base.update(overrides)
    return Sequence(**base)


def _write_newick(path: Path, body: str) -> None:
    path.write_text(body)


# ---------------------------------------------------------------------------
# Confidence helpers (unit tests, no I/O)
# ---------------------------------------------------------------------------

def test_normalize_confidence_sh_like_rescales():
    assert _normalize_confidence(0.95, "sh_like") == 95
    assert _normalize_confidence(0.5, "sh_like") == 50
    assert _normalize_confidence(1.0, "sh_like") == 100


def test_normalize_confidence_ufboot_passthrough():
    assert _normalize_confidence(95.0, "ufboot") == 95
    assert _normalize_confidence(50.6, "ufboot") == 51  # rounds


def test_normalize_confidence_clamps_to_0_100():
    assert _normalize_confidence(120.0, "ufboot") == 100
    assert _normalize_confidence(-5.0, "ufboot") == 0


def test_normalize_confidence_none_returns_none():
    assert _normalize_confidence(None, "ufboot") is None


def test_confidence_type_for_default_iqtree_is_ufboot():
    assert _confidence_type_for("iqtree", None) == "ufboot"
    assert _confidence_type_for("IQ-TREE", None) == "ufboot"


def test_confidence_type_for_default_fasttree_is_sh_like():
    assert _confidence_type_for("fasttree", None) == "sh_like"
    assert _confidence_type_for("FastTree", None) == "sh_like"


def test_confidence_type_for_override():
    assert _confidence_type_for("fasttree", "bootstrap") == "bootstrap"
    assert _confidence_type_for("iqtree", "sh_alrt") == "sh_alrt"


def test_confidence_type_for_auto_falls_through():
    assert _confidence_type_for("iqtree", "auto") == "ufboot"


def test_confidence_type_for_unknown_tool_defaults_to_bootstrap():
    assert _confidence_type_for("raxml", None) == "bootstrap"


# ---------------------------------------------------------------------------
# Full write: end-to-end through a tiny Newick file
# ---------------------------------------------------------------------------

def _run_write(tmp_path, reps, *, cfg=None, tree_tool="FastTree",
               model="GTR", ufboot=None, alphabet="nucleotide"):
    newick = tmp_path / "tree.nwk"
    # Build a Newick where leaf labels are the short ids that id_map
    # will resolve back. Three taxa is the minimum the orchestrator
    # would let through, so use that.
    _write_newick(newick, "((S0001:0.1,S0002:0.1)0.95:0.2,S0003:0.3);")
    id_map = {f"S{i + 1:04d}": rep.id for i, rep in enumerate(reps)}
    out = tmp_path / "tree.xml"
    write_phyloxml(
        newick, out, reps, id_map,
        cfg=cfg or {},
        prefix="test",
        alphabet=alphabet,
        msa_tool="MAFFT",
        msa_version="v7.520",
        tree_tool=tree_tool,
        tree_version="2.1.11",
        model=model,
        ufboot=ufboot,
    )
    return out


def test_writes_well_formed_xml_with_required_elements(tmp_path):
    reps = [_make_seq("A"), _make_seq("B"), _make_seq("C")]
    out = _run_write(tmp_path, reps)
    root = ET.parse(out).getroot()
    assert root.tag == _ns("phyloxml")
    phylogeny = root.find(_ns("phylogeny"))
    assert phylogeny is not None
    assert phylogeny.find(_ns("name")) is not None
    assert phylogeny.find(_ns("description")) is not None


def test_leaf_carries_taxonomy_block(tmp_path):
    reps = [_make_seq("A"), _make_seq("B"), _make_seq("C")]
    out = _run_write(tmp_path, reps)
    root = ET.parse(out).getroot()
    # Find any leaf-level taxonomy. We expect 3.
    taxonomies = root.findall(
        f".//{_ns('clade')}/{_ns('taxonomy')}"
    )
    assert len(taxonomies) == 3
    for tax in taxonomies:
        id_el = tax.find(_ns("id"))
        sci_name = tax.find(_ns("scientific_name"))
        assert id_el is not None
        assert id_el.get("provider") == "ncbi"
        assert id_el.text == "1980456"
        assert sci_name is not None
        assert sci_name.text == "Hantaan orthohantavirus"


def test_leaf_carries_one_dna_sequence_when_no_proteins(tmp_path):
    """Bare non-segmented input (no .proteins): one <sequence type="dna">
    per leaf, accession + name as before."""
    reps = [_make_seq("ACC1"), _make_seq("ACC2"), _make_seq("ACC3")]
    out = _run_write(tmp_path, reps)
    root = ET.parse(out).getroot()
    seqs = root.findall(f".//{_ns('clade')}/{_ns('sequence')}")
    assert len(seqs) == 3
    for s in seqs:
        assert s.get("type") == "dna"
        acc = s.find(_ns("accession"))
        assert acc is not None
        assert acc.get("source") == "ncbi"
        assert acc.text in ("ACC1", "ACC2", "ACC3")
        name = s.find(_ns("name"))
        assert name is not None
        assert "fake header" in name.text


def test_leaf_with_proteins_emits_multiple_sequences_markers_first(tmp_path):
    """Non-segmented input with .proteins populated: writer emits
    one <sequence type="protein"> per CDS (marker first), then one
    <sequence type="dna"> for the nuc accession."""
    reps = [_make_seq("ACC1"), _make_seq("ACC2"), _make_seq("ACC3")]
    reps[0].proteins = [
        {"protein_id": "P_aux", "product": "auxiliary protein",
         "length": 50, "sequence": "M" * 50},
        {"protein_id": "P_main", "product": "polymerase",
         "length": 100, "sequence": "M" * 100},
    ]
    reps[0].marker_protein_ids = ["P_main"]
    out = _run_write(tmp_path, reps)
    root = ET.parse(out).getroot()
    # Find the clade that carries ACC1 (via the dna <sequence>).
    target = None
    for clade in root.iter(_ns("clade")):
        for s in clade.findall(_ns("sequence")):
            acc = s.find(_ns("accession"))
            if acc is not None and acc.text == "ACC1":
                target = clade
                break
        if target is not None:
            break
    assert target is not None
    seqs = target.findall(_ns("sequence"))
    # 2 proteins + 1 nuc = 3 sequences on ACC1's leaf.
    assert len(seqs) == 3
    # Order: marker protein first (P_main), then other protein (P_aux),
    # then the nuc accession.
    types = [s.get("type") for s in seqs]
    accs = [s.find(_ns("accession")).text for s in seqs]
    assert types == ["protein", "protein", "dna"]
    assert accs == ["P_main", "P_aux", "ACC1"]
    # And the <name> elements track the product / header strings.
    names = [s.find(_ns("name")).text for s in seqs]
    assert names[0] == "polymerase"
    assert names[1] == "auxiliary protein"
    assert "fake header" in names[2]


def test_leaf_emits_three_summary_property_lists(tmp_path):
    """repseq:nuc_acc / protein_acc / protein_names list the same
    data as the <sequence> elements, comma-joined — for renderers
    that only display the first <sequence>."""
    reps = [_make_seq("ACC1"), _make_seq("ACC2"), _make_seq("ACC3")]
    reps[0].proteins = [
        {"protein_id": "P_aux", "product": "auxiliary protein",
         "length": 50, "sequence": "M" * 50},
        {"protein_id": "P_main", "product": "polymerase",
         "length": 100, "sequence": "M" * 100},
    ]
    reps[0].marker_protein_ids = ["P_main"]
    out = _run_write(tmp_path, reps)
    root = ET.parse(out).getroot()
    target = None
    for clade in root.iter(_ns("clade")):
        accs = [a.text for a in clade.findall(f"{_ns('sequence')}/{_ns('accession')}")]
        if "ACC1" in accs:
            target = clade
            break
    assert target is not None
    props = {
        p.get("ref"): p.text
        for p in target.findall(_ns("property"))
    }
    assert props["repseq:nuc_acc"] == "ACC1"
    # Marker-first ordering replicated in the summary lists.
    assert props["repseq:protein_acc"] == "P_main, P_aux"
    assert props["repseq:protein_names"] == "polymerase, auxiliary protein"


def test_segmented_concat_emits_one_sequence_per_segment_and_protein(tmp_path):
    """A CONCAT leaf with three underlying segments (L/M/S) and four
    proteins (one each on L and M, two on S) emits 4 + 3 = 7
    <sequence> elements — markers first."""
    from repseq.models import Sequence, SequenceType
    seg_L = Sequence(
        id="seg_L", header="seg_L", sequence="A" * 100,
        seq_type=SequenceType.NUCLEOTIDE, accession="L_ACC",
        segment="L", proteins=[
            {"protein_id": "L_pol", "product": "L polymerase",
             "length": 50, "sequence": "M" * 50},
        ],
    )
    seg_M = Sequence(
        id="seg_M", header="seg_M", sequence="C" * 100,
        seq_type=SequenceType.NUCLEOTIDE, accession="M_ACC",
        segment="M", proteins=[
            {"protein_id": "M_gly", "product": "glycoprotein",
             "length": 40, "sequence": "M" * 40},
        ],
    )
    seg_S = Sequence(
        id="seg_S", header="seg_S", sequence="G" * 100,
        seq_type=SequenceType.NUCLEOTIDE, accession="S_ACC",
        segment="S", proteins=[
            {"protein_id": "S_N", "product": "nucleoprotein",
             "length": 30, "sequence": "M" * 30},
            {"protein_id": "S_NSs", "product": "NSs",
             "length": 20, "sequence": "M" * 20},
        ],
    )
    concat = _make_seq("CONCAT|iso1")
    concat.id = "CONCAT|iso1"
    concat.concat_segments = [seg_L, seg_M, seg_S]
    concat.marker_protein_ids = ["L_pol", "M_gly", "S_N"]  # one per segment
    reps = [concat, _make_seq("B"), _make_seq("C")]
    out = _run_write(tmp_path, reps)
    root = ET.parse(out).getroot()
    # Find the CONCAT clade.
    target = None
    for clade in root.iter(_ns("clade")):
        accs = [a.text for a in clade.findall(f"{_ns('sequence')}/{_ns('accession')}")]
        if "L_ACC" in accs:
            target = clade
            break
    assert target is not None
    seqs = target.findall(_ns("sequence"))
    assert len(seqs) == 7  # 4 proteins + 3 nucs
    types = [s.get("type") for s in seqs]
    accs = [s.find(_ns("accession")).text for s in seqs]
    # Markers first (L_pol, M_gly, S_N), then non-marker protein (S_NSs),
    # then nucs in segment order.
    assert types == ["protein"] * 4 + ["dna"] * 3
    assert accs[:4] == ["L_pol", "M_gly", "S_N", "S_NSs"]
    assert accs[4:] == ["L_ACC", "M_ACC", "S_ACC"]
    # Summary properties.
    props = {p.get("ref"): p.text for p in target.findall(_ns("property"))}
    assert props["repseq:nuc_acc"] == "L_ACC, M_ACC, S_ACC"
    assert props["repseq:protein_acc"] == "L_pol, M_gly, S_N, S_NSs"
    assert props["repseq:protein_names"] == (
        "L polymerase, glycoprotein, nucleoprotein, NSs"
    )


def test_sequence_name_strips_ncbi_virus_leading_pipe(tmp_path):
    """NCBI Virus headers carry pipe-separated metadata after the
    accession (``NC_078889.1 |species|host|...|segment``); the parsed
    ``description`` then starts with ``|``. The writer must strip that
    so leaf <name> doesn't show as |Turlock orthobunyavirus segment L..."""
    reps = [
        _make_seq(
            "ACC1",
            description="|Turlock orthobunyavirus segment L|Turlock virus|Orthobunyavirus turlockense||L",
        ),
        _make_seq("ACC2"),
        _make_seq("ACC3"),
    ]
    out = _run_write(tmp_path, reps)
    root = ET.parse(out).getroot()
    # Find ACC1's <sequence type="dna"><name> — it should NOT start with '|'.
    for clade in root.iter(_ns("clade")):
        for s in clade.findall(_ns("sequence")):
            acc = s.find(_ns("accession"))
            if acc is None or acc.text != "ACC1":
                continue
            name = s.find(_ns("name")).text
            assert not name.startswith("|"), f"leading pipe leaked: {name!r}"
            assert not name.endswith("|"), f"trailing pipe leaked: {name!r}"
            # The empty NCBI-Virus field (``||``) should have been
            # collapsed, not left as a double-pipe.
            assert "||" not in name


def test_summary_properties_omitted_when_no_proteins(tmp_path):
    """A bare non-segmented input shouldn't emit empty
    repseq:protein_acc or repseq:protein_names properties."""
    reps = [_make_seq("ACC1"), _make_seq("ACC2"), _make_seq("ACC3")]
    out = _run_write(tmp_path, reps)
    root = ET.parse(out).getroot()
    refs = {p.get("ref") for p in root.findall(f".//{_ns('property')}")}
    # nuc_acc is present (every leaf has an accession), but the
    # protein lists must not be emitted as empty.
    assert "repseq:nuc_acc" in refs
    assert "repseq:protein_acc" not in refs
    assert "repseq:protein_names" not in refs


def test_leaf_property_elements_use_repseq_namespace(tmp_path):
    reps = [_make_seq("A"), _make_seq("B"), _make_seq("C")]
    out = _run_write(tmp_path, reps)
    root = ET.parse(out).getroot()
    props = root.findall(f".//{_ns('clade')}/{_ns('property')}")
    assert len(props) > 0
    # Every property must be repseq:... — no foreign namespaces leaked.
    for p in props:
        assert p.get("ref", "").startswith("repseq:")
        assert p.get("datatype") == "xsd:string"
        assert p.get("applies_to") == "clade"
    refs = {p.get("ref") for p in props}
    expected = {
        "repseq:host", "repseq:collection_date", "repseq:country",
        "repseq:strain", "repseq:isolate_id", "repseq:year",
        "repseq:species", "repseq:subgenus", "repseq:genus",
        "repseq:subfamily", "repseq:family",
    }
    assert expected.issubset(refs)


def test_subgenus_property_emitted_on_leaves(tmp_path):
    """The full _TAX_RANKS ladder reaches leaf properties — subgenus
    (lineage-only, commonly meaningful for coronaviruses) included."""
    reps = [_make_seq("A"), _make_seq("B"), _make_seq("C")]
    out = _run_write(tmp_path, reps)
    root = ET.parse(out).getroot()
    subgenus = [
        p for p in root.findall(f".//{_ns('clade')}/{_ns('property')}")
        if p.get("ref") == "repseq:subgenus"
    ]
    assert len(subgenus) == 3  # one per leaf
    assert all(p.text == "Hantaanvirus" for p in subgenus)


def test_absent_subranks_omitted_from_leaf_properties(tmp_path):
    """suborder/order/subclass/class aren't in the fixture lineage, so
    they must not appear as empty stubs."""
    reps = [_make_seq("A"), _make_seq("B"), _make_seq("C")]
    out = _run_write(tmp_path, reps)
    root = ET.parse(out).getroot()
    refs = {
        p.get("ref")
        for p in root.findall(f".//{_ns('clade')}/{_ns('property')}")
    }
    for absent in ("repseq:suborder", "repseq:order",
                   "repseq:subclass", "repseq:class"):
        assert absent not in refs


def test_empty_metadata_omitted_from_properties(tmp_path):
    """A leaf with no host/country still validates; properties are dropped
    rather than emitted as empty stubs."""
    reps = [
        _make_seq("A", host=None, country=None, strain=None),
        _make_seq("B"),
        _make_seq("C"),
    ]
    out = _run_write(tmp_path, reps)
    root = ET.parse(out).getroot()
    # The first leaf should have no host/country/strain properties.
    clades = root.findall(f".//{_ns('clade')}/{_ns('sequence')}/..")
    # Find which clade carries ACC=A.
    target = None
    for clade in clades:
        seq = clade.find(_ns("sequence"))
        acc = seq.find(_ns("accession")) if seq is not None else None
        if acc is not None and acc.text == "A":
            target = clade
    assert target is not None
    refs = {p.get("ref") for p in target.findall(_ns("property"))}
    assert "repseq:host" not in refs
    assert "repseq:country" not in refs
    assert "repseq:strain" not in refs
    # But the ones that *are* set on this leaf should still be there.
    assert "repseq:species" in refs  # taxonomy still populated


def test_year_property_derived_from_collection_date(tmp_path):
    reps = [
        _make_seq("A", collection_date="04-Apr-1976"),
        _make_seq("B"),
        _make_seq("C"),
    ]
    out = _run_write(tmp_path, reps)
    root = ET.parse(out).getroot()
    years = [
        p.text for p in root.findall(f".//{_ns('property')}")
        if p.get("ref") == "repseq:year"
    ]
    assert "1976" in years


def test_phylogeny_name_contains_alphabet_and_tools(tmp_path):
    reps = [_make_seq("A"), _make_seq("B"), _make_seq("C")]
    out = _run_write(
        tmp_path, reps,
        tree_tool="IQ-TREE", model="LG+G4", alphabet="protein",
    )
    root = ET.parse(out).getroot()
    name = root.find(f"{_ns('phylogeny')}/{_ns('name')}").text
    assert "test" in name
    assert "protein" in name
    assert "MAFFT" in name
    assert "IQ-TREE" in name
    assert "LG+G4" in name


def test_phylogeny_description_records_versions_and_model(tmp_path):
    reps = [_make_seq("A"), _make_seq("B"), _make_seq("C")]
    out = _run_write(
        tmp_path, reps,
        tree_tool="IQ-TREE", model="MFP", ufboot=1000, alphabet="protein",
    )
    root = ET.parse(out).getroot()
    desc = root.find(f"{_ns('phylogeny')}/{_ns('description')}").text
    assert "repseq" in desc
    assert "MAFFT v7.520" in desc
    assert "IQ-TREE 2.1.11" in desc
    assert "model=MFP" in desc
    assert "bootstrap=1000" in desc


def test_phylogeny_description_records_no_bootstrap(tmp_path):
    reps = [_make_seq("A"), _make_seq("B"), _make_seq("C")]
    out = _run_write(tmp_path, reps, ufboot=0)
    root = ET.parse(out).getroot()
    desc = root.find(f"{_ns('phylogeny')}/{_ns('description')}").text
    assert "bootstrap=none" in desc


def _make_segmented_rep(iso_id: str, marker_product_by_seg: dict[str, str]):
    """Build a CONCAT rep with one segment per entry in
    ``marker_product_by_seg`` and exactly one protein per segment
    matching the requested product."""
    from repseq.models import Sequence, SequenceType
    segs = []
    marker_ids: list[str] = []
    for seg_name, product in marker_product_by_seg.items():
        pid = f"P_{iso_id}_{seg_name}"
        seg = Sequence(
            id=f"{iso_id}_{seg_name}", header=f"{iso_id}_{seg_name}",
            sequence="A" * 100, seq_type=SequenceType.NUCLEOTIDE,
            accession=f"{iso_id}_{seg_name}_acc", segment=seg_name,
            proteins=[
                {"protein_id": pid, "product": product,
                 "length": 50, "sequence": "M" * 50},
            ],
        )
        segs.append(seg)
        marker_ids.append(pid)
    concat = _make_seq(f"CONCAT|{iso_id}")
    concat.id = f"CONCAT|{iso_id}"
    concat.concat_segments = segs
    concat.marker_protein_ids = marker_ids
    return concat


def test_description_lists_segmented_markers_per_segment(tmp_path):
    """All three reps picked the same per-segment marker → one
    product per segment, no pipe-join."""
    reps = [
        _make_segmented_rep("iso1", {"L": "polymerase", "M": "glycoprotein",
                                      "S": "nucleoprotein"}),
        _make_segmented_rep("iso2", {"L": "polymerase", "M": "glycoprotein",
                                      "S": "nucleoprotein"}),
        _make_segmented_rep("iso3", {"L": "polymerase", "M": "glycoprotein",
                                      "S": "nucleoprotein"}),
    ]
    out = _run_write(tmp_path, reps, alphabet="protein")
    root = ET.parse(out).getroot()
    desc = root.find(f"{_ns('phylogeny')}/{_ns('description')}").text
    assert "markers=L:polymerase, M:glycoprotein, S:nucleoprotein" in desc


def test_description_pipe_joins_mixed_segmented_markers(tmp_path):
    """When different reps picked different products on the same
    segment, they're pipe-joined within that segment (alphabetised)."""
    reps = [
        _make_segmented_rep("iso1", {"L": "polymerase", "M": "glycoprotein",
                                      "S": "nucleoprotein"}),
        # iso2 has "RdRp" on segment L instead of "polymerase"
        _make_segmented_rep("iso2", {"L": "RdRp", "M": "glycoprotein",
                                      "S": "nucleoprotein"}),
        _make_segmented_rep("iso3", {"L": "polymerase", "M": "glycoprotein",
                                      "S": "nucleoprotein"}),
    ]
    out = _run_write(tmp_path, reps, alphabet="protein")
    root = ET.parse(out).getroot()
    desc = root.find(f"{_ns('phylogeny')}/{_ns('description')}").text
    # Alphabetical inside the segment: 'RdRp' before 'polymerase'? No —
    # default sort is ASCII so uppercase letters come first.
    assert "markers=L:RdRp|polymerase, M:glycoprotein, S:nucleoprotein" in desc


def test_description_lists_non_segmented_marker_products(tmp_path):
    """Non-segmented: flat unique-product list across all reps."""
    reps = [_make_seq("A"), _make_seq("B"), _make_seq("C")]
    reps[0].proteins = [
        {"protein_id": "P_pol_a", "product": "polymerase",
         "length": 100, "sequence": "M" * 100},
    ]
    reps[0].marker_protein_ids = ["P_pol_a"]
    reps[1].proteins = [
        {"protein_id": "P_pol_b", "product": "polymerase",
         "length": 100, "sequence": "M" * 100},
    ]
    reps[1].marker_protein_ids = ["P_pol_b"]
    reps[2].proteins = [
        {"protein_id": "P_ns5_c", "product": "NS5",
         "length": 80, "sequence": "M" * 80},
    ]
    reps[2].marker_protein_ids = ["P_ns5_c"]
    out = _run_write(tmp_path, reps, alphabet="protein")
    root = ET.parse(out).getroot()
    desc = root.find(f"{_ns('phylogeny')}/{_ns('description')}").text
    assert "markers=NS5, polymerase" in desc


def test_description_omits_markers_for_nt_runs(tmp_path):
    """No rep has marker_protein_ids set → no markers field."""
    reps = [_make_seq("A"), _make_seq("B"), _make_seq("C")]
    out = _run_write(tmp_path, reps, alphabet="nucleotide")
    root = ET.parse(out).getroot()
    desc = root.find(f"{_ns('phylogeny')}/{_ns('description')}").text
    assert "markers=" not in desc


def test_internal_confidence_normalized_to_integer(tmp_path):
    """The Newick has '0.95' on an internal node; FastTree default →
    sh_like → rescale to 95."""
    reps = [_make_seq("A"), _make_seq("B"), _make_seq("C")]
    out = _run_write(tmp_path, reps, tree_tool="FastTree")
    root = ET.parse(out).getroot()
    confs = root.findall(f".//{_ns('clade')}/{_ns('confidence')}")
    # At least one internal confidence.
    assert any(c.text == "95" for c in confs)
    # All confidence elements declare a type.
    for c in confs:
        assert c.get("type") == "sh_like"


def test_iqtree_confidence_type_attribute(tmp_path):
    reps = [_make_seq("A"), _make_seq("B"), _make_seq("C")]
    out = _run_write(tmp_path, reps, tree_tool="IQ-TREE", model="MFP")
    root = ET.parse(out).getroot()
    confs = root.findall(f".//{_ns('confidence')}")
    assert all(c.get("type") == "ufboot" for c in confs)


def test_confidence_type_override(tmp_path):
    cfg = {"phylo": {"phyloxml": {"confidence_type": "bootstrap"}}}
    reps = [_make_seq("A"), _make_seq("B"), _make_seq("C")]
    out = _run_write(tmp_path, reps, cfg=cfg, tree_tool="FastTree")
    root = ET.parse(out).getroot()
    confs = root.findall(f".//{_ns('confidence')}")
    assert all(c.get("type") == "bootstrap" for c in confs)


def test_schema_element_order_on_leaf(tmp_path):
    """PhyloXML schema requires: name → branch_length → confidence →
    taxonomy → sequence → property → child clades."""
    reps = [_make_seq("A"), _make_seq("B"), _make_seq("C")]
    out = _run_write(tmp_path, reps)
    root = ET.parse(out).getroot()
    # Pick the first terminal clade (S0001 → leaf, no child clades).
    leaves = [c for c in root.iter(_ns("clade")) if c.find(_ns("taxonomy")) is not None]
    assert leaves
    leaf = leaves[0]
    tag_order = [
        child.tag for child in leaf
        if child.tag in {
            _ns("name"), _ns("branch_length"), _ns("confidence"),
            _ns("taxonomy"), _ns("sequence"), _ns("property"),
        }
    ]
    expected_order = [
        _ns("name"), _ns("branch_length"), _ns("confidence"),
        _ns("taxonomy"), _ns("sequence"), _ns("property"),
    ]
    # Compress duplicates to validate ORDER, not COUNT.
    seen = []
    for t in tag_order:
        if not seen or seen[-1] != t:
            seen.append(t)
    # Every seen tag should appear in expected_order, in order.
    indices = [expected_order.index(t) for t in seen]
    assert indices == sorted(indices), f"order broken: {seen}"


def test_label_uses_configured_format(tmp_path):
    cfg = {"phylo": {"labeling": {"format": "{species}|{id}"}}}
    reps = [_make_seq("A"), _make_seq("B"), _make_seq("C")]
    out = _run_write(tmp_path, reps, cfg=cfg)
    root = ET.parse(out).getroot()
    names = [
        n.text for n in root.findall(f".//{_ns('clade')}/{_ns('name')}")
        if n.text and "|" in n.text
    ]
    assert any("Hantaan_orthohantavirus|A" in n for n in names)


def test_segmented_label_uses_strain(tmp_path):
    """When segmented mode is on, the segmented_format wins."""
    cfg = {
        "segmented": {"enabled": True},
        "phylo": {
            "labeling": {
                "format": "{species}|{id}",
                "segmented_format": "{species}|{strain}",
            },
        },
    }
    reps = [_make_seq("A"), _make_seq("B"), _make_seq("C")]
    out = _run_write(tmp_path, reps, cfg=cfg)
    root = ET.parse(out).getroot()
    names = [
        n.text for n in root.findall(f".//{_ns('clade')}/{_ns('name')}")
        if n.text and "|" in n.text
    ]
    # Strain "76-118" should appear, not the accession.
    assert any("76-118" in n for n in names)
    assert not any("|A|" in n for n in names)


def test_internal_clade_emits_name_and_taxonomy_when_lca_set(tmp_path):
    """When the pipeline pre-sets _lca_name / _lca_rank on internal
    clades, the writer emits <name> and <taxonomy><scientific_name>
    + <taxonomy><rank> on those clades."""
    from io import StringIO
    from Bio import Phylo
    reps = [_make_seq("A"), _make_seq("B"), _make_seq("C")]
    id_map = {f"S{i + 1:04d}": rep.id for i, rep in enumerate(reps)}
    tree = Phylo.read(StringIO("((S0001:0.1,S0002:0.1):0.2,S0003:0.3);"), "newick")
    # Label the (A,B) internal as the family Hantaviridae.
    for node in tree.get_nonterminals():
        if {t.name for t in node.get_terminals()} == {"S0001", "S0002"}:
            node._lca_name = "Hantaviridae"
            node._lca_rank = "family"
    out = tmp_path / "tree.xml"
    write_phyloxml(
        None, out, reps, id_map,
        cfg={}, prefix="t", alphabet="nucleotide",
        msa_tool="MAFFT", msa_version="v7",
        tree_tool="FastTree", tree_version="2.1",
        model="GTR", ufboot=None, tree=tree,
    )
    root = ET.parse(out).getroot()
    # Find the labelled internal clade.
    labelled = [
        c for c in root.iter(_ns("clade"))
        if (c.find(_ns("name")) is not None
            and c.find(_ns("name")).text == "Hantaviridae")
    ]
    assert len(labelled) == 1
    tax = labelled[0].find(_ns("taxonomy"))
    assert tax is not None
    assert tax.find(_ns("scientific_name")).text == "Hantaviridae"
    assert tax.find(_ns("rank")).text == "family"


def test_internal_rank_falls_back_to_other_for_unknown_rank(tmp_path):
    """A rank that isn't in PhyloXML's enumeration (e.g. NCBI's
    "no rank") must be written as "other" — otherwise the file
    fails to validate."""
    from io import StringIO
    from Bio import Phylo
    reps = [_make_seq("A"), _make_seq("B"), _make_seq("C")]
    id_map = {f"S{i + 1:04d}": rep.id for i, rep in enumerate(reps)}
    tree = Phylo.read(StringIO("((S0001:0.1,S0002:0.1):0.2,S0003:0.3);"), "newick")
    for node in tree.get_nonterminals():
        if {t.name for t in node.get_terminals()} == {"S0001", "S0002"}:
            node._lca_name = "Some clade"
            node._lca_rank = "no rank"
    out = tmp_path / "tree.xml"
    write_phyloxml(
        None, out, reps, id_map,
        cfg={}, prefix="t", alphabet="nucleotide",
        msa_tool="MAFFT", msa_version="v7",
        tree_tool="FastTree", tree_version="2.1",
        model="GTR", ufboot=None, tree=tree,
    )
    root = ET.parse(out).getroot()
    ranks = [r.text for r in root.iter(_ns("rank"))]
    assert "other" in ranks


def test_internal_without_lca_emits_no_name_or_taxonomy(tmp_path):
    """An internal node the LCA annotator skipped (coverage too low,
    etc.) should not have a <name> or <taxonomy> block — just
    branch_length + confidence, like before Pass B."""
    reps = [_make_seq("A"), _make_seq("B"), _make_seq("C")]
    out = _run_write(tmp_path, reps)  # no LCA pre-set
    root = ET.parse(out).getroot()
    # Internal clades have child <clade> elements; leaf clades don't.
    internals = [
        c for c in root.iter(_ns("clade"))
        if c.findall(_ns("clade"))  # has at least one nested clade
    ]
    for internal in internals:
        assert internal.find(_ns("taxonomy")) is None
        # The leaf-formatted <name> elements (e.g. "Hantaan_..|A|...")
        # are children of leaf clades, not internal clades. An
        # internal's direct <name> should be absent.
        direct_name = internal.find(_ns("name"))
        assert direct_name is None or direct_name.text not in (
            "Hantaviridae", "Orthohantavirus",
        )


def test_rooting_method_in_description(tmp_path):
    """The <phylogeny><description> should record which rooting method
    actually fired, so the user can tell auto-chain outcomes apart."""
    reps = [_make_seq("A"), _make_seq("B"), _make_seq("C")]
    newick = tmp_path / "tree.nwk"
    _write_newick(newick, "((S0001:0.1,S0002:0.1):0.2,S0003:0.3);")
    id_map = {f"S{i + 1:04d}": rep.id for i, rep in enumerate(reps)}
    out = tmp_path / "tree.xml"
    write_phyloxml(
        newick, out, reps, id_map,
        cfg={}, prefix="t", alphabet="nucleotide",
        msa_tool="MAFFT", msa_version="v7",
        tree_tool="FastTree", tree_version="2.1",
        model="GTR", ufboot=None,
        rooting_method="taxonomy",
    )
    root = ET.parse(out).getroot()
    desc = root.find(f"{_ns('phylogeny')}/{_ns('description')}").text
    assert "rooting=taxonomy" in desc


def test_tree_is_ladderized(tmp_path):
    """Subtrees are reordered largest-first; the asymmetric Newick
    ((A,B),C) should land C as the first child of the root after
    ladderize(reverse=True)."""
    reps = [_make_seq("A"), _make_seq("B"), _make_seq("C")]
    newick = tmp_path / "tree.nwk"
    _write_newick(newick, "((S0001:0.1,S0002:0.1):0.2,S0003:0.3);")
    id_map = {f"S{i + 1:04d}": rep.id for i, rep in enumerate(reps)}
    out = tmp_path / "tree.xml"
    write_phyloxml(
        newick, out, reps, id_map,
        cfg={}, prefix="test", alphabet="nucleotide",
        msa_tool="MAFFT", msa_version="v7",
        tree_tool="FastTree", tree_version="2.1",
        model="GTR", ufboot=None,
    )
    root = ET.parse(out).getroot()
    # The root phylogeny's first child clade should be the larger subtree
    # (the one containing A and B) after ladderize(reverse=True).
    phylogeny = root.find(_ns("phylogeny"))
    root_clade = phylogeny.find(_ns("clade"))
    first_subtree = root_clade.findall(_ns("clade"))[0]
    # The larger subtree has children; a leaf does not (leaf has no nested <clade>).
    assert first_subtree.findall(_ns("clade")), "ladderize(reverse=True) put smaller subtree first"
