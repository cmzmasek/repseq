"""Plain-English analysis flags synthesised from the conflict tables."""
from __future__ import annotations

from repseq.output.flags import collect_flags, write_flags_report

_MONO = [
    "tree", "rank", "taxon", "n_leaves", "status",
    "n_clusters", "n_intruders", "intruder_clusters", "intruder_taxa",
]
_INC = ["tree_a", "tree_b", "rf", "norm_rf", "n_common_taxa"]
_REV = [
    "accession", "organism", "rank", "current_value", "suggested_value",
    "action", "clade_support", "clade_purity", "n_agreeing",
    "anchor_refseq", "confidence",
]


def _write(path, header, rows):
    lines = ["\t".join(header)] + ["\t".join(map(str, r)) for r in rows]
    path.write_text("\n".join(lines) + "\n")


def test_polyphyletic_genome_is_flagged_mono_is_not(tmp_path):
    _write(tmp_path / "x_monophyly.tsv", _MONO, [
        ["x_tree.xml", "genus", "Alphacoronavirus", "199", "polyphyletic",
         "70", "285", "67", "Betacoronavirus"],
        ["x_tree.xml", "genus", "Deltacoronavirus", "24", "monophyletic",
         "1", "0", "0", ""],
    ])
    msgs = [f.message for f in collect_flags(tmp_path, "x")]
    assert any("Alphacoronavirus" in m and "polyphyletic" in m for m in msgs)
    assert not any("Deltacoronavirus" in m for m in msgs)


def test_reassortment_marker_disagreement(tmp_path):
    _write(tmp_path / "x_monophyly.tsv", _MONO, [
        ["x_tree.xml", "genus", "Betacoronavirus", "100", "monophyletic",
         "1", "0", "0", ""],
        ["x_per_protein/Spike_tree.xml", "genus", "Betacoronavirus", "100",
         "polyphyletic", "5", "10", "3", "Alphacoronavirus"],
    ])
    flags = collect_flags(tmp_path, "x")
    assert any(
        f.category == "reassortment"
        and "Betacoronavirus" in f.message and "Spike" in f.message
        for f in flags
    )


def test_incongruence_threshold(tmp_path):
    (tmp_path / "x_per_protein").mkdir()
    _write(tmp_path / "x_per_protein" / "x_incongruence.tsv", _INC, [
        ["Spike", "RdRP", "20", "0.45", "30"],
        ["Spike", "Nucleocapsid", "2", "0.05", "30"],  # below threshold
        ["RdRP", "Helicase", "1", "NA", "2"],          # NA ignored
    ])
    msgs = [f.message for f in collect_flags(tmp_path, "x")]
    assert any("Spike" in m and "RdRP" in m and "0.45" in m for m in msgs)
    assert not any("Nucleocapsid" in m for m in msgs)
    assert not any("Helicase" in m for m in msgs)


def test_taxonomy_conflict_flags_only_genuine_disagreement(tmp_path):
    _write(tmp_path / "x_taxonomy_review.tsv", _REV, [
        ["AB1", "Some virus", "genus", "Alphacoronavirus", "Betacoronavirus",
         "flag", "95", "0.95", "5", "TRUE", "high"],
        ["AB2", "Other", "genus", "", "Betacoronavirus",  # blank-fill impute
         "impute", "90", "0.9", "4", "TRUE", "high"],
    ])
    msgs = [f.message for f in collect_flags(tmp_path, "x")]
    assert any("AB1" in m for m in msgs)
    assert not any("AB2" in m for m in msgs)


def test_write_flags_report_clean_run(tmp_path):
    _write(tmp_path / "x_monophyly.tsv", _MONO, [
        ["x_tree.xml", "genus", "Deltacoronavirus", "24", "monophyletic",
         "1", "0", "0", ""],
    ])
    out = write_flags_report(tmp_path, "x")
    assert out is not None and out.name == "x_flags.txt"
    assert "No flags raised" in out.read_text()


def test_write_flags_report_none_without_sources(tmp_path):
    assert write_flags_report(tmp_path, "x") is None
