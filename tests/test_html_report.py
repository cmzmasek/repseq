"""Single-file HTML run report."""
from __future__ import annotations

from repseq.output.html_report import write_html_report

_MONO_HDR = (
    "tree\trank\ttaxon\tn_leaves\tstatus\t"
    "n_clusters\tn_intruders\tintruder_clusters\tintruder_taxa\n"
)


def _mono(tmp_path, taxon="Alphacoronavirus", status="polyphyletic"):
    (tmp_path / "x_monophyly.tsv").write_text(
        _MONO_HDR
        + f"x_tree.xml\tgenus\t{taxon}\t10\t{status}\t3\t5\t2\tBetacoronavirus\n"
    )


def test_html_report_bundles_flags_and_trees(tmp_path):
    _mono(tmp_path)
    (tmp_path / "x_tree.png").write_bytes(b"\x89PNG\r\n\x1a\n fake")
    (tmp_path / "x_clusters.tsv").write_text("a\tb\n1\t2\n")
    out = write_html_report(tmp_path, "x")
    assert out is not None and out.name == "x_report.html"
    txt = out.read_text()
    assert "repseq run report" in txt
    assert "Analysis flags" in txt
    assert "Alphacoronavirus" in txt                 # flag embedded
    assert "data:image/png;base64," in txt           # tree figure embedded
    assert "Output files" in txt
    assert "x_clusters.tsv" in txt                    # file index
    assert "x_report.html" not in txt.split("Output files")[1]  # not self-listed


def test_html_report_escapes_content(tmp_path):
    _mono(tmp_path, taxon="<script>evil</script>")
    out = write_html_report(tmp_path, "x")
    txt = out.read_text()
    assert "<script>evil</script>" not in txt
    assert "&lt;script&gt;" in txt


def test_html_report_clean_run_says_so(tmp_path):
    _mono(tmp_path, status="monophyletic")
    out = write_html_report(tmp_path, "x")
    assert "No taxonomy / tree conflicts flagged" in out.read_text()


def test_html_report_none_when_empty(tmp_path):
    assert write_html_report(tmp_path, "x") is None
