"""Parsing the ModelFinder pick(s) out of IQ-TREE's ``.iqtree`` report.

The .iqtree report is buried-but-stable text. These tests pin the two
shapes the parser must handle (single-model and per-partition) and the
soft-fail contracts (missing file / malformed text → ``{}`` so the
caller can render a degraded but non-crashing summary).
"""
from __future__ import annotations

from pathlib import Path

from repseq.phylo.iqtree_parse import (
    format_models_for_description,
    parse_chosen_models,
    write_iqtree_model_file,
)


# ---------------------------------------------------------------------------
# parse_chosen_models — non-partitioned
# ---------------------------------------------------------------------------

_NON_PARTITIONED_BIC = """
ModelFinder
-----------

Best-fit model according to BIC: LG+I+G4

Tree in newick format: (..);
"""

_NON_PARTITIONED_LEGACY = """
ModelFinder
-----------

Best-fit model: WAG+G4 chosen according to BIC.

Tree in newick format: (..);
"""


def test_parse_non_partitioned_modern_phrasing(tmp_path):
    path = tmp_path / "run.iqtree"
    path.write_text(_NON_PARTITIONED_BIC)
    assert parse_chosen_models(path) == {"GENOME": "LG+I+G4"}


def test_parse_non_partitioned_legacy_phrasing(tmp_path):
    path = tmp_path / "run.iqtree"
    path.write_text(_NON_PARTITIONED_LEGACY)
    assert parse_chosen_models(path) == {"GENOME": "WAG+G4"}


def test_parse_non_partitioned_missing_returns_empty(tmp_path):
    """A report that doesn't carry a best-fit line (e.g. truncated /
    unexpected IQ-TREE format) returns ``{}`` — the caller falls back
    to the configured model string."""
    path = tmp_path / "run.iqtree"
    path.write_text("Something else entirely.")
    assert parse_chosen_models(path) == {}


def test_parse_missing_file_returns_empty(tmp_path):
    """Soft-fail when the report file is absent."""
    assert parse_chosen_models(tmp_path / "nope.iqtree") == {}


# ---------------------------------------------------------------------------
# parse_chosen_models — partitioned
# ---------------------------------------------------------------------------

_PARTITIONED = """
ModelFinder
-----------

List of best-fit models per partition:

  ID  Model           Speed  Parameters
   1  LG+I+G4         1.000  LG+I{0.123}+G4{0.456}
   2  JTT+G4          0.500  JTT+G4{0.789}
   3  WAG+I+G4        0.300  WAG+I{0.111}+G4{0.222}

Best-fit partition model selected by ModelFinder.

PARTITION       MODEL
1               LG+I+G4
"""


def test_parse_partitioned_zips_to_labels(tmp_path):
    path = tmp_path / "run.iqtree"
    path.write_text(_PARTITIONED)
    out = parse_chosen_models(
        path, partition_labels=["CoV_S1", "CoV_M", "CoV_N"],
    )
    assert out == {
        "CoV_S1": "LG+I+G4",
        "CoV_M": "JTT+G4",
        "CoV_N": "WAG+I+G4",
    }


def test_parse_partitioned_fewer_labels_than_rows(tmp_path):
    """Extra rows beyond the caller-supplied labels are silently dropped —
    a partial annotation is better than none."""
    path = tmp_path / "run.iqtree"
    path.write_text(_PARTITIONED)
    out = parse_chosen_models(path, partition_labels=["CoV_S1", "CoV_M"])
    assert out == {"CoV_S1": "LG+I+G4", "CoV_M": "JTT+G4"}


def test_parse_partitioned_more_labels_than_rows(tmp_path):
    """Extra labels with no matching row are dropped — same soft-fail rule."""
    path = tmp_path / "run.iqtree"
    path.write_text(_PARTITIONED)
    out = parse_chosen_models(
        path, partition_labels=["A", "B", "C", "D", "E"],
    )
    # Only 3 rows in the synthetic table.
    assert out == {"A": "LG+I+G4", "B": "JTT+G4", "C": "WAG+I+G4"}


def test_parse_partitioned_no_table_returns_empty(tmp_path):
    """A 'partitioned' call against a non-partitioned report falls back
    to empty (the table header is absent)."""
    path = tmp_path / "run.iqtree"
    path.write_text(_NON_PARTITIONED_BIC)
    assert parse_chosen_models(path, partition_labels=["X", "Y"]) == {}


# ---------------------------------------------------------------------------
# format_models_for_description
# ---------------------------------------------------------------------------

def test_format_empty_is_none():
    assert format_models_for_description({}) is None


def test_format_single_returns_bare_model():
    assert format_models_for_description({"GENOME": "LG+I+G4"}) == "LG+I+G4"


def test_format_multi_pipes_label_model_pairs():
    out = format_models_for_description(
        {"CoV_S1": "LG+I+G4", "CoV_M": "JTT+G4"},
    )
    assert out == "CoV_S1:LG+I+G4|CoV_M:JTT+G4"


# ---------------------------------------------------------------------------
# write_iqtree_model_file
# ---------------------------------------------------------------------------

def test_write_model_file_writes_one_line_per_entry(tmp_path):
    path = tmp_path / "out" / "x_iqtree_model.txt"
    write_iqtree_model_file(
        {"CoV_S1": "LG+I+G4", "CoV_M": "JTT+G4"}, path,
    )
    assert path.read_text() == "CoV_S1: LG+I+G4\nCoV_M: JTT+G4\n"


def test_write_model_file_skips_empty_dict(tmp_path):
    path = tmp_path / "x_iqtree_model.txt"
    write_iqtree_model_file({}, path)
    assert not path.exists()
