"""Unit tests for pairwise unrooted Robinson-Foulds incongruence
(``repseq.phylo.incongruence``).

RF is computed dependency-free from Bio.Phylo clades, so these build
tiny Newick files (+ id maps) and check the bipartition arithmetic
directly.
"""
from __future__ import annotations

from pathlib import Path

from repseq.phylo.incongruence import (
    compute_incongruence,
    load_tree_clusters,
    rf_distance,
    write_incongruence_tsv,
)


def _write(tmp_path: Path, name: str, newick: str, mapping: dict[str, str]) -> tuple[Path, Path]:
    nwk = tmp_path / f"{name}_tree.nwk"
    idm = tmp_path / f"{name}_tree_id_map.tsv"
    nwk.write_text(newick + "\n")
    lines = ["short_id\taccession"] + [f"{k}\t{v}" for k, v in mapping.items()]
    idm.write_text("\n".join(lines) + "\n")
    return nwk, idm


_IDMAP4 = {"S1": "A", "S2": "B", "S3": "C", "S4": "D"}


# ---------------------------------------------------------------------------
# load_tree_clusters
# ---------------------------------------------------------------------------

def test_load_relabels_short_ids_and_collects_clusters(tmp_path):
    nwk, idm = _write(tmp_path, "t", "((S1,S2),(S3,S4));", _IDMAP4)
    taxa, clusters = load_tree_clusters(nwk, idm)
    assert taxa == frozenset({"A", "B", "C", "D"})
    # Two non-trivial internal clusters: {A,B} and {C,D} (root dropped).
    assert frozenset({"A", "B"}) in clusters
    assert frozenset({"C", "D"}) in clusters
    assert frozenset({"A", "B", "C", "D"}) not in clusters


def test_load_returns_none_for_unparseable(tmp_path):
    nwk = tmp_path / "bad_tree.nwk"
    nwk.write_text("((((;")
    idm = tmp_path / "bad_tree_id_map.tsv"
    idm.write_text("short_id\taccession\n")
    assert load_tree_clusters(nwk, idm) is None


# ---------------------------------------------------------------------------
# rf_distance
# ---------------------------------------------------------------------------

def test_identical_topologies_rf_zero():
    a = [frozenset({"A", "B"}), frozenset({"C", "D"})]
    rf, norm, n = rf_distance(
        frozenset("ABCD"), a, frozenset("ABCD"), list(a)
    )
    assert rf == 0
    assert norm == 0.0
    assert n == 4


def test_maximally_different_4taxon_trees_rf_two():
    # AB|CD  vs  AC|BD : the two distinct unrooted 4-taxon topologies.
    ab = [frozenset({"A", "B"}), frozenset({"C", "D"})]
    ac = [frozenset({"A", "C"}), frozenset({"B", "D"})]
    rf, norm, n = rf_distance(frozenset("ABCD"), ab, frozenset("ABCD"), ac)
    assert rf == 2
    assert norm == 1.0  # 2 / (2*(4-3))


def test_mirror_bipartition_is_same():
    # A cluster and its complement describe the same unrooted split.
    a = [frozenset({"A", "B"})]            # AB | CD
    b = [frozenset({"C", "D"})]            # CD | AB  (same split)
    rf, _, _ = rf_distance(frozenset("ABCD"), a, frozenset("ABCD"), b)
    assert rf == 0


def test_fewer_than_four_common_taxa_norm_is_na():
    a = [frozenset({"A", "B"})]
    rf, norm, n = rf_distance(
        frozenset("ABC"), a, frozenset("ABC"), list(a)
    )
    assert n == 3
    assert norm is None  # undefined below 4 shared taxa


def test_partial_overlap_scored_on_common_taxa():
    # Tree1 on {A,B,C,D,E}, Tree2 on {A,B,C,D,F}; common = {A,B,C,D}.
    t1 = [frozenset({"A", "B"}), frozenset({"A", "B", "E"}), frozenset({"C", "D"})]
    t2 = [frozenset({"A", "B"}), frozenset({"C", "D", "F"}), frozenset({"C", "D"})]
    rf, norm, n = rf_distance(
        frozenset("ABCDE"), t1, frozenset("ABCDF"), t2
    )
    assert n == 4
    # Restricted to {A,B,C,D} both reduce to AB|CD → identical.
    assert rf == 0


# ---------------------------------------------------------------------------
# compute_incongruence + writer
# ---------------------------------------------------------------------------

def test_compute_incongruence_pairs_and_tsv(tmp_path):
    n1 = _write(tmp_path, "Spike", "((S1,S2),(S3,S4));", _IDMAP4)
    n2 = _write(tmp_path, "N", "((S1,S3),(S2,S4));", _IDMAP4)
    n3 = _write(tmp_path, "GENOME", "((S1,S2),(S3,S4));", _IDMAP4)
    rows = compute_incongruence([
        ("Spike", *n1), ("N", *n2), ("GENOME", *n3),
    ])
    # 3 trees → 3 pairs, in input order.
    assert [(r["tree_a"], r["tree_b"]) for r in rows] == [
        ("Spike", "N"), ("Spike", "GENOME"), ("N", "GENOME"),
    ]
    by_pair = {(r["tree_a"], r["tree_b"]): r for r in rows}
    assert by_pair[("Spike", "N")]["rf"] == 2
    assert by_pair[("Spike", "GENOME")]["rf"] == 0      # identical topo
    assert by_pair[("N", "GENOME")]["rf"] == 2

    out = tmp_path / "inc.tsv"
    write_incongruence_tsv(rows, out)
    lines = out.read_text().splitlines()
    assert lines[0] == "tree_a\ttree_b\trf\tnorm_rf\tn_common_taxa"
    assert "Spike\tN\t2\t1.0000\t4" in lines
    assert "Spike\tGENOME\t0\t0.0000\t4" in lines


def test_unparseable_tree_dropped_from_pairs(tmp_path):
    good1 = _write(tmp_path, "A", "((S1,S2),(S3,S4));", _IDMAP4)
    good2 = _write(tmp_path, "B", "((S1,S3),(S2,S4));", _IDMAP4)
    bad_nwk = tmp_path / "C_tree.nwk"
    bad_nwk.write_text("(((;")
    bad_map = tmp_path / "C_tree_id_map.tsv"
    bad_map.write_text("short_id\taccession\n")
    rows = compute_incongruence([
        ("A", *good1), ("B", *good2), ("C", bad_nwk, bad_map),
    ])
    # Only the A×B pair survives.
    assert [(r["tree_a"], r["tree_b"]) for r in rows] == [("A", "B")]
