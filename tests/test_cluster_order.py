"""Clusters come back in a fixed order, whatever order the backend used.

A multi-threaded backend does not return clusters in a stable order:
clustering one 215-sequence file six times with ``mmseqs easy-cluster
--threads 10`` gave the same clusters in four different orders (2026-11-25),
so two identical repseq runs wrote files that differed line by line.
"""
from __future__ import annotations

import itertools
from unittest.mock import patch

from repseq.clustering import canonical_order, run_clustering
from repseq.models import Cluster
from repseq.representative.selector import apply_representative_selection

_CFG = {
    "clustering": {"backend": "mmseqs2"},
    "representative": {"priority": ["refseq", "reviewed_uniprot", "longest"]},
}


def _clusters(make_seq, layout):
    """``layout`` = [(rep_id, [member_ids])] -> fresh Cluster objects."""
    return [
        Cluster(
            cluster_id=f"c{i}",
            representative=make_seq(rep, "ACGT" * 10),
            members=[make_seq(m, "ACGT" * 10) for m in members],
        )
        for i, (rep, members) in enumerate(layout)
    ]


def _shape(clusters):
    return [
        (c.representative.id, [m.id for m in c.members]) for c in clusters
    ]


def test_canonical_order_sorts_clusters_and_members(make_seq):
    clusters = _clusters(make_seq, [("C", ["z", "a"]), ("A", []), ("B", ["m", "b"])])
    assert _shape(canonical_order(clusters)) == [
        ("A", []), ("B", ["b", "m"]), ("C", ["a", "z"]),
    ]


def test_every_backend_order_gives_the_same_result(make_seq):
    """The property that matters: any order in, one order out."""
    layout = [("C", ["z", "a"]), ("A", ["q"]), ("B", ["m", "b"])]
    seen = set()
    for perm in itertools.permutations(layout):
        shuffled = [(rep, list(reversed(members))) for rep, members in perm]
        with patch(
            "repseq.clustering.mmseqs2.run_clustering",
            return_value=_clusters(make_seq, shuffled),
        ):
            seen.add(str(_shape(run_clustering([], 0.9, _CFG))))
    assert len(seen) == 1


def test_canonical_order_keeps_the_same_sequence_objects(make_seq):
    """force_select matches by object identity — sorting must not copy."""
    clusters = _clusters(make_seq, [("B", ["y", "x"]), ("A", [])])
    before = {id(s) for c in clusters for s in [c.representative, *c.members]}
    ordered = canonical_order(clusters)
    after = {id(s) for c in ordered for s in [c.representative, *c.members]}
    assert before == after
    assert ordered is clusters


def test_selection_orders_clusters_by_their_final_representative(make_seq):
    """Re-election changes who leads a cluster; the output order follows."""
    clusters = _clusters(make_seq, [("B", []), ("Z", ["A_refseq"])])
    clusters[1].members[0].is_refseq = True       # wins cluster 2's slot
    out = apply_representative_selection(clusters, _CFG)
    assert [c.representative.id for c in out] == ["A_refseq", "B"]


def test_representative_wins_a_tie_against_an_equal_member(make_seq):
    """Unchanged behaviour: an equally good member does not displace the
    backend's representative, in any member order."""
    for members in (["a", "b"], ["b", "a"]):
        clusters = canonical_order(_clusters(make_seq, [("m", members)]))
        out = apply_representative_selection(clusters, _CFG)
        assert out[0].representative.id == "m"


def test_cluster_ids_follow_the_fixed_order(make_seq):
    """Backends number clusters in the order they returned them, so the ids
    in _clusters.tsv (global mode) moved between identical runs."""
    clusters = _clusters(make_seq, [("C", []), ("A", []), ("B", [])])
    out = canonical_order(clusters)
    assert [(c.cluster_id, c.representative.id) for c in out] == [
        ("cluster_000001", "A"), ("cluster_000002", "B"), ("cluster_000003", "C"),
    ]


def test_sets_are_written_as_sorted_lists():
    """Resolved override id sets reached the lockfile and the run log in
    Python's per-process set order."""
    from repseq.config import sanitize_config, sets_to_sorted_lists

    cfg = {"_overrides_runtime": {"ids": frozenset({"b.1", "a.1", "c.1"}),
                                  "stages": {"hmm", "annotation"}}}
    want = {"_overrides_runtime": {"ids": ["a.1", "b.1", "c.1"],
                                   "stages": ["annotation", "hmm"]}}
    assert sanitize_config(cfg, drop_private=False) == want
    assert sets_to_sorted_lists(cfg) == want
    # The live config is untouched (the pipeline still needs real sets).
    assert isinstance(cfg["_overrides_runtime"]["ids"], frozenset)
