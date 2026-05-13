"""Representative selection: RefSeq > reviewed UniProt > longest tiebreaker."""
from __future__ import annotations

from repseq.models import Cluster
from repseq.representative.selector import (
    apply_representative_selection,
    select_representative,
)


PRIORITY = ["refseq", "reviewed_uniprot", "longest"]


def test_select_refseq_wins_over_longer_non_refseq(make_seq):
    refseq = make_seq("rs", "A" * 100, is_refseq=True)
    longer = make_seq("ln", "A" * 200)
    chosen = select_representative([refseq, longer], PRIORITY)
    assert chosen.id == "rs"


def test_select_reviewed_uniprot_beats_unreviewed_when_no_refseq(make_seq):
    reviewed = make_seq("rev", "A" * 100, is_reviewed=True)
    unreviewed = make_seq("un", "A" * 200)
    chosen = select_representative([reviewed, unreviewed], PRIORITY)
    assert chosen.id == "rev"


def test_select_longest_when_tied(make_seq):
    a = make_seq("short", "A" * 50)
    b = make_seq("long", "A" * 500)
    assert select_representative([a, b], PRIORITY).id == "long"


def test_apply_representative_selection_swaps_in_cluster(make_seq):
    """If the cluster's current rep is dominated, it should be swapped."""
    original_rep = make_seq("rep", "A" * 300)  # longest but not RefSeq
    member_refseq = make_seq("m1", "A" * 50, is_refseq=True)
    member_short = make_seq("m2", "A" * 10)
    cluster = Cluster(
        cluster_id="c1",
        representative=original_rep,
        members=[member_refseq, member_short],
    )
    apply_representative_selection([cluster], {"representative": {"priority": PRIORITY}})

    assert cluster.representative.id == "m1"
    assert original_rep in cluster.members
    assert member_refseq not in cluster.members


def test_apply_representative_keeps_rep_when_already_best(make_seq):
    rep = make_seq("rep", "A" * 100, is_refseq=True)
    member = make_seq("m1", "A" * 200)
    cluster = Cluster(cluster_id="c1", representative=rep, members=[member])
    apply_representative_selection([cluster], {"representative": {"priority": PRIORITY}})
    assert cluster.representative.id == "rep"
    assert cluster.members == [member]
