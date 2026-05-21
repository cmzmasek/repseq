"""2F per-protein trees: family collection, CDS selection, orchestration.

MAFFT / FastTree subprocesses are mocked (patched on
``repseq.phylo.pipeline`` where ``_build_tree`` calls them) so the
token→family grouping, segment scoping, min-taxa floor, and soft-fail
behaviour can be locked without real binaries.
"""
from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from repseq.models import Sequence, SequenceType
from repseq.phylo.pipeline import PhyloError
from repseq.phylo.per_protein import (
    _best_satisfying_cds,
    _segment_proteins,
    collect_family_specs,
    run_per_protein_phylogeny,
)


# ---------------------------------------------------------------------------
# Fixtures / builders
# ---------------------------------------------------------------------------

def _hit(target, e=1e-20, ali_from=1, ali_to=100, passing=True):
    return {
        "target": target, "dom_evalue": e,
        "ali_from": ali_from, "ali_to": ali_to, "passing": passing,
    }


def _prot(pid, product, seq, hits):
    return {
        "protein_id": pid, "product": product,
        "length": len(seq), "sequence": seq, "hmm_hits": hits,
    }


def _segseq(segment, proteins, accession):
    return Sequence(
        id=accession, header=accession, sequence="ACGT" * 30,
        seq_type=SequenceType.NUCLEOTIDE, accession=accession,
        segment=segment, proteins=proteins,
    )


def _concat_rep(iso, seg_to_proteins):
    """A CONCAT representative carrying per-segment proteins."""
    segs = [_segseq(s, ps, f"{iso}_{s}") for s, ps in seg_to_proteins.items()]
    return Sequence(
        id=f"CONCAT|{iso}", header=iso, sequence="ACGT" * 90,
        seq_type=SequenceType.NUCLEOTIDE, isolate_id=iso, concat_segments=segs,
    )


# A nucleocapsid hit (single-HMM token) and an ordered G1+G2 pair
# (multidomain token "Bunya_G1--Bunya_G2": G1 must sit C-terminal to G2,
# i.e. G1.ali_from > G2.ali_to).
def _S_proteins():
    return [_prot("S_N", "nucleoprotein", "M" * 200, [_hit("Bunya_nucleocap")])]


def _M_proteins():
    return [_prot(
        "M_GP", "glycoprotein polyprotein", "M" * 400,
        [_hit("Bunya_G2", ali_from=10, ali_to=100),
         _hit("Bunya_G1", ali_from=150, ali_to=250)],
    )]


def _seg_cfg(segment_markers):
    return {
        "phylo": {"tool": "fasttree", "per_protein": {"min_taxa": 3}},
        "segmented": {
            "enabled": True, "virus": "Bunya",
            "viruses": {"Bunya": {
                "segments": ["L", "M", "S"],
                "segment_markers": segment_markers,
            }},
        },
        "_hmm_runtime": {"active": True},
    }


# ---------------------------------------------------------------------------
# collect_family_specs
# ---------------------------------------------------------------------------

def test_collect_family_specs_segmented_scopes_and_prefixes():
    cfg = _seg_cfg({
        "S": {"hmms": ["Bunya_nucleocap"]},
        "M": {"hmms": ["Bunya_G1--Bunya_G2"]},
    })
    specs = collect_family_specs(cfg)
    assert specs == [
        ("M_Bunya_G1--Bunya_G2", "Bunya_G1--Bunya_G2", "M"),
        ("S_Bunya_nucleocap", "Bunya_nucleocap", "S"),
    ]


def test_collect_family_specs_non_segmented_from_cluster_protein():
    cfg = {"segmented": {"enabled": False}, "clustering": {"cluster_protein": [
        {"name": "RdRp", "hmms": ["RdRp_4"]},
        {"hmms": ["Nucleocap"]},
        ["polymerase"],  # alias-only legacy entry: no hmms → ignored
    ]}}
    specs = collect_family_specs(cfg)
    assert specs == [("RdRp_4", "RdRp_4", None), ("Nucleocap", "Nucleocap", None)]


def test_collect_family_specs_empty_when_no_tokens():
    cfg = _seg_cfg({})  # no hmms anywhere
    assert collect_family_specs(cfg) == []


# ---------------------------------------------------------------------------
# CDS selection + segment scoping
# ---------------------------------------------------------------------------

def test_best_satisfying_cds_picks_longest():
    short = _prot("p1", "nuc", "M" * 100, [_hit("Bunya_nucleocap")])
    longer = _prot("p2", "nuc", "M" * 300, [_hit("Bunya_nucleocap")])
    miss = _prot("p3", "other", "M" * 999, [_hit("Other_HMM")])
    chosen = _best_satisfying_cds([short, miss, longer], ["Bunya_nucleocap"])
    assert chosen["protein_id"] == "p2"


def test_best_satisfying_cds_none_when_unsatisfied():
    p = _prot("p1", "x", "M" * 100, [_hit("Bunya_nucleocap", passing=False)])
    assert _best_satisfying_cds([p], ["Bunya_nucleocap"]) is None


def test_segment_proteins_scopes_to_named_segment():
    rep = _concat_rep("iso1", {"S": _S_proteins(), "M": _M_proteins()})
    s_prots = _segment_proteins(rep, "S")
    assert [p["protein_id"] for p in s_prots] == ["S_N"]
    assert _segment_proteins(rep, "L") == []  # no L segment present


# ---------------------------------------------------------------------------
# Orchestration (mocked binaries)
# ---------------------------------------------------------------------------

def _stub_mafft(input_fasta: Path, output_fasta: Path, cfg):
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    output_fasta.write_text(input_fasta.read_text())


def _stub_fasttree(msa_fasta: Path, output_newick: Path, cfg, is_protein):
    # Build a ladder over whatever short ids the MSA carries — works for
    # any family size.
    short_ids = [
        line[1:].split()[0]
        for line in msa_fasta.read_text().splitlines()
        if line.startswith(">")
    ]
    body = short_ids[0]
    for sid in short_ids[1:]:
        body = f"({body}:0.1,{sid}:0.1)"
    output_newick.parent.mkdir(parents=True, exist_ok=True)
    output_newick.write_text(body + ";\n")


def _run(reps, cfg, tmp_path, prefix="bunya"):
    with patch("repseq.phylo.pipeline.run_mafft", side_effect=_stub_mafft), \
         patch("repseq.phylo.pipeline.run_fasttree", side_effect=_stub_fasttree):
        return run_per_protein_phylogeny(reps, cfg, tmp_path, prefix)


def test_builds_one_tree_per_token(tmp_path):
    cfg = _seg_cfg({
        "S": {"hmms": ["Bunya_nucleocap"]},
        "M": {"hmms": ["Bunya_G1--Bunya_G2"]},
    })
    reps = [
        _concat_rep(f"iso{i}", {"S": _S_proteins(), "M": _M_proteins()})
        for i in range(3)
    ]
    files = _run(reps, cfg, tmp_path)
    names = sorted(f.name for f in files)
    sub = tmp_path / "bunya_per_protein"
    # Both families clear min_taxa=3.
    for fam in ("M_Bunya_G1--Bunya_G2", "S_Bunya_nucleocap"):
        for ext in ("_msa.fasta", "_tree.nwk", "_tree.xml", "_tree_id_map.tsv"):
            assert (sub / f"{fam}{ext}").exists()
    assert all(str(f).startswith(str(sub)) for f in files)
    assert any("M_Bunya_G1--Bunya_G2_tree.xml" in n for n in names)


def test_family_below_min_taxa_is_skipped(tmp_path):
    cfg = _seg_cfg({
        "S": {"hmms": ["Bunya_nucleocap"]},
        "M": {"hmms": ["Bunya_G1--Bunya_G2"]},
    })
    # 3 isolates carry S; only 2 carry the full M architecture (the third
    # M CDS is missing the G1 domain).
    reps = [
        _concat_rep("iso0", {"S": _S_proteins(), "M": _M_proteins()}),
        _concat_rep("iso1", {"S": _S_proteins(), "M": _M_proteins()}),
        _concat_rep("iso2", {
            "S": _S_proteins(),
            "M": [_prot("M_partial", "gp", "M" * 400, [_hit("Bunya_G2")])],
        }),
    ]
    files = _run(reps, cfg, tmp_path)
    sub = tmp_path / "bunya_per_protein"
    assert (sub / "S_Bunya_nucleocap_tree.xml").exists()
    assert not (sub / "M_Bunya_G1--Bunya_G2_tree.xml").exists()


def test_incongruence_table_written_for_two_families(tmp_path):
    cfg = _seg_cfg({
        "S": {"hmms": ["Bunya_nucleocap"]},
        "M": {"hmms": ["Bunya_G1--Bunya_G2"]},
    })
    reps = [
        _concat_rep(f"iso{i}", {"S": _S_proteins(), "M": _M_proteins()})
        for i in range(3)
    ]
    files = _run(reps, cfg, tmp_path)
    inc = tmp_path / "bunya_per_protein" / "bunya_incongruence.tsv"
    assert inc.exists()
    assert inc in files
    lines = inc.read_text().splitlines()
    assert lines[0] == "tree_a\ttree_b\trf\tnorm_rf\tn_common_taxa"
    # One pair (S × M); all 3 isolates carry both → 3 common taxa.
    body = [ln for ln in lines[1:] if ln.strip()]
    assert len(body) == 1
    assert body[0].endswith("\t3")  # n_common_taxa


def test_incongruence_table_skipped_when_disabled(tmp_path):
    cfg = _seg_cfg({
        "S": {"hmms": ["Bunya_nucleocap"]},
        "M": {"hmms": ["Bunya_G1--Bunya_G2"]},
    })
    cfg["phylo"]["per_protein"]["incongruence"] = False
    reps = [
        _concat_rep(f"iso{i}", {"S": _S_proteins(), "M": _M_proteins()})
        for i in range(3)
    ]
    _run(reps, cfg, tmp_path)
    assert not (tmp_path / "bunya_per_protein" / "bunya_incongruence.tsv").exists()


def test_incongruence_skipped_for_single_family(tmp_path):
    cfg = _seg_cfg({"S": {"hmms": ["Bunya_nucleocap"]}})
    reps = [_concat_rep(f"iso{i}", {"S": _S_proteins()}) for i in range(3)]
    _run(reps, cfg, tmp_path)
    # One family, no --phylo genome tree → no pair to compare.
    assert not (tmp_path / "bunya_per_protein" / "bunya_incongruence.tsv").exists()


def test_raises_when_no_tokens_configured(tmp_path):
    cfg = _seg_cfg({})
    reps = [_concat_rep(f"iso{i}", {"S": _S_proteins()}) for i in range(3)]
    with pytest.raises(PhyloError, match="no HMM marker tokens"):
        _run(reps, cfg, tmp_path)


def test_raises_when_hmm_tier_did_not_run(tmp_path):
    cfg = _seg_cfg({"S": {"hmms": ["Bunya_nucleocap"]}})
    cfg["_hmm_runtime"] = {"active": False}
    # Strip hmm_hits so the fallback scan also reports the tier inactive.
    def _bare_S():
        return [_prot("S_N", "nucleoprotein", "M" * 200, [])]
    reps = [_concat_rep(f"iso{i}", {"S": _bare_S()}) for i in range(3)]
    with pytest.raises(PhyloError, match="HMM tier did not run"):
        _run(reps, cfg, tmp_path)


def test_raises_when_no_family_clears_min_taxa(tmp_path):
    cfg = _seg_cfg({"S": {"hmms": ["Bunya_nucleocap"]}})
    reps = [_concat_rep("iso0", {"S": _S_proteins()})]  # only 1 rep
    with pytest.raises(PhyloError, match="nothing built"):
        _run(reps, cfg, tmp_path)
