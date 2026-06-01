"""Phylogeny-based taxonomy review: clade-walk impute / conflict logic."""
from __future__ import annotations

from io import StringIO

from Bio import Phylo

from repseq.models import Sequence, SequenceType, TaxonomyInfo
from repseq.phylo.taxonomy_review import apply_imputations, run_taxonomy_review


def _tree(newick: str):
    return Phylo.read(StringIO(newick), "newick")


def _seq(sid, *, family="Coronaviridae", genus=None, subgenus=None,
         refseq=False, organism=None):
    lineage = {}
    if subgenus:
        lineage["subgenus"] = subgenus
    tax = TaxonomyInfo(family=family, genus=genus, lineage=lineage)
    return Sequence(
        id=sid, header=sid, sequence="ACGT", seq_type=SequenceType.NUCLEOTIDE,
        accession=sid, organism=organism or sid, taxonomy=tax,
        is_refseq=refseq,
    )


def _cfg(**over):
    base = {
        "enabled": True,
        "ranks": ["family", "genus", "subgenus"],
        "min_support": 90, "min_purity": 0.9, "min_agreeing": 3,
        "require_refseq_anchor": True,
    }
    base.update(over)
    return {"phylo": {"taxonomy_review": base}}


# A clade of four Betacoronaviruses (one RefSeq anchor) + a blank-genus
# leaf, a sister clade of three Alphacoronaviruses + a mislabelled leaf,
# and a lone outgroup. Support 95/92 on the two clades, 50 at the root.
_NEWICK = (
    "((X1:0.1,X2:0.1,X3:0.1,X4:0.1,Xb:0.1)95:0.2,"
    "(Y1:0.1,Y2:0.1,Y3:0.1,Ym:0.1)92:0.2,O:0.3)50:0.0;"
)


def _reps(**extra):
    reps = {
        "X1": _seq("X1", genus="Betacoronavirus", subgenus="Sarbecovirus", refseq=True),
        "X2": _seq("X2", genus="Betacoronavirus", subgenus="Sarbecovirus"),
        "X3": _seq("X3", genus="Betacoronavirus", subgenus="Sarbecovirus"),
        "X4": _seq("X4", genus="Betacoronavirus", subgenus="Sarbecovirus"),
        "Xb": _seq("Xb", genus=None, subgenus=None),          # blank → impute
        "Y1": _seq("Y1", genus="Alphacoronavirus", refseq=True),
        "Y2": _seq("Y2", genus="Alphacoronavirus"),
        "Y3": _seq("Y3", genus="Alphacoronavirus"),
        "Ym": _seq("Ym", genus="Betacoronavirus"),            # wrong → conflict
        "O": _seq("O", genus="Gammacoronavirus"),
    }
    reps.update(extra)
    return reps


def test_impute_missing_high_confidence(tmp_path):
    out = run_taxonomy_review(
        _tree(_NEWICK), _reps(), tree_tool="iqtree",
        cfg=_cfg(), out_dir=tmp_path, file_prefix="t",
    )
    imp = out["imputations"]
    # Xb's blank genus filled from the pure, anchored, well-supported clade.
    assert imp["Xb"]["genus"] == "Betacoronavirus"
    # And its subgenus filled too (hierarchy-consistent with the genus).
    assert imp["Xb"]["subgenus"] == "Sarbecovirus"
    # Written to the TSV.
    text = (tmp_path / "t_taxonomy_review.tsv").read_text()
    assert "impute_missing" in text
    assert "Betacoronavirus" in text


def test_conflict_flagged_not_imputed(tmp_path):
    out = run_taxonomy_review(
        _tree(_NEWICK), _reps(), tree_tool="iqtree",
        cfg=_cfg(), out_dir=tmp_path, file_prefix="t",
    )
    # Ym claims Beta but sits among Alphas → flagged, never auto-applied.
    assert "Ym" not in out["imputations"]
    ym = [v for v in out["verdicts"] if v.seq_id == "Ym" and v.rank == "genus"]
    assert ym and ym[0].action == "conflict_flag"
    assert ym[0].current_value == "Betacoronavirus"
    assert ym[0].suggested_value == "Alphacoronavirus"


def test_hierarchy_constraint_blocks_cross_genus_subgenus(tmp_path):
    # Zb sits in the Beta clade but is populated genus=Alpha. Its subgenus
    # must NOT be imputed Sarbecovirus (that would contradict its genus).
    reps = _reps(Xb=_seq("Xb", genus="Alphacoronavirus", subgenus=None))
    # rename the blank leaf role: reuse Xb slot as the Alpha-in-Beta-clade leaf
    out = run_taxonomy_review(
        _tree(_NEWICK), reps, tree_tool="iqtree",
        cfg=_cfg(), out_dir=tmp_path, file_prefix="t",
    )
    assert "subgenus" not in out["imputations"].get("Xb", {})


def test_anchor_required_demotes_to_medium(tmp_path):
    # Remove the RefSeq anchor from the Beta clade → no high-confidence call,
    # so nothing is written to the corrected-copy imputation map.
    reps = _reps(X1=_seq("X1", genus="Betacoronavirus", subgenus="Sarbecovirus"))
    out = run_taxonomy_review(
        _tree(_NEWICK), reps, tree_tool="iqtree",
        cfg=_cfg(), out_dir=tmp_path, file_prefix="t",
    )
    assert "Xb" not in out["imputations"]
    # But it still appears in the review TSV as a medium-confidence suggestion.
    xb = [v for v in out["verdicts"] if v.seq_id == "Xb" and v.rank == "genus"]
    assert xb and xb[0].confidence == "medium"


def test_low_support_clade_yields_no_call(tmp_path):
    # Drop the Beta clade's support below the medium floor (70).
    nwk = _NEWICK.replace(")95:", ")40:")
    out = run_taxonomy_review(
        _tree(nwk), _reps(), tree_tool="iqtree",
        cfg=_cfg(), out_dir=tmp_path, file_prefix="t",
    )
    assert "Xb" not in out["imputations"]
    assert not [v for v in out["verdicts"] if v.seq_id == "Xb"]


def test_disabled_returns_empty(tmp_path):
    out = run_taxonomy_review(
        _tree(_NEWICK), _reps(), tree_tool="iqtree",
        cfg=_cfg(enabled=False), out_dir=tmp_path, file_prefix="t",
    )
    assert out == {}
    assert not (tmp_path / "t_taxonomy_review.tsv").exists()


def test_apply_imputations_fills_clean_values_without_mutating_originals():
    reps = [
        _seq("Xb", genus=None, subgenus=None),
        _seq("X1", genus="Betacoronavirus", subgenus="Sarbecovirus"),
    ]
    imputations = {"Xb": {"genus": "Betacoronavirus", "subgenus": "Sarbecovirus"}}
    corrected, ci = apply_imputations(reps, None, imputations)
    # Corrected copy carries the filled values...
    assert corrected[0].taxonomy.get_rank("genus") == "Betacoronavirus"
    assert corrected[0].taxonomy.get_rank("subgenus") == "Sarbecovirus"
    # ...the values are clean (no ";imputed" pollution — the review TSV is
    # the ledger)...
    assert ";" not in corrected[0].taxonomy.get_rank("genus")
    # ...and the ORIGINAL is untouched.
    assert reps[0].taxonomy.get_rank("genus") is None
    # Reps with no imputation pass through unchanged (same object).
    assert corrected[1] is reps[1]
    assert ci is None


def test_apply_imputations_propagates_to_segments():
    # Segmented: the imputation keyed on the CONCAT rep propagates to every
    # segment of that isolate in complete_isolates.
    rep = _seq("CONCAT|iso1", genus=None)
    rep.id = "CONCAT|iso1"
    rep.isolate_id = "iso1"
    seg_l = _seq("L1", genus=None)
    seg_s = _seq("S1", genus=None)
    corrected, ci = apply_imputations(
        [rep], {"iso1": [seg_l, seg_s]}, {"CONCAT|iso1": {"genus": "Orthohantavirus"}},
    )
    assert ci["iso1"][0].taxonomy.get_rank("genus") == "Orthohantavirus"
    assert ci["iso1"][1].taxonomy.get_rank("genus") == "Orthohantavirus"
    # originals untouched
    assert seg_l.taxonomy.get_rank("genus") is None


def test_corrected_rep_tsv_carries_imputed_value(tmp_path):
    """End-to-end of the corrected-copy glue: imputed rep copies fed to the
    real rep-TSV writer produce a TSV whose genus cell is filled (clean,
    no ';imputed' pollution)."""
    from repseq.output.report import write_representative_sequences_tsv
    reps = [_seq("Xb", genus=None), _seq("X1", genus="Betacoronavirus", refseq=True)]
    corrected, _ = apply_imputations(
        reps, None, {"Xb": {"genus": "Betacoronavirus"}}
    )
    path = tmp_path / "reps_corrected.tsv"
    write_representative_sequences_tsv(corrected, path)
    lines = path.read_text().splitlines()
    header = lines[0].split("\t")
    gi = header.index("genus")
    ai = header.index("accessions")
    xb_row = next(
        cells for r in lines[1:]
        if (cells := r.split("\t"))[ai] == "Xb"
    )
    assert xb_row[gi] == "Betacoronavirus"
    assert ";" not in xb_row[gi]


def test_no_verdicts_writes_no_file(tmp_path):
    # Every leaf fully and consistently labelled → no rows, no file.
    reps = {
        sid: _seq(sid, genus="Betacoronavirus", subgenus="Sarbecovirus",
                  refseq=(sid == "X1"))
        for sid in ["X1", "X2", "X3", "X4", "Xb"]
    }
    nwk = "((X1:0.1,X2:0.1,X3:0.1,X4:0.1,Xb:0.1)95:0.2,O:0.3)50:0.0;"
    reps["O"] = _seq("O", genus="Betacoronavirus", subgenus="Sarbecovirus")
    out = run_taxonomy_review(
        _tree(nwk), reps, tree_tool="iqtree",
        cfg=_cfg(), out_dir=tmp_path, file_prefix="t",
    )
    assert out["verdicts"] == []
    assert out["path"] is None
    assert not (tmp_path / "t_taxonomy_review.tsv").exists()
