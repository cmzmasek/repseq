"""Partitioned-supermatrix tree (2E): per-family alignment, column-wise
concatenation, NEXUS charsets, IQ-TREE partition dispatch, and the
soft-fallback contract.

MAFFT / IQ-TREE subprocesses are mocked on ``repseq.phylo.partition`` (where
``build_partitioned_phylogeny`` calls them) so the supermatrix assembly,
partition-file generation, and run_phylogeny dispatch can be locked without
real binaries.
"""
from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

from repseq.models import Sequence, SequenceType
from repseq.phylo.partition import (
    _nexus_safe,
    build_partitioned_phylogeny,
    build_supermatrix,
    read_msa,
    write_partition_nexus,
)
from repseq.phylo.pipeline import run_phylogeny


# ---------------------------------------------------------------------------
# Builders (mirror tests/test_per_protein.py)
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
    segs = [_segseq(s, ps, f"{iso}_{s}") for s, ps in seg_to_proteins.items()]
    rep = Sequence(
        id=f"CONCAT|{iso}", header=iso, sequence="ACGT" * 90,
        seq_type=SequenceType.NUCLEOTIDE, isolate_id=iso, concat_segments=segs,
    )
    # alphabet=protein gate needs a populated protein_sequence on every rep.
    rep.protein_sequence = "M" * 120
    return rep


def _S_prot(n):
    return [_prot("S_N", "nucleoprotein", "M" * n, [_hit("Bunya_nucleocap")])]


def _M_prot(n):
    return [_prot(
        "M_GP", "glycoprotein", "K" * n,
        [_hit("Bunya_G2", ali_from=10, ali_to=100),
         _hit("Bunya_G1", ali_from=150, ali_to=250)],
    )]


def _cfg(segment_markers, *, tool="iqtree", enabled=True, linkage="proportional"):
    return {
        "clustering": {"alphabet_for_clustering": "protein"},
        "phylo": {
            "tool": tool,
            "partition": {"enabled": enabled, "linkage": linkage, "models": {}},
        },
        "segmented": {
            "enabled": True, "virus": "Bunya",
            "viruses": {"Bunya": {
                "segments": ["M", "S"],
                "segment_markers": segment_markers,
            }},
        },
        "_hmm_runtime": {"active": True},
    }


_TWO_FAMILIES = {
    "S": {"hmms": ["Bunya_nucleocap"]},
    "M": {"hmms": ["Bunya_G2--Bunya_G1"]},
}


def _three_reps():
    return [
        _concat_rep(f"iso{i}", {"S": _S_prot(200 + i), "M": _M_prot(400 + i)})
        for i in range(3)
    ]


# ---------------------------------------------------------------------------
# Pure helpers
# ---------------------------------------------------------------------------

def test_nexus_safe_collapses_multidomain_separator_and_specials():
    assert _nexus_safe("M_Bunya_G2--Bunya_G1") == "M_Bunya_G2_Bunya_G1"
    assert _nexus_safe("PF00937.24") == "PF00937_24"
    assert _nexus_safe("") == "part"


def test_build_supermatrix_concatenates_in_order_with_column_offsets():
    fam = [
        ("FamA", "FamA", {"S0001": "AAA", "S0002": "A-A"}),
        ("FamB", "FamB", {"S0001": "CCCC", "S0002": "CC--"}),
    ]
    sup, blocks = build_supermatrix(fam, ["S0001", "S0002"])
    assert sup["S0001"] == "AAACCCC"
    assert sup["S0002"] == "A-ACC--"
    # 1-based inclusive ranges, contiguous.
    assert blocks == [("FamA", "FamA", 1, 3), ("FamB", "FamB", 4, 7)]


def test_build_supermatrix_gap_fills_missing_family():
    """A row absent from a family gets that block filled with '-' × width."""
    fam = [
        ("FamA", "FamA", {"S0001": "MMM", "S0002": "MM-"}),
        ("FamB", "FamB", {"S0001": "KKKK"}),  # S0002 lacks FamB
    ]
    sup, blocks = build_supermatrix(fam, ["S0001", "S0002"])
    assert sup["S0001"] == "MMMKKKK"
    assert sup["S0002"] == "MM-" + "----"   # 4-wide gap block
    assert blocks[1] == ("FamB", "FamB", 4, 7)


def test_write_partition_nexus_modelfinder_by_default(tmp_path):
    p = tmp_path / "part.nex"
    write_partition_nexus(
        [("L_RdRP_4", "L_RdRP_4", 1, 741), ("S_Nuc", "S_Nuc", 742, 941)],
        {}, p,
    )
    text = p.read_text()
    assert "charset L_RdRP_4 = 1-741;" in text
    assert "charset S_Nuc = 742-941;" in text
    # Unpinned → MFP per partition.
    assert "charpartition repseq = MFP:L_RdRP_4, MFP:S_Nuc;" in text


def test_write_partition_nexus_pins_model_by_family_label(tmp_path):
    p = tmp_path / "part.nex"
    write_partition_nexus(
        [("L_RdRP_4", "L_RdRP_4", 1, 741), ("S_Nuc", "S_Nuc", 742, 941)],
        {"L_RdRP_4": "LG+G4"}, p,
    )
    text = p.read_text()
    # Pinned family uses its model; the other still falls to MFP.
    assert "charpartition repseq = LG+G4:L_RdRP_4, MFP:S_Nuc;" in text


def test_read_msa_keys_on_short_id_first_token(tmp_path):
    p = tmp_path / "msa.fasta"
    p.write_text(">S0001 Foo virus|host\nMK\nLP\n>S0002 Bar\nMKLP\n")
    assert read_msa(p) == {"S0001": "MKLP", "S0002": "MKLP"}


# ---------------------------------------------------------------------------
# Mocks for the integration tests
# ---------------------------------------------------------------------------

def _stub_mafft(input_fasta: Path, output_fasta: Path, cfg, **kwargs):
    # Stand-in alignment: right-pad every record with '-' to the longest
    # body so the family MSA has equal-width rows, like real MAFFT.
    records: list[tuple[str, str]] = []
    header, buf = None, []
    for line in input_fasta.read_text().splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(buf)))
            header, buf = line, []
        else:
            buf.append(line.strip())
    if header is not None:
        records.append((header, "".join(buf)))
    width = max((len(b) for _, b in records), default=0)
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    with open(output_fasta, "w") as fh:
        for hdr, body in records:
            fh.write(hdr + "\n" + body + ("-" * (width - len(body))) + "\n")


def _stub_iqtree(msa_fasta, output_newick, cfg, is_protein,
                 summary_path=None, partition_file=None,
                 partition_linkage="proportional"):
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
    if summary_path is not None:
        summary_path.write_text("partitioned run\n")


# ---------------------------------------------------------------------------
# build_partitioned_phylogeny
# ---------------------------------------------------------------------------

def test_build_partitioned_happy_path_writes_supermatrix_and_nexus(tmp_path):
    reps = _three_reps()
    captured: dict = {}

    def _iq(*a, **kw):
        captured.update(kw)
        _stub_iqtree(*a, **kw)

    with patch("repseq.phylo.partition.run_mafft", side_effect=_stub_mafft), \
         patch("repseq.phylo.partition.run_iqtree", side_effect=_iq):
        files = build_partitioned_phylogeny(reps, _cfg(_TWO_FAMILIES), tmp_path, "test")

    assert files is not None
    names = {f.name for f in files}
    # Canonical outputs + the partition extras.
    assert {"test_msa.fasta", "test_tree.nwk", "test_tree.xml",
            "test_tree_id_map.tsv", "test_partition.nex"} <= names
    # One per-family MSA per family (label carries the segment prefix).
    assert "test_msa_M_Bunya_G2--Bunya_G1.fasta" in names
    assert "test_msa_S_Bunya_nucleocap.fasta" in names

    # IQ-TREE was invoked in partition mode with the requested linkage.
    assert captured["partition_file"].name == "test_partition.nex"
    assert captured["partition_linkage"] == "proportional"

    # NEXUS declares two charsets and per-partition ModelFinder.
    nex = (tmp_path / "test_partition.nex").read_text()
    assert nex.count("charset ") == 2
    assert "MFP:" in nex

    # Supermatrix width = sum of the two family widths; every rep is a row.
    sup = read_msa(tmp_path / "test_msa.fasta")
    assert len(sup) == 3
    assert len({len(v) for v in sup.values()}) == 1  # all rows equal length

    # Temp per-family inputs are cleaned up.
    assert not list(tmp_path.glob("test_partin_*.fasta"))


def test_build_partitioned_returns_none_without_hmm_tokens(tmp_path):
    reps = _three_reps()
    with patch("repseq.phylo.partition.run_mafft", side_effect=_stub_mafft), \
         patch("repseq.phylo.partition.run_iqtree", side_effect=_stub_iqtree):
        out = build_partitioned_phylogeny(reps, _cfg({}), tmp_path, "test")
    assert out is None


def test_build_partitioned_returns_none_with_single_family(tmp_path):
    reps = _three_reps()
    cfg = _cfg({"S": {"hmms": ["Bunya_nucleocap"]}})  # only one family
    with patch("repseq.phylo.partition.run_mafft", side_effect=_stub_mafft), \
         patch("repseq.phylo.partition.run_iqtree", side_effect=_stub_iqtree):
        out = build_partitioned_phylogeny(reps, cfg, tmp_path, "test")
    assert out is None


def test_build_partitioned_returns_none_when_hmm_tier_inactive(tmp_path):
    # Families configured, but no hmm_hits anywhere and runtime inactive.
    reps = [
        _concat_rep(f"iso{i}", {
            "S": [_prot("S_N", "nuc", "M" * 100, [])],
            "M": [_prot("M_GP", "gp", "K" * 100, [])],
        })
        for i in range(3)
    ]
    cfg = _cfg(_TWO_FAMILIES)
    cfg["_hmm_runtime"] = {"active": False}
    with patch("repseq.phylo.partition.run_mafft", side_effect=_stub_mafft), \
         patch("repseq.phylo.partition.run_iqtree", side_effect=_stub_iqtree):
        out = build_partitioned_phylogeny(reps, cfg, tmp_path, "test")
    assert out is None


# ---------------------------------------------------------------------------
# run_phylogeny dispatch
# ---------------------------------------------------------------------------

def test_run_phylogeny_dispatches_to_partition_for_protein_iqtree(tmp_path):
    reps = _three_reps()
    with patch("repseq.phylo.partition.run_mafft", side_effect=_stub_mafft), \
         patch("repseq.phylo.partition.run_iqtree", side_effect=_stub_iqtree):
        files = run_phylogeny(reps, _cfg(_TWO_FAMILIES), tmp_path, "test")

    # The partition NEXUS is the tell that the partitioned path ran.
    assert (tmp_path / "test_partition.nex").exists()
    assert any(f.name == "test_partition.nex" for f in files)


def test_run_phylogeny_falls_back_to_concat_when_partition_disabled(tmp_path):
    reps = _three_reps()
    cfg = _cfg(_TWO_FAMILIES, enabled=False)
    # Concat path goes through pipeline.run_iqtree / run_mafft, not partition.
    with patch("repseq.phylo.pipeline.run_mafft", side_effect=_stub_mafft), \
         patch("repseq.phylo.pipeline.run_iqtree", side_effect=_stub_iqtree):
        files = run_phylogeny(reps, cfg, tmp_path, "test")

    assert not (tmp_path / "test_partition.nex").exists()
    assert files  # concat-then-align still produced a tree


def test_build_partitioned_trims_per_family_and_retains_untrimmed(tmp_path):
    """phylo.trimal on → each family trimmed before concat; per-family +
    supermatrix untrimmed companions retained; charsets reflect trimmed widths."""
    reps = _three_reps()
    cfg = _cfg(_TWO_FAMILIES)
    cfg["phylo"]["trimal"] = {"enabled": True, "mode": "automated1"}

    def _trim_half(input_fasta, output_fasta, cfg, settings, *, label=""):
        # Drop the back half of every row (a deterministic stand-in for
        # column trimming), keeping equal width across the family.
        recs, hdr, buf = [], None, []
        for line in Path(input_fasta).read_text().splitlines():
            if line.startswith(">"):
                if hdr is not None:
                    recs.append((hdr, "".join(buf)))
                hdr, buf = line, []
            else:
                buf.append(line.strip())
        if hdr is not None:
            recs.append((hdr, "".join(buf)))
        keep = max(1, (max((len(b) for _, b in recs), default=0)) // 2)
        with open(output_fasta, "w") as fh:
            for h, b in recs:
                fh.write(h + "\n" + b[:keep] + "\n")
        return True

    with patch("repseq.phylo.partition.run_mafft", side_effect=_stub_mafft), \
         patch("repseq.phylo.partition.run_iqtree", side_effect=_stub_iqtree), \
         patch("repseq.phylo.partition.maybe_trim", side_effect=_trim_half):
        files = build_partitioned_phylogeny(reps, cfg, tmp_path, "test")

    names = {f.name for f in files}
    # Per-family + supermatrix untrimmed companions retained.
    assert "test_msa_untrimmed.fasta" in names
    assert "test_msa_S_Bunya_nucleocap_untrimmed.fasta" in names
    # Charsets reflect the TRIMMED widths: supermatrix width == sum of
    # trimmed family widths, and is shorter than the untrimmed supermatrix.
    trimmed = read_msa(tmp_path / "test_msa.fasta")
    untrimmed = read_msa(tmp_path / "test_msa_untrimmed.fasta")
    w_trim = len(next(iter(trimmed.values())))
    w_raw = len(next(iter(untrimmed.values())))
    assert w_trim < w_raw
    nex = (tmp_path / "test_partition.nex").read_text()
    last_end = int(nex.rsplit("-", 1)[1].split(";")[0])
    assert last_end == w_trim   # final charset end == trimmed supermatrix width
    # Provenance.
    assert "trim=" in (tmp_path / "test_tree.xml").read_text()


def test_run_phylogeny_skips_partition_for_fasttree(tmp_path):
    reps = _three_reps()
    cfg = _cfg(_TWO_FAMILIES, tool="fasttree")

    def _stub_fasttree(msa_fasta, output_newick, cfg, is_protein):
        short_ids = [
            line[1:].split()[0]
            for line in msa_fasta.read_text().splitlines()
            if line.startswith(">")
        ]
        body = short_ids[0]
        for sid in short_ids[1:]:
            body = f"({body}:0.1,{sid}:0.1)"
        output_newick.write_text(body + ";\n")

    with patch("repseq.phylo.pipeline.run_mafft", side_effect=_stub_mafft), \
         patch("repseq.phylo.pipeline.run_fasttree", side_effect=_stub_fasttree):
        run_phylogeny(reps, cfg, tmp_path, "test")

    # FastTree can't partition → concat path, no NEXUS.
    assert not (tmp_path / "test_partition.nex").exists()
