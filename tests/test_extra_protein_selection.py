"""extra_protein CDS selection (v0.75.0): "HMM specs choose first" + length_aa.

The motivating case is real (Elliovirales, 2026-10-16): a tospovirus isolate
whose GenBank M segment labels its 1,129 aa glycoprotein precursor "NSm" and
its real 303 aa NSm "Gn/Gc". The alias spec ``["NSm", "movement protein"]``
kept the LONGEST match, so the mislabelled glycoprotein landed in the NSm
file. The fixtures below reproduce that shape, plus the two regression probes
that shaped the design:

* influenza PA-X / PA — only the CDS an HMM spec PICKS is claimed, so a
  shorter CDS sharing the marker's domain stays available to an alias spec;
* Paramyxoviridae N / P / M / F — alias-only specs target proteins that DO
  carry passing HMM hits (to profiles no spec declares); those must still be
  selected, i.e. "has an HMM hit" is NOT the claiming criterion.

Protein dicts mirror the production schema (``protein_id``, ``product``,
``length``, ``sequence``, ``hmm_hits`` with ``target`` / ``dom_evalue`` /
``ali_from`` / ``ali_to`` / ``passing``) — same builders as
``tests/test_per_protein.py``.
"""
from __future__ import annotations

import csv
from pathlib import Path
from unittest.mock import patch

from repseq.config import load_config, validate_config
from repseq.models import RunResult, Sequence, SequenceType, TaxonomyInfo
from repseq.output.report import (
    write_extra_protein_fastas,
    write_protein_taxonomic_report,
    write_protein_taxonomic_report_tsv,
)
from repseq.phylo import per_protein
from repseq.phylo.per_protein import (
    build_extra_pick_context,
    collect_extra_length_bounds,
    collect_extra_specs,
    pick_extra_cds,
    pick_marker_cds,
)


# ---------------------------------------------------------------------------
# Builders (production schema)
# ---------------------------------------------------------------------------

def _hit(target, e=1e-20, ali_from=1, ali_to=100, passing=True):
    return {
        "target": target, "dom_evalue": e,
        "ali_from": ali_from, "ali_to": ali_to, "passing": passing,
    }


def _prot(pid, product, n_aa, hits=()):
    seq = "M" * n_aa
    return {
        "protein_id": pid, "product": product,
        "length": n_aa, "sequence": seq, "hmm_hits": list(hits),
    }


def _segseq(segment, proteins, accession):
    return Sequence(
        id=accession, header=accession, sequence="ACGT" * 30,
        seq_type=SequenceType.NUCLEOTIDE, accession=accession,
        segment=segment, proteins=proteins,
    )


def _concat_rep(iso, seg_to_proteins, genus="Orthotospovirus"):
    segs = [_segseq(s, ps, f"{iso}_{s}") for s, ps in seg_to_proteins.items()]
    return Sequence(
        id=f"CONCAT|{iso}", header=iso, sequence="ACGT" * 90,
        seq_type=SequenceType.NUCLEOTIDE, isolate_id=iso,
        concat_segments=segs,
        taxonomy=TaxonomyInfo(genus=genus, family="Tospoviridae"),
    )


def _l_proteins():
    return [_prot("L1", "RNA-dependent RNA polymerase", 2900, [_hit("Bunya_RdRp")])]


def _gpc(pid, product):
    """Tospovirus glycoprotein precursor — carries a passing Bunya_G1 hit."""
    return _prot(pid, product, 1129, [_hit("Bunya_G1", ali_from=40, ali_to=520)])


def _nsm(pid, product, n_aa=303):
    """Tospovirus NSm movement protein — no profile, so no HMM hit."""
    return _prot(pid, product, n_aa)


def _seg_cfg(m_extra, *, hmm_active=True, l_marker=True):
    cfg = {
        "phylo": {"tool": "fasttree", "per_protein": {"min_taxa": 3}},
        "segmented": {
            "enabled": True, "virus": "Tospo",
            "viruses": {"Tospo": {
                "expected_segments": 2, "segments": ["L", "M"],
                "isolate_regex": r"(?P<isolate>.+)",
                "segment_markers": ({"L": {"hmms": ["Bunya_RdRp"]}} if l_marker else {}),
                "extra_protein": {"M": m_extra},
            }},
        },
        "hmm": {"multidomain_overlap_tolerance": 30},
    }
    if hmm_active:
        cfg["_hmm_runtime"] = {"active": True}
    return cfg


GLYCO = {"name": "M_glycoprotein", "hmms": ["Bunya_G1"]}
NSM = {"name": "M_NSm", "aliases": ["NSm", "movement protein"]}


def _pick(cfg, proteins, label, segment="M", hmm_active=True):
    specs = {s[0]: s for s in collect_extra_specs(cfg)}
    _lab, tokens, aliases, seg = specs[label]
    ctx = build_extra_pick_context(cfg, hmm_active)
    return pick_extra_cds(proteins, label, tokens, aliases, seg, ctx)


# ---------------------------------------------------------------------------
# Claiming
# ---------------------------------------------------------------------------

def test_mislabelled_glycoprotein_is_not_picked_as_nsm():
    """The m24 shape: GPC labelled "NSm", real NSm labelled "Gn/Gc"."""
    cfg = _seg_cfg([GLYCO, NSM])
    proteins = [_gpc("GPC", "NSm"), _nsm("REAL", "Gn/Gc")]

    # Pre-v0.75 behaviour, kept on the plain picker: the longest alias match.
    assert pick_marker_cds(proteins, [], NSM["aliases"], True)["protein_id"] == "GPC"

    cds, set_aside = _pick(cfg, proteins, "M_NSm")
    assert cds is None  # better no record than the wrong protein
    assert [(p["protein_id"], why) for p, why in set_aside] == [
        ("GPC", "already selected by HMM spec M_glycoprotein"),
    ]
    # The glycoprotein spec still gets it.
    gcds, _ = _pick(cfg, proteins, "M_glycoprotein")
    assert gcds["protein_id"] == "GPC"


def test_real_nsm_still_picked_beside_a_mislabelled_glycoprotein():
    cfg = _seg_cfg([GLYCO, NSM])
    proteins = [_gpc("GPC", "movement protein"), _nsm("REAL", "NSm protein")]
    cds, set_aside = _pick(cfg, proteins, "M_NSm")
    assert cds["protein_id"] == "REAL"
    assert [p["protein_id"] for p, _ in set_aside] == ["GPC"]


def test_cluster_marker_claims_too():
    """A segment_markers HMM spec claims its pick, same as an extra HMM spec."""
    cfg = _seg_cfg([{"name": "L_alias", "aliases": ["polymerase"]}])
    cfg["segmented"]["viruses"]["Tospo"]["extra_protein"] = {
        "L": [{"name": "L_alias", "aliases": ["polymerase"]}],
    }
    cds, set_aside = _pick(cfg, _l_proteins(), "L_alias", segment="L")
    assert cds is None
    assert set_aside[0][1] == "already selected by HMM spec L_Bunya_RdRp"


def test_only_the_picked_cds_is_claimed_pa_x_shape():
    """PA-X shares PA's N-terminal domain. Even if PA-X carried a passing
    Flu_PA hit, the marker PICKS the longer PA, so PA-X stays available."""
    cfg = {
        "segmented": {
            "enabled": True, "virus": "flu",
            "viruses": {"flu": {
                "expected_segments": 1, "segments": ["3"],
                "isolate_regex": r"(?P<isolate>.+)",
                "segment_markers": {"3": {"hmms": ["Flu_PA"]}},
                "extra_protein": {"3": [{"name": "PA-X", "aliases": ["PA-X"]}]},
            }},
        },
        "_hmm_runtime": {"active": True},
    }
    pa = _prot("PA", "polymerase PA", 716, [_hit("Flu_PA", ali_to=700)])
    pax = _prot("PAX", "PA-X protein", 252, [_hit("Flu_PA", ali_to=190)])
    cds, set_aside = _pick(cfg, [pa, pax], "PA-X", segment="3")
    assert cds["protein_id"] == "PAX"
    assert set_aside == []


def test_hit_to_an_undeclared_profile_does_not_claim_paramyxo_shape():
    """Non-segmented; the alias target carries a passing HMM hit to a
    profile no spec declares — it must still be selected."""
    cfg = {
        "clustering": {
            "cluster_protein": [{"name": "L", "hmms": ["Mononeg_RNA_pol"]}],
            "extra_protein": [{"name": "N", "aliases": ["nucleoprotein"]}],
        },
        "_hmm_runtime": {"active": True},
    }
    proteins = [
        _prot("L", "large protein", 2200, [_hit("Mononeg_RNA_pol")]),
        _prot("N", "nucleoprotein", 532, [_hit("Paramyxo_ncap")]),
    ]
    cds, set_aside = _pick(cfg, proteins, "N", segment=None)
    assert cds["protein_id"] == "N"
    assert set_aside == []


def test_no_claiming_without_the_hmm_tier():
    """HMM tier off (--no-resolve / hmm disabled): nothing can be claimed,
    so alias picking is exactly the pre-v0.75 behaviour."""
    cfg = _seg_cfg([GLYCO, NSM], hmm_active=False)
    proteins = [_gpc("GPC", "NSm"), _nsm("REAL", "Gn/Gc")]
    cds, set_aside = _pick(cfg, proteins, "M_NSm", hmm_active=False)
    assert cds["protein_id"] == "GPC"
    assert set_aside == []


def test_claims_are_scoped_to_the_segment():
    """A claim on segment L never removes a CDS from an M-segment spec."""
    cfg = _seg_cfg([{"name": "M_any", "aliases": ["polymerase"]}])
    m_proteins = [_prot("MPOL", "polymerase-like M protein", 400)]
    cds, set_aside = _pick(cfg, m_proteins, "M_any")
    assert cds["protein_id"] == "MPOL"
    assert set_aside == []


# ---------------------------------------------------------------------------
# length_aa
# ---------------------------------------------------------------------------

def test_length_aa_filters_before_the_pick():
    """The longer out-of-range match must not hide the in-range one."""
    nsm = dict(NSM, aliases=["NSm", "nonstructural protein"], length_aa=[250, 400])
    cfg = _seg_cfg([GLYCO, nsm])
    proteins = [
        _prot("FUNGAL", "nonstructural protein", 549),
        _nsm("REAL", "nonstructural protein", 308),
    ]
    cds, set_aside = _pick(cfg, proteins, "M_NSm")
    assert cds["protein_id"] == "REAL"
    assert [(p["protein_id"], why) for p, why in set_aside] == [
        ("FUNGAL", "outside length_aa [250, 400]"),
    ]


def test_length_aa_alone_can_leave_a_spec_empty():
    nsm = dict(NSM, aliases=["nonstructural protein"], length_aa=[250, 400])
    cfg = _seg_cfg([nsm])
    cds, set_aside = _pick(cfg, [_prot("FUNGAL", "nonstructural protein", 549)], "M_NSm")
    assert cds is None
    assert set_aside[0][1] == "outside length_aa [250, 400]"


def test_length_aa_open_end():
    nsm = dict(NSM, length_aa=[None, 400])
    cfg = _seg_cfg([nsm])
    proteins = [_nsm("SHORT", "NSm", 120), _nsm("LONG", "NSm", 900)]
    cds, set_aside = _pick(cfg, proteins, "M_NSm")
    assert cds["protein_id"] == "SHORT"
    assert set_aside[0][1] == "outside length_aa [null, 400]"
    assert collect_extra_length_bounds(cfg) == {("M_NSm", "M"): (None, 400)}


def test_length_aa_on_an_hmm_extra_spec():
    glyco = dict(GLYCO, length_aa=[1000, 1200])
    cfg = _seg_cfg([glyco])
    long_gpc = _prot("READTHROUGH", "glycoprotein", 1500, [_hit("Bunya_G1")])
    cds, set_aside = _pick(cfg, [long_gpc, _gpc("GPC", "glycoprotein")], "M_glycoprotein")
    assert cds["protein_id"] == "GPC"
    assert [p["protein_id"] for p, _ in set_aside] == ["READTHROUGH"]


def test_a_bounded_hmm_claimant_claims_its_bounded_pick():
    """The claimant's claim is the CDS it actually outputs."""
    glyco = dict(GLYCO, length_aa=[1000, 1200])
    cfg = _seg_cfg([glyco, {"name": "M_any", "aliases": ["glycoprotein"]}])
    long_gpc = _prot("READTHROUGH", "glycoprotein", 1500, [_hit("Bunya_G1")])
    gpc = _gpc("GPC", "glycoprotein")
    cds, _ = _pick(cfg, [long_gpc, gpc], "M_any")
    assert cds["protein_id"] == "READTHROUGH"  # GPC was claimed; this was not


# ---------------------------------------------------------------------------
# All three consumers agree: FASTA writer, extra trees, protein report
# ---------------------------------------------------------------------------

def _three_reps():
    return [
        _concat_rep("m24", {"L": _l_proteins(), "M": [_gpc("GPC24", "NSm"), _nsm("N24", "Gn/Gc")]}),
        _concat_rep("ok1", {"L": _l_proteins(), "M": [_gpc("GPC1", "glycoprotein"), _nsm("N1", "NSm")]}),
        _concat_rep("ok2", {"L": _l_proteins(), "M": [_gpc("GPC2", "glycoprotein"), _nsm("N2", "NSm")]}),
    ]


def _fasta_ids(path: Path) -> list[str]:
    return [ln[1:].split()[0] for ln in path.read_text().splitlines() if ln.startswith(">")]


def test_fasta_writer_uses_the_rule_and_reports_set_aside(tmp_path, capsys):
    cfg = _seg_cfg([GLYCO, NSM])
    result = RunResult(mode="global", representatives=_three_reps())
    written = write_extra_protein_fastas(result, cfg, None, tmp_path, "pfx")
    by_name = {p.name: p for p in written}
    assert _fasta_ids(by_name["pfx_M_NSm.fasta"]) == ["N1", "N2"]
    assert _fasta_ids(by_name["pfx_M_glycoprotein.fasta"]) == ["GPC24", "GPC1", "GPC2"]
    err = capsys.readouterr().err
    assert "[extra-protein] M_NSm: left out 1 CDS" in err
    assert 'GPC24 (1129 aa, "NSm", isolate m24) — already selected by HMM spec M_glycoprotein' in err


def test_extra_trees_use_the_same_picks(tmp_path):
    cfg = _seg_cfg([GLYCO, NSM])
    calls: dict[str, dict] = {}

    def _fake_build_tree(leaf_reps, bodies, **kw):
        calls[kw["file_prefix"]] = dict(bodies)
        return []

    with patch.object(per_protein, "_build_tree", _fake_build_tree):
        per_protein.run_per_protein_phylogeny(_three_reps(), cfg, tmp_path, "pfx")
    # m24 has no NSm leaf (only 2 carriers < min_taxa 3, so no NSm tree at
    # all); the glycoprotein tree has all three.
    assert "M_NSm" not in calls
    assert set(calls["M_glycoprotein"]) == {"CONCAT|m24", "CONCAT|ok1", "CONCAT|ok2"}


def test_protein_report_counts_match_the_fasta(tmp_path):
    cfg = _seg_cfg([GLYCO, dict(NSM, length_aa=[250, 400])])
    reps = _three_reps()
    tsv = tmp_path / "r.tsv"
    assert write_protein_taxonomic_report_tsv(reps, reps, cfg, tsv)
    rows = list(csv.DictReader(tsv.open(), delimiter="\t"))
    nsm = {
        r["pool"]: r["value"] for r in rows
        if r["rank"] == "genus" and r["spec"] == "M_NSm" and r["metric"] == "coverage_count"
    }
    assert nsm == {"post_qc": "2", "reps": "2"}  # m24 excluded, as in the FASTA

    txt = tmp_path / "r.txt"
    assert write_protein_taxonomic_report(reps, reps, cfg, True, txt)
    text = txt.read_text()
    assert "== Accessory-protein (extra_protein) selection ==" in text
    assert 'M_NSm (segment M): aliases "NSm" | "movement protein"; length_aa [250, 400]' in text


# ---------------------------------------------------------------------------
# Config validation
# ---------------------------------------------------------------------------

def _seg_validation_cfg(extra_entry):
    cfg = load_config(None)
    cfg["segmented"] = {
        "enabled": True, "virus": "T",
        "viruses": {"T": {
            "expected_segments": 1, "segments": ["M"],
            "isolate_regex": r"(?P<isolate>X)",
            "extra_protein": {"M": [extra_entry]},
        }},
    }
    return cfg


def test_validate_accepts_length_aa_on_extra_protein():
    for bounds in ([250, 400], [None, 400], [250, None], None):
        entry = {"name": "NSm", "aliases": ["NSm"], "length_aa": bounds}
        assert validate_config(_seg_validation_cfg(entry)) == [], bounds
    cfg = load_config(None)
    cfg["clustering"]["extra_protein"] = [
        {"name": "NSm", "aliases": ["NSm"], "length_aa": [250, 400]},
    ]
    assert validate_config(cfg) == []


def test_validate_rejects_malformed_length_aa():
    cases = {
        "two-item list": ([250], "250-400", [250, 400, 500], [0, 400], [250.5, 400], [True, 400]),
        "at least one of min / max": ([None, None],),
        "min (400) is larger than max (250)": ([400, 250],),
    }
    for needle, values in cases.items():
        for bounds in values:
            entry = {"name": "NSm", "aliases": ["NSm"], "length_aa": bounds}
            errs = validate_config(_seg_validation_cfg(entry))
            assert any(needle in e for e in errs), (bounds, errs)


def test_validate_rejects_length_aa_on_a_cluster_marker():
    cfg = load_config(None)
    cfg["clustering"]["cluster_protein"] = [
        {"name": "L", "hmms": ["RdRP_1"], "length_aa": [100, 200]},
    ]
    errs = validate_config(cfg)
    assert any("'length_aa' is only supported on extra_protein entries" in e for e in errs)
    # Tailored message replaces the generic one — no duplicate line.
    assert not any("unknown key 'length_aa'" in e for e in errs)
