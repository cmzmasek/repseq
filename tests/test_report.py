"""Report writers: run-log secret redaction, per-group counts TSV,
per-mode representative TSVs."""
from __future__ import annotations

from repseq.models import (
    GroupStat,
    QCReport,
    RunResult,
    Sequence,
    SequenceType,
    TaxonomyInfo,
)
from pathlib import Path

from repseq.models import Cluster
from repseq.output.report import (
    write_all_reports,
    write_cluster_tsv,
    write_group_counts_tsv,
    write_hmm_diagnostic_tsv,
    write_nucleotide_taxonomic_report,
    write_representative_isolates_tsv,
    write_representative_sequences_tsv,
    write_run_log,
    write_taxonomic_report,
)


def test_write_run_log_redacts_ncbi_api_key(tmp_path):
    # The full config is dumped into the plaintext run log; a configured
    # NCBI API key must be redacted so it isn't leaked into the log file.
    cfg = {
        "output": {"dir": str(tmp_path), "prefix": "x"},
        "taxonomy": {"ncbi_email": "me@example.org", "ncbi_api_key": "SECRET123"},
    }
    log_path = tmp_path / "x_run.log"
    write_run_log(
        RunResult(mode="test"),
        QCReport(),
        cfg,
        input_paths=["in.fasta"],
        output_files=[],
        log_path=log_path,
    )
    text = log_path.read_text()
    assert "SECRET123" not in text
    assert "***redacted***" in text
    # Non-secret config is still recorded, and the caller's dict is untouched.
    assert "me@example.org" in text
    assert cfg["taxonomy"]["ncbi_api_key"] == "SECRET123"


# ---------------------------------------------------------------------------
# Per-group selection counts TSV
# ---------------------------------------------------------------------------

def test_write_group_counts_tsv(tmp_path):
    result = RunResult(
        mode="taxonomic1",
        group_stats=[
            GroupStat(grouping="genus", group="Alphainfluenzavirus",
                      n_before=1487, n_after=5, clustered=True, cutoff=0.8342),
            GroupStat(grouping="genus", group="Betainfluenzavirus",
                      n_before=3, n_after=3, clustered=False),
        ],
    )
    path = tmp_path / "x_group_counts.tsv"
    assert write_group_counts_tsv(result, path) is True
    lines = path.read_text().splitlines()
    assert lines[0] == "stratified_by\tstratum\tstratum_size_before\tstratum_size_after\tclustered\tcutoff"
    assert lines[1] == "genus\tAlphainfluenzavirus\t1487\t5\tTRUE\t0.8342"
    # Group kept without clustering: cutoff column left blank.
    assert lines[2] == "genus\tBetainfluenzavirus\t3\t3\tFALSE\t"


def test_write_group_counts_tsv_skips_when_no_stats(tmp_path):
    # A mode that recorded no group stats writes no file.
    path = tmp_path / "x_group_counts.tsv"
    assert write_group_counts_tsv(RunResult(mode="test"), path) is False
    assert not path.exists()


def test_write_group_counts_tsv_emits_diversity_curve_columns(tmp_path):
    """When any row has cutoff_counts populated, the TSV gains
    n_clusters_<c> columns sorted from most-to-least stringent.
    Rows without cutoff_counts (unclustered groups) leave those cells
    empty; below-floor cutoffs are emitted as NA."""
    result = RunResult(
        mode="taxonomic1",
        group_stats=[
            GroupStat(
                grouping="genus", group="A",
                n_before=100, n_after=20, clustered=True, cutoff=0.91,
                cutoff_counts={0.99: 87, 0.95: 42, 0.9: 18, 0.8: 5, 0.7: None},
            ),
            # Below-target group with no clustering: no curve cells.
            GroupStat(
                grouping="genus", group="B",
                n_before=3, n_after=3, clustered=False,
            ),
        ],
    )
    path = tmp_path / "g.tsv"
    assert write_group_counts_tsv(result, path) is True
    lines = path.read_text().splitlines()
    # Header gains 5 trailing curve columns, descending cutoff order.
    assert lines[0] == (
        "stratified_by\tstratum\tstratum_size_before\tstratum_size_after\t"
        "clustered\tcutoff\t"
        "n_clusters_0.99\tn_clusters_0.95\tn_clusters_0.9\t"
        "n_clusters_0.8\tn_clusters_0.7"
    )
    # Row A: cutoff_counts populated, 0.7 below cd-hit-est floor → NA.
    assert lines[1] == (
        "genus\tA\t100\t20\tTRUE\t0.9100\t87\t42\t18\t5\tNA"
    )
    # Row B: unclustered, curve cells empty (not NA — different meaning).
    assert lines[2] == "genus\tB\t3\t3\tFALSE\t\t\t\t\t\t"


def test_write_group_counts_tsv_no_curve_columns_when_feature_off(tmp_path):
    """When no row has cutoff_counts, the schema reverts to the original
    6-column layout. Backward-compat for runs that set
    diversity_curve_cutoffs: []."""
    result = RunResult(
        mode="taxonomic1",
        group_stats=[
            GroupStat(grouping="genus", group="A",
                      n_before=10, n_after=5, clustered=True, cutoff=0.9),
        ],
    )
    path = tmp_path / "g.tsv"
    write_group_counts_tsv(result, path)
    lines = path.read_text().splitlines()
    assert "n_clusters" not in lines[0]
    assert lines[1] == "genus\tA\t10\t5\tTRUE\t0.9000"


def test_write_representative_isolates_tsv_emits_isolate_columns(tmp_path):
    """A CONCAT|<isolate_id> rep with concat_segments populated must
    produce: n_segments, comma-joined segments, comma-joined accessions
    (in concat order), and total_length_nt = sum of segment lengths.
    Per-sequence noise columns (accession, segment, description,
    molecule_type) must be absent."""
    tax = TaxonomyInfo(species="Schmallenberg virus", genus="Orthobunyavirus",
                       family="Peribunyaviridae")
    s_l = Sequence(id="L_acc", header=">L_acc", sequence="A" * 6800,
                   accession="L_acc.1", segment="L",
                   seq_type=SequenceType.NUCLEOTIDE)
    s_m = Sequence(id="M_acc", header=">M_acc", sequence="A" * 4400,
                   accession="M_acc.1", segment="M",
                   seq_type=SequenceType.NUCLEOTIDE)
    s_s = Sequence(id="S_acc", header=">S_acc", sequence="A" * 800,
                   accession="S_acc.1", segment="S",
                   seq_type=SequenceType.NUCLEOTIDE)
    concat = Sequence(
        id="CONCAT|ISO_X", header="CONCAT|ISO_X", sequence="A" * 12000,
        accession=None, organism="Schmallenberg virus",
        strain="X-1", host="bovine", collection_date="2024",
        country="DE", isolate_id="ISO_X",
        is_refseq=False, is_reviewed=False, taxonomy=tax,
        seq_type=SequenceType.NUCLEOTIDE,
        concat_segments=[s_l, s_m, s_s],
    )
    path = tmp_path / "x_representative_isolates.tsv"
    write_representative_isolates_tsv([concat], path)

    lines = path.read_text().splitlines()
    header = lines[0].split("\t")
    # Per-sequence columns must NOT appear.
    for col in ("accession", "segment", "description", "molecule_type",
                "length_nt"):
        assert col not in header, f"unexpected column {col!r} in isolate TSV"
    # Per-isolate columns must appear in the documented order.
    expected_prefix = [
        "isolate_id", "isolate_id_source", "organism", "strain", "host",
        "collection_date", "country", "n_segments", "segments", "accessions",
        "total_length_nt",
    ]
    assert header[: len(expected_prefix)] == expected_prefix

    row = dict(zip(header, lines[1].split("\t")))
    assert row["isolate_id"] == "ISO_X"
    assert row["n_segments"] == "3"
    assert row["segments"] == "L,M,S"
    assert row["accessions"] == "L_acc.1,M_acc.1,S_acc.1"
    # Sum of the three segment NT lengths.
    assert row["total_length_nt"] == str(6800 + 4400 + 800)
    assert row["is_refseq"] == "FALSE"


def test_write_representative_sequences_tsv_matches_isolates_schema(tmp_path):
    """Non-segmented rep table is column-identical to the segmented
    isolates table: isolate-style columns present, per-sequence-only
    columns (accession, segment, description, molecule_type, length_nt)
    absent. Per-sequence values map onto the isolate schema: accessions =
    the single accession, total_length_nt = the sequence NT length, and
    the isolate-only cells are blank."""
    rep = Sequence(
        id="ACC.1", header=">ACC.1", sequence="ACGT" * 100,
        accession="ACC.1", organism="Some virus",
        description="hypothetical", seq_type=SequenceType.NUCLEOTIDE,
    )
    path = tmp_path / "x_representative_sequences.tsv"
    write_representative_sequences_tsv([rep], path)

    lines = path.read_text().splitlines()
    header = lines[0].split("\t")
    # Schema-identical to write_representative_isolates_tsv.
    for col in ("isolate_id", "isolate_id_source", "n_segments", "segments",
                "accessions", "total_length_nt"):
        assert col in header
    # Per-sequence-only columns must NOT appear (no slot in the shared schema).
    for col in ("accession", "segment", "description", "molecule_type",
                "length_nt"):
        assert col not in header

    row = dict(zip(header, lines[1].split("\t")))
    assert row["isolate_id"] == ""
    assert row["isolate_id_source"] == ""
    assert row["n_segments"] == ""
    assert row["segments"] == ""
    assert row["accessions"] == "ACC.1"
    assert row["total_length_nt"] == str(len("ACGT" * 100))
    assert row["organism"] == "Some virus"


def test_write_sequence_proteins_tsv_non_segmented_schema(tmp_path):
    """The non-segmented per-CDS protein TSV must (a) share the exact
    column schema of _isolate_proteins.tsv, (b) blank the columns with no
    non-segmented meaning, (c) populate segment_length_nt with the parent
    sequence NT length and accession with the parent accession, and (d)
    set `representative` from the supplied id set."""
    from repseq.output.report import write_sequence_proteins_tsv

    rep = Sequence(
        id="R1", header=">R1", sequence="ACGT" * 50,
        accession="R1.1", organism="Some virus",
        seq_type=SequenceType.NUCLEOTIDE,
    )
    rep.proteins = [
        {"protein_id": "YP_1", "product": "polymerase", "length": 660},
    ]
    nonrep = Sequence(
        id="N1", header=">N1", sequence="ACGT" * 40,
        accession="N1.1", seq_type=SequenceType.NUCLEOTIDE,
    )
    nonrep.proteins = [
        {"protein_id": "YP_2", "product": "polymerase", "length": 500},
    ]
    path = tmp_path / "x_sequence_proteins.tsv"
    written = write_sequence_proteins_tsv(
        [rep, nonrep], path, representative_ids={"R1"}
    )
    assert written is True

    lines = path.read_text().splitlines()
    header = lines[0].split("\t")
    # Identical to _isolate_proteins.tsv column schema.
    assert header[:10] == [
        "protein_id", "product", "length_aa", "isolate_id",
        "isolate_id_source", "segment", "segment_length_nt", "accession",
        "representative", "hmmscan",
    ]
    rows = {ln.split("\t")[0]: dict(zip(header, ln.split("\t")))
            for ln in lines[1:]}
    r = rows["YP_1"]
    assert r["isolate_id"] == "" and r["isolate_id_source"] == ""
    assert r["segment"] == ""
    assert r["segment_length_nt"] == str(len("ACGT" * 50))
    assert r["accession"] == "R1.1"
    assert r["representative"] == "TRUE"
    assert rows["YP_2"]["representative"] == "FALSE"


def test_write_sequence_proteins_tsv_skips_when_no_proteins(tmp_path):
    """No file when nothing carries proteins (returns False)."""
    from repseq.output.report import write_sequence_proteins_tsv

    seq = Sequence(
        id="R1", header=">R1", sequence="ACGT", accession="R1.1",
        seq_type=SequenceType.NUCLEOTIDE,
    )
    path = tmp_path / "x_sequence_proteins.tsv"
    assert write_sequence_proteins_tsv([seq], path) is False
    assert not path.exists()


def test_write_all_reports_emits_non_segmented_protein_tsvs(tmp_path):
    """Non-segmented run emits both _sequence_proteins.tsv (all post-QC,
    representative TRUE/FALSE) and _representative_sequence_proteins.tsv
    (reps only, all TRUE), sharing one header."""
    rep = Sequence(
        id="R1", header=">R1", sequence="ACGT" * 50, accession="R1.1",
        seq_type=SequenceType.NUCLEOTIDE,
    )
    rep.proteins = [{"protein_id": "YP_1", "product": "pol", "length": 660}]
    dropped = Sequence(
        id="N1", header=">N1", sequence="ACGT" * 40, accession="N1.1",
        seq_type=SequenceType.NUCLEOTIDE,
    )
    dropped.proteins = [{"protein_id": "YP_2", "product": "pol", "length": 500}]
    result = RunResult(mode="global", representatives=[rep], clusters=[])

    cfg = {"output": {"dir": str(tmp_path), "prefix": "x"}}
    write_all_reports(
        result, QCReport(), cfg, [], [],
        pre_clustering_sequences=[rep, dropped],
    )

    full = tmp_path / "x_sequence_proteins.tsv"
    subset = tmp_path / "x_representative_sequence_proteins.tsv"
    assert full.exists() and subset.exists()

    full_lines = full.read_text().splitlines()
    subset_lines = subset.read_text().splitlines()
    assert full_lines[0] == subset_lines[0]  # same schema
    assert len(full_lines) == 3   # both proteins
    assert len(subset_lines) == 2  # rep only

    header = full_lines[0].split("\t")
    rep_col = header.index("representative")
    full_reps = {ln.split("\t")[0]: ln.split("\t")[rep_col]
                 for ln in full_lines[1:]}
    assert full_reps == {"YP_1": "TRUE", "YP_2": "FALSE"}
    assert subset_lines[1].split("\t")[rep_col] == "TRUE"


def test_write_all_reports_emits_representative_isolate_proteins_subset(
    tmp_path: Path,
):
    """The rep-only TSV must (a) exist alongside _isolate_proteins.tsv,
    (b) share the exact same header, and (c) contain only rows of
    representative isolates (`representative` column all TRUE)."""
    s_picked = Sequence(
        id="P_acc", header=">P_acc", sequence="ACGT",
        accession="P_acc.1", segment="L",
        seq_type=SequenceType.NUCLEOTIDE,
    )
    s_picked.proteins = [
        {"protein_id": "YP_PICK_1", "product": "polymerase", "length": 2200},
    ]
    s_dropped = Sequence(
        id="D_acc", header=">D_acc", sequence="ACGT",
        accession="D_acc.1", segment="L",
        seq_type=SequenceType.NUCLEOTIDE,
    )
    s_dropped.proteins = [
        {"protein_id": "YP_DROP_1", "product": "polymerase", "length": 2200},
    ]
    complete_isolates = {"ISO_PICKED": [s_picked], "ISO_DROPPED": [s_dropped]}

    concat = Sequence(
        id="CONCAT|ISO_PICKED", header="CONCAT|ISO_PICKED", sequence="ACGT",
        accession=None, isolate_id="ISO_PICKED",
        seq_type=SequenceType.NUCLEOTIDE,
    )
    result = RunResult(mode="segmented", representatives=[concat], clusters=[])

    cfg = {"output": {"dir": str(tmp_path), "prefix": "x"}}
    write_all_reports(
        result, QCReport(), cfg, [], [], complete_isolates=complete_isolates,
    )

    full = tmp_path / "x_isolate_proteins.tsv"
    subset = tmp_path / "x_representative_isolate_proteins.tsv"
    assert full.exists() and subset.exists()

    full_lines = full.read_text().splitlines()
    subset_lines = subset.read_text().splitlines()
    # Same header.
    assert full_lines[0] == subset_lines[0]
    # Full file has both isolates' rows; subset has only ISO_PICKED.
    assert len(full_lines) == 3
    assert len(subset_lines) == 2

    header = subset_lines[0].split("\t")
    rep_col = header.index("representative")
    iso_col = header.index("isolate_id")
    row = subset_lines[1].split("\t")
    assert row[iso_col] == "ISO_PICKED"
    assert row[rep_col] == "TRUE"


# ---------------------------------------------------------------------------
# write_taxonomic_report
# ---------------------------------------------------------------------------

def _diversity_rows(text):
    """Parse Section-1 rows (rank + two integer columns) from the report.

    Skips the Section-2 breakdown headers like ``species (N distinct):``
    by accepting only lines whose two trailing tokens are both integers.
    """
    rows = {}
    for ln in text.splitlines():
        parts = ln.split()
        if len(parts) == 3 and parts[1].isdigit() and parts[2].isdigit():
            rows[parts[0]] = parts[1:]
    return rows


def _tax_seq(seq_id, species=None, genus=None, family=None, subgenus=None):
    lineage = {}
    if subgenus:
        lineage["subgenus"] = subgenus
    tax = TaxonomyInfo(species=species, genus=genus, family=family, lineage=lineage)
    return Sequence(
        id=seq_id, header=f">{seq_id}", sequence="ACGT",
        seq_type=SequenceType.NUCLEOTIDE, taxonomy=tax,
    )


def test_taxonomic_report_prepends_provenance_header_when_given(tmp_path):
    before = [_tax_seq("a", species="Sp A", genus="Gen X", family="Fam")]
    after = list(before)
    path = tmp_path / "x_taxonomic_report.txt"
    write_taxonomic_report(
        before, after, segmented=False, path=path,
        provenance="# repseq 9.9.9 | global selection | non-segmented",
    )
    lines = path.read_text().splitlines()
    # The provenance comment is the very first line.
    assert lines[0] == "# repseq 9.9.9 | global selection | non-segmented"
    assert lines[1] == "Taxonomic report"


def test_taxonomic_report_omits_header_when_no_provenance(tmp_path):
    before = [_tax_seq("a", species="Sp A", genus="Gen X", family="Fam")]
    path = tmp_path / "x_taxonomic_report.txt"
    write_taxonomic_report(before, list(before), segmented=False, path=path)
    assert path.read_text().splitlines()[0] == "Taxonomic report"


def test_taxonomic_report_rank_diversity_before_after(tmp_path):
    """Section 1 counts distinct non-empty taxa per rank before vs after."""
    before = [
        _tax_seq("a", species="Sp A", genus="Gen X", family="Fam"),
        _tax_seq("b", species="Sp B", genus="Gen X", family="Fam"),
        _tax_seq("c", species="Sp C", genus="Gen Y", family="Fam"),
    ]
    after = [
        _tax_seq("a", species="Sp A", genus="Gen X", family="Fam"),
        _tax_seq("c", species="Sp C", genus="Gen Y", family="Fam"),
    ]
    path = tmp_path / "x_taxonomic_report.txt"
    write_taxonomic_report(before, after, segmented=False, path=path)
    text = path.read_text()
    assert "Counting unit: sequences" in text
    # 3 species before, 2 after; 2 genera before, 2 after; 1 family both.
    diversity = _diversity_rows(text)
    assert diversity["species"] == ["3", "2"]
    assert diversity["genus"] == ["2", "2"]
    assert diversity["family"] == ["1", "1"]


def test_taxonomic_report_breakdown_lists_low_diversity_ranks(tmp_path):
    """Ranks with <=15 distinct taxa get a per-taxon before/after table."""
    before = [
        _tax_seq("a", species="SpA", genus="GenX"),
        _tax_seq("b", species="SpB", genus="GenX"),
        _tax_seq("c", species="SpC", genus="GenY"),
    ]
    after = [_tax_seq("a", species="SpA", genus="GenX")]
    path = tmp_path / "x_taxonomic_report.txt"
    write_taxonomic_report(before, after, segmented=False, path=path)
    text = path.read_text()
    assert "genus (2 distinct):" in text
    rows = {ln.split()[0]: ln.split()[1:] for ln in text.splitlines()
            if ln.strip().startswith(("GenX", "GenY"))}
    # GenX: 2 before, 1 after; GenY: 1 before, 0 after.
    assert rows["GenX"] == ["2", "1"]
    assert rows["GenY"] == ["1", "0"]


def test_taxonomic_report_truncates_breakdown_above_threshold(tmp_path):
    """A rank with more than max_breakdown distinct taxa shows the top
    max_breakdown by member count, with a note in the rank label."""
    # 20 species with descending member counts: Sp00 has 20 seqs, ... Sp19 has 1.
    before = []
    for i in range(20):
        for j in range(20 - i):
            before.append(_tax_seq(f"s{i}_{j}", species=f"Sp{i:02d}", genus="GenX"))
    after = []
    path = tmp_path / "x_taxonomic_report.txt"
    write_taxonomic_report(before, after, segmented=False, path=path, max_breakdown=15)
    text = path.read_text()
    assert "species (20 distinct, top 15 by member count shown):" in text
    # Top 15 by member count are Sp00..Sp14; Sp15..Sp19 (smallest) are dropped.
    assert "Sp00" in text and "Sp14" in text
    assert "Sp15" not in text and "Sp19" not in text
    # Genus has 1 distinct so it gets a full (un-noted) breakdown.
    assert "genus (1 distinct):" in text


def test_taxonomic_report_segmented_unit_is_isolates(tmp_path):
    """In segmented mode the unit label is 'isolates'."""
    before = [_tax_seq("iso1", species="Sp A"), _tax_seq("iso2", species="Sp B")]
    after = [_tax_seq("iso1", species="Sp A")]
    path = tmp_path / "x_taxonomic_report.txt"
    write_taxonomic_report(before, after, segmented=True, path=path)
    text = path.read_text()
    assert "Counting unit: isolates" in text


def test_taxonomic_report_excludes_blank_taxa(tmp_path):
    """Sequences with no taxonomy at a rank are not counted as a taxon."""
    before = [
        _tax_seq("a", species="Sp A"),
        _tax_seq("b", species=None),   # no species
        _tax_seq("c", species=""),     # empty species
    ]
    after = [_tax_seq("a", species="Sp A")]
    path = tmp_path / "x_taxonomic_report.txt"
    write_taxonomic_report(before, after, segmented=False, path=path)
    diversity = _diversity_rows(path.read_text())
    # Only "Sp A" counts; blanks excluded.
    assert diversity["species"] == ["1", "1"]


def test_taxonomic_report_handles_no_taxonomy(tmp_path):
    """A run with no taxonomy at any rank still writes a valid report."""
    before = [
        Sequence(id="a", header=">a", sequence="ACGT", seq_type=SequenceType.NUCLEOTIDE),
    ]
    path = tmp_path / "x_taxonomic_report.txt"
    write_taxonomic_report(before, [], segmented=False, path=path)
    text = path.read_text()
    assert "(no taxonomy available at any rank)" in text


# ---------------------------------------------------------------------------
# write_nucleotide_taxonomic_report
# ---------------------------------------------------------------------------

def _nt_seq(seq_id, *, length, genus=None, species=None, family=None):
    tax = TaxonomyInfo(species=species, genus=genus, family=family)
    return Sequence(
        id=seq_id, header=f">{seq_id}", sequence="A" * length,
        seq_type=SequenceType.NUCLEOTIDE, taxonomy=tax,
    )


def test_nucleotide_report_non_segmented_emits_genome_column(tmp_path):
    """Non-segmented mode: one `genome` column with NT-length stats per rank."""
    before = [
        _nt_seq("a", length=27000, genus="Alphacoronavirus", family="Coronaviridae"),
        _nt_seq("b", length=28000, genus="Alphacoronavirus", family="Coronaviridae"),
        _nt_seq("c", length=29000, genus="Alphacoronavirus", family="Coronaviridae"),
        _nt_seq("d", length=30000, genus="Betacoronavirus",  family="Coronaviridae"),
    ]
    after = [
        _nt_seq("a", length=27000, genus="Alphacoronavirus", family="Coronaviridae"),
        _nt_seq("d", length=30000, genus="Betacoronavirus",  family="Coronaviridae"),
    ]
    path = tmp_path / "x_nucleotide_taxonomic_report.txt"
    cfg = {"output": {"protein_report": {"max_breakdown": 20}}}
    assert write_nucleotide_taxonomic_report(
        before, after, cfg, segmented=False, path=path,
    )
    text = path.read_text()
    # Header + units.
    assert "Nucleotide taxonomic report" in text
    assert "Counting unit: sequences" in text
    # `genus` rank table exists with a `genome` column.
    assert "genus (2 distinct):" in text
    assert "genome" in text
    # Post-QC pool stats for Alphacoronavirus genome lengths
    # 27000/28000/29000 → min 27000, max 29000, median 28000, n=3.
    assert "27000, 29000, 28000, " in text  # min, max, median, ...
    assert ", 3" in text  # n=3 contributors


def test_nucleotide_report_segmented_emits_per_segment_and_total(tmp_path):
    """Segmented mode: columns S/M/L + a trailing `total` column."""
    # Build two CONCAT isolates of the same hanta genus, each with three
    # per-segment Sequence objects feeding `concat_segments`.
    def isolate(iso_id, s_len, m_len, l_len, *, genus="Orthohantavirus"):
        tax = TaxonomyInfo(genus=genus, family="Hantaviridae")
        s = Sequence(id=f"{iso_id}.S", header="", sequence="A" * s_len,
                     seq_type=SequenceType.NUCLEOTIDE, segment="S")
        m = Sequence(id=f"{iso_id}.M", header="", sequence="A" * m_len,
                     seq_type=SequenceType.NUCLEOTIDE, segment="M")
        l_ = Sequence(id=f"{iso_id}.L", header="", sequence="A" * l_len,
                      seq_type=SequenceType.NUCLEOTIDE, segment="L")
        return Sequence(
            id=f"CONCAT|{iso_id}",
            header=f"CONCAT|{iso_id}",
            sequence=("A" * (s_len + m_len + l_len)),
            seq_type=SequenceType.NUCLEOTIDE,
            isolate_id=iso_id,
            taxonomy=tax,
            concat_segments=[s, m, l_],
        )

    before = [
        isolate("iso1", 1850, 3650, 6500),
        isolate("iso2", 1900, 3700, 6600),
    ]
    after = [isolate("iso1", 1850, 3650, 6500)]
    cfg = {
        "output": {"protein_report": {"max_breakdown": 20}},
        "segmented": {
            "enabled": True,
            "virus": "hantaviridae",
            "viruses": {"hantaviridae": {"segments": ["S", "M", "L"]}},
        },
    }
    path = tmp_path / "x_nucleotide_taxonomic_report.txt"
    assert write_nucleotide_taxonomic_report(
        before, after, cfg, segmented=True, path=path,
    )
    text = path.read_text()
    assert "Counting unit: isolates" in text
    # Column order S, M, L, total.
    header_line = next(
        ln for ln in text.splitlines()
        if "taxon" in ln and "S" in ln and "M" in ln and "L" in ln
    )
    s_idx = header_line.index(" S")
    m_idx = header_line.index(" M")
    l_idx = header_line.index(" L")
    t_idx = header_line.index("total")
    assert s_idx < m_idx < l_idx < t_idx
    # The L column should carry n=2 contributors in the post-QC table.
    # (Min 6500, max 6600, median 6550, IQR 0 for two values? statistics
    # gives Q3-Q1 = 50 for [6500,6600] under exclusive quantiles —
    # check `n=2` is what's reported regardless of IQR rounding.)
    assert "6500" in text and "6600" in text
    # Total column for Orthohantavirus, n=2 isolates.
    total_min = 1850 + 3650 + 6500  # 12000
    total_max = 1900 + 3700 + 6600  # 12200
    assert str(total_min) in text and str(total_max) in text


def test_nucleotide_report_skips_species_rank(tmp_path):
    """Species rank is intentionally skipped (matches the protein report)."""
    before = [
        _nt_seq("a", length=1000, species="Sp A", genus="Gen X"),
        _nt_seq("b", length=1100, species="Sp A", genus="Gen X"),
    ]
    cfg = {"output": {"protein_report": {"max_breakdown": 20}}}
    path = tmp_path / "x_nucleotide_taxonomic_report.txt"
    write_nucleotide_taxonomic_report(before, [], cfg, segmented=False, path=path)
    text = path.read_text()
    # No `species (...)` rank header; genus is present.
    assert "species (" not in text
    assert "genus (1 distinct):" in text


def test_nucleotide_report_truncates_rank_above_max_breakdown(tmp_path):
    """A rank with > max_breakdown distinct taxa shows top-N + a 'more taxa' line."""
    before = [
        _nt_seq(f"s{i}", length=1000 + i, genus=f"Genus{i:02d}")
        for i in range(25)
    ]
    cfg = {"output": {"protein_report": {"max_breakdown": 5}}}
    path = tmp_path / "x_nucleotide_taxonomic_report.txt"
    write_nucleotide_taxonomic_report(before, [], cfg, segmented=False, path=path)
    text = path.read_text()
    assert "genus (25 distinct, top 5 by member count shown):" in text
    assert "... +20 more taxa not shown" in text


def test_nucleotide_report_handles_no_taxonomy(tmp_path):
    """No taxonomy at any rank still writes a valid report."""
    before = [
        Sequence(id="a", header=">a", sequence="ACGT", seq_type=SequenceType.NUCLEOTIDE),
    ]
    cfg = {"output": {"protein_report": {"max_breakdown": 20}}}
    path = tmp_path / "x_nucleotide_taxonomic_report.txt"
    assert write_nucleotide_taxonomic_report(before, [], cfg, segmented=False, path=path)
    text = path.read_text()
    assert "(no taxonomy available at any rank from subgenus to class)" in text


# ---------------------------------------------------------------------------
# Cluster purity diagnostic columns
# ---------------------------------------------------------------------------


def _cluster_seq(sid, *, length=1000, species=None, genus=None, host=None):
    tax = None
    if species or genus:
        tax = TaxonomyInfo(species=species, genus=genus)
    return Sequence(
        id=sid, header=sid, sequence="A" * length,
        seq_type=SequenceType.NUCLEOTIDE, accession=sid,
        host=host, taxonomy=tax,
    )


def test_write_cluster_tsv_emits_purity_columns(tmp_path):
    rep = _cluster_seq("acc1", length=1000, species="Foo virus", genus="Foovirus", host="human")
    members = [
        _cluster_seq("acc2", length=900, species="Foo virus", genus="Foovirus", host="human"),
        _cluster_seq("acc3", length=1100, species="Foo virus", genus="Foovirus", host="human"),
    ]
    result = RunResult(
        mode="global", representatives=[rep],
        clusters=[Cluster(cluster_id="c1", representative=rep, members=members)],
    )
    path = tmp_path / "c.tsv"
    write_cluster_tsv(result, path)
    lines = path.read_text().splitlines()
    header = lines[0].split("\t")
    assert "species_purity" in header
    assert "genus_purity" in header
    assert "host_purity" in header
    assert "member_length_min" in header
    assert "member_length_max" in header
    assert "member_length_median" in header

    row = dict(zip(header, lines[1].split("\t")))
    # Every member agrees → purity 1.0, distinct 1.
    assert row["species_purity"] == "1.000"
    assert row["species_distinct"] == "1"
    assert row["genus_purity"] == "1.000"
    assert row["host_purity"] == "1.000"
    assert row["member_length_min"] == "900"
    assert row["member_length_max"] == "1100"
    assert row["member_length_median"] == "1000"


def test_write_cluster_tsv_flags_mixed_genus_cluster(tmp_path):
    """A cluster spanning two genera should report < 1.0 purity and
    distinct=2 — the actionable diagnostic for "threshold too lax"."""
    rep = _cluster_seq("acc1", genus="Foovirus")
    members = [
        _cluster_seq("acc2", genus="Foovirus"),
        _cluster_seq("acc3", genus="Barvirus"),  # the odd one out
    ]
    result = RunResult(
        mode="global", representatives=[rep],
        clusters=[Cluster(cluster_id="c1", representative=rep, members=members)],
    )
    path = tmp_path / "c.tsv"
    write_cluster_tsv(result, path)
    header, row_str = path.read_text().splitlines()
    row = dict(zip(header.split("\t"), row_str.split("\t")))
    # 2 of 3 share the same genus.
    assert row["genus_purity"] == "0.667"
    assert row["genus_distinct"] == "2"


def test_write_cluster_tsv_blank_purity_when_no_taxonomy(tmp_path):
    """Cluster with zero populated values at a rank → blank purity, distinct=0."""
    rep = _cluster_seq("acc1")  # no taxonomy
    result = RunResult(
        mode="global", representatives=[rep],
        clusters=[Cluster(cluster_id="c1", representative=rep, members=[])],
    )
    path = tmp_path / "c.tsv"
    write_cluster_tsv(result, path)
    header, row_str = path.read_text().splitlines()
    row = dict(zip(header.split("\t"), row_str.split("\t")))
    assert row["genus_purity"] == ""
    assert row["genus_distinct"] == "0"
    assert row["host_purity"] == ""


def test_write_cluster_tsv_ignores_empty_values_in_purity(tmp_path):
    """One member with no host shouldn't drag down host purity — empty
    values are excluded from both numerator and denominator."""
    rep = _cluster_seq("acc1", host="human")
    members = [
        _cluster_seq("acc2", host="human"),
        _cluster_seq("acc3", host=None),  # no host annotation
    ]
    result = RunResult(
        mode="global", representatives=[rep],
        clusters=[Cluster(cluster_id="c1", representative=rep, members=members)],
    )
    path = tmp_path / "c.tsv"
    write_cluster_tsv(result, path)
    header, row_str = path.read_text().splitlines()
    row = dict(zip(header.split("\t"), row_str.split("\t")))
    # 2 of 2 populated values agree → purity 1.0.
    assert row["host_purity"] == "1.000"
    assert row["host_distinct"] == "1"


# ---------------------------------------------------------------------------
# HMM diagnostic TSV (#4)
# ---------------------------------------------------------------------------

def _hmm_hit(target, *, e, score=100.0, ali_from=1, ali_to=200, hmm_from=1,
             hmm_to=200, hmm_len=250, passing=True):
    return {
        "target": target, "dom_evalue": e, "dom_score": score,
        "ali_from": ali_from, "ali_to": ali_to,
        "hmm_from": hmm_from, "hmm_to": hmm_to, "hmm_len": hmm_len,
        "passing": passing,
    }


def _rep_with_hmm(acc, *, hits):
    s = Sequence(
        id=acc, header=acc, sequence="ACGT" * 100,
        seq_type=SequenceType.NUCLEOTIDE, accession=acc,
    )
    s.proteins = [{
        "protein_id": f"p_{acc}", "product": "marker",
        "sequence": "M" * 300, "length": 300, "hmm_hits": list(hits),
    }]
    return s


def test_hmm_diagnostic_emits_row_per_rep_per_profile(tmp_path):
    """One declared HMM-bearing spec, two reps → 2 rows (single profile)."""
    cfg = {
        "segmented": {"enabled": False},
        "clustering": {"cluster_protein": [
            {"name": "RdRp", "hmms": ["RdRP_4"]},
        ]},
        "_hmm_runtime": {
            "active": True,
            "ga_cutoffs": {"RdRP_4": 25.0},
            "hmm_cfg": {
                "use_ga_when_available": True,
                "default_evalue": 1e-5,
                "relative_length_cutoff": 0.5,
            },
        },
    }
    reps = [
        _rep_with_hmm("acc1", hits=[_hmm_hit("RdRP_4", e=1e-30, score=80.0)]),
        _rep_with_hmm("acc2", hits=[_hmm_hit("RdRP_4", e=1e-3, score=15.0,
                                              passing=False)]),
    ]
    result = RunResult(mode="global", representatives=reps)
    path = tmp_path / "hmm.tsv"
    assert write_hmm_diagnostic_tsv(result, cfg, None, path)
    lines = path.read_text().splitlines()
    assert len(lines) == 3  # header + 2 reps
    header = lines[0].split("\t")
    for col in ("spec_name", "profile", "hit", "best_dom_evalue",
                "ga_cutoff", "cutoff_used", "passing"):
        assert col in header
    row1 = dict(zip(header, lines[1].split("\t")))
    row2 = dict(zip(header, lines[2].split("\t")))
    assert row1["accession"] == "acc1"
    assert row1["profile"] == "RdRP_4"
    assert row1["hit"] == "TRUE"
    assert row1["passing"] == "TRUE"
    assert row1["cutoff_used"] == "GA=25"
    assert row1["ga_cutoff"] == "25.0"
    assert row2["passing"] == "FALSE"


def test_hmm_diagnostic_emits_hit_false_for_missing_profile(tmp_path):
    """When a CDS has no hit on a declared profile, the row still appears
    with hit=FALSE and the stats columns blank — this is the silent-
    failure visibility the report exists to provide."""
    cfg = {
        "segmented": {"enabled": False},
        "clustering": {"cluster_protein": [
            {"name": "RdRp", "hmms": ["RdRP_4"]},
            {"name": "Cap", "hmms": ["Capsid_X"]},
        ]},
        "_hmm_runtime": {
            "active": True,
            "ga_cutoffs": {},
            "hmm_cfg": {"default_evalue": 1e-5},
        },
    }
    reps = [_rep_with_hmm("acc1", hits=[_hmm_hit("RdRP_4", e=1e-30)])]
    result = RunResult(mode="global", representatives=reps)
    path = tmp_path / "hmm.tsv"
    assert write_hmm_diagnostic_tsv(result, cfg, None, path)
    lines = path.read_text().splitlines()
    # 1 rep × 2 specs × 1 profile each = 2 rows.
    assert len(lines) == 3
    header = lines[0].split("\t")
    rows = [dict(zip(header, l.split("\t"))) for l in lines[1:]]
    by_profile = {r["profile"]: r for r in rows}
    assert by_profile["RdRP_4"]["hit"] == "TRUE"
    assert by_profile["Capsid_X"]["hit"] == "FALSE"
    assert by_profile["Capsid_X"]["best_dom_evalue"] == ""
    assert by_profile["Capsid_X"]["passing"] == ""


def test_hmm_diagnostic_expands_multidomain_token_to_components(tmp_path):
    """A multidomain token "A--B" must produce one row per HMM
    (A, B) — diagnostic granularity is per profile, not per token."""
    cfg = {
        "segmented": {"enabled": False},
        "clustering": {"cluster_protein": [
            {"name": "Spike", "hmms": ["CoV_S1--CoV_S2"]},
        ]},
        "_hmm_runtime": {"active": True, "ga_cutoffs": {}, "hmm_cfg": {}},
    }
    reps = [_rep_with_hmm("acc1", hits=[
        _hmm_hit("CoV_S1", e=1e-50, ali_from=10, ali_to=600),
        _hmm_hit("CoV_S2", e=1e-40, ali_from=610, ali_to=1200),
    ])]
    result = RunResult(mode="global", representatives=reps)
    path = tmp_path / "hmm.tsv"
    write_hmm_diagnostic_tsv(result, cfg, None, path)
    lines = path.read_text().splitlines()
    profiles = {dict(zip(lines[0].split("\t"), l.split("\t")))["profile"]
                for l in lines[1:]}
    assert profiles == {"CoV_S1", "CoV_S2"}


def test_hmm_diagnostic_skips_when_hmm_tier_inactive(tmp_path):
    cfg = {
        "segmented": {"enabled": False},
        "clustering": {"cluster_protein": [
            {"name": "RdRp", "hmms": ["RdRP_4"]},
        ]},
        # No _hmm_runtime.active flag, no populated hmm_hits.
    }
    reps = [Sequence(id="acc1", header="acc1", sequence="ACGT",
                     seq_type=SequenceType.NUCLEOTIDE, accession="acc1")]
    result = RunResult(mode="global", representatives=reps)
    path = tmp_path / "hmm.tsv"
    assert write_hmm_diagnostic_tsv(result, cfg, None, path) is False
    assert not path.exists()


def test_hmm_diagnostic_skips_when_no_specs_declare_hmms(tmp_path):
    cfg = {
        "segmented": {"enabled": False},
        # Alias-only spec — no hmms declared.
        "clustering": {"cluster_protein": [{"name": "RdRp", "aliases": ["polymerase"]}]},
        "_hmm_runtime": {"active": True, "ga_cutoffs": {}, "hmm_cfg": {}},
    }
    reps = [_rep_with_hmm("acc1", hits=[])]
    result = RunResult(mode="global", representatives=reps)
    path = tmp_path / "hmm.tsv"
    assert write_hmm_diagnostic_tsv(result, cfg, None, path) is False
