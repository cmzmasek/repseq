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

from repseq.output.report import (
    write_all_reports,
    write_group_counts_tsv,
    write_representative_isolates_tsv,
    write_representative_sequences_tsv,
    write_run_log,
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
        "isolate_id", "organism", "strain", "host", "collection_date",
        "country", "n_segments", "segments", "accessions",
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


def test_write_representative_sequences_tsv_keeps_per_sequence_schema(tmp_path):
    """Non-segmented writer keeps the original per-sequence columns
    including `accession`, `segment`, `description`, `molecule_type`,
    `length_nt`."""
    rep = Sequence(
        id="ACC.1", header=">ACC.1", sequence="ACGT" * 100,
        accession="ACC.1", organism="Some virus",
        description="hypothetical", seq_type=SequenceType.NUCLEOTIDE,
    )
    path = tmp_path / "x_representative_sequences.tsv"
    write_representative_sequences_tsv([rep], path)

    header = path.read_text().splitlines()[0].split("\t")
    for col in ("accession", "segment", "description", "molecule_type",
                "length_nt"):
        assert col in header
    # Isolate-only columns must NOT appear.
    for col in ("n_segments", "segments", "accessions", "total_length_nt"):
        assert col not in header


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
