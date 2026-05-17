"""Brittleness probes: unusual characters in sequence names and metadata.

Biological FASTA headers carry strain names, host descriptors, country
strings, and isolate identifiers that — in practice — contain characters
the pipeline did not originally plan for: whitespace beyond the regular
space, the pipe used as our own internal separator, tabs that would
break a TSV row, embedded newlines from poorly-cleaned data, brackets,
slashes, etc.

These tests document and lock in safe behavior at every surface where
those characters used to (or could) corrupt outputs:

  * FASTA writer headers must not produce phantom records (newline / CR)
  * MMseqs2 cluster round-trip must not silently drop sequences whose
    seq.id contains whitespace (MMseqs2 truncates header tokens)
  * CONCAT|<isolate_id> parsing must round-trip even when the isolate
    name contains a pipe
  * TSV writers must keep every row to its expected column count even
    when a field value contains tabs or newlines
"""
from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from repseq.clustering.mmseqs2 import (
    MMseqs2Error,
    _parse_cluster_tsv,
    _write_id_fasta,
    run_clustering,
)
from repseq.io.fasta import read_fasta, write_fasta
from repseq.models import (
    Cluster,
    GroupStat,
    QCReport,
    RunResult,
    Sequence,
    SequenceType,
)
from repseq.output.report import (
    write_cluster_tsv,
    write_group_counts_tsv,
    write_qc_tsv,
    write_representative_sequences_tsv,
)
from repseq.segmented.completeness import (
    _normalise_isolate_id,
    concatenate_isolate,
)


# Characters that have caused or could plausibly cause issues. Grouped
# by failure mode so the parametrized IDs read well in pytest output.
WHITESPACE_CHARS = [
    ("space", " "),
    ("tab", "\t"),
    ("newline", "\n"),
    ("carriage", "\r"),
    ("nbsp", "\xa0"),
]

PUNCTUATION_CHARS = [
    ("pipe", "|"),
    ("forward", "/"),
    ("backslash", "\\"),
    ("open-paren", "("),
    ("close-paren", ")"),
    ("open-brack", "["),
    ("close-brack", "]"),
    ("single-q", "'"),
    ("double-q", '"'),
    ("comma", ","),
    ("semicolon", ";"),
    ("colon", ":"),
    ("equals", "="),
    ("less-than", "<"),
    ("greater-than", ">"),
    ("ampersand", "&"),
    ("hash", "#"),
    ("question", "?"),
    ("star", "*"),
    ("at", "@"),
]

ALL_CHARS = WHITESPACE_CHARS + PUNCTUATION_CHARS


def _seq(sid, header=None, accession=None):
    return Sequence(
        id=sid,
        header=header if header is not None else sid,
        sequence="ACGTACGTACGT" * 20,
        seq_type=SequenceType.NUCLEOTIDE,
        accession=accession or sid,
    )


# ---------------------------------------------------------------------------
# FASTA writer must not produce phantom records when a header carries an
# embedded newline / carriage return. Other unusual characters should
# pass through unharmed.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("name,ch", ALL_CHARS)
def test_fasta_write_read_roundtrip_single_record(name, ch, tmp_path):
    sequence = "ACGTACGTACGT" * 20
    header = f"acc1 description with {ch} embedded"
    seq = Sequence(
        id="acc1", header=header, sequence=sequence,
        seq_type=SequenceType.NUCLEOTIDE, accession="acc1",
    )
    p = tmp_path / "x.fa"
    write_fasta([seq], p)

    records = list(read_fasta(p))
    # The critical invariants: one input sequence yields one output record,
    # AND the original biological sequence comes back uncorrupted.
    # Embedded newlines in a header would otherwise leak the post-newline
    # description fragment into the sequence body, since the reader treats
    # any non-`>` line as sequence content.
    assert len(records) == 1, (
        f"Header containing {name!r} produced {len(records)} records "
        f"(expected 1). File contents:\n{p.read_text()!r}"
    )
    assert records[0].sequence == sequence, (
        f"Header containing {name!r} corrupted the sequence body. "
        f"Got {records[0].sequence!r}, expected {sequence!r}."
    )


# ---------------------------------------------------------------------------
# MMseqs2 cluster round-trip: every input sequence must be accounted for
# in the parsed cluster list. The round-trip check inside run_clustering
# enforces this — but it needs the underlying _write_id_fasta to keep
# seq.id intact as the FASTA header token (MMseqs2 splits on any
# whitespace and takes the first token as the identifier).
# ---------------------------------------------------------------------------

def _fake_mmseqs(headers_to_tsv):
    """Return a subprocess.run replacement that drops `headers_to_tsv`
    into the cluster TSV the wrapper expects to read back.
    """
    def _runner(cmd, **kwargs):
        result_prefix = cmd[3]
        tsv = Path(result_prefix + "_cluster.tsv")
        tsv.write_text("".join(f"{h}\t{h}\n" for h in headers_to_tsv))
        class _R:
            stderr = ""
            stdout = ""
        return _R()
    return _runner


@pytest.mark.parametrize("name,ch", PUNCTUATION_CHARS)
def test_run_clustering_handles_punctuation_in_id(name, ch, tmp_path):
    # Non-whitespace punctuation in seq.id should pass straight through
    # MMseqs2 (which splits on whitespace) — every sequence must come back.
    sid = f"acc{ch}1"
    seqs = [_seq(sid), _seq(f"other{ch}id")]

    # Simulate MMseqs2 emitting both IDs as their own clusters.
    runner = _fake_mmseqs([sid, f"other{ch}id"])
    with patch("repseq.clustering.mmseqs2._check_mmseqs2", return_value="mmseqs"), \
         patch("repseq.clustering.mmseqs2.subprocess.run", side_effect=runner):
        clusters = run_clustering(seqs, 0.9, {"temp_dir": str(tmp_path)})

    assert len(clusters) == 2
    assert {c.representative.id for c in clusters} == {sid, f"other{ch}id"}


@pytest.mark.parametrize("name,ch", WHITESPACE_CHARS)
def test_run_clustering_raises_on_whitespace_in_id(name, ch, tmp_path):
    # Whitespace in seq.id is fatal for the round-trip: MMseqs2 keeps
    # only the leading token, so the parser cannot reconnect cluster
    # rows to input sequences. The wrapper's sanity check must catch
    # this and raise rather than silently returning zero clusters.
    sid = f"acc{ch}1"
    seqs = [_seq(sid)]

    # Simulate MMseqs2 emitting the truncated token.
    truncated = sid.split()[0] if ch.strip() == "" else sid.split(ch, 1)[0]
    runner = _fake_mmseqs([truncated])

    with patch("repseq.clustering.mmseqs2._check_mmseqs2", return_value="mmseqs"), \
         patch("repseq.clustering.mmseqs2.subprocess.run", side_effect=runner):
        with pytest.raises(MMseqs2Error, match="round-trip"):
            run_clustering(seqs, 0.9, {"temp_dir": str(tmp_path)})


# ---------------------------------------------------------------------------
# CONCAT|<isolate_id> identity must survive an isolate name that contains
# any of the characters the regex might capture. The normalisation step
# is what makes this safe — both whitespace and our own internal pipe
# separator must be neutralised.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("name,ch", ALL_CHARS)
def test_normalise_isolate_id_yields_safe_id(name, ch):
    # The normalised isolate_id is propagated into seq.id as
    # CONCAT|<isolate_id>. It must contain neither whitespace
    # (MMseqs2 truncation) nor a pipe (split-based CONCAT parsers).
    raw = f"abc{ch}def"
    out = _normalise_isolate_id(raw)
    assert " " not in out and "\t" not in out and "\n" not in out and "\r" not in out
    assert "\xa0" not in out
    assert "|" not in out


@pytest.mark.parametrize("name,ch", ALL_CHARS)
def test_concat_id_split_recovers_isolate(name, ch):
    # Both the segmented FASTA writer and the protein-fasta report
    # extract the isolate id from CONCAT seq.id via
    # `seq.id.split("|")[1]`. The normalisation step keeps that
    # operation honest no matter what character the captured isolate
    # name contained.
    raw_iso = f"strain{ch}variant"
    normalised = _normalise_isolate_id(raw_iso)

    concat = concatenate_isolate(
        [_seq("acc1"), _seq("acc2")],
        normalised,
    )
    parts = concat.id.split("|")
    assert len(parts) >= 2
    assert parts[0] == "CONCAT"
    # The slot the writer reads must equal the normalised id we passed in.
    assert parts[1] == normalised
    # And it must be the *whole* id segment — no truncation by '|'.
    assert concat.id == f"CONCAT|{normalised}"


# ---------------------------------------------------------------------------
# TSV writers: every row must keep its expected column count even when
# field values contain tabs or embedded newlines. We don't dictate how
# the writer escapes such characters (replacement, quoting) — we only
# verify that the row structure survives.
# ---------------------------------------------------------------------------

def _tsv_row_shape(text: str) -> tuple[int, list[int]]:
    """Return (n_data_rows, tab_counts_per_row). Trailing blank lines are dropped."""
    lines = [ln for ln in text.splitlines() if ln != ""]
    # The first line is the header row — count tabs from every line.
    return len(lines), [ln.count("\t") for ln in lines]


@pytest.mark.parametrize("name,ch", ALL_CHARS)
def test_representative_tsv_row_shape_survives(name, ch, tmp_path):
    seq = _seq("acc1")
    seq.organism = f"genus{ch}species"
    seq.description = f"desc with {ch}"
    seq.host = f"host{ch}"
    seq.country = f"country{ch}"
    p = tmp_path / "reps.tsv"
    write_representative_sequences_tsv([seq], p)

    n_rows, tab_counts = _tsv_row_shape(p.read_text())
    # 1 header + 1 data row, every row must have the same tab count.
    assert n_rows == 2, (
        f"Field with {name!r} produced {n_rows} rows (expected 2). "
        f"Contents:\n{p.read_text()!r}"
    )
    assert tab_counts[0] == tab_counts[1], (
        f"Field with {name!r} broke column alignment: "
        f"header={tab_counts[0]} tabs, row={tab_counts[1]} tabs."
    )


@pytest.mark.parametrize("name,ch", ALL_CHARS)
def test_cluster_tsv_row_shape_survives(name, ch, tmp_path):
    rep = _seq("acc1")
    rep.organism = f"genus{ch}species"
    result = RunResult(
        mode="global",
        representatives=[rep],
        clusters=[Cluster(cluster_id=f"c|{rep.id}", representative=rep, members=[])],
    )
    p = tmp_path / "clusters.tsv"
    write_cluster_tsv(result, p)

    n_rows, tab_counts = _tsv_row_shape(p.read_text())
    assert n_rows == 2
    assert tab_counts[0] == tab_counts[1]


@pytest.mark.parametrize("name,ch", ALL_CHARS)
def test_group_counts_tsv_row_shape_survives(name, ch, tmp_path):
    result = RunResult(
        mode="taxonomic1",
        representatives=[],
        clusters=[],
        group_stats=[
            GroupStat(
                grouping="genus",
                group=f"Some{ch}Genus",
                n_before=10,
                n_after=5,
                clustered=True,
                cutoff=0.9,
            )
        ],
    )
    p = tmp_path / "gc.tsv"
    assert write_group_counts_tsv(result, p)

    n_rows, tab_counts = _tsv_row_shape(p.read_text())
    assert n_rows == 2
    assert tab_counts[0] == tab_counts[1]


@pytest.mark.parametrize("name,ch", ALL_CHARS)
def test_qc_removed_tsv_row_shape_survives(name, ch, tmp_path):
    report = QCReport()
    report.add_removed(f"acc{ch}1", f"failed because {ch} reason")
    p = tmp_path / "qc.tsv"
    write_qc_tsv(report, p)

    n_rows, tab_counts = _tsv_row_shape(p.read_text())
    assert n_rows >= 2
    # Every row must have the same number of tabs as the header.
    assert all(t == tab_counts[0] for t in tab_counts), (
        f"qc_removed.tsv lost column alignment with {name!r}: {tab_counts}"
    )
