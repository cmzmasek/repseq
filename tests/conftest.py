"""Shared fixtures for repseq regression tests."""
from __future__ import annotations

from pathlib import Path

import pytest

from repseq.models import Sequence, SequenceSource, SequenceType, TaxonomyInfo


def _seq(
    sid: str,
    sequence: str,
    *,
    header: str | None = None,
    source: SequenceSource = SequenceSource.NCBI,
    seq_type: SequenceType = SequenceType.PROTEIN,
    accession: str | None = None,
    is_refseq: bool = False,
    is_reviewed: bool = False,
    organism: str | None = None,
    host: str | None = None,
    country: str | None = None,
    collection_date: str | None = None,
    segment: str | None = None,
    isolate_id: str | None = None,
    taxonomy: TaxonomyInfo | None = None,
) -> Sequence:
    return Sequence(
        id=sid,
        header=header or sid,
        sequence=sequence,
        seq_type=seq_type,
        source=source,
        accession=accession or sid,
        organism=organism,
        host=host,
        country=country,
        collection_date=collection_date,
        segment=segment,
        isolate_id=isolate_id,
        is_refseq=is_refseq,
        is_reviewed=is_reviewed,
        taxonomy=taxonomy,
    )


@pytest.fixture
def make_seq():
    """Factory for building a Sequence with sensible defaults."""
    return _seq


@pytest.fixture
def tmp_cache_dir(tmp_path: Path) -> Path:
    d = tmp_path / "cache"
    d.mkdir()
    return d
