"""Per-isolate species-consistency filter for segmented viruses."""
from __future__ import annotations

from repseq.models import TaxonomyInfo
from repseq.segmented.taxonomy_consistency import (
    filter_taxonomy_consistent_isolates,
)


def _tax(species: str | None = None, genus: str | None = None) -> TaxonomyInfo:
    return TaxonomyInfo(species=species, genus=genus)


def test_consistent_species_kept(make_seq):
    seqs = [
        make_seq("acc1", "AAA", isolate_id="iso1", segment="L",
                 taxonomy=_tax(species="Bunyamwera virus")),
        make_seq("acc2", "CCC", isolate_id="iso1", segment="M",
                 taxonomy=_tax(species="Bunyamwera virus")),
        make_seq("acc3", "GGG", isolate_id="iso1", segment="S",
                 taxonomy=_tax(species="Bunyamwera virus")),
    ]
    kept, removed = filter_taxonomy_consistent_isolates(seqs, rank="species")
    assert len(kept) == 3
    assert removed == []


def test_two_distinct_species_drops_whole_isolate(make_seq):
    seqs = [
        make_seq("acc1", "AAA", isolate_id="iso1", segment="L",
                 taxonomy=_tax(species="Bunyamwera virus")),
        make_seq("acc2", "CCC", isolate_id="iso1", segment="M",
                 taxonomy=_tax(species="Cache Valley virus")),
        make_seq("acc3", "GGG", isolate_id="iso1", segment="S",
                 taxonomy=_tax(species="Bunyamwera virus")),
    ]
    kept, removed = filter_taxonomy_consistent_isolates(seqs, rank="species")
    assert kept == []
    accs = [pair[0] for pair in removed]
    reasons = {pair[1] for pair in removed}
    assert set(accs) == {"acc1", "acc2", "acc3"}
    assert reasons == {"taxonomy_mismatch:species"}


def test_mixed_populated_and_missing_populated_agree_kept(make_seq):
    seqs = [
        make_seq("acc1", "AAA", isolate_id="iso1", segment="L",
                 taxonomy=_tax(species="Bunyamwera virus")),
        make_seq("acc2", "CCC", isolate_id="iso1", segment="M",
                 taxonomy=_tax(species=None)),
        make_seq("acc3", "GGG", isolate_id="iso1", segment="S",
                 taxonomy=_tax(species="Bunyamwera virus")),
    ]
    kept, removed = filter_taxonomy_consistent_isolates(seqs, rank="species")
    assert len(kept) == 3
    assert removed == []


def test_all_missing_species_kept(make_seq):
    seqs = [
        make_seq("acc1", "AAA", isolate_id="iso1", segment="L",
                 taxonomy=_tax()),
        make_seq("acc2", "CCC", isolate_id="iso1", segment="M",
                 taxonomy=_tax()),
    ]
    kept, removed = filter_taxonomy_consistent_isolates(seqs, rank="species")
    assert len(kept) == 2
    assert removed == []


def test_single_segment_isolate_kept(make_seq):
    seqs = [
        make_seq("acc1", "AAA", isolate_id="iso1", segment="L",
                 taxonomy=_tax(species="Bunyamwera virus")),
    ]
    kept, removed = filter_taxonomy_consistent_isolates(seqs, rank="species")
    assert len(kept) == 1
    assert removed == []


def test_no_isolate_id_passes_through_untouched(make_seq):
    # Floaters (no isolate_id yet — e.g. UniProt input or accessions
    # that GenBank lookup couldn't fill) bypass the filter. The regex
    # fallback in filter_complete_isolates will catch them later.
    seqs = [
        make_seq("acc1", "AAA", segment="L",
                 taxonomy=_tax(species="Bunyamwera virus")),
        make_seq("acc2", "CCC", segment="M",
                 taxonomy=_tax(species="Cache Valley virus")),
    ]
    kept, removed = filter_taxonomy_consistent_isolates(seqs, rank="species")
    assert len(kept) == 2
    assert removed == []


def test_case_and_whitespace_insensitive(make_seq):
    # "bunyamwera  virus " (extra spaces, trailing space) and
    # "Bunyamwera virus" must compare equal — same species, just
    # noisy formatting.
    seqs = [
        make_seq("acc1", "AAA", isolate_id="iso1", segment="L",
                 taxonomy=_tax(species="Bunyamwera virus")),
        make_seq("acc2", "CCC", isolate_id="iso1", segment="M",
                 taxonomy=_tax(species="bunyamwera  virus ")),
    ]
    kept, removed = filter_taxonomy_consistent_isolates(seqs, rank="species")
    assert len(kept) == 2
    assert removed == []


def test_only_offending_isolate_dropped_other_kept(make_seq):
    seqs = [
        # iso1: consistent
        make_seq("acc1", "AAA", isolate_id="iso1", segment="L",
                 taxonomy=_tax(species="Bunyamwera virus")),
        make_seq("acc2", "CCC", isolate_id="iso1", segment="M",
                 taxonomy=_tax(species="Bunyamwera virus")),
        # iso2: mismatch — only this isolate should be dropped
        make_seq("acc3", "GGG", isolate_id="iso2", segment="L",
                 taxonomy=_tax(species="Bunyamwera virus")),
        make_seq("acc4", "TTT", isolate_id="iso2", segment="M",
                 taxonomy=_tax(species="Cache Valley virus")),
    ]
    kept, removed = filter_taxonomy_consistent_isolates(seqs, rank="species")
    kept_accs = {s.accession for s in kept}
    assert kept_accs == {"acc1", "acc2"}
    assert {pair[0] for pair in removed} == {"acc3", "acc4"}


def test_rank_genus_compares_genus(make_seq):
    # When the rank knob is set to genus, species mismatches that
    # share a genus must not trigger a drop.
    seqs = [
        make_seq("acc1", "AAA", isolate_id="iso1", segment="L",
                 taxonomy=_tax(species="Bunyamwera virus",
                               genus="Orthobunyavirus")),
        make_seq("acc2", "CCC", isolate_id="iso1", segment="M",
                 taxonomy=_tax(species="Cache Valley virus",
                               genus="Orthobunyavirus")),
    ]
    kept, removed = filter_taxonomy_consistent_isolates(seqs, rank="genus")
    assert len(kept) == 2
    assert removed == []


def test_input_order_preserved(make_seq):
    # Downstream steps (concat segment order, deterministic output)
    # rely on input order being stable.
    seqs = [
        make_seq("acc1", "AAA", isolate_id="iso1", segment="L",
                 taxonomy=_tax(species="Bunyamwera virus")),
        make_seq("acc2", "CCC", isolate_id="iso2", segment="L",
                 taxonomy=_tax(species="Cache Valley virus")),
        make_seq("acc3", "GGG", isolate_id="iso1", segment="M",
                 taxonomy=_tax(species="Bunyamwera virus")),
    ]
    kept, _ = filter_taxonomy_consistent_isolates(seqs, rank="species")
    assert [s.accession for s in kept] == ["acc1", "acc2", "acc3"]
