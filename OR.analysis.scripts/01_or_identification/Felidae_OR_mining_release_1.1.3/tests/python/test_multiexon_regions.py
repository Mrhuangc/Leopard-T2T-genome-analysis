from or_mining.blast import BlastHit
from or_mining.multiexon_regions import KnownLocus, build_multiexon_regions


def hit(query: str, start: int, end: int, strand: int = 1, slen: int = 1000) -> BlastHit:
    sstart, send = (start, end) if strand > 0 else (end, start)
    return BlastHit(
        query, "chr1", 90.0, 20, 0, 0, 1, 20, sstart, send,
        1e-20, 80.0, 0, strand, 20, slen,
    )


def test_excludes_confirmed_and_merges_extended_regions_across_strands() -> None:
    hits = [
        hit("known", 101, 150),
        hit("left", 301, 330),
        hit("right", 440, 470, -1),
    ]
    regions, counts = build_multiexon_regions(
        hits,
        [KnownLocus("chr1", 100, 150)],
        maximum_intron_length=200,
        gene_clearance=50,
    )
    assert counts["regions_overlapping_confirmed_or"] == 1
    assert counts["left_clipped_by_confirmed_or"] == 1
    assert len(regions) == 1
    assert (regions[0].start_0, regions[0].end_0) == (199, 670)
    assert regions[0].source_interval_count == 2


def test_clips_to_subject_boundaries() -> None:
    regions, counts = build_multiexon_regions(
        [hit("edge", 10, 30, slen=100)], [], maximum_intron_length=80
    )
    assert (regions[0].start_0, regions[0].end_0) == (0, 100)
    assert counts["left_clipped_by_sequence_boundary"] == 1
    assert counts["right_clipped_by_sequence_boundary"] == 1
