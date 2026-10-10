"""Tests for NC-compatible candidate-region construction."""

from or_mining.blast import parse_blast_row
from or_mining.regions import build_candidate_regions


def _hit(query: str, start: int, end: int, frame: int, slen: int = 1000) -> str:
    return (
        f"{query}\tChr1\t50\t40\t0\t0\t1\t40\t{start}\t{end}"
        f"\t1e-10\t80\t0\t{frame}\t100\t{slen}"
    )


def test_reduces_adjacent_hits_on_same_strand_only() -> None:
    hits = [
        parse_blast_row(_hit("q1", 100, 199, 1)),
        parse_blast_row(_hit("q2", 200, 299, 2)),
        parse_blast_row(_hit("q3", 299, 200, -1)),
    ]

    regions, stats = build_candidate_regions(hits, min_width=100, flank=50)

    assert len(regions) == 2
    assert (regions[0].strand, regions[0].hit_start_0, regions[0].hit_end_0) == (
        "+",
        99,
        299,
    )
    assert regions[0].hit_count == 2
    assert regions[0].query_count == 2
    assert regions[1].strand == "-"
    assert stats["candidate_regions"] == 2


def test_filters_before_flanking_and_clips_scaffold_boundaries() -> None:
    hits = [
        parse_blast_row(_hit("short", 1, 50, 1)),
        parse_blast_row(_hit("kept", 901, 1000, 1)),
    ]

    regions, stats = build_candidate_regions(hits, min_width=100, flank=1000)

    assert len(regions) == 1
    assert (regions[0].expanded_start_0, regions[0].expanded_end_0) == (0, 1000)
    assert stats["filtered_short_regions"] == 1
    assert stats["right_boundary_clipped"] == 1

