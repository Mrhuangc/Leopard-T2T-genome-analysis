"""Construct deterministic candidate regions from validated TBLASTN hits."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable
import json

from or_mining.blast import BlastHit, _open_text, parse_blast_row


@dataclass(frozen=True, slots=True)
class CandidateRegion:
    """A strand-aware reduced hit interval using zero-based half-open coordinates."""

    sseqid: str
    strand: str
    hit_start_0: int
    hit_end_0: int
    expanded_start_0: int
    expanded_end_0: int
    slen: int
    hit_count: int
    query_count: int
    best_evalue: float
    max_bitscore: float

    @property
    def hit_width(self) -> int:
        return self.hit_end_0 - self.hit_start_0

    @property
    def expanded_width(self) -> int:
        return self.expanded_end_0 - self.expanded_start_0


@dataclass(slots=True)
class _RegionAccumulator:
    sseqid: str
    strand: str
    start_0: int
    end_0: int
    slen: int
    hit_count: int = 0
    query_ids: set[str] = field(default_factory=set)
    best_evalue: float = float("inf")
    max_bitscore: float = 0.0

    def add(self, hit: BlastHit) -> None:
        self.end_0 = max(self.end_0, hit.genome_end_0)
        self.hit_count += 1
        self.query_ids.add(hit.qseqid)
        self.best_evalue = min(self.best_evalue, hit.evalue)
        self.max_bitscore = max(self.max_bitscore, hit.bitscore)


def _new_accumulator(hit: BlastHit) -> _RegionAccumulator:
    region = _RegionAccumulator(
        sseqid=hit.sseqid,
        strand=hit.strand,
        start_0=hit.genome_start_0,
        end_0=hit.genome_end_0,
        slen=hit.slen,
    )
    region.add(hit)
    return region


def _finalize_region(region: _RegionAccumulator, flank: int) -> CandidateRegion:
    return CandidateRegion(
        sseqid=region.sseqid,
        strand=region.strand,
        hit_start_0=region.start_0,
        hit_end_0=region.end_0,
        expanded_start_0=max(0, region.start_0 - flank),
        expanded_end_0=min(region.slen, region.end_0 + flank),
        slen=region.slen,
        hit_count=region.hit_count,
        query_count=len(region.query_ids),
        best_evalue=region.best_evalue,
        max_bitscore=region.max_bitscore,
    )


def build_candidate_regions(
    hits: Iterable[BlastHit], min_width: int = 100, flank: int = 1000
) -> tuple[list[CandidateRegion], dict[str, int]]:
    """Reduce overlapping/adjacent hits by scaffold and strand, matching GenomicRanges."""

    if min_width < 1:
        raise ValueError("min_width must be at least 1")
    if flank < 0:
        raise ValueError("flank must be non-negative")

    sorted_hits = sorted(
        hits,
        key=lambda hit: (
            hit.sseqid,
            hit.strand,
            hit.genome_start_0,
            hit.genome_end_0,
            hit.qseqid,
        ),
    )
    candidates: list[CandidateRegion] = []
    reduced_count = 0
    filtered_short = 0
    current: _RegionAccumulator | None = None

    def keep_or_count(region: _RegionAccumulator) -> None:
        nonlocal reduced_count, filtered_short
        reduced_count += 1
        if region.end_0 - region.start_0 < min_width:
            filtered_short += 1
        else:
            candidates.append(_finalize_region(region, flank))

    for hit in sorted_hits:
        if current is None:
            current = _new_accumulator(hit)
            continue

        same_group = hit.sseqid == current.sseqid and hit.strand == current.strand
        overlaps_or_adjacent = hit.genome_start_0 <= current.end_0
        if same_group and overlaps_or_adjacent:
            current.add(hit)
        else:
            keep_or_count(current)
            current = _new_accumulator(hit)

    if current is not None:
        keep_or_count(current)

    statistics = {
        "input_hits": len(sorted_hits),
        "reduced_regions": reduced_count,
        "filtered_short_regions": filtered_short,
        "candidate_regions": len(candidates),
        "left_boundary_clipped": sum(region.expanded_start_0 == 0 for region in candidates),
        "right_boundary_clipped": sum(
            region.expanded_end_0 == region.slen for region in candidates
        ),
        "candidate_expanded_bp": sum(region.expanded_width for region in candidates),
    }
    return candidates, statistics


def read_blast_hits(path: Path) -> list[BlastHit]:
    """Read a previously validated table into compact slotted hit objects."""

    hits: list[BlastHit] = []
    with _open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            try:
                hits.append(parse_blast_row(line))
            except ValueError as exc:
                raise ValueError(f"{path}: line {line_number}: {exc}") from exc
    return hits


def write_candidate_region_outputs(
    candidates: list[CandidateRegion],
    statistics: dict[str, int],
    output_prefix: Path,
    species: str,
    min_width: int,
    flank: int,
) -> dict[str, Path]:
    """Write explicit TSV/BED/faidx contracts and a JSON audit summary."""

    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    tsv_path = Path(f"{output_prefix}.candidate_regions.tsv")
    bed_path = Path(f"{output_prefix}.candidate_regions.bed")
    faidx_path = Path(f"{output_prefix}.candidate_regions.faidx.txt")
    json_path = Path(f"{output_prefix}.candidate_regions.summary.json")
    header = (
        "region_id\tsseqid\tstrand\thit_start_1\thit_end_1\t"
        "expanded_start_1\texpanded_end_1\thit_width_bp\texpanded_width_bp\t"
        "hit_count\tquery_count\tbest_evalue\tmax_bitscore\tslen\tfaidx_region\n"
    )

    with (
        tsv_path.open("wt", encoding="utf-8", newline="\n") as tsv,
        bed_path.open("wt", encoding="utf-8", newline="\n") as bed,
        faidx_path.open("wt", encoding="utf-8", newline="\n") as faidx,
    ):
        tsv.write(header)
        for index, region in enumerate(candidates, start=1):
            region_id = f"{species}_ORC{index:06d}"
            faidx_region = (
                f"{region.sseqid}:{region.expanded_start_0 + 1}-{region.expanded_end_0}"
            )
            tsv.write(
                "\t".join(
                    (
                        region_id,
                        region.sseqid,
                        region.strand,
                        str(region.hit_start_0 + 1),
                        str(region.hit_end_0),
                        str(region.expanded_start_0 + 1),
                        str(region.expanded_end_0),
                        str(region.hit_width),
                        str(region.expanded_width),
                        str(region.hit_count),
                        str(region.query_count),
                        f"{region.best_evalue:.6g}",
                        f"{region.max_bitscore:.6g}",
                        str(region.slen),
                        faidx_region,
                    )
                )
                + "\n"
            )
            bed.write(
                f"{region.sseqid}\t{region.expanded_start_0}\t{region.expanded_end_0}"
                f"\t{region_id}\t0\t{region.strand}\n"
            )
            faidx.write(f"{faidx_region}\n")

    summary: dict[str, Any] = {
        "schema_version": "1.0",
        "status": "pass",
        "species": species,
        "method": "strand-aware GenomicRanges-compatible reduce",
        "coordinates": {
            "tsv": "1-based inclusive",
            "bed": "0-based half-open",
            "faidx": "1-based inclusive",
        },
        "parameters": {"min_hit_region_width_bp": min_width, "flank_bp": flank},
        "counts": statistics,
    }
    with json_path.open("wt", encoding="utf-8", newline="\n") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")

    return {"tsv": tsv_path, "bed": bed_path, "faidx": faidx_path, "summary": json_path}
