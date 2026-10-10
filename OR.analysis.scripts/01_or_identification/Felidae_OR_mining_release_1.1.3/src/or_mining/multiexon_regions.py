"""Build NC-compatible candidate regions for multi-exon OR discovery."""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

from or_mining.blast import BlastHit
from or_mining.regions import read_blast_hits


@dataclass(frozen=True, slots=True)
class KnownLocus:
    sseqid: str
    start_0: int
    end_0: int


@dataclass(slots=True)
class _ReducedHit:
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


@dataclass(frozen=True, slots=True)
class MultiExonRegion:
    sseqid: str
    start_0: int
    end_0: int
    slen: int
    source_interval_count: int
    hit_count: int
    query_count: int
    best_evalue: float
    max_bitscore: float

    @property
    def width(self) -> int:
        return self.end_0 - self.start_0


def read_known_loci(path: Path) -> list[KnownLocus]:
    """Read the tree-filtered coordinate table (1-based inclusive)."""

    loci: list[KnownLocus] = []
    with path.open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sseqid", "start_1", "end_1"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path}: required columns are {sorted(required)}")
        for line_number, row in enumerate(reader, start=2):
            try:
                start_1 = int(row["start_1"])
                end_1 = int(row["end_1"])
            except (TypeError, ValueError) as exc:
                raise ValueError(f"{path}: line {line_number}: invalid coordinates") from exc
            if start_1 < 1 or end_1 < start_1:
                raise ValueError(f"{path}: line {line_number}: invalid 1-based interval")
            loci.append(KnownLocus(row["sseqid"], start_1 - 1, end_1))
    return loci


def _reduce_hits(hits: Iterable[BlastHit]) -> list[_ReducedHit]:
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
    reduced: list[_ReducedHit] = []
    current: _ReducedHit | None = None
    for hit in sorted_hits:
        if current is None or not (
            hit.sseqid == current.sseqid
            and hit.strand == current.strand
            and hit.genome_start_0 <= current.end_0
        ):
            if current is not None:
                reduced.append(current)
            current = _ReducedHit(
                hit.sseqid,
                hit.strand,
                hit.genome_start_0,
                hit.genome_end_0,
                hit.slen,
            )
        elif hit.slen != current.slen:
            raise ValueError(f"inconsistent slen for {hit.sseqid!r}")
        current.add(hit)
    if current is not None:
        reduced.append(current)
    return reduced


def _overlaps_locus(region: _ReducedHit, loci: list[KnownLocus]) -> bool:
    return any(
        locus.sseqid == region.sseqid
        and locus.start_0 < region.end_0
        and region.start_0 < locus.end_0
        for locus in loci
    )


def build_multiexon_regions(
    hits: Iterable[BlastHit],
    known_loci: list[KnownLocus],
    maximum_intron_length: int = 10_000,
    gene_clearance: int = 50,
) -> tuple[list[MultiExonRegion], dict[str, int]]:
    """Reproduce the NC R merge logic with explicit coordinate safeguards."""

    if maximum_intron_length < 1:
        raise ValueError("maximum_intron_length must be at least 1")
    if gene_clearance < 0:
        raise ValueError("gene_clearance must be non-negative")

    hit_list = list(hits)
    initial = _reduce_hits(hit_list)
    remaining = [region for region in initial if not _overlaps_locus(region, known_loci)]
    loci_by_subject: dict[str, list[KnownLocus]] = {}
    for locus in known_loci:
        loci_by_subject.setdefault(locus.sseqid, []).append(locus)
    for loci in loci_by_subject.values():
        loci.sort(key=lambda locus: (locus.start_0, locus.end_0))

    left_gene_clips = 0
    right_gene_clips = 0
    left_sequence_clips = 0
    right_sequence_clips = 0
    expanded: list[tuple[str, int, int, int, _ReducedHit]] = []
    for region in remaining:
        desired_start = region.start_0 - maximum_intron_length
        desired_end = region.end_0 + maximum_intron_length
        start_0 = max(0, desired_start)
        end_0 = min(region.slen, desired_end)
        left_sequence_clips += desired_start < 0
        right_sequence_clips += desired_end > region.slen

        left_loci = [
            locus
            for locus in loci_by_subject.get(region.sseqid, [])
            if locus.end_0 <= region.start_0 and locus.end_0 >= desired_start + 1
        ]
        if left_loci:
            nearest = max(left_loci, key=lambda locus: locus.end_0)
            # Matches the published R expression: known_end_1 + 50.
            start_0 = max(start_0, nearest.end_0 + gene_clearance - 1)
            left_gene_clips += 1

        right_loci = [
            locus
            for locus in loci_by_subject.get(region.sseqid, [])
            if locus.start_0 >= region.end_0 and locus.start_0 < desired_end
        ]
        if right_loci:
            nearest = min(right_loci, key=lambda locus: locus.start_0)
            # Matches the published R expression: known_start_1 - 50.
            end_0 = min(end_0, nearest.start_0 - gene_clearance + 1)
            right_gene_clips += 1

        if start_0 < end_0:
            expanded.append((region.sseqid, start_0, end_0, region.slen, region))

    expanded.sort(key=lambda item: (item[0], item[1], item[2]))
    final: list[MultiExonRegion] = []
    group: list[tuple[str, int, int, int, _ReducedHit]] = []

    def finish_group(items: list[tuple[str, int, int, int, _ReducedHit]]) -> None:
        if not items:
            return
        query_ids: set[str] = set()
        for item in items:
            query_ids.update(item[4].query_ids)
        final.append(
            MultiExonRegion(
                sseqid=items[0][0],
                start_0=min(item[1] for item in items),
                end_0=max(item[2] for item in items),
                slen=items[0][3],
                source_interval_count=len(items),
                hit_count=sum(item[4].hit_count for item in items),
                query_count=len(query_ids),
                best_evalue=min(item[4].best_evalue for item in items),
                max_bitscore=max(item[4].max_bitscore for item in items),
            )
        )

    for item in expanded:
        if group and (item[0] != group[0][0] or item[1] > max(x[2] for x in group)):
            finish_group(group)
            group = []
        group.append(item)
    finish_group(group)

    counts = {
        "input_hits": len(hit_list),
        "initial_strand_aware_regions": len(initial),
        "regions_overlapping_confirmed_or": len(initial) - len(remaining),
        "regions_after_confirmed_or_removal": len(remaining),
        "left_clipped_by_confirmed_or": left_gene_clips,
        "right_clipped_by_confirmed_or": right_gene_clips,
        "left_clipped_by_sequence_boundary": left_sequence_clips,
        "right_clipped_by_sequence_boundary": right_sequence_clips,
        "final_multiexon_candidate_regions": len(final),
        "final_candidate_bp": sum(region.width for region in final),
    }
    return final, counts


def write_multiexon_outputs(
    regions: list[MultiExonRegion],
    counts: dict[str, int],
    output_prefix: Path,
    species: str,
    maximum_intron_length: int,
    gene_clearance: int,
) -> dict[str, Path]:
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    tsv_path = Path(f"{output_prefix}.multiexon_regions.tsv")
    bed_path = Path(f"{output_prefix}.multiexon_regions.bed")
    faidx_path = Path(f"{output_prefix}.multiexon_regions.faidx.txt")
    summary_path = Path(f"{output_prefix}.multiexon_regions.summary.json")
    header = (
        "region_id\tsseqid\tstrand\thit_start_1\thit_end_1\t"
        "expanded_start_1\texpanded_end_1\thit_width_bp\texpanded_width_bp\t"
        "hit_count\tquery_count\tbest_evalue\tmax_bitscore\tslen\tfaidx_region\t"
        "source_interval_count\n"
    )
    with (
        tsv_path.open("wt", encoding="utf-8", newline="\n") as tsv,
        bed_path.open("wt", encoding="utf-8", newline="\n") as bed,
        faidx_path.open("wt", encoding="utf-8", newline="\n") as faidx,
    ):
        tsv.write(header)
        for index, region in enumerate(regions, start=1):
            region_id = f"{species}_ORME{index:06d}"
            faidx_region = f"{region.sseqid}:{region.start_0 + 1}-{region.end_0}"
            values = (
                region_id,
                region.sseqid,
                ".",
                str(region.start_0 + 1),
                str(region.end_0),
                str(region.start_0 + 1),
                str(region.end_0),
                str(region.width),
                str(region.width),
                str(region.hit_count),
                str(region.query_count),
                f"{region.best_evalue:.6g}",
                f"{region.max_bitscore:.6g}",
                str(region.slen),
                faidx_region,
                str(region.source_interval_count),
            )
            tsv.write("\t".join(values) + "\n")
            bed.write(f"{region.sseqid}\t{region.start_0}\t{region.end_0}\t{region_id}\t0\t.\n")
            faidx.write(f"{faidx_region}\n")

    summary: dict[str, Any] = {
        "schema_version": "1.0",
        "status": "pass",
        "species": species,
        "method": "NC R_sript_merge_multiexon.R-compatible two-stage reduction",
        "coordinates": {
            "input_known_or": "1-based inclusive",
            "tsv": "1-based inclusive",
            "bed": "0-based half-open",
            "faidx": "1-based inclusive",
        },
        "parameters": {
            "maximum_intron_length_bp": maximum_intron_length,
            "confirmed_or_clearance_bp": gene_clearance,
        },
        "counts": counts,
    }
    with summary_path.open("wt", encoding="utf-8", newline="\n") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
    return {"tsv": tsv_path, "bed": bed_path, "faidx": faidx_path, "summary": summary_path}


def run_multiexon_region_builder(
    merged_hits: Path,
    known_coordinates: Path,
    output_prefix: Path,
    species: str,
    maximum_intron_length: int,
    gene_clearance: int,
) -> dict[str, Path]:
    hits = read_blast_hits(merged_hits)
    known_loci = read_known_loci(known_coordinates)
    regions, counts = build_multiexon_regions(
        hits, known_loci, maximum_intron_length, gene_clearance
    )
    return write_multiexon_outputs(
        regions, counts, output_prefix, species, maximum_intron_length, gene_clearance
    )
