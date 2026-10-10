"""Strict parsing and genome-coordinate normalization of EMBOSS getorf output."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
import csv
import json
import re

from or_mining.fasta import read_fasta, write_fasta_record


GETORF_HEADER = re.compile(
    r"^(?P<region_id>.+)_(?P<ordinal>[0-9]+)\s+"
    r"\[(?P<start>[0-9]+)\s*-\s*(?P<end>[0-9]+)\]"
    r"(?P<reverse>\s+\(REVERSE SENSE\))?(?:\s+.*)?$"
)


@dataclass(frozen=True, slots=True)
class ParsedGetorfHeader:
    region_id: str
    ordinal: int
    local_start_1: int
    local_end_1: int
    strand: str


@dataclass(slots=True)
class NormalizedOrf:
    sseqid: str
    start_1: int
    end_1: int
    strand: str
    sequence: str
    source_region_ids: list[str] = field(default_factory=list)
    source_getorf_ids: list[str] = field(default_factory=list)
    local_coordinates: list[str] = field(default_factory=list)
    candidate_hit_strands: set[str] = field(default_factory=set)


def parse_getorf_header(header: str) -> ParsedGetorfHeader:
    """Parse the documented EMBOSS 6.6 getorf FASTA description."""

    match = GETORF_HEADER.fullmatch(header)
    if match is None:
        raise ValueError(f"unrecognized EMBOSS getorf header: {header!r}")
    local_start = int(match.group("start"))
    local_end = int(match.group("end"))
    reverse = match.group("reverse") is not None
    if reverse and local_start <= local_end:
        raise ValueError(f"reverse getorf coordinates must descend: {header!r}")
    if not reverse and local_start > local_end:
        raise ValueError(f"forward getorf coordinates must ascend: {header!r}")
    return ParsedGetorfHeader(
        region_id=match.group("region_id"),
        ordinal=int(match.group("ordinal")),
        local_start_1=local_start,
        local_end_1=local_end,
        strand="-" if reverse else "+",
    )


def normalize_getorf(
    raw_orfs: Path,
    regions_tsv: Path,
    output_fasta: Path,
    output_tsv: Path,
    summary_path: Path,
    species: str,
) -> dict[str, Any]:
    """Map local getorf coordinates to the genome and remove duplicate loci."""

    with regions_tsv.open("rt", encoding="utf-8", newline="") as handle:
        region_rows = list(csv.DictReader(handle, delimiter="\t"))
    regions = {row["region_id"]: row for row in region_rows}
    if len(regions) != len(region_rows):
        raise ValueError("candidate region IDs are not unique")

    unique: dict[tuple[str, int, int, str, str], NormalizedOrf] = {}
    raw_count = 0
    for header, raw_sequence in read_fasta(raw_orfs):
        raw_count += 1
        parsed = parse_getorf_header(header)
        region = regions.get(parsed.region_id)
        if region is None:
            raise ValueError(f"getorf references unknown region {parsed.region_id!r}")

        region_length = int(region["expanded_width_bp"])
        local_left = min(parsed.local_start_1, parsed.local_end_1)
        local_right = max(parsed.local_start_1, parsed.local_end_1)
        if not (1 <= local_left <= local_right <= region_length):
            raise ValueError(f"getorf coordinates outside {parsed.region_id}: {header!r}")

        sequence = raw_sequence.upper()
        coordinate_length = local_right - local_left + 1
        if len(sequence) != coordinate_length:
            raise ValueError(
                f"{header!r}: coordinate length {coordinate_length} != sequence length {len(sequence)}"
            )
        if len(sequence) % 3 != 0:
            raise ValueError(f"{header!r}: ORF nucleotide length is not divisible by three")

        expanded_start_1 = int(region["expanded_start_1"])
        genome_start_1 = expanded_start_1 + local_left - 1
        genome_end_1 = expanded_start_1 + local_right - 1
        getorf_id = f"{parsed.region_id}_{parsed.ordinal}"
        key = (region["sseqid"], genome_start_1, genome_end_1, parsed.strand, sequence)
        normalized = unique.get(key)
        if normalized is None:
            normalized = NormalizedOrf(
                sseqid=region["sseqid"],
                start_1=genome_start_1,
                end_1=genome_end_1,
                strand=parsed.strand,
                sequence=sequence,
            )
            unique[key] = normalized
        normalized.source_region_ids.append(parsed.region_id)
        normalized.source_getorf_ids.append(getorf_id)
        normalized.local_coordinates.append(
            f"{parsed.local_start_1}-{parsed.local_end_1}"
        )
        normalized.candidate_hit_strands.add(region["strand"])

    ordered = sorted(
        unique.values(), key=lambda orf: (orf.sseqid, orf.start_1, orf.end_1, orf.strand, orf.sequence)
    )
    header = (
        "orf_id\tsseqid\tstart_1\tend_1\tstrand\tlength_nt\t"
        "source_region_ids\tsource_getorf_ids\tlocal_coordinates\t"
        "candidate_hit_strands\tduplicate_source_count\n"
    )
    with (
        output_fasta.open("wt", encoding="utf-8", newline="\n") as fasta,
        output_tsv.open("wt", encoding="utf-8", newline="\n") as table,
    ):
        table.write(header)
        for index, orf in enumerate(ordered, start=1):
            orf_id = f"{species}_ORF{index:06d}"
            source_regions = ";".join(orf.source_region_ids)
            source_getorf = ";".join(orf.source_getorf_ids)
            local_coordinates = ";".join(orf.local_coordinates)
            hit_strands = ";".join(sorted(orf.candidate_hit_strands))
            table.write(
                f"{orf_id}\t{orf.sseqid}\t{orf.start_1}\t{orf.end_1}\t{orf.strand}\t"
                f"{len(orf.sequence)}\t{source_regions}\t{source_getorf}\t"
                f"{local_coordinates}\t{hit_strands}\t{len(orf.source_getorf_ids)}\n"
            )
            fasta_header = (
                f"{orf_id} {orf.sseqid}:{orf.start_1}-{orf.end_1} "
                f"strand={orf.strand} source_regions={source_regions}"
            )
            write_fasta_record(fasta, fasta_header, orf.sequence)

    summary: dict[str, Any] = {
        "schema_version": "1.0",
        "status": "pass",
        "species": species,
        "coordinates": "genome 1-based inclusive",
        "counts": {
            "raw_getorf_records": raw_count,
            "unique_orfs": len(ordered),
            "duplicate_records_removed": raw_count - len(ordered),
            "plus_strand_orfs": sum(orf.strand == "+" for orf in ordered),
            "minus_strand_orfs": sum(orf.strand == "-" for orf in ordered),
            "orfs_with_ambiguous_bases": sum("N" in orf.sequence for orf in ordered),
        },
    }
    with summary_path.open("wt", encoding="utf-8", newline="\n") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
    return summary

