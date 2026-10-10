"""Rename samtools-faidx output to stable candidate-region identifiers."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from or_mining.fasta import read_fasta, write_fasta_record


def rename_extracted_regions(
    raw_fasta: Path, regions_tsv: Path, output_fasta: Path, summary_path: Path
) -> dict[str, int | str]:
    with regions_tsv.open("rt", encoding="utf-8", newline="") as handle:
        regions = list(csv.DictReader(handle, delimiter="\t"))
    records = list(read_fasta(raw_fasta))
    if len(records) != len(regions):
        raise ValueError(
            f"samtools returned {len(records)} records for {len(regions)} candidate regions"
        )

    total_bases = 0
    with output_fasta.open("wt", encoding="utf-8", newline="\n") as output:
        for region, (raw_header, sequence) in zip(regions, records, strict=True):
            expected_header = region["faidx_region"]
            if raw_header.split()[0] != expected_header:
                raise ValueError(
                    f"faidx order/header mismatch: expected {expected_header!r}, got {raw_header!r}"
                )
            expected_length = int(region["expanded_width_bp"])
            if len(sequence) != expected_length:
                raise ValueError(
                    f"{region['region_id']}: expected {expected_length} bp, got {len(sequence)}"
                )
            sequence = sequence.upper()
            total_bases += len(sequence)
            write_fasta_record(output, region["region_id"], sequence)

    summary: dict[str, int | str] = {
        "schema_version": "1.0",
        "status": "pass",
        "candidate_regions": len(regions),
        "extracted_sequences": len(records),
        "extracted_bases": total_bases,
    }
    with summary_path.open("wt", encoding="utf-8", newline="\n") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description="Rename and verify samtools candidate FASTA")
    parser.add_argument("raw_fasta", type=Path)
    parser.add_argument("regions_tsv", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    args = parser.parse_args()
    summary = rename_extracted_regions(
        args.raw_fasta, args.regions_tsv, args.output, args.summary
    )
    print(
        f"status=pass sequences={summary['extracted_sequences']} "
        f"bases={summary['extracted_bases']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

