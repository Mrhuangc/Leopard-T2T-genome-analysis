"""Standalone CLI for strict getorf coordinate normalization."""

from __future__ import annotations

import argparse
from pathlib import Path

from or_mining.orfs import normalize_getorf


def main() -> int:
    parser = argparse.ArgumentParser(description="Normalize EMBOSS getorf coordinates")
    parser.add_argument("raw_orfs", type=Path)
    parser.add_argument("regions_tsv", type=Path)
    parser.add_argument("--species", required=True)
    parser.add_argument("--output-fasta", required=True, type=Path)
    parser.add_argument("--output-tsv", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    args = parser.parse_args()
    summary = normalize_getorf(
        args.raw_orfs,
        args.regions_tsv,
        args.output_fasta,
        args.output_tsv,
        args.summary,
        args.species,
    )
    counts = summary["counts"]
    print(
        f"status=pass raw_orfs={counts['raw_getorf_records']} "
        f"unique_orfs={counts['unique_orfs']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

