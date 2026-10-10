"""CLI for TBLASTN-seeded targeted Exonerate pair preparation."""

from __future__ import annotations

import argparse
from pathlib import Path

from or_mining.targeted_exonerate import prepare_targeted_pairs


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference_or", type=Path)
    parser.add_argument("species_or", type=Path)
    parser.add_argument("target_fasta", type=Path)
    parser.add_argument("regions_tsv", type=Path)
    parser.add_argument("merged_hits", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--top-queries-per-region", type=int, default=3)
    args = parser.parse_args()
    prepare_targeted_pairs(
        args.reference_or,
        args.species_or,
        args.target_fasta,
        args.regions_tsv,
        args.merged_hits,
        args.output_dir,
        args.manifest,
        args.summary,
        args.top_queries_per_region,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
