"""Standalone candidate-region builder used by Nextflow."""

from __future__ import annotations

import argparse
from pathlib import Path

from or_mining.regions import (
    build_candidate_regions,
    read_blast_hits,
    write_candidate_region_outputs,
)


def main() -> int:
    parser = argparse.ArgumentParser(description="Build candidate regions from merged TBLASTN hits")
    parser.add_argument("input", type=Path)
    parser.add_argument("--species", required=True)
    parser.add_argument("--output-prefix", required=True, type=Path)
    parser.add_argument("--min-width", type=int, default=100)
    parser.add_argument("--flank", type=int, default=1000)
    args = parser.parse_args()

    hits = read_blast_hits(args.input)
    candidates, statistics = build_candidate_regions(hits, args.min_width, args.flank)
    write_candidate_region_outputs(
        candidates,
        statistics,
        args.output_prefix,
        args.species,
        args.min_width,
        args.flank,
    )
    print(
        f"status=pass hits={statistics['input_hits']} "
        f"candidate_regions={statistics['candidate_regions']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

