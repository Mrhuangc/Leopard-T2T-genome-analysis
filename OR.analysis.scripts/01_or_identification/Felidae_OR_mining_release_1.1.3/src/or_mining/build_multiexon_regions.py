"""CLI for the NC-compatible multi-exon candidate-region builder."""

from __future__ import annotations

import argparse
from pathlib import Path

from or_mining.multiexon_regions import run_multiexon_region_builder


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("merged_hits", type=Path)
    parser.add_argument("known_coordinates", type=Path)
    parser.add_argument("--output-prefix", required=True, type=Path)
    parser.add_argument("--species", required=True)
    parser.add_argument("--maximum-intron-length", type=int, default=10_000)
    parser.add_argument("--gene-clearance", type=int, default=50)
    args = parser.parse_args()
    outputs = run_multiexon_region_builder(
        args.merged_hits,
        args.known_coordinates,
        args.output_prefix,
        args.species,
        args.maximum_intron_length,
        args.gene_clearance,
    )
    print(f"status=pass summary={outputs['summary']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
