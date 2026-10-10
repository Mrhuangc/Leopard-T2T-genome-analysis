"""Command-line adapter for the species-specific TBLASTN merge."""

from __future__ import annotations

from pathlib import Path
import argparse

from or_mining.expansion_hits import merge_species_expansion_hits


def main() -> int:
    parser = argparse.ArgumentParser(description="Filter and merge species OR TBLASTN hits")
    parser.add_argument("initial_hits", type=Path)
    parser.add_argument("species_hits", type=Path)
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument("--max-hits-per-query", type=int, default=10_000)
    args = parser.parse_args()
    merge_species_expansion_hits(
        args.initial_hits,
        args.species_hits,
        args.output_prefix,
        max_hits_per_query=args.max_hits_per_query,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
