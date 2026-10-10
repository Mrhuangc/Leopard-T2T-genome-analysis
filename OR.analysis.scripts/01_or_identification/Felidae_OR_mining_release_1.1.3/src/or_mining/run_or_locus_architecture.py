"""CLI for single-species OR genomic-architecture summaries."""

from __future__ import annotations

from argparse import ArgumentParser
from pathlib import Path

from or_mining.or_locus_architecture import analyze_or_locus_architecture


def main() -> int:
    parser = ArgumentParser()
    parser.add_argument("catalog", type=Path)
    parser.add_argument("genome_fai", type=Path)
    parser.add_argument("--species", required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument("--cluster-gap-bp", type=int, default=1_000_000)
    args = parser.parse_args()
    analyze_or_locus_architecture(
        args.species, args.catalog, args.genome_fai, args.output_prefix, args.cluster_gap_bp
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
