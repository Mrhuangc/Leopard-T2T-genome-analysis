"""Command-line adapter for species-specific OR query preparation."""

from __future__ import annotations

from pathlib import Path
import argparse

from or_mining.species_queries import prepare_species_queries


def main() -> int:
    parser = argparse.ArgumentParser(description="Prepare species-specific OR TBLASTN queries")
    parser.add_argument("protein_fasta", type=Path)
    parser.add_argument("nucleotide_fasta", type=Path)
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument("--cluster-threshold", type=int, default=500)
    parser.add_argument("--cluster-identity", type=float, default=0.90)
    parser.add_argument("--cpus", type=int, default=1)
    args = parser.parse_args()
    prepare_species_queries(
        args.protein_fasta,
        args.nucleotide_fasta,
        args.output_prefix,
        cluster_threshold=args.cluster_threshold,
        cluster_identity=args.cluster_identity,
        cpus=args.cpus,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
