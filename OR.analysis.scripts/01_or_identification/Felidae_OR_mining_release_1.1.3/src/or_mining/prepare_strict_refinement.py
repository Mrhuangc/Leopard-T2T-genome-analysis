"""CLI for preparing strict local Exonerate refinement pairs."""

from __future__ import annotations

import argparse
from pathlib import Path

from or_mining.strict_refinement import prepare_strict_pairs


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference_or", type=Path)
    parser.add_argument("species_or", type=Path)
    parser.add_argument("region_fasta", type=Path)
    parser.add_argument("regions_tsv", type=Path)
    parser.add_argument("candidates_tsv", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    args = parser.parse_args()
    prepare_strict_pairs(
        args.reference_or,
        args.species_or,
        args.region_fasta,
        args.regions_tsv,
        args.candidates_tsv,
        args.output_dir,
        args.manifest,
        args.summary,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
