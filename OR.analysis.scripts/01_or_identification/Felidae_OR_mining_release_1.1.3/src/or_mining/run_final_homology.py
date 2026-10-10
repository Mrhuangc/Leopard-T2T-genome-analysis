"""CLI for the NC-compatible final OR homology filter."""

from __future__ import annotations

import argparse
from pathlib import Path

from or_mining.final_homology import build_nc_final_catalog


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("catalog", type=Path)
    parser.add_argument("gff3", type=Path)
    parser.add_argument("functional_nucleotide", type=Path)
    parser.add_argument("functional_protein", type=Path)
    parser.add_argument("pseudogene_nucleotide", type=Path)
    parser.add_argument("truncated_nucleotide", type=Path)
    parser.add_argument("ambiguous_nucleotide", type=Path)
    parser.add_argument("functional_blastp", type=Path)
    parser.add_argument("pseudogene_blastx", type=Path)
    parser.add_argument("truncated_blastx", type=Path)
    parser.add_argument("ambiguous_blastx", type=Path)
    parser.add_argument("--species", required=True)
    parser.add_argument("--blast-evalue", required=True)
    parser.add_argument("--output-prefix", required=True, type=Path)
    args = parser.parse_args()
    build_nc_final_catalog(
        args.species,
        args.catalog,
        args.gff3,
        args.functional_nucleotide,
        args.functional_protein,
        args.pseudogene_nucleotide,
        args.truncated_nucleotide,
        args.ambiguous_nucleotide,
        args.functional_blastp,
        args.pseudogene_blastx,
        args.truncated_blastx,
        args.ambiguous_blastx,
        args.output_prefix,
        args.blast_evalue,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
