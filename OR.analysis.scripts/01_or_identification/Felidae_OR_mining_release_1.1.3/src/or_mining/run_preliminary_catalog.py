"""CLI for pre-7TM OR catalog consolidation."""

from __future__ import annotations

import argparse
from pathlib import Path

from or_mining.preliminary_catalog import build_preliminary_catalog


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("genome_fai", type=Path)
    parser.add_argument("single_nucleotide", type=Path)
    parser.add_argument("single_protein", type=Path)
    parser.add_argument("single_coordinates", type=Path)
    parser.add_argument("multiexon_calls", type=Path)
    parser.add_argument("multiexon_gff3", type=Path)
    parser.add_argument("multiexon_functional_nucleotide", type=Path)
    parser.add_argument("multiexon_functional_protein", type=Path)
    parser.add_argument("multiexon_pseudogene_nucleotide", type=Path)
    parser.add_argument("multiexon_truncated_nucleotide", type=Path)
    parser.add_argument("multiexon_ambiguous_nucleotide", type=Path)
    parser.add_argument("--species", required=True)
    parser.add_argument("--output-prefix", required=True, type=Path)
    args = parser.parse_args()
    build_preliminary_catalog(
        args.species, args.genome_fai, args.single_nucleotide, args.single_protein,
        args.single_coordinates, args.multiexon_calls, args.multiexon_gff3,
        args.multiexon_functional_nucleotide, args.multiexon_functional_protein,
        args.multiexon_pseudogene_nucleotide, args.multiexon_truncated_nucleotide,
        args.multiexon_ambiguous_nucleotide, args.output_prefix,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
