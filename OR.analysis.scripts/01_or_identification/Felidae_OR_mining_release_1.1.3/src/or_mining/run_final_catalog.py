"""CLI for finalizing the OR catalog with InterProScan IPR017452 evidence."""

from __future__ import annotations

import argparse
from pathlib import Path

from or_mining.final_catalog import build_final_catalog


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pre7tm_catalog", type=Path)
    parser.add_argument("pre7tm_gff3", type=Path)
    parser.add_argument("functional_nucleotide", type=Path)
    parser.add_argument("functional_protein", type=Path)
    parser.add_argument("pseudogene_nucleotide", type=Path)
    parser.add_argument("truncated_nucleotide", type=Path)
    parser.add_argument("ambiguous_nucleotide", type=Path)
    parser.add_argument("interpro_tsv", type=Path)
    parser.add_argument("version_report", type=Path)
    parser.add_argument("--species", required=True)
    parser.add_argument("--expected-version", required=True)
    parser.add_argument("--expected-ipr", default="IPR017452")
    parser.add_argument("--output-prefix", required=True, type=Path)
    args = parser.parse_args()
    build_final_catalog(
        args.species,
        args.pre7tm_catalog,
        args.pre7tm_gff3,
        args.functional_nucleotide,
        args.functional_protein,
        args.pseudogene_nucleotide,
        args.truncated_nucleotide,
        args.ambiguous_nucleotide,
        args.interpro_tsv,
        args.version_report,
        args.expected_version,
        args.output_prefix,
        args.expected_ipr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
