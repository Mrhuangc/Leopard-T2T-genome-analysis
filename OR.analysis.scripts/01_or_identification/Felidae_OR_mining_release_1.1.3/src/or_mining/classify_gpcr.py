"""Standalone GPCR best-hit classifier."""

from __future__ import annotations

import argparse
from pathlib import Path

from or_mining.classification import classify_best_gpcr_hits


def main() -> int:
    parser = argparse.ArgumentParser(description="Classify ORFs by best GPCR BLASTP hit")
    parser.add_argument("blastp_tsv", type=Path)
    parser.add_argument("protein_fasta", type=Path)
    parser.add_argument("nucleotide_fasta", type=Path)
    parser.add_argument("orf_table", type=Path)
    parser.add_argument("--output-prefix", required=True, type=Path)
    args = parser.parse_args()
    summary = classify_best_gpcr_hits(
        args.blastp_tsv,
        args.protein_fasta,
        args.nucleotide_fasta,
        args.orf_table,
        args.output_prefix,
    )
    counts = summary["counts"]
    print(f"status=pass retained_or={counts['retained_or']} input_orfs={counts['input_orfs']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

