"""Standalone entry point for translation normalization."""

from __future__ import annotations

import argparse
from pathlib import Path

from or_mining.translations import normalize_transeq


def main() -> int:
    parser = argparse.ArgumentParser(description="Normalize EMBOSS transeq output")
    parser.add_argument("nucleotide_fasta", type=Path)
    parser.add_argument("raw_protein_fasta", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    args = parser.parse_args()
    summary = normalize_transeq(
        args.nucleotide_fasta, args.raw_protein_fasta, args.output, args.summary
    )
    counts = summary["counts"]
    print(f"status=pass proteins={counts['protein_orfs']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

