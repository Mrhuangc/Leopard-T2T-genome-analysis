"""Command-line adapters used by the initial TBLASTN Nextflow modules."""

from __future__ import annotations

import argparse
from pathlib import Path

from or_mining.initial_search import merge_initial_tblastn, split_initial_queries


def main() -> int:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    split = subparsers.add_parser("split")
    split.add_argument("query_fasta", type=Path)
    split.add_argument("--output-dir", type=Path, required=True)
    split.add_argument("--output-prefix", required=True)
    split.add_argument("--shards", type=int, required=True)

    merge = subparsers.add_parser("merge")
    merge.add_argument("parts", nargs="+", type=Path)
    merge.add_argument("--manifest", type=Path, required=True)
    merge.add_argument("--output-merged", type=Path, required=True)
    merge.add_argument("--validation", type=Path, required=True)
    merge.add_argument("--summary", type=Path, required=True)
    merge.add_argument("--species", required=True)
    merge.add_argument("--evalue", required=True)

    args = parser.parse_args()
    if args.command == "split":
        split_initial_queries(
            args.query_fasta,
            args.output_dir,
            args.output_prefix,
            args.shards,
        )
    else:
        merge_initial_tblastn(
            args.parts,
            args.manifest,
            args.output_merged,
            args.validation,
            args.summary,
            species=args.species,
            evalue=args.evalue,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
