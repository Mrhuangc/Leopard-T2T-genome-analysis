"""CLI adapters for NC-compatible single-exon nonfunctional OR recovery."""

from __future__ import annotations

import argparse
from pathlib import Path

from or_mining.single_exon_recovery import (
    build_complete_catalog,
    build_single_exon_candidates,
    classify_recovered_candidates,
    normalize_extracted_targets,
    run_targeted_tblastn,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare = subparsers.add_parser("prepare")
    prepare.add_argument("merged_hits", type=Path)
    prepare.add_argument("catalog", type=Path)
    prepare.add_argument("reference_or", type=Path)
    prepare.add_argument("species_or", type=Path)
    prepare.add_argument("--species", required=True)
    prepare.add_argument("--output-prefix", required=True, type=Path)
    prepare.add_argument("--clearance", type=int, default=100)
    prepare.add_argument("--flank", type=int, default=100)

    extract = subparsers.add_parser("normalize-targets")
    extract.add_argument("raw_fasta", type=Path)
    extract.add_argument("candidates", type=Path)
    extract.add_argument("--output", required=True, type=Path)
    extract.add_argument("--summary", required=True, type=Path)

    search = subparsers.add_parser("search")
    search.add_argument("candidates", type=Path)
    search.add_argument("queries", type=Path)
    search.add_argument("targets", type=Path)
    search.add_argument("--standard-output", required=True, type=Path)
    search.add_argument("--gapless-output", required=True, type=Path)
    search.add_argument("--summary", required=True, type=Path)
    search.add_argument("--gate", type=Path)
    search.add_argument("--evalue", required=True)
    search.add_argument("--workers", required=True, type=int)
    search.add_argument("--timeout-seconds", required=True, type=int)
    search.add_argument("--limit", type=int)

    classify = subparsers.add_parser("classify")
    classify.add_argument("candidates", type=Path)
    classify.add_argument("targets", type=Path)
    classify.add_argument("standard_tsv", type=Path)
    classify.add_argument("gapless_tsv", type=Path)
    classify.add_argument("--output-prefix", required=True, type=Path)

    complete = subparsers.add_parser("complete")
    for name in (
        "existing_catalog", "existing_gff3", "existing_functional_nucleotide",
        "existing_functional_protein", "existing_pseudogene", "existing_truncated",
        "existing_ambiguous", "recovery_calls", "recovered_pseudogene",
        "recovered_truncated", "recovered_ambiguous", "recovery_blastx",
    ):
        complete.add_argument(name, type=Path)
    complete.add_argument("--species", required=True)
    complete.add_argument("--output-prefix", required=True, type=Path)

    args = parser.parse_args()
    if args.command == "prepare":
        build_single_exon_candidates(
            args.merged_hits, args.catalog, args.reference_or, args.species_or,
            args.species, args.output_prefix, args.clearance, args.flank,
        )
    elif args.command == "normalize-targets":
        normalize_extracted_targets(args.raw_fasta, args.candidates, args.output, args.summary)
    elif args.command == "search":
        run_targeted_tblastn(
            args.candidates, args.queries, args.targets, args.standard_output,
            args.gapless_output, args.summary, args.gate, args.evalue, args.workers,
            args.timeout_seconds, args.limit,
        )
    elif args.command == "classify":
        classify_recovered_candidates(
            args.candidates, args.targets, args.standard_tsv, args.gapless_tsv,
            args.output_prefix,
        )
    else:
        build_complete_catalog(
            args.species, args.existing_catalog, args.existing_gff3,
            args.existing_functional_nucleotide, args.existing_functional_protein,
            args.existing_pseudogene, args.existing_truncated, args.existing_ambiguous,
            args.recovery_calls, args.recovered_pseudogene, args.recovered_truncated,
            args.recovered_ambiguous, args.recovery_blastx, args.output_prefix,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
