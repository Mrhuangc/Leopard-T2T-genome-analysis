"""CLI for preparing and validating the InterProScan known-OR smoke test."""

from __future__ import annotations

import argparse
from pathlib import Path

from or_mining.interproscan_smoke import (
    prepare_interproscan_smoke,
    validate_interproscan_smoke,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare = subparsers.add_parser("prepare")
    prepare.add_argument("reference_or", type=Path)
    prepare.add_argument("--records", type=int, default=5)
    prepare.add_argument("--output", required=True, type=Path)
    prepare.add_argument("--summary", required=True, type=Path)

    validate = subparsers.add_parser("validate")
    validate.add_argument("input_fasta", type=Path)
    validate.add_argument("result_tsv", type=Path)
    validate.add_argument("version_report", type=Path)
    validate.add_argument("--expected-version", required=True)
    validate.add_argument("--output", required=True, type=Path)
    validate.add_argument("--gate", required=True, type=Path)

    args = parser.parse_args()
    if args.command == "prepare":
        prepare_interproscan_smoke(
            args.reference_or,
            args.output,
            args.summary,
            args.records,
        )
    else:
        validate_interproscan_smoke(
            args.input_fasta,
            args.result_tsv,
            args.version_report,
            args.output,
            args.gate,
            args.expected_version,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
