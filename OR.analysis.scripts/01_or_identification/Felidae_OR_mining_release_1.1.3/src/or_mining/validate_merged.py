"""Standalone entry point used inside the Nextflow validation process."""

from __future__ import annotations

import argparse
from pathlib import Path

from or_mining.blast import validate_merged_file, write_validation_report


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate a canonical 16-column TBLASTN file")
    parser.add_argument("input", type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--max-errors", type=int, default=20)
    args = parser.parse_args()

    report = validate_merged_file(args.input, args.max_errors)
    write_validation_report(report, args.report)
    counts = report["counts"]
    print(
        f"status={report['status']} total={counts['total_rows']} "
        f"valid={counts['valid_rows']} invalid={counts['invalid_rows']}"
    )
    return 0 if report["status"] == "pass" else 2


if __name__ == "__main__":
    raise SystemExit(main())

