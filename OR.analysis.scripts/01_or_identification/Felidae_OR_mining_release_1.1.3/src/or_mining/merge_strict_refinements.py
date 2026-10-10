"""CLI for validating and merging strict local Exonerate refinements."""

from __future__ import annotations

import argparse
from pathlib import Path

from or_mining.strict_refinement import merge_strict_results


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", nargs="+", required=True, type=Path)
    parser.add_argument("--statuses", nargs="+", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    args = parser.parse_args()
    merge_strict_results(
        args.results,
        args.statuses,
        args.manifest,
        args.output,
        args.summary,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
