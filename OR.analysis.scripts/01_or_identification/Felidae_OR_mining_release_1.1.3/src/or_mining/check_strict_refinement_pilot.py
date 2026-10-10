"""CLI for the strict local-refinement pilot gate."""

from __future__ import annotations

import argparse
from pathlib import Path

from or_mining.strict_refinement import validate_strict_pilot


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", nargs="+", required=True, type=Path)
    parser.add_argument("--statuses", nargs="+", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--gate", required=True, type=Path)
    parser.add_argument("--expected-pairs", required=True, type=int)
    args = parser.parse_args()
    validate_strict_pilot(
        args.results, args.statuses, args.output, args.gate, args.expected_pairs
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
