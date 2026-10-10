"""CLI for the real-data Exonerate pilot gate."""

from __future__ import annotations

import argparse
from pathlib import Path

from or_mining.exonerate_pilot import validate_pilot


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", nargs="+", required=True, type=Path)
    parser.add_argument("--statuses", nargs="+", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--gate", required=True, type=Path)
    parser.add_argument("--expected-pairs", required=True, type=int)
    parser.add_argument("--minimum-positive-pairs", type=int, default=1)
    args = parser.parse_args()
    validate_pilot(
        args.results,
        args.statuses,
        args.output,
        args.gate,
        args.expected_pairs,
        args.minimum_positive_pairs,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
