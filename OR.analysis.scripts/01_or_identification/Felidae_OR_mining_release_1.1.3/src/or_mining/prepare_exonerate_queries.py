"""CLI adapter for Exonerate query preparation."""

from __future__ import annotations

import argparse
from pathlib import Path

from or_mining.exonerate_queries import prepare_exonerate_queries


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference_or", type=Path)
    parser.add_argument("species_or", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--chunks", type=int, default=32)
    args = parser.parse_args()
    prepare_exonerate_queries(
        args.reference_or,
        args.species_or,
        args.output_dir,
        args.summary,
        args.chunks,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
