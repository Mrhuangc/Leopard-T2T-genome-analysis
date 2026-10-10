"""CLI for deterministic Exonerate target-region sharding."""

from __future__ import annotations

import argparse
from pathlib import Path

from or_mining.exonerate_targets import prepare_exonerate_targets


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target_fasta", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--chunks", type=int, default=8)
    args = parser.parse_args()
    prepare_exonerate_targets(
        args.target_fasta,
        args.output_dir,
        args.summary,
        requested_chunks=args.chunks,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
