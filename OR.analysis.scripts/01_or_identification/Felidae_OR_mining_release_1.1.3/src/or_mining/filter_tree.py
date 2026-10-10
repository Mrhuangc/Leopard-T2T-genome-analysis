"""Standalone CLI for audited OR tree filtering."""

from __future__ import annotations

import argparse
from pathlib import Path

from or_mining.tree_filter import filter_or_clade


def main() -> int:
    parser = argparse.ArgumentParser(description="Filter OR candidates by the known-OR clade")
    parser.add_argument("tree", type=Path)
    parser.add_argument("template_alignment", type=Path)
    parser.add_argument("type1_ids", type=Path)
    parser.add_argument("type2_ids", type=Path)
    parser.add_argument("candidate_proteins", type=Path)
    parser.add_argument("candidate_nucleotides", type=Path)
    parser.add_argument("candidate_coordinates", type=Path)
    parser.add_argument("--output-prefix", required=True, type=Path)
    args = parser.parse_args()
    summary = filter_or_clade(
        args.tree,
        args.template_alignment,
        args.type1_ids,
        args.type2_ids,
        args.candidate_proteins,
        args.candidate_nucleotides,
        args.candidate_coordinates,
        args.output_prefix,
    )
    counts = summary["counts"]
    print(
        f"status=pass retained={counts['retained_candidates']} "
        f"rejected={counts['rejected_candidates']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

