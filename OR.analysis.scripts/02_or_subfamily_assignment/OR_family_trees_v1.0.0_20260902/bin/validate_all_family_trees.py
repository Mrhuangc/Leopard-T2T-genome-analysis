#!/usr/bin/env python3
"""Aggregate the 17 guarded family-tree completion markers."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


FAMILIES = [
    "OR1", "OR2", "OR3", "OR4", "OR5", "OR6", "OR7", "OR8", "OR9",
    "OR10", "OR11", "OR12", "OR13", "OR14", "OR51", "OR52", "OR56",
]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-root", required=True, type=Path)
    args = parser.parse_args()
    summaries = {}
    missing = []
    failed = []
    for family in FAMILIES:
        marker = args.output_root / "families" / family / "COMPLETE.json"
        if not marker.is_file() or marker.stat().st_size == 0:
            missing.append(family)
            continue
        payload = json.loads(marker.read_text(encoding="utf-8"))
        summaries[family] = {
            "status": payload.get("status"),
            "input_records": payload.get("counts", {}).get("input_records"),
            "trimmed_sites": payload.get("counts", {}).get("trimmed_sites"),
            "tree_leaves": payload.get("counts", {}).get("tree_leaves"),
            "tree": payload.get("outputs", {}).get("tree"),
        }
        if payload.get("status") != "pass":
            failed.append(family)

    checks = {
        "all_17_family_markers_present": not missing and len(summaries) == 17,
        "all_family_markers_pass": not failed and len(summaries) == 17,
        "all_tree_counts_match_inputs": len(summaries) == 17
        and all(x["input_records"] == x["tree_leaves"] for x in summaries.values()),
    }
    status = "pass" if all(checks.values()) else "fail"
    payload = {
        "schema_version": "1.0",
        "status": status,
        "stage": "formal_OR_family_trees",
        "counts": {
            "expected_families": 17,
            "completed_families": len(summaries),
            "total_family_records": sum(
                x["input_records"] or 0 for x in summaries.values()
            ),
        },
        "checks": checks,
        "missing_families": missing,
        "failed_families": failed,
        "families": summaries,
    }
    marker = args.output_root / "state" / "FORMAL_FAMILY_TREES_COMPLETE.json"
    marker.parent.mkdir(parents=True, exist_ok=True)
    marker.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(marker)
    return 0 if status == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())

