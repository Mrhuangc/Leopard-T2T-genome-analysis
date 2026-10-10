#!/usr/bin/env python3
"""Validate frozen macro-synteny inputs without accessing server genome FASTAs."""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path


UNITS = ("Fca", "Lge", "Pbe", "Ple", "Ppa", "Ppa_or", "Pti")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    manifest = read_tsv(Path(os.environ.get("REGIONS_TSV", root / "input" / "regions.tsv")))
    architecture = Path(os.environ.get("ARCHITECTURE_DIR", root / "input" / "architecture"))
    problems: list[str] = []
    if tuple(row["analysis_unit"] for row in manifest) != UNITS:
        problems.append("regions.tsv must contain the seven expected units in fixed order")

    for row in manifest:
        unit = row["analysis_unit"]
        cluster_path = architecture / f"{unit}.OR_operational_clusters.tsv"
        locus_path = architecture / f"{unit}.OR_locus_architecture.tsv"
        if not cluster_path.is_file() or cluster_path.stat().st_size == 0:
            problems.append(f"missing cluster table for {unit}")
            continue
        if not locus_path.is_file() or locus_path.stat().st_size == 0:
            problems.append(f"missing locus table for {unit}")
            continue
        clusters = read_tsv(cluster_path)
        top2 = sorted(clusters, key=lambda x: int(x["member_count"]), reverse=True)[:2]
        top2_ids = {x["cluster_id"] for x in top2}
        manifest_ids = {row["left_cluster_id"], row["right_cluster_id"]}
        if top2_ids != manifest_ids:
            problems.append(f"manifest cluster IDs are not the two largest clusters for {unit}")
        expected_start = max(1, min(int(x["start_1"]) for x in top2) - 1_000_000)
        expected_end = max(int(x["end_1"]) for x in top2) + 1_000_000
        if (int(row["start_1"]), int(row["end_1"])) != (expected_start, expected_end):
            problems.append(f"manifest interval is not top-two union plus 1 Mb flanks for {unit}")
        loci = read_tsv(locus_path)
        if not loci or any(int(x["start_1"]) < 1 or int(x["end_1"]) < int(x["start_1"]) for x in loci):
            problems.append(f"invalid locus coordinates for {unit}")

    payload = {
        "schema_version": "1.0",
        "status": "pass" if not problems else "fail",
        "counts": {"analysis_units": len(manifest), "problems": len(problems)},
        "checks": {
            "seven_expected_units": len(manifest) == 7,
            "top_two_cluster_ids_verified": not any("two largest" in x for x in problems),
            "interval_rule_verified": not any("top-two union" in x for x in problems),
            "locus_coordinates_valid": not any("locus coordinates" in x for x in problems),
        },
        "problems": problems,
    }
    out = root / "PACKAGE_VALIDATION.json"
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"status={payload['status']} report={out}")
    return 0 if not problems else 1


if __name__ == "__main__":
    raise SystemExit(main())
