#!/usr/bin/env python3
"""Validate reviewed, mutually exclusive ChrY regions (0-based half-open)."""
from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


def validate(path: Path, length: int, chrom: str = "ChrY") -> dict[str, int]:
    if length <= 0:
        raise ValueError("chromosome length must be positive")
    rows: list[tuple[int, int, str, str]] = []
    with path.open(encoding="utf-8", newline="") as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if fields[0].lower() == "chrom":
                continue
            if len(fields) < 5 or fields[0] != chrom:
                raise ValueError(f"{path}: expected chrom,start0,end0,region_id,class on {chrom}")
            rows.append((int(fields[1]), int(fields[2]), fields[3], fields[4]))
    if not rows:
        raise ValueError("empty partition")
    totals: Counter[str] = Counter()
    expected = 0
    ids: set[str] = set()
    for start, end, region_id, region_class in rows:
        if start != expected:
            raise ValueError(f"partition gap or overlap before {region_id}: {expected} vs {start}")
        if end <= start or end > length:
            raise ValueError(f"invalid region {region_id}: {start}-{end}")
        if region_id in ids:
            raise ValueError(f"duplicate region_id {region_id}")
        ids.add(region_id)
        totals[region_class] += end - start
        expected = end
    if expected != length:
        raise ValueError(f"partition ends at {expected}, expected {length}")
    return dict(totals)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bed", type=Path, required=True)
    parser.add_argument("--length", type=int, required=True)
    parser.add_argument("--chrom", default="ChrY")
    args = parser.parse_args()
    totals = validate(args.bed, args.length, args.chrom)
    writer = csv.writer(__import__("sys").stdout, delimiter="\t", lineterminator="\n")
    writer.writerow(("class", "bp"))
    for region_class, bp in sorted(totals.items()):
        writer.writerow((region_class, bp))
    writer.writerow(("TOTAL", sum(totals.values())))


if __name__ == "__main__":
    main()
