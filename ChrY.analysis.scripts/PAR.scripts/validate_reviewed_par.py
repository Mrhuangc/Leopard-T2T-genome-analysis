#!/usr/bin/env python3
"""Validate a reviewed PAR BED against the chromosome and terminal repeat cap.

This checks coordinates; it does not infer the pseudoautosomal boundary from
alignment alone. Inputs are 0-based half-open BED intervals.
"""
from __future__ import annotations

import argparse
from pathlib import Path


def interval(path: Path, chrom: str = "ChrY") -> tuple[int, int]:
    records = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if fields[0].lower() in {"chrom", "#chrom"}:
            continue
        if len(fields) < 3 or fields[0] != chrom:
            raise ValueError(f"{path}: expected a BED interval on {chrom}")
        records.append((int(fields[1]), int(fields[2])))
    if len(records) != 1:
        raise ValueError(f"{path}: expected exactly one interval, got {len(records)}")
    return records[0]


def validate(par_path: Path, cap_path: Path | None, length: int,
             chrom: str = "ChrY") -> tuple[tuple[int, int], tuple[int, int] | None]:
    if length <= 0:
        raise ValueError("chromosome length must be positive")
    par = interval(par_path, chrom)
    cap = interval(cap_path, chrom) if cap_path is not None else None
    for label, (start, end) in [("PAR", par)] + ([("cap", cap)] if cap else []):
        if start < 0 or end <= start or end > length:
            raise ValueError(f"{label} interval is outside ChrY or empty: {start}-{end}")
    if cap and max(par[0], cap[0]) < min(par[1], cap[1]):
        raise ValueError("PAR overlaps the terminal repeat cap")
    return par, cap


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--par", type=Path, required=True)
    parser.add_argument("--cap", type=Path)
    parser.add_argument("--length", type=int, required=True)
    parser.add_argument("--chrom", default="ChrY")
    args = parser.parse_args()
    par, cap = validate(args.par, args.cap, args.length, args.chrom)
    print(f"PASS\t{args.chrom}\tPAR\t{par[0]}\t{par[1]}\t{par[1]-par[0]}")
    if cap:
        print(f"PASS\t{args.chrom}\tterminal_repeat_cap\t{cap[0]}\t{cap[1]}")


if __name__ == "__main__":
    main()
