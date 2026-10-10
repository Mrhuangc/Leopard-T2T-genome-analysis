#!/usr/bin/env python3
"""Summarize primary long regional alignments from minimap2 PAF output."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def merge_covered(intervals: list[tuple[int, int]]) -> int:
    if not intervals:
        return 0
    merged = 0
    left, right = sorted(intervals)[0]
    for start, end in sorted(intervals)[1:]:
        if start <= right:
            right = max(right, end)
        else:
            merged += right - left
            left, right = start, end
    return merged + right - left


def parse_paf(path: Path, min_block: int, min_identity: float, min_mapq: int) -> dict:
    kept: list[list[str]] = []
    all_rows = 0
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            all_rows += 1
            fields = line.rstrip("\n").split("\t")
            tags = {tag[:2]: tag[5:] for tag in fields[12:] if len(tag) >= 5}
            aln_len = int(fields[10])
            identity = int(fields[9]) / aln_len if aln_len else 0.0
            if (
                tags.get("tp") == "P"
                and aln_len >= min_block
                and identity >= min_identity
                and int(fields[11]) >= min_mapq
            ):
                kept.append(fields)
    if not kept:
        raise ValueError(f"No primary alignment passed filters for {path}")
    query_len = int(kept[0][1])
    target_len = int(kept[0][6])
    query_cov = merge_covered([(int(r[2]), int(r[3])) for r in kept])
    target_cov = merge_covered([(int(r[7]), int(r[8])) for r in kept])
    aligned_bp = sum(int(r[10]) for r in kept)
    matched_bp = sum(int(r[9]) for r in kept)
    plus_bp = sum(int(r[10]) for r in kept if r[4] == "+")
    minus_bp = aligned_bp - plus_bp
    return {
        "raw_paf_rows": all_rows,
        "retained_primary_blocks": len(kept),
        "query_length_bp": query_len,
        "target_length_bp": target_len,
        "query_covered_bp": query_cov,
        "target_covered_bp": target_cov,
        "query_coverage": query_cov / query_len,
        "target_coverage": target_cov / target_len,
        "weighted_identity": matched_bp / aligned_bp,
        "plus_aligned_bp": plus_bp,
        "minus_aligned_bp": minus_bp,
        "dominant_orientation": "+" if plus_bp >= minus_bp else "-",
        "dominant_orientation_fraction": max(plus_bp, minus_bp) / aligned_bp,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--paf-dir", type=Path, required=True)
    parser.add_argument("--output-tsv", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--min-block", type=int, default=10_000)
    parser.add_argument("--min-identity", type=float, default=0.90)
    parser.add_argument("--min-mapq", type=int, default=20)
    args = parser.parse_args()

    rows = []
    for paf in sorted(args.paf_dir.glob("Ppa_vs_*.asm5.paf")):
        target_unit = paf.name.removeprefix("Ppa_vs_").removesuffix(".asm5.paf")
        stats = parse_paf(paf, args.min_block, args.min_identity, args.min_mapq)
        rows.append({"reference_unit": "Ppa", "comparison_unit": target_unit, **stats})
    if len(rows) != 6:
        raise ValueError(f"Expected 6 PAF comparisons, found {len(rows)}")

    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    payload = {
        "schema_version": "1.0",
        "status": "pass",
        "method": "minimap2 asm5 regional alignment; primary blocks filtered by length, identity and MAPQ",
        "parameters": {
            "minimum_alignment_block_bp": args.min_block,
            "minimum_identity": args.min_identity,
            "minimum_mapq": args.min_mapq,
        },
        "counts": {"comparisons": len(rows)},
        "comparisons": rows,
        "interpretation_boundary": "Continuous alignment metrics support regional homology and orientation; they do not by themselves prove locus-level orthology or causal rearrangement.",
    }
    args.output_json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
