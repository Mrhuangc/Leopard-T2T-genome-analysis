#!/usr/bin/env python3
"""Summarize RepeatMasker annotations in fixed ChrY windows."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

from common import add_interval_coverage, covered_bases, is_real_path, load_config, read_fai


CLASSES = ("LINE", "SINE", "LTR", "DNA", "Satellite", "Tandem", "Other")


def normalize_class(raw: str) -> str:
    text = raw.lower()
    if text.startswith("line"):
        return "LINE"
    if text.startswith("sine"):
        return "SINE"
    if text.startswith("ltr") or "erv" in text:
        return "LTR"
    if text.startswith("dna") or text.startswith("rc/"):
        return "DNA"
    if "satellite" in text:
        return "Satellite"
    if "simple_repeat" in text or "low_complexity" in text:
        return "Tandem"
    return "Other"


def read_repeatmasker(path: str, chrom: str) -> list[tuple[int, int, str]]:
    records: list[tuple[int, int, str]] = []
    with open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.strip().split()
            if len(fields) < 11 or not fields[0].isdigit():
                continue
            if fields[4] != chrom:
                continue
            start = int(fields[5]) - 1
            end = int(fields[6])
            if end > start:
                records.append((start, end, normalize_class(fields[10])))
    return records


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--fai", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    cfg = load_config(args.config)
    chrom = cfg["reference"]["chr_y"]
    length = read_fai(args.fai)[chrom]
    window = int(cfg.get("parameters", {}).get("repeat_window", 50_000))
    n_windows = (length + window - 1) // window
    per_class = {name: [[] for _ in range(n_windows)] for name in CLASSES}
    total = [[] for _ in range(n_windows)]

    rm_path = cfg.get("annotations", {}).get("repeatmasker_out", "")
    available = is_real_path(rm_path)
    if available:
        for start, end, repeat_class in read_repeatmasker(str(rm_path), chrom):
            add_interval_coverage(total, start, end, window)
            add_interval_coverage(per_class[repeat_class], start, end, window)

    fields = [
        "chrom", "start", "end", "repeat_fraction",
        *[f"{name.lower()}_fraction" for name in CLASSES],
        "repeat_annotation_available",
    ]
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for index in range(n_windows):
            start = index * window
            end = min(length, start + window)
            width = end - start
            row: dict[str, object] = {
                "chrom": chrom,
                "start": start,
                "end": end,
                "repeat_fraction": f"{min(1.0, covered_bases(total[index]) / width):.6f}" if available else "NA",
                "repeat_annotation_available": int(available),
            }
            for repeat_class in CLASSES:
                key = f"{repeat_class.lower()}_fraction"
                row[key] = f"{min(1.0, covered_bases(per_class[repeat_class][index]) / width):.6f}" if available else "NA"
            writer.writerow(row)


if __name__ == "__main__":
    main()
