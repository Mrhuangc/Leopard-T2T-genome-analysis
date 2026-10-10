#!/usr/bin/env python3
"""Create mutually exclusive, evidence-labelled draft ChrY sequence classes."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from common import (
    add_interval_coverage,
    covered_bases,
    interval_overlap,
    is_real_path,
    load_config,
    merge_intervals,
    parse_paf,
    read_fai,
)


def read_best_par(path: str) -> tuple[int, int] | None:
    with open(path, "rt", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        return None
    row = next((item for item in rows if item.get("status") == "best_terminal_candidate"), rows[0])
    return int(row["y_start"]), int(row["y_end"])


def read_repeat_windows(path: str) -> dict[tuple[int, int], tuple[float | None, float | None]]:
    data: dict[tuple[int, int], tuple[float | None, float | None]] = {}
    with open(path, "rt", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            repeat = None if row["repeat_fraction"] == "NA" else float(row["repeat_fraction"])
            satellite = None if row["satellite_fraction"] == "NA" else float(row["satellite_fraction"])
            data[(int(row["start"]), int(row["end"]))] = (repeat, satellite)
    return data


def values_for_window(data: dict[tuple[int, int], tuple[float | None, float | None]], start: int, end: int) -> tuple[float | None, float | None]:
    overlaps = []
    for (r_start, r_end), values in data.items():
        overlap = interval_overlap(start, end, r_start, r_end)
        if overlap:
            overlaps.append((overlap, values))
    if not overlaps:
        return None, None
    total = sum(item[0] for item in overlaps)
    repeat_values = [(bp, vals[0]) for bp, vals in overlaps if vals[0] is not None]
    sat_values = [(bp, vals[1]) for bp, vals in overlaps if vals[1] is not None]
    repeat = sum(bp * value for bp, value in repeat_values) / total if repeat_values else None
    satellite = sum(bp * value for bp, value in sat_values) / total if sat_values else None
    return repeat, satellite


def read_bed_intervals(path: str, chrom: str) -> list[tuple[int, int]]:
    intervals = []
    with open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3 and fields[0] == chrom:
                intervals.append((int(fields[1]), int(fields[2])))
    return merge_intervals(intervals)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--fai", required=True)
    parser.add_argument("--xy-paf", required=True)
    parser.add_argument("--self-paf", required=True)
    parser.add_argument("--par", required=True)
    parser.add_argument("--repeats", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    cfg = load_config(args.config)
    params = cfg["parameters"]
    chrom = cfg["reference"]["chr_y"]
    length = read_fai(args.fai)[chrom]
    window = int(params.get("class_window", 5_000))
    n_windows = (length + window - 1) // window
    par = read_best_par(args.par)
    repeat_data = read_repeat_windows(args.repeats)

    self_buckets = [[] for _ in range(n_windows)]
    min_amp_id = float(params.get("ampliconic_min_identity", 0.99))
    for rec in parse_paf(args.self_paf):
        if rec.qname != chrom or rec.tname != chrom or rec.identity < min_amp_id:
            continue
        if rec.aln_len < max(1_000, window // 2):
            continue
        if rec.strand == "+" and abs(rec.qstart - rec.tstart) < window:
            continue
        add_interval_coverage(self_buckets, rec.qstart, rec.qend, window)
        add_interval_coverage(self_buckets, rec.tstart, rec.tend, window)

    x_intervals: list[tuple[int, int]] = []
    x_min_id = float(params.get("x_transposed_min_identity", 0.90))
    x_merge_gap = int(params.get("x_transposed_merge_gap", 15_000))
    x_min_block = int(params.get("x_transposed_min_block", 1_000_000))
    for rec in parse_paf(args.xy_paf):
        if rec.qname != chrom or rec.identity < x_min_id or rec.aln_len < 10_000:
            continue
        if par and interval_overlap(rec.qstart, rec.qend, par[0], par[1]) / rec.aln_len >= 0.5:
            continue
        x_intervals.append((rec.qstart, rec.qend))
    x_blocks = [item for item in merge_intervals(x_intervals, x_merge_gap) if item[1] - item[0] >= x_min_block]

    centromere_path = cfg.get("annotations", {}).get("putative_centromere_bed", "")
    centromere = read_bed_intervals(str(centromere_path), chrom) if is_real_path(centromere_path) else []
    amp_cov_cutoff = float(params.get("ampliconic_min_coverage", 0.50))
    repeat_cutoff = float(params.get("repeat_rich_fraction", 0.80))
    sat_cutoff = float(params.get("satellite_rich_fraction", 0.80))

    fields = ["chrom", "start", "end", "class", "ampliconic_coverage", "repeat_fraction", "satellite_fraction", "evidence"]
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for index in range(n_windows):
            start = index * window
            end = min(length, start + window)
            width = end - start
            amp_cov = min(1.0, covered_bases(self_buckets[index]) / width)
            repeat, satellite = values_for_window(repeat_data, start, end)
            if par and interval_overlap(start, end, par[0], par[1]) / width >= 0.5:
                sequence_class, evidence = "PAR_candidate", "terminal_XY_chain"
            elif any(interval_overlap(start, end, c_start, c_end) > 0 for c_start, c_end in centromere):
                sequence_class, evidence = "putative_centromere", "user_supplied_centromere_bed"
            elif any(interval_overlap(start, end, x_start, x_end) / width >= 0.5 for x_start, x_end in x_blocks):
                sequence_class, evidence = "X_transposed_candidate", "XY_identity_and_block_length"
            elif amp_cov >= amp_cov_cutoff:
                sequence_class, evidence = "high_identity_ampliconic", "non_diagonal_self_alignment"
            elif satellite is not None and satellite >= sat_cutoff:
                sequence_class, evidence = "satellite_rich", "RepeatMasker_satellite_fraction"
            elif repeat is not None and repeat >= repeat_cutoff:
                sequence_class, evidence = "repeat_rich", "RepeatMasker_total_fraction"
            else:
                sequence_class, evidence = "other_MSY", "requires_gametolog_curation_for_X_degenerate"
            writer.writerow({
                "chrom": chrom,
                "start": start,
                "end": end,
                "class": sequence_class,
                "ampliconic_coverage": f"{amp_cov:.6f}",
                "repeat_fraction": "NA" if repeat is None else f"{repeat:.6f}",
                "satellite_fraction": "NA" if satellite is None else f"{satellite:.6f}",
                "evidence": evidence,
            })


if __name__ == "__main__":
    main()
