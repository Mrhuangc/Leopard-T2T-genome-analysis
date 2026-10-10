#!/usr/bin/env python3
"""Map OR locus midpoints through the detailed Ppa/Ppa_or CIGAR chain."""

from __future__ import annotations

import argparse
import bisect
import csv
import json
import re
from pathlib import Path


CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("paf", type=Path)
    parser.add_argument("intervals", type=Path)
    parser.add_argument("locus_map", type=Path)
    parser.add_argument("ppa_catalog", type=Path)
    parser.add_argument("ppa_or_catalog", type=Path)
    parser.add_argument("output_tsv", type=Path)
    parser.add_argument("summary_json", type=Path)
    parser.add_argument("--maximum-nearest-distance-bp", type=int, default=50000)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def tags(fields: list[str]) -> dict[str, str]:
    result = {}
    for value in fields[12:]:
        if len(value) >= 6 and value[2:5] in {":A:", ":Z:", ":i:", ":f:"}:
            result[value[:2]] = value[5:]
    return result


def cigar_ops(cigar: str) -> list[tuple[int, str]]:
    parsed = [(int(length), operation) for length, operation in CIGAR_RE.findall(cigar)]
    if "".join(f"{length}{operation}" for length, operation in parsed) != cigar:
        raise ValueError("unsupported or malformed CIGAR")
    return parsed


def consume(ops: list[tuple[int, str]]) -> tuple[int, int]:
    q = sum(length for length, op in ops if op in {"M", "I", "=", "X", "S"})
    t = sum(length for length, op in ops if op in {"M", "D", "N", "=", "X"})
    return q, t


def map_query(segment: dict[str, object], coordinate: int) -> tuple[int | None, str]:
    if not int(segment["qstart"]) <= coordinate < int(segment["qend"]):
        return None, "outside_segment"
    qpos, tpos = int(segment["qstart"]), int(segment["tstart"])
    for length, operation in segment["ops"]:
        if operation in {"M", "=", "X"}:
            if qpos <= coordinate < qpos + length:
                return tpos + coordinate - qpos, "aligned"
            qpos += length
            tpos += length
        elif operation in {"I", "S"}:
            if qpos <= coordinate < qpos + length:
                return None, "Ppa_insertion_or_unaligned_query"
            qpos += length
        elif operation in {"D", "N"}:
            tpos += length
        elif operation in {"H", "P"}:
            continue
    return None, "unmapped_within_segment"


def map_target(segment: dict[str, object], coordinate: int) -> tuple[int | None, str]:
    if not int(segment["tstart"]) <= coordinate < int(segment["tend"]):
        return None, "outside_segment"
    qpos, tpos = int(segment["qstart"]), int(segment["tstart"])
    for length, operation in segment["ops"]:
        if operation in {"M", "=", "X"}:
            if tpos <= coordinate < tpos + length:
                return qpos + coordinate - tpos, "aligned"
            qpos += length
            tpos += length
        elif operation in {"I", "S"}:
            qpos += length
        elif operation in {"D", "N"}:
            if tpos <= coordinate < tpos + length:
                return None, "Ppa_or_insertion_or_unaligned_target"
            tpos += length
        elif operation in {"H", "P"}:
            continue
    return None, "unmapped_within_segment"


def nearest(rows: list[dict[str, object]], positions: list[int], coordinate: int) -> tuple[int, int]:
    insertion = bisect.bisect_left(positions, coordinate)
    candidates = [index for index in (insertion - 1, insertion) if 0 <= index < len(rows)]
    best = min(candidates, key=lambda index: abs(positions[index] - coordinate))
    return best, abs(positions[best] - coordinate)


def main() -> int:
    config = args()
    interval_rows = read_tsv(config.intervals)
    intervals = {row["analysis_unit"]: row for row in interval_rows}
    if set(intervals) != {"Ppa", "Ppa_or"}:
        raise ValueError("interval table must contain only Ppa and Ppa_or")

    subfamily = {
        (row["analysis_unit"], row["catalog_id"]): row
        for row in read_tsv(config.locus_map)
        if row["analysis_unit"] in {"Ppa", "Ppa_or"}
    }

    def catalog(path: Path, unit: str) -> list[dict[str, object]]:
        interval_row = intervals[unit]
        contig = interval_row["contig"]
        region_start = int(interval_row["extracted_start_1"])
        region_end = int(interval_row["extracted_end_1"])
        output = []
        for row in read_tsv(path):
            start, end = int(row["start_1"]), int(row["end_1"])
            if row["sseqid"] != contig or start > region_end or end < region_start:
                continue
            annotation = subfamily.get((unit, row["catalog_id"]), {})
            row = dict(row)
            row.update(
                {
                    "unit": unit,
                    "start": start,
                    "end": end,
                    "midpoint_abs": (start + end) // 2,
                    "midpoint_rel0": (start + end) // 2 - region_start,
                    "final_subfamily": annotation.get("final_subfamily", ""),
                    "confidence": annotation.get("assignment_confidence", ""),
                    "representative_id": annotation.get("representative_id", ""),
                }
            )
            output.append(row)
        output.sort(key=lambda row: int(row["midpoint_abs"]))
        return output

    ppa = catalog(config.ppa_catalog, "Ppa")
    ppa_or = catalog(config.ppa_or_catalog, "Ppa_or")
    if len(ppa) != 65 or len(ppa_or) != 61:
        raise ValueError(f"unexpected extracted-region catalog counts: Ppa={len(ppa)}, Ppa_or={len(ppa_or)}")

    segments = []
    with config.paf.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            value_tags = tags(fields)
            if fields[4] != "+" or int(fields[11]) < 30 or "cg" not in value_tags:
                continue
            operations = cigar_ops(value_tags["cg"])
            q_consumed, t_consumed = consume(operations)
            qstart, qend = int(fields[2]), int(fields[3])
            tstart, tend = int(fields[7]), int(fields[8])
            if q_consumed != qend - qstart or t_consumed != tend - tstart:
                raise ValueError(f"CIGAR consumption mismatch on PAF line {line_number}")
            segments.append(
                {
                    "line": line_number,
                    "qstart": qstart,
                    "qend": qend,
                    "tstart": tstart,
                    "tend": tend,
                    "ops": operations,
                }
            )
    segments.sort(key=lambda segment: int(segment["qstart"]))

    def map_ppa(row: dict[str, object]) -> tuple[int | None, str]:
        coordinate = int(row["midpoint_rel0"])
        for segment in segments:
            if int(segment["qstart"]) <= coordinate < int(segment["qend"]):
                return map_query(segment, coordinate)
        return None, "between_forward_segments"

    def map_ppa_or(row: dict[str, object]) -> tuple[int | None, str]:
        coordinate = int(row["midpoint_rel0"])
        for segment in segments:
            if int(segment["tstart"]) <= coordinate < int(segment["tend"]):
                return map_target(segment, coordinate)
        return None, "between_forward_segments"

    ppa_or_positions = [int(row["midpoint_rel0"]) for row in ppa_or]
    ppa_positions = [int(row["midpoint_rel0"]) for row in ppa]
    forward_matches: dict[int, int] = {}
    reverse_matches: dict[int, int] = {}
    forward_meta: dict[int, tuple[int | None, str, int | None]] = {}

    for index, row in enumerate(ppa):
        mapped, reason = map_ppa(row)
        distance = None
        if mapped is not None:
            target_index, distance = nearest(ppa_or, ppa_or_positions, mapped)
            if distance <= config.maximum_nearest_distance_bp:
                forward_matches[index] = target_index
        forward_meta[index] = (mapped, reason, distance)

    for index, row in enumerate(ppa_or):
        mapped, _ = map_ppa_or(row)
        if mapped is not None:
            query_index, distance = nearest(ppa, ppa_positions, mapped)
            if distance <= config.maximum_nearest_distance_bp:
                reverse_matches[index] = query_index

    output_rows = []
    for index, row in enumerate(ppa):
        mapped, reason, distance = forward_meta[index]
        target_index = forward_matches.get(index)
        target = ppa_or[target_index] if target_index is not None else None
        reciprocal = target_index is not None and reverse_matches.get(target_index) == index
        exact_protein = bool(
            target
            and row["representative_id"]
            and row["representative_id"] == target["representative_id"]
        )
        output_rows.append(
            {
                "Ppa_catalog_id": row["catalog_id"],
                "Ppa_classification": row["classification"],
                "Ppa_subfamily": row["final_subfamily"],
                "Ppa_start_1": row["start"],
                "Ppa_end_1": row["end"],
                "mapping_status": reason,
                "mapped_Ppa_or_position_1": ""
                if mapped is None
                else int(intervals["Ppa_or"]["extracted_start_1"]) + mapped,
                "nearest_distance_bp": "" if distance is None else distance,
                "Ppa_or_catalog_id": "" if target is None else target["catalog_id"],
                "Ppa_or_classification": "" if target is None else target["classification"],
                "Ppa_or_subfamily": "" if target is None else target["final_subfamily"],
                "Ppa_or_start_1": "" if target is None else target["start"],
                "Ppa_or_end_1": "" if target is None else target["end"],
                "reciprocal_nearest_locus": str(reciprocal).lower(),
                "identical_protein_representative": str(exact_protein).lower(),
            }
        )

    config.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with config.output_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)

    or6c = [row for row in output_rows if row["Ppa_subfamily"] == "OR6C"]
    target_states: dict[str, int] = {}
    for row in or6c:
        if not row["Ppa_or_catalog_id"]:
            state = f"unmapped:{row['mapping_status']}"
        else:
            target_label = row["Ppa_or_subfamily"] or row["Ppa_or_classification"]
            state = f"nearest:{target_label}"
        target_states[state] = target_states.get(state, 0) + 1
    summary = {
        "schema_version": "1.0",
        "status": "pass",
        "method": "CIGAR midpoint liftover plus reciprocal nearest OR-locus audit",
        "maximum_nearest_distance_bp": config.maximum_nearest_distance_bp,
        "counts": {
            "Ppa_catalog_loci": len(ppa),
            "Ppa_or_catalog_loci": len(ppa_or),
            "Ppa_OR6C_loci": len(or6c),
            "Ppa_OR6C_reciprocal_nearest": sum(row["reciprocal_nearest_locus"] == "true" for row in or6c),
            "Ppa_OR6C_identical_protein_pair": sum(
                row["identical_protein_representative"] == "true" for row in or6c
            ),
        },
        "Ppa_OR6C_target_states": target_states,
    }
    config.summary_json.parent.mkdir(parents=True, exist_ok=True)
    config.summary_json.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
