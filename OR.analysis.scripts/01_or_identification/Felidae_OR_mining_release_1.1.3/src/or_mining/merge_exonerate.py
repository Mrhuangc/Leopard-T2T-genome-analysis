"""Merge and audit completed Exonerate shard outputs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


def merge_outputs(
    inputs: list[Path],
    statuses: list[Path],
    output: Path,
    summary_path: Path,
) -> dict[str, Any]:
    if not inputs:
        raise ValueError("no Exonerate outputs supplied")
    if len(inputs) != len(statuses):
        raise ValueError("Exonerate output and status counts differ")

    status_records: dict[str, dict[str, Any]] = {}
    for path in statuses:
        record = json.loads(path.read_text(encoding="utf-8"))
        if record.get("status") != "pass":
            raise ValueError(f"{path}: chunk status is not pass")
        search_parameters = record.get("search_parameters")
        required_parameters = {
            "model",
            "exhaustive",
            "suboptimal",
            "bestn",
            "timeout_is_fatal",
        }
        if not isinstance(search_parameters, dict) or not required_parameters.issubset(
            search_parameters
        ):
            raise ValueError(f"{path}: missing auditable Exonerate search parameters")
        output_file = str(record["output_file"])
        if output_file in status_records:
            raise ValueError(f"duplicate chunk status for {output_file}")
        status_records[output_file] = record

    input_names = {path.name for path in inputs}
    if input_names != set(status_records):
        raise ValueError("Exonerate outputs and status manifests do not match")

    vulgar_lines = 0
    gff_dumps = 0
    total_bytes = 0
    ordered = sorted(inputs, key=lambda path: path.name)
    with output.open("wt", encoding="utf-8", newline="\n") as merged:
        for path in ordered:
            with path.open("rt", encoding="utf-8", errors="replace") as source:
                for line in source:
                    vulgar_lines += line.startswith("vulgar:")
                    gff_dumps += line.startswith("# --- START OF GFF DUMP ---")
                    merged.write(line)
            total_bytes += path.stat().st_size

    records = list(status_records.values())
    search_parameter_sets = {
        json.dumps(record.get("search_parameters", {}), sort_keys=True)
        for record in records
    }
    if len(search_parameter_sets) != 1:
        raise ValueError("Exonerate chunks used inconsistent search parameters")
    search_parameters = records[0].get("search_parameters", {})
    completed_chunks = sum(record["termination"] == "completed" for record in records)
    timed_out_chunks = sum(record["termination"] == "soft_timeout" for record in records)
    retained_blocks = sum(
        int(record["counts"]["retained_complete_alignment_blocks"]) for record in records
    )
    dropped_blocks = sum(
        int(record["counts"]["dropped_incomplete_alignment_blocks"]) for record in records
    )
    elapsed_values = [
        int(record["elapsed_seconds"])
        for record in records
        if record.get("elapsed_seconds") is not None
    ]
    if timed_out_chunks:
        raise ValueError("refusing partial Exonerate result: one or more tasks timed out")
    if completed_chunks != len(ordered):
        raise ValueError("not every Exonerate task completed")
    if retained_blocks == 0 or vulgar_lines == 0 or gff_dumps == 0:
        raise ValueError("Exonerate completed but produced no biological alignment output")

    summary: dict[str, Any] = {
        "schema_version": "1.2",
        "status": "pass",
        "method": "ordered merge of complete blocks from time-capped Exonerate chunks",
        "search_completeness": "complete",
        "search_parameters": search_parameters,
        "counts": {
            "input_chunks": len(ordered),
            "completed_chunks": completed_chunks,
            "soft_timed_out_chunks": timed_out_chunks,
            "retained_complete_alignment_blocks": retained_blocks,
            "dropped_incomplete_alignment_blocks": dropped_blocks,
            "vulgar_alignments": vulgar_lines,
            "gff_dumps": gff_dumps,
            "merged_bytes": total_bytes,
            "total_task_elapsed_seconds": sum(elapsed_values),
            "maximum_task_elapsed_seconds": max(elapsed_values, default=0),
        },
    }
    with summary_path.open("wt", encoding="utf-8", newline="\n") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", nargs="+", required=True, type=Path)
    parser.add_argument("--statuses", nargs="+", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    args = parser.parse_args()
    merge_outputs(args.inputs, args.statuses, args.output, args.summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
