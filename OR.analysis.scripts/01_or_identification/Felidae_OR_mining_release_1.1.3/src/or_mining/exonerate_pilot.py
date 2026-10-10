"""Validate a bounded real-data Exonerate pilot before releasing production tasks."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def validate_pilot(
    results: list[Path],
    statuses: list[Path],
    output_path: Path,
    gate_path: Path,
    expected_pairs: int,
    minimum_positive_pairs: int = 1,
) -> dict[str, Any]:
    if expected_pairs < 1:
        raise ValueError("expected_pairs must be positive")
    if not 1 <= minimum_positive_pairs <= expected_pairs:
        raise ValueError("minimum_positive_pairs must be within 1..expected_pairs")
    if len(results) != len(statuses) or len(statuses) != expected_pairs:
        raise ValueError("pilot result count does not match expected pair count")

    records = [json.loads(path.read_text(encoding="utf-8")) for path in statuses]
    if any(record.get("status") != "pass" for record in records):
        raise ValueError("pilot contains a failed chunk status")
    if {path.name for path in results} != {
        str(record.get("output_file")) for record in records
    }:
        raise ValueError("pilot result files and status manifests do not match")
    parameter_sets = {
        json.dumps(record.get("search_parameters", {}), sort_keys=True)
        for record in records
    }
    if len(parameter_sets) != 1:
        raise ValueError("pilot tasks used inconsistent search parameters")
    search_parameters = records[0].get("search_parameters", {})
    expected_parameters = {
        "model": "protein2genome",
        "exhaustive": False,
        "suboptimal": False,
        "bestn": None,
        "timeout_is_fatal": True,
    }
    if search_parameters != expected_parameters:
        raise ValueError("pilot did not use the bounded targeted search parameters")
    timed_out = sum(record.get("termination") != "completed" for record in records)
    alignment_counts = [
        int(record.get("counts", {}).get("retained_complete_alignment_blocks", 0))
        for record in records
    ]
    alignments = sum(alignment_counts)
    positive_pairs = sum(count > 0 for count in alignment_counts)
    vulgar = sum(
        int(record.get("counts", {}).get("retained_vulgar_alignments", 0))
        for record in records
    )
    gff = sum(
        int(record.get("counts", {}).get("retained_gff_dumps", 0))
        for record in records
    )
    elapsed_values = [
        int(record["elapsed_seconds"])
        for record in records
        if record.get("elapsed_seconds") is not None
    ]
    if timed_out:
        raise ValueError(f"pilot timed out for {timed_out} of {expected_pairs} pairs")
    if positive_pairs < minimum_positive_pairs or vulgar == 0 or gff == 0:
        raise ValueError(
            f"pilot positive pairs {positive_pairs} below required {minimum_positive_pairs}"
        )

    summary: dict[str, Any] = {
        "schema_version": "1.0",
        "status": "pass",
        "method": "bounded highest-confidence real-data Exonerate pilot gate",
        "parameters": {
            "expected_pairs": expected_pairs,
            "minimum_positive_pairs": minimum_positive_pairs,
        },
        "search_parameters": search_parameters,
        "counts": {
            "completed_pairs": expected_pairs,
            "timed_out_pairs": timed_out,
            "positive_pairs": positive_pairs,
            "retained_complete_alignment_blocks": alignments,
            "vulgar_alignments": vulgar,
            "gff_dumps": gff,
            "total_elapsed_seconds": sum(elapsed_values),
            "maximum_pair_elapsed_seconds": max(elapsed_values, default=0),
        },
    }
    output_path.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    gate_path.write_text("status=pass\n", encoding="utf-8")
    return summary
