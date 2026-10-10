"""Prepare, gate and merge bounded strict Exonerate refinements."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from or_mining.fasta import read_fasta, write_fasta_record
from or_mining.merge_exonerate import merge_outputs


STRICT_SEARCH_PARAMETERS = {
    "model": "protein2genome",
    "exhaustive": False,
    "suboptimal": False,
    "bestn": 1,
    "timeout_is_fatal": True,
}


def _unique_fasta(paths: list[Path]) -> dict[str, tuple[str, str]]:
    records: dict[str, tuple[str, str]] = {}
    for path in paths:
        for header, sequence in read_fasta(path):
            identifier = header.split()[0]
            if identifier in records:
                raise ValueError(f"duplicate FASTA identifier: {identifier}")
            if not sequence:
                raise ValueError(f"empty FASTA sequence: {identifier}")
            records[identifier] = (header, sequence.upper())
    return records


def prepare_strict_pairs(
    reference_or: Path,
    species_or: Path,
    region_fasta: Path,
    regions_tsv: Path,
    candidates_tsv: Path,
    output_dir: Path,
    manifest_path: Path,
    summary_path: Path,
) -> dict[str, Any]:
    queries = _unique_fasta([reference_or, species_or])
    targets = _unique_fasta([region_fasta])

    with regions_tsv.open("rt", encoding="utf-8", newline="") as handle:
        region_rows = list(csv.DictReader(handle, delimiter="\t"))
    regions = {row["region_id"]: row for row in region_rows}
    if len(regions) != len(region_rows):
        raise ValueError("duplicate region IDs in region table")
    if set(targets) != set(regions):
        raise ValueError("region FASTA and region table IDs differ")

    with candidates_tsv.open("rt", encoding="utf-8", newline="") as handle:
        candidates = list(csv.DictReader(handle, delimiter="\t"))
    required = {
        "candidate_id",
        "region_id",
        "query_id",
        "sseqid",
        "strand",
        "genomic_start_1",
        "genomic_end_1",
        "refinement_start_1",
        "refinement_end_1",
        "query_span_aa",
        "exonerate_score",
    }
    if not candidates or not required.issubset(candidates[0]):
        raise ValueError(f"candidate table is empty or missing columns {sorted(required)}")
    if len({row["candidate_id"] for row in candidates}) != len(candidates):
        raise ValueError("duplicate candidate IDs")

    ranked = sorted(
        candidates,
        key=lambda row: (
            -int(row["query_span_aa"]),
            -float(row["exonerate_score"]),
            row["candidate_id"],
        ),
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_fields = (
        "strict_pair_id",
        "refinement_rank",
        "candidate_id",
        "region_id",
        "query_id",
        "sseqid",
        "strand",
        "model_start_1",
        "model_end_1",
        "refinement_start_1",
        "refinement_end_1",
        "target_bases",
        "query_residues",
    )
    total_target_bases = 0
    with manifest_path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=manifest_fields, delimiter="\t")
        writer.writeheader()
        for rank, candidate in enumerate(ranked, start=1):
            candidate_id = candidate["candidate_id"]
            region_id = candidate["region_id"]
            query_id = candidate["query_id"]
            region = regions.get(region_id)
            target = targets.get(region_id)
            query = queries.get(query_id)
            if region is None or target is None:
                raise ValueError(f"candidate {candidate_id}: unknown region {region_id}")
            if query is None:
                raise ValueError(f"candidate {candidate_id}: unknown query {query_id}")
            if candidate["sseqid"] != region["sseqid"]:
                raise ValueError(f"candidate {candidate_id}: scaffold disagrees with region table")
            region_start_1 = int(region["expanded_start_1"])
            region_end_1 = int(region["expanded_end_1"])
            refinement_start_1 = max(
                region_start_1, int(candidate["refinement_start_1"])
            )
            refinement_end_1 = min(region_end_1, int(candidate["refinement_end_1"]))
            if refinement_start_1 > refinement_end_1:
                raise ValueError(f"candidate {candidate_id}: empty refinement interval")
            local_start_0 = refinement_start_1 - region_start_1
            local_end_0 = refinement_end_1 - region_start_1 + 1
            target_sequence = target[1][local_start_0:local_end_0]
            expected_bases = refinement_end_1 - refinement_start_1 + 1
            if len(target_sequence) != expected_bases:
                raise ValueError(f"candidate {candidate_id}: refinement target extraction failed")

            strict_pair_id = f"strict_{rank:06d}__{candidate_id}"
            query_path = output_dir / f"{strict_pair_id}.query.prot.fa"
            target_path = output_dir / f"{strict_pair_id}.target.nucl.fa"
            with query_path.open("wt", encoding="utf-8", newline="\n") as query_handle:
                write_fasta_record(query_handle, query[0], query[1])
            with target_path.open("wt", encoding="utf-8", newline="\n") as target_handle:
                write_fasta_record(target_handle, candidate_id, target_sequence)
            total_target_bases += expected_bases
            writer.writerow(
                {
                    "strict_pair_id": strict_pair_id,
                    "refinement_rank": rank,
                    "candidate_id": candidate_id,
                    "region_id": region_id,
                    "query_id": query_id,
                    "sseqid": candidate["sseqid"],
                    "strand": candidate["strand"],
                    "model_start_1": candidate["genomic_start_1"],
                    "model_end_1": candidate["genomic_end_1"],
                    "refinement_start_1": refinement_start_1,
                    "refinement_end_1": refinement_end_1,
                    "target_bases": expected_bases,
                    "query_residues": len(query[1]),
                }
            )

    summary: dict[str, Any] = {
        "schema_version": "1.0",
        "status": "pass",
        "method": "one evidence-selected OR query and one local target per screened locus",
        "search_mode_rationale": (
            "candidate loci already passed complete targeted discovery and "
            "NC-compatible model selection; bounded confirmation avoids repeating "
            "the infeasible exhaustive dynamic-programming search"
        ),
        "ranking": "query span, Exonerate score, candidate ID",
        "counts": {
            "input_candidates": len(candidates),
            "strict_pairs": len(candidates),
            "combined_query_records": len(queries),
            "target_regions": len(targets),
            "total_refinement_target_bases": total_target_bases,
        },
    }
    summary_path.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return summary


def _load_statuses(results: list[Path], statuses: list[Path]) -> list[dict[str, Any]]:
    if len(results) != len(statuses):
        raise ValueError("strict result and status counts differ")
    records = [json.loads(path.read_text(encoding="utf-8")) for path in statuses]
    if {path.name for path in results} != {
        str(record.get("output_file")) for record in records
    }:
        raise ValueError("strict result files and status manifests do not match")
    for record in records:
        if record.get("status") != "pass" or record.get("termination") != "completed":
            raise ValueError("strict refinement contains a failed or timed-out task")
        if record.get("search_parameters") != STRICT_SEARCH_PARAMETERS:
            raise ValueError("strict refinement search parameters disagree with the protocol")
        count = int(record.get("counts", {}).get("retained_complete_alignment_blocks", 0))
        if count != 1:
            raise ValueError(f"strict bestn=1 task retained {count} alignment blocks")
    return records


def validate_strict_pilot(
    results: list[Path],
    statuses: list[Path],
    output_path: Path,
    gate_path: Path,
    expected_pairs: int,
) -> dict[str, Any]:
    if expected_pairs < 1 or len(results) != expected_pairs:
        raise ValueError("strict pilot result count does not match expected pairs")
    records = _load_statuses(results, statuses)
    elapsed = [int(record.get("elapsed_seconds") or 0) for record in records]
    summary: dict[str, Any] = {
        "schema_version": "1.0",
        "status": "pass",
        "method": "bounded seeded bestn=1 real-candidate pilot gate",
        "search_parameters": STRICT_SEARCH_PARAMETERS,
        "counts": {
            "completed_pairs": len(records),
            "positive_pairs": len(records),
            "timed_out_pairs": 0,
            "retained_complete_alignment_blocks": len(records),
            "total_elapsed_seconds": sum(elapsed),
            "maximum_pair_elapsed_seconds": max(elapsed, default=0),
        },
    }
    output_path.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    gate_path.write_text("status=pass\n", encoding="utf-8")
    return summary


def merge_strict_results(
    results: list[Path],
    statuses: list[Path],
    manifest_path: Path,
    output_path: Path,
    summary_path: Path,
) -> dict[str, Any]:
    with manifest_path.open("rt", encoding="utf-8", newline="") as handle:
        expected_pairs = sum(1 for _row in csv.DictReader(handle, delimiter="\t"))
    if expected_pairs < 1 or len(results) != expected_pairs:
        raise ValueError("strict result count does not match expected pairs")
    records = _load_statuses(results, statuses)
    summary = merge_outputs(results, statuses, output_path, summary_path)
    summary["schema_version"] = "1.0"
    summary["stage"] = "strict local multi-exon OR refinement"
    summary["method"] = "ordered merge of bounded seeded protein2genome bestn=1 loci"
    summary["protocol_relation"] = {
        "candidate_selection": "NC-compatible",
        "published_final_replay": "exhaustive protein2genome bestn=1",
        "implemented_confirmation": "seeded protein2genome bestn=1",
        "reason": "real-candidate exhaustive pilot exceeded the bounded runtime",
    }
    summary["counts"]["expected_candidates"] = expected_pairs
    summary["counts"]["positive_candidates"] = len(records)
    summary_path.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return summary
