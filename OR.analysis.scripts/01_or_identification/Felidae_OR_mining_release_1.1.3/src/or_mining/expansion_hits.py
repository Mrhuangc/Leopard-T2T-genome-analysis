"""Filter species-query TBLASTN outliers and append valid hits to the initial table."""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Any, TextIO
import gzip
import json

from or_mining.blast import BlastFormatError, parse_blast_row


def _open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("rt", encoding="utf-8", newline="")


def _validated_lines(path: Path):
    with _open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                raise BlastFormatError(f"{path}: line {line_number}: blank line")
            try:
                hit = parse_blast_row(line)
            except BlastFormatError as exc:
                raise BlastFormatError(f"{path}: line {line_number}: {exc}") from exc
            yield line.rstrip("\r\n"), hit


def merge_species_expansion_hits(
    initial_hits: Path,
    species_hits: Path,
    output_prefix: Path,
    max_hits_per_query: int = 10_000,
) -> dict[str, Any]:
    """Apply the NC <10,000-hit query rule, then concatenate valid hit tables."""

    if max_hits_per_query < 1:
        raise ValueError("max_hits_per_query must be at least 1")

    query_counts: Counter[str] = Counter()
    species_rows = 0
    for _, hit in _validated_lines(species_hits):
        query_counts[hit.qseqid] += 1
        species_rows += 1

    kept_queries = sorted(
        query for query, count in query_counts.items() if count < max_hits_per_query
    )
    outlier_queries = sorted(set(query_counts) - set(kept_queries))
    kept_set = set(kept_queries)

    filtered_path = Path(f"{output_prefix}.species.filtered.tblastn.tsv")
    merged_path = Path(f"{output_prefix}.expanded.merged.tblastn.tsv")
    kept_path = Path(f"{output_prefix}.species.kept_query_ids.txt")
    outlier_path = Path(f"{output_prefix}.species.outlier_query_ids.txt")
    summary_path = Path(f"{output_prefix}.species_expansion.summary.json")

    kept_path.write_text("".join(f"{query}\n" for query in kept_queries), encoding="utf-8")
    outlier_path.write_text(
        "".join(f"{query}\n" for query in outlier_queries), encoding="utf-8"
    )

    filtered_rows = 0
    with filtered_path.open("wt", encoding="utf-8", newline="\n") as output:
        for line, hit in _validated_lines(species_hits):
            if hit.qseqid in kept_set:
                output.write(line + "\n")
                filtered_rows += 1

    initial_rows = 0
    with merged_path.open("wt", encoding="utf-8", newline="\n") as output:
        for line, _ in _validated_lines(initial_hits):
            initial_rows += 1
            output.write(line + "\n")
        for line, _ in _validated_lines(filtered_path):
            output.write(line + "\n")

    summary: dict[str, Any] = {
        "schema_version": "1.0",
        "status": "pass",
        "rule": "retain species queries with fewer than max_hits_per_query hits; concatenate",
        "parameters": {"max_hits_per_query": max_hits_per_query},
        "counts": {
            "initial_rows": initial_rows,
            "species_rows": species_rows,
            "species_queries": len(query_counts),
            "kept_species_queries": len(kept_queries),
            "outlier_species_queries": len(outlier_queries),
            "filtered_species_rows": filtered_rows,
            "merged_rows": initial_rows + filtered_rows,
        },
        "duplicate_policy": "retain rows exactly as the published concatenation step",
    }
    with summary_path.open("wt", encoding="utf-8", newline="\n") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
    return summary
