"""Deterministically shard Exonerate target regions by total nucleotide length."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from or_mining.fasta import read_fasta, write_fasta_record


def prepare_exonerate_targets(
    target_fasta: Path,
    output_dir: Path,
    summary_path: Path,
    requested_chunks: int = 8,
) -> dict[str, Any]:
    if requested_chunks < 1:
        raise ValueError("requested_chunks must be at least 1")

    records = list(read_fasta(target_fasta))
    if not records:
        raise ValueError("Exonerate target FASTA is empty")
    if any(not sequence for _, sequence in records):
        raise ValueError("empty nucleotide sequence in Exonerate target FASTA")

    identifiers = [header.split()[0] for header, _ in records]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("duplicate FASTA identifiers in Exonerate target FASTA")

    chunk_count = min(requested_chunks, len(records))
    output_dir.mkdir(parents=True, exist_ok=True)
    chunk_paths = [
        output_dir / f"target_{index:03d}.nucl.fa"
        for index in range(1, chunk_count + 1)
    ]
    handles = [path.open("wt", encoding="utf-8", newline="\n") for path in chunk_paths]
    chunk_records = [0] * chunk_count
    chunk_bases = [0] * chunk_count
    try:
        # Longest-first greedy assignment minimizes imbalance without changing records.
        for header, sequence in sorted(records, key=lambda record: len(record[1]), reverse=True):
            index = min(range(chunk_count), key=lambda item: chunk_bases[item])
            write_fasta_record(handles[index], header, sequence)
            chunk_records[index] += 1
            chunk_bases[index] += len(sequence)
    finally:
        for handle in handles:
            handle.close()

    summary: dict[str, Any] = {
        "schema_version": "1.0",
        "status": "pass",
        "method": "deterministic longest-first nucleotide-balanced target sharding",
        "counts": {
            "target_records": len(records),
            "target_bases": sum(len(sequence) for _, sequence in records),
            "target_chunks": chunk_count,
            "minimum_chunk_records": min(chunk_records),
            "maximum_chunk_records": max(chunk_records),
            "minimum_chunk_bases": min(chunk_bases),
            "maximum_chunk_bases": max(chunk_bases),
        },
        "parameters": {"requested_chunks": requested_chunks},
    }
    with summary_path.open("wt", encoding="utf-8", newline="\n") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
    return summary
