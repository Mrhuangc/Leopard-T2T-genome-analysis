"""Prepare and deterministically shard the NC Exonerate protein query set."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from or_mining.fasta import read_fasta, write_fasta_record


def prepare_exonerate_queries(
    reference_or: Path,
    species_or: Path,
    output_dir: Path,
    summary_path: Path,
    requested_chunks: int = 32,
) -> dict[str, Any]:
    if requested_chunks < 1:
        raise ValueError("requested_chunks must be at least 1")

    reference_records = list(read_fasta(reference_or))
    species_records = list(read_fasta(species_or))
    records = reference_records + species_records
    if not records:
        raise ValueError("combined Exonerate query set is empty")

    identifiers = [header.split()[0] for header, _ in records]
    duplicate_ids = sorted(
        identifier for identifier in set(identifiers) if identifiers.count(identifier) > 1
    )
    if duplicate_ids:
        preview = ", ".join(duplicate_ids[:10])
        raise ValueError(f"duplicate FASTA identifiers in combined query set: {preview}")
    if any(not sequence for _, sequence in records):
        raise ValueError("empty protein sequence in combined Exonerate query set")

    chunk_count = min(requested_chunks, len(records))
    output_dir.mkdir(parents=True, exist_ok=True)
    chunk_paths = [output_dir / f"chunk_{index:03d}.prot.fa" for index in range(1, chunk_count + 1)]
    handles = [path.open("wt", encoding="utf-8", newline="\n") for path in chunk_paths]
    chunk_records = [0] * chunk_count
    chunk_residues = [0] * chunk_count
    try:
        # Greedy length balancing keeps the expensive exhaustive tasks similar in size.
        for header, sequence in sorted(records, key=lambda record: len(record[1]), reverse=True):
            index = min(range(chunk_count), key=lambda item: chunk_residues[item])
            write_fasta_record(handles[index], header, sequence)
            chunk_records[index] += 1
            chunk_residues[index] += len(sequence)
    finally:
        for handle in handles:
            handle.close()

    summary: dict[str, Any] = {
        "schema_version": "1.0",
        "status": "pass",
        "method": "NC Complete_OR_db.prot with deterministic length-balanced sharding",
        "counts": {
            "reference_or_records": len(reference_records),
            "species_or_records": len(species_records),
            "combined_query_records": len(records),
            "query_chunks": chunk_count,
            "combined_residues": sum(len(sequence) for _, sequence in records),
            "minimum_chunk_records": min(chunk_records),
            "maximum_chunk_records": max(chunk_records),
            "minimum_chunk_residues": min(chunk_residues),
            "maximum_chunk_residues": max(chunk_residues),
        },
        "parameters": {"requested_chunks": requested_chunks},
    }
    with summary_path.open("wt", encoding="utf-8", newline="\n") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
    return summary
