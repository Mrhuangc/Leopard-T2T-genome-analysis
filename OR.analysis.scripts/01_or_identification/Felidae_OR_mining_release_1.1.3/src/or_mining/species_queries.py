"""Prepare clear, optionally CD-HIT-reduced species-specific OR queries."""

from __future__ import annotations

from pathlib import Path
from typing import Any
import json
import shutil
import subprocess

from or_mining.fasta import read_fasta, write_fasta_record


def _records(path: Path) -> tuple[list[str], dict[str, tuple[str, str]]]:
    order: list[str] = []
    records: dict[str, tuple[str, str]] = {}
    for header, sequence in read_fasta(path):
        identifier = header.split()[0]
        if identifier in records:
            raise ValueError(f"duplicate FASTA identifier {identifier!r} in {path}")
        order.append(identifier)
        records[identifier] = (header, sequence)
    return order, records


def prepare_species_queries(
    protein_fasta: Path,
    nucleotide_fasta: Path,
    output_prefix: Path,
    cluster_threshold: int = 500,
    cluster_identity: float = 0.90,
    cpus: int = 1,
    cd_hit_executable: str = "cd-hit",
) -> dict[str, Any]:
    """Remove N-containing records and apply the published >=500 CD-HIT rule."""

    if cluster_threshold < 1:
        raise ValueError("cluster_threshold must be at least 1")
    if not 0.0 < cluster_identity <= 1.0:
        raise ValueError("cluster_identity must be in (0, 1]")
    if cpus < 1:
        raise ValueError("cpus must be at least 1")

    order, proteins = _records(protein_fasta)
    nucleotide_order, nucleotides = _records(nucleotide_fasta)
    if order != nucleotide_order or set(proteins) != set(nucleotides):
        raise ValueError("protein and nucleotide FASTA identifiers or order differ")

    clear_ids = [identifier for identifier in order if "N" not in nucleotides[identifier][1].upper()]
    clear_set = set(clear_ids)
    ambiguous_ids = [identifier for identifier in order if identifier not in clear_set]
    clear_path = Path(f"{output_prefix}.clear_OR.prot.fa")
    query_path = Path(f"{output_prefix}.species_OR_queries.prot.fa")
    cluster_path = Path(f"{query_path}.clstr")
    clear_ids_path = Path(f"{output_prefix}.clear_OR.ids.txt")
    ambiguous_ids_path = Path(f"{output_prefix}.ambiguous_OR.ids.txt")
    summary_path = Path(f"{output_prefix}.species_OR_queries.summary.json")

    with (
        clear_path.open("wt", encoding="utf-8", newline="\n") as clear_out,
        clear_ids_path.open("wt", encoding="utf-8", newline="\n") as ids_out,
    ):
        for identifier in clear_ids:
            header, sequence = proteins[identifier]
            write_fasta_record(clear_out, header, sequence)
            ids_out.write(f"{identifier}\n")
    ambiguous_ids_path.write_text(
        "".join(f"{identifier}\n" for identifier in ambiguous_ids), encoding="utf-8"
    )

    clustered = len(clear_ids) >= cluster_threshold
    if clustered:
        subprocess.run(
            [
                cd_hit_executable,
                "-i",
                str(clear_path),
                "-o",
                str(query_path),
                "-c",
                str(cluster_identity),
                "-T",
                str(cpus),
                "-M",
                "0",
            ],
            check=True,
        )
    else:
        shutil.copyfile(clear_path, query_path)
        cluster_path.write_text("", encoding="utf-8")

    query_order, _ = _records(query_path)
    if not set(query_order).issubset(clear_set):
        raise ValueError("CD-HIT output contains identifiers absent from the clear input")

    summary: dict[str, Any] = {
        "schema_version": "1.0",
        "status": "pass",
        "rule": "remove nucleotide records containing N; CD-HIT at >= threshold",
        "parameters": {
            "cluster_threshold": cluster_threshold,
            "cluster_identity": cluster_identity,
        },
        "counts": {
            "input_candidates": len(order),
            "clear_candidates": len(clear_ids),
            "ambiguous_candidates": len(ambiguous_ids),
            "query_representatives": len(query_order),
        },
        "cd_hit_applied": clustered,
    }
    with summary_path.open("wt", encoding="utf-8", newline="\n") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
    return summary
