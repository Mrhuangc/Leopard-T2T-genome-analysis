"""Validate and normalize EMBOSS transeq output."""

from __future__ import annotations

from pathlib import Path
from typing import Any
import json

from or_mining.fasta import read_fasta, write_fasta_record


def normalize_transeq(
    nucleotide_fasta: Path,
    raw_protein_fasta: Path,
    output_fasta: Path,
    summary_path: Path,
) -> dict[str, Any]:
    """Restore stable ORF IDs and remove only the expected terminal stop symbol."""

    nucleotides = {header.split()[0]: sequence.upper() for header, sequence in read_fasta(nucleotide_fasta)}
    if len(nucleotides) == 0:
        raise ValueError("nucleotide ORF FASTA is empty")

    proteins: dict[str, str] = {}
    terminal_stops = 0
    translations_without_terminal_stop = 0
    with output_fasta.open("wt", encoding="utf-8", newline="\n") as output:
        for raw_header, raw_protein in read_fasta(raw_protein_fasta):
            raw_id = raw_header.split()[0]
            if not raw_id.endswith("_1"):
                raise ValueError(f"unexpected transeq identifier: {raw_id!r}")
            orf_id = raw_id[:-2]
            if orf_id not in nucleotides:
                raise ValueError(f"transeq references unknown ORF {orf_id!r}")
            if orf_id in proteins:
                raise ValueError(f"duplicate transeq output for {orf_id!r}")

            protein = raw_protein.upper()
            if protein.endswith("*"):
                protein = protein[:-1]
                terminal_stops += 1
                expected_aa = len(nucleotides[orf_id]) // 3 - 1
            else:
                translations_without_terminal_stop += 1
                expected_aa = len(nucleotides[orf_id]) // 3
            if "*" in protein:
                raise ValueError(f"internal stop codon in translated ORF {orf_id!r}")
            if len(protein) != expected_aa:
                raise ValueError(
                    f"{orf_id}: expected {expected_aa} aa after terminal stop removal, got {len(protein)}"
                )
            proteins[orf_id] = protein
            write_fasta_record(output, orf_id, protein)

    missing = sorted(set(nucleotides) - set(proteins))
    if missing:
        raise ValueError(f"transeq omitted {len(missing)} ORFs; first missing ID: {missing[0]}")

    summary: dict[str, Any] = {
        "schema_version": "1.0",
        "status": "pass",
        "counts": {
            "nucleotide_orfs": len(nucleotides),
            "protein_orfs": len(proteins),
            "terminal_stops_removed": terminal_stops,
            "translations_without_terminal_stop": translations_without_terminal_stop,
            "internal_stop_orfs": 0,
        },
    }
    with summary_path.open("wt", encoding="utf-8", newline="\n") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
    return summary
