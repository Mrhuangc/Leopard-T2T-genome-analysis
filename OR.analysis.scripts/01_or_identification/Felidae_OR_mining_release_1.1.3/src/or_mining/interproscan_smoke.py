"""Prepare and validate a small, real-OR InterProScan positive control."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from or_mining.fasta import read_fasta, write_fasta_record


EXPECTED_IPR = "IPR017452"
PROTEIN_ALPHABET = set("ABCDEFGHIKLMNPQRSTVWXYZUO")


def _normalise_protein(header: str, sequence: str) -> tuple[str, str]:
    identifier = header.split()[0]
    protein = sequence.upper().replace("-", "").rstrip("*")
    if not identifier:
        raise ValueError("empty protein identifier")
    if not protein:
        raise ValueError(f"{identifier}: empty protein sequence")
    invalid = sorted(set(protein) - PROTEIN_ALPHABET)
    if invalid:
        raise ValueError(f"{identifier}: invalid protein residues: {''.join(invalid)}")
    return identifier, protein


def prepare_interproscan_smoke(
    reference_or: Path,
    output_fasta: Path,
    summary_path: Path,
    record_count: int = 5,
) -> dict[str, Any]:
    """Select deterministic, evenly distributed records from the known-OR database."""

    if record_count < 1:
        raise ValueError("record_count must be positive")

    records = [_normalise_protein(header, sequence) for header, sequence in read_fasta(reference_or)]
    if len(records) < record_count:
        raise ValueError(
            f"known-OR database has {len(records)} records; {record_count} required for smoke test"
        )

    ids = [identifier for identifier, _ in records]
    if len(ids) != len(set(ids)):
        raise ValueError("known-OR database contains duplicate first-token identifiers")

    if record_count == 1:
        indices = [len(records) // 2]
    else:
        indices = [round(i * (len(records) - 1) / (record_count - 1)) for i in range(record_count)]
    selected = [records[index] for index in indices]

    with output_fasta.open("wt", encoding="utf-8", newline="\n") as handle:
        for identifier, protein in selected:
            write_fasta_record(handle, identifier, protein)

    summary: dict[str, Any] = {
        "schema_version": "1.0",
        "status": "pass",
        "test_type": "known_or_interproscan_positive_control",
        "source_database": str(reference_or),
        "selection_method": "deterministic evenly spaced records",
        "selected_ids": [identifier for identifier, _ in selected],
        "counts": {
            "source_records": len(records),
            "selected_records": len(selected),
            "selected_residues": sum(len(protein) for _, protein in selected),
        },
    }
    with summary_path.open("wt", encoding="utf-8", newline="\n") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
    return summary


def validate_interproscan_smoke(
    input_fasta: Path,
    result_tsv: Path,
    version_report: Path,
    output_summary: Path,
    gate_path: Path,
    expected_version: str,
    expected_ipr: str = EXPECTED_IPR,
) -> dict[str, Any]:
    """Require a completed local scan with a rhodopsin-like 7TM InterPro match."""

    input_ids = {header.split()[0] for header, _ in read_fasta(input_fasta)}
    matched_ids: set[str] = set()
    ipr_ids: set[str] = set()
    analyses: set[str] = set()
    row_count = 0

    with result_tsv.open("rt", encoding="utf-8", newline="") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n\r")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 11:
                raise ValueError(
                    f"{result_tsv}: expected at least 11 TSV columns at line {line_number}"
                )
            row_count += 1
            identifier = fields[0]
            matched_ids.add(identifier)
            analyses.add(fields[3])
            if len(fields) > 11 and fields[11] == expected_ipr:
                ipr_ids.add(identifier)

    unexpected_ids = matched_ids - input_ids
    version_text = version_report.read_text(encoding="utf-8", errors="replace")
    checks = {
        "input_records_present": bool(input_ids),
        "result_rows_present": row_count > 0,
        "member_database_matches_present": bool(analyses),
        "all_result_ids_are_inputs": not unexpected_ids,
        "expected_interpro_accession_present": bool(ipr_ids),
        "expected_version_present": expected_version in version_text,
    }
    failed = [name for name, passed in checks.items() if not passed]
    if failed:
        raise ValueError("InterProScan functional smoke test failed: " + ", ".join(failed))

    summary: dict[str, Any] = {
        "schema_version": "1.0",
        "status": "pass",
        "test_type": "interproscan_known_or_functional_smoke",
        "expected_interpro_accession": expected_ipr,
        "expected_interproscan_version": expected_version,
        "checks": checks,
        "counts": {
            "input_records": len(input_ids),
            "result_rows": row_count,
            "sequences_with_matches": len(matched_ids),
            "sequences_with_expected_interpro_accession": len(ipr_ids),
            "member_databases_with_matches": len(analyses),
        },
        "expected_interpro_accession_ids": sorted(ipr_ids),
        "member_databases": sorted(analyses),
        "unexpected_result_ids": sorted(unexpected_ids),
    }
    with output_summary.open("wt", encoding="utf-8", newline="\n") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
    gate_path.write_text(
        f"status=pass\nexpected_interpro_accession={expected_ipr}\n",
        encoding="utf-8",
    )
    return summary
