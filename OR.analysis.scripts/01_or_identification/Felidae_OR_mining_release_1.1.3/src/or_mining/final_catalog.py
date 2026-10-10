"""Apply the NC-compatible IPR017452 gate to the pre-7TM OR catalog."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Any
import csv
import json

from or_mining.fasta import read_fasta, write_fasta_record


EXPECTED_PRE7TM_CLASSES = {
    "functional_candidate",
    "pseudogene",
    "truncated",
    "ambiguous",
}


def _read_fasta_map(path: Path) -> dict[str, tuple[str, str]]:
    records: dict[str, tuple[str, str]] = {}
    for header, sequence in read_fasta(path):
        identifier = header.split()[0]
        if identifier in records:
            raise ValueError(f"duplicate FASTA identifier {identifier!r} in {path}")
        if not sequence:
            raise ValueError(f"empty FASTA sequence for {identifier!r} in {path}")
        records[identifier] = (header, sequence.upper())
    return records


def _read_catalog(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "catalog_id", "source_id", "source_stage", "classification", "evidence",
            "sseqid", "start_1", "end_1", "strand", "exon_count",
            "requires_7tm_gate",
        }
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"pre-7TM catalog lacks fields: {sorted(required)}")
        rows = list(reader)
        fields = list(reader.fieldnames)
    if not rows:
        raise ValueError("pre-7TM catalog is empty")
    identifiers = [row["catalog_id"] for row in rows]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("pre-7TM catalog contains duplicate catalog IDs")
    invalid = sorted({row["classification"] for row in rows} - EXPECTED_PRE7TM_CLASSES)
    if invalid:
        raise ValueError(f"pre-7TM catalog contains invalid classifications: {invalid}")
    for row in rows:
        expected_gate = row["classification"] == "functional_candidate"
        observed_gate = row["requires_7tm_gate"].lower() == "true"
        if expected_gate != observed_gate:
            raise ValueError(
                f"{row['catalog_id']}: classification and requires_7tm_gate disagree"
            )
    return fields, rows


def _parse_interpro(
    path: Path,
    candidate_ids: set[str],
    expected_ipr: str,
) -> tuple[dict[str, dict[str, set[str]]], int]:
    evidence: dict[str, dict[str, set[str]]] = defaultdict(
        lambda: {"analyses": set(), "signatures": set(), "interpro": set()}
    )
    row_count = 0
    unexpected: set[str] = set()
    with path.open("rt", encoding="utf-8", newline="") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 12:
                raise ValueError(
                    f"{path}: expected at least 12 TSV columns at line {line_number}"
                )
            row_count += 1
            identifier = fields[0]
            if identifier not in candidate_ids:
                unexpected.add(identifier)
                continue
            evidence[identifier]["analyses"].add(fields[3])
            evidence[identifier]["signatures"].add(fields[4])
            if fields[11] not in {"", "-"}:
                evidence[identifier]["interpro"].add(fields[11])
    if unexpected:
        raise ValueError(
            "InterProScan output contains IDs absent from the candidate FASTA: "
            + ", ".join(sorted(unexpected)[:10])
        )
    passing = {
        identifier for identifier, values in evidence.items()
        if expected_ipr in values["interpro"]
    }
    if not passing:
        raise ValueError(
            f"formal InterProScan output contains no {expected_ipr} matches; "
            "refusing to reclassify the entire candidate set"
        )
    return evidence, row_count


def _format_header(row: dict[str, str], classification: str) -> str:
    return (
        f"{row['catalog_id']} source_id={row['source_id']} "
        f"{row['sseqid']}:{row['start_1']}-{row['end_1']} "
        f"strand={row['strand']} class={classification} exons={row['exon_count']}"
    )


def _write_gff3(
    source: Path,
    destination: Path,
    classifications: dict[str, str],
    gate_status: dict[str, str],
) -> None:
    observed_genes: set[str] = set()
    with source.open("rt", encoding="utf-8") as input_handle, destination.open(
        "wt", encoding="utf-8", newline="\n"
    ) as output_handle:
        for line_number, raw_line in enumerate(input_handle, start=1):
            line = raw_line.rstrip("\r\n")
            if not line or line.startswith("#"):
                output_handle.write(line + "\n")
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                raise ValueError(f"{source}: malformed GFF3 line {line_number}")
            if fields[2] == "gene":
                items = fields[8].split(";")
                attributes: dict[str, str] = {}
                order: list[str] = []
                for item in items:
                    if "=" not in item:
                        raise ValueError(f"{source}: malformed attributes at line {line_number}")
                    key, value = item.split("=", 1)
                    attributes[key] = value
                    order.append(key)
                identifier = attributes.get("ID", "")
                if identifier not in classifications:
                    raise ValueError(f"{source}: unknown gene ID {identifier!r}")
                observed_genes.add(identifier)
                attributes["pre7tm_classification"] = attributes.get("classification", "")
                attributes["classification"] = classifications[identifier]
                attributes["interpro_7tm"] = gate_status[identifier]
                for key in ("pre7tm_classification", "interpro_7tm"):
                    if key not in order:
                        order.append(key)
                fields[8] = ";".join(f"{key}={attributes[key]}" for key in order)
            output_handle.write("\t".join(fields) + "\n")
    if observed_genes != set(classifications):
        missing = sorted(set(classifications) - observed_genes)
        raise ValueError(f"pre-7TM GFF3 omitted catalog genes: {missing[:10]}")


def build_final_catalog(
    species: str,
    pre7tm_catalog: Path,
    pre7tm_gff3: Path,
    functional_nucleotide: Path,
    functional_protein: Path,
    pseudogene_nucleotide: Path,
    truncated_nucleotide: Path,
    ambiguous_nucleotide: Path,
    interpro_tsv: Path,
    version_report: Path,
    expected_version: str,
    output_prefix: Path,
    expected_ipr: str = "IPR017452",
) -> dict[str, Any]:
    """Build the final locus catalog after the rhodopsin-like 7TM gate."""

    if not species or any(character.isspace() for character in species):
        raise ValueError("species must be a non-empty identifier without whitespace")
    version_text = version_report.read_text(encoding="utf-8", errors="replace")
    if expected_version not in version_text:
        raise ValueError(
            f"expected InterProScan version {expected_version!r} is absent from version report"
        )

    original_fields, rows = _read_catalog(pre7tm_catalog)
    rows_by_id = {row["catalog_id"]: row for row in rows}
    class_ids = {
        classification: {
            row["catalog_id"] for row in rows if row["classification"] == classification
        }
        for classification in EXPECTED_PRE7TM_CLASSES
    }
    fasta_records = {
        "functional_nucleotide": _read_fasta_map(functional_nucleotide),
        "functional_protein": _read_fasta_map(functional_protein),
        "pseudogene": _read_fasta_map(pseudogene_nucleotide),
        "truncated": _read_fasta_map(truncated_nucleotide),
        "ambiguous": _read_fasta_map(ambiguous_nucleotide),
    }
    candidate_ids = class_ids["functional_candidate"]
    if set(fasta_records["functional_nucleotide"]) != candidate_ids:
        raise ValueError("functional nucleotide FASTA and pre-7TM catalog IDs differ")
    if set(fasta_records["functional_protein"]) != candidate_ids:
        raise ValueError("functional protein FASTA and pre-7TM catalog IDs differ")
    for classification in ("pseudogene", "truncated", "ambiguous"):
        if set(fasta_records[classification]) != class_ids[classification]:
            raise ValueError(f"{classification} FASTA and pre-7TM catalog IDs differ")

    evidence, result_rows = _parse_interpro(interpro_tsv, candidate_ids, expected_ipr)
    passing_ids = {
        identifier for identifier in candidate_ids
        if expected_ipr in evidence[identifier]["interpro"]
    }
    failing_ids = candidate_ids - passing_ids

    final_classification: dict[str, str] = {}
    gate_status: dict[str, str] = {}
    for row in rows:
        identifier = row["catalog_id"]
        if identifier in passing_ids:
            final_classification[identifier] = "functional"
            gate_status[identifier] = "pass"
        elif identifier in failing_ids:
            final_classification[identifier] = "truncated"
            gate_status[identifier] = "fail"
        else:
            final_classification[identifier] = row["classification"]
            gate_status[identifier] = "not_applicable"

    output_paths = {
        "catalog": Path(f"{output_prefix}.final_OR_catalog.tsv"),
        "gff3": Path(f"{output_prefix}.final_OR_catalog.gff3"),
        "audit": Path(f"{output_prefix}.interpro_7tm_audit.tsv"),
        "functional_nucleotide": Path(f"{output_prefix}.final_functional_OR.nucl.fa"),
        "functional_protein": Path(f"{output_prefix}.final_functional_OR.prot.fa"),
        "pseudogene": Path(f"{output_prefix}.final_pseudogene_OR.nucl.fa"),
        "truncated": Path(f"{output_prefix}.final_truncated_OR.nucl.fa"),
        "ambiguous": Path(f"{output_prefix}.final_ambiguous_OR.nucl.fa"),
        "rejected_nucleotide": Path(f"{output_prefix}.non7tm_rejected.nucl.fa"),
        "rejected_protein": Path(f"{output_prefix}.non7tm_rejected.prot.fa"),
    }

    added_fields = (
        "pre7tm_classification", "interpro_7tm_status",
        "interpro_expected_accession", "interpro_member_databases",
        "interpro_signature_accessions", "interpro_accessions",
    )
    final_fields = original_fields + [field for field in added_fields if field not in original_fields]
    with output_paths["catalog"].open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=final_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for original in rows:
            identifier = original["catalog_id"]
            row = dict(original)
            row["pre7tm_classification"] = original["classification"]
            row["classification"] = final_classification[identifier]
            row["requires_7tm_gate"] = "false"
            row["interpro_7tm_status"] = gate_status[identifier]
            row["interpro_expected_accession"] = (
                expected_ipr if identifier in passing_ids else ""
            )
            row["interpro_member_databases"] = ";".join(sorted(evidence[identifier]["analyses"]))
            row["interpro_signature_accessions"] = ";".join(
                sorted(evidence[identifier]["signatures"])
            )
            row["interpro_accessions"] = ";".join(sorted(evidence[identifier]["interpro"]))
            if identifier in passing_ids:
                row["evidence"] = original["evidence"] + f"; InterProScan {expected_ipr}"
            elif identifier in failing_ids:
                row["evidence"] = original["evidence"] + f"; failed InterProScan {expected_ipr} gate"
            writer.writerow(row)

    with output_paths["audit"].open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow((
            "catalog_id", "status", "member_database_matches", "member_databases",
            "signature_accessions", "interpro_accessions", "expected_interpro_accession",
        ))
        for identifier in sorted(candidate_ids):
            values = evidence[identifier]
            writer.writerow((
                identifier, gate_status[identifier], len(values["signatures"]),
                ";".join(sorted(values["analyses"])),
                ";".join(sorted(values["signatures"])),
                ";".join(sorted(values["interpro"])), expected_ipr,
            ))

    _write_gff3(
        pre7tm_gff3, output_paths["gff3"], final_classification, gate_status
    )

    fasta_handles = {
        key: output_paths[key].open("wt", encoding="utf-8", newline="\n")
        for key in (
            "functional_nucleotide", "functional_protein", "pseudogene", "truncated",
            "ambiguous", "rejected_nucleotide", "rejected_protein",
        )
    }
    try:
        for row in rows:
            identifier = row["catalog_id"]
            classification = final_classification[identifier]
            header = _format_header(row, classification)
            if identifier in passing_ids:
                write_fasta_record(
                    fasta_handles["functional_nucleotide"], header,
                    fasta_records["functional_nucleotide"][identifier][1],
                )
                write_fasta_record(
                    fasta_handles["functional_protein"], header,
                    fasta_records["functional_protein"][identifier][1],
                )
            elif identifier in failing_ids:
                nucleotide = fasta_records["functional_nucleotide"][identifier][1]
                protein = fasta_records["functional_protein"][identifier][1]
                write_fasta_record(fasta_handles["truncated"], header, nucleotide)
                write_fasta_record(fasta_handles["rejected_nucleotide"], header, nucleotide)
                write_fasta_record(fasta_handles["rejected_protein"], header, protein)
            else:
                sequence = fasta_records[classification][identifier][1]
                write_fasta_record(fasta_handles[classification], header, sequence)
    finally:
        for handle in fasta_handles.values():
            handle.close()

    final_counts = {
        "functional": len(passing_ids),
        "pseudogene": len(class_ids["pseudogene"]),
        "truncated": len(class_ids["truncated"]) + len(failing_ids),
        "ambiguous": len(class_ids["ambiguous"]),
    }
    if sum(final_counts.values()) != len(rows):
        raise ValueError("internal error: final classifications do not sum to catalog loci")

    summary: dict[str, Any] = {
        "schema_version": "1.0",
        "status": "pass",
        "stage": "final OR locus catalog after 7TM validation",
        "reportable_final_catalog": True,
        "species": species,
        "method": (
            "NC OR_Finder.bash-compatible InterProScan IPR017452 gate; intact ORF "
            "candidates without the expected rhodopsin-like 7TM accession are reclassified "
            "as truncated and retained in an explicit rejection audit"
        ),
        "expected_interpro_accession": expected_ipr,
        "interproscan_version": expected_version,
        "counts": {
            "catalog_loci": len(rows),
            "input_functional_candidates": len(candidate_ids),
            "interpro_result_rows": result_rows,
            "candidates_with_any_member_match": sum(bool(evidence[i]["signatures"]) for i in candidate_ids),
            "candidates_with_no_member_match": sum(not evidence[i]["signatures"] for i in candidate_ids),
            "candidates_passing_7tm_gate": len(passing_ids),
            "candidates_failing_7tm_gate": len(failing_ids),
            **final_counts,
        },
        "checks": {
            "all_candidate_proteins_evaluated": True,
            "all_interpro_result_ids_are_candidates": True,
            "expected_interpro_accession_present": True,
            "expected_interproscan_version_present": True,
            "catalog_partition_is_complete": True,
        },
        "policies": {
            "functional": f"intact ORF candidate with InterPro accession {expected_ipr}",
            "failed_7tm_gate": "reclassify as truncated; never discard silently",
            "locus_identity": "unchanged from the coordinate-aware pre-7TM catalog",
        },
        "outputs": {key: path.name for key, path in output_paths.items()},
    }
    summary_path = Path(f"{output_prefix}.final_OR_catalog.summary.json")
    summary["outputs"]["summary"] = summary_path.name
    summary_path.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return summary
