"""NC-compatible final homology audit for the coordinate-aware OR catalog."""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Any
import csv
import json

from or_mining.fasta import read_fasta, write_fasta_record


BLAST_COLUMNS = (
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
    "qlen", "slen", "stitle",
)
CATALOG_CLASSES = ("functional", "pseudogene", "truncated", "ambiguous")


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
        required = {"catalog_id", "classification", "evidence"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"final catalog lacks required fields: {sorted(required)}")
        fields = list(reader.fieldnames)
        rows = list(reader)
    if not rows:
        raise ValueError("final catalog is empty")
    identifiers = [row["catalog_id"] for row in rows]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("final catalog contains duplicate catalog IDs")
    invalid = sorted({row["classification"] for row in rows} - set(CATALOG_CLASSES))
    if invalid:
        raise ValueError(f"final catalog contains invalid classifications: {invalid}")
    return fields, rows


def _parse_best_hits(
    path: Path,
    expected_ids: set[str],
) -> tuple[dict[str, dict[str, str]], int, int]:
    best_hits: dict[str, dict[str, str]] = {}
    hit_rows: Counter[str] = Counter()
    row_count = 0
    with path.open("rt", encoding="utf-8", newline="") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\r\n")
            if not line:
                continue
            fields = line.split("\t", maxsplit=len(BLAST_COLUMNS) - 1)
            if len(fields) != len(BLAST_COLUMNS):
                raise ValueError(
                    f"{path}: line {line_number}: expected 15 columns, found {len(fields)}"
                )
            row_count += 1
            row = dict(zip(BLAST_COLUMNS, fields, strict=True))
            query = row["qseqid"]
            if query not in expected_ids:
                raise ValueError(f"{path}: BLAST returned unknown query ID {query!r}")
            try:
                rank = (
                    float(row["evalue"]), -float(row["bitscore"]),
                    -int(row["length"]), row["sseqid"],
                )
            except ValueError as error:
                raise ValueError(f"{path}: invalid numeric BLAST field at line {line_number}") from error
            hit_rows[query] += 1
            current = best_hits.get(query)
            if current is None:
                best_hits[query] = row
                continue
            current_rank = (
                float(current["evalue"]), -float(current["bitscore"]),
                -int(current["length"]), current["sseqid"],
            )
            if rank < current_rank:
                best_hits[query] = row
    return best_hits, row_count, sum(count > 1 for count in hit_rows.values())


def _decision(hit: dict[str, str] | None) -> tuple[str, str]:
    if hit is None:
        return "exclude", "no_hit"
    annotation = f"{hit['sseqid']} {hit['stitle']}".lower()
    if "growth-factor-receptor" in annotation:
        return "exclude", "growth_factor_receptor"
    if "or-receptor" in annotation:
        return "retain", "or_receptor_best_hit"
    return "exclude", "non_or_best_hit"


def _write_filtered_gff3(source: Path, destination: Path, retained: set[str]) -> None:
    observed: set[str] = set()
    current_gene: str | None = None
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
            attributes: dict[str, str] = {}
            for item in fields[8].split(";"):
                if "=" in item:
                    key, value = item.split("=", 1)
                    attributes[key] = value
            if fields[2] == "gene":
                current_gene = attributes.get("ID")
                if not current_gene:
                    raise ValueError(f"{source}: gene without ID at line {line_number}")
                observed.add(current_gene)
            else:
                parents = [parent for parent in attributes.get("Parent", "").split(",") if parent]
                if parents:
                    current_gene = parents[0]
            if current_gene in retained:
                if fields[2] == "gene" and "final_homology=pass" not in fields[8]:
                    fields[8] += ";final_homology=pass"
                output_handle.write("\t".join(fields) + "\n")
    missing = retained - observed
    if missing:
        raise ValueError(f"final GFF3 omitted retained genes: {sorted(missing)[:10]}")


def build_nc_final_catalog(
    species: str,
    catalog: Path,
    gff3: Path,
    functional_nucleotide: Path,
    functional_protein: Path,
    pseudogene_nucleotide: Path,
    truncated_nucleotide: Path,
    ambiguous_nucleotide: Path,
    functional_blastp: Path,
    pseudogene_blastx: Path,
    truncated_blastx: Path,
    ambiguous_blastx: Path,
    output_prefix: Path,
    blast_evalue: str,
) -> dict[str, Any]:
    """Apply the published final OR-Receptor best-hit rule to every catalog locus."""

    if not species or any(character.isspace() for character in species):
        raise ValueError("species must be a non-empty identifier without whitespace")
    original_fields, rows = _read_catalog(catalog)
    class_ids = {
        classification: {
            row["catalog_id"] for row in rows if row["classification"] == classification
        }
        for classification in CATALOG_CLASSES
    }
    nucleotide_records = {
        "functional": _read_fasta_map(functional_nucleotide),
        "pseudogene": _read_fasta_map(pseudogene_nucleotide),
        "truncated": _read_fasta_map(truncated_nucleotide),
        "ambiguous": _read_fasta_map(ambiguous_nucleotide),
    }
    protein_records = _read_fasta_map(functional_protein)
    for classification in CATALOG_CLASSES:
        if set(nucleotide_records[classification]) != class_ids[classification]:
            raise ValueError(f"{classification} nucleotide FASTA and catalog IDs differ")
    if set(protein_records) != class_ids["functional"]:
        raise ValueError("functional protein FASTA and catalog IDs differ")

    blast_inputs = {
        "functional": functional_blastp,
        "pseudogene": pseudogene_blastx,
        "truncated": truncated_blastx,
        "ambiguous": ambiguous_blastx,
    }
    best_hits: dict[str, dict[str, str]] = {}
    blast_rows: dict[str, int] = {}
    multiple_hsps: dict[str, int] = {}
    for classification in CATALOG_CLASSES:
        parsed, row_count, multi_count = _parse_best_hits(
            blast_inputs[classification], class_ids[classification]
        )
        best_hits.update(parsed)
        blast_rows[classification] = row_count
        multiple_hsps[classification] = multi_count

    decisions: dict[str, tuple[str, str]] = {
        row["catalog_id"]: _decision(best_hits.get(row["catalog_id"])) for row in rows
    }
    retained_ids = {identifier for identifier, value in decisions.items() if value[0] == "retain"}
    excluded_ids = set(decisions) - retained_ids
    retained_functional = retained_ids & class_ids["functional"]
    if not retained_functional:
        raise ValueError("final homology audit retained no functional ORs")

    output_paths = {
        "catalog": Path(f"{output_prefix}.NC_final_OR_catalog.tsv"),
        "gff3": Path(f"{output_prefix}.NC_final_OR_catalog.gff3"),
        "audit": Path(f"{output_prefix}.NC_final_OR_homology_audit.tsv"),
        "excluded_catalog": Path(f"{output_prefix}.NC_excluded_nonOR_catalog.tsv"),
        "functional_nucleotide": Path(f"{output_prefix}.NC_final_functional_OR.nucl.fa"),
        "functional_protein": Path(f"{output_prefix}.NC_final_functional_OR.prot.fa"),
        "pseudogene": Path(f"{output_prefix}.NC_final_pseudogene_OR.nucl.fa"),
        "truncated": Path(f"{output_prefix}.NC_final_truncated_OR.nucl.fa"),
        "ambiguous": Path(f"{output_prefix}.NC_final_ambiguous_OR.nucl.fa"),
        "excluded_nucleotide": Path(f"{output_prefix}.NC_excluded_nonOR.nucl.fa"),
        "excluded_protein": Path(f"{output_prefix}.NC_excluded_nonOR.prot.fa"),
    }
    added_fields = (
        "final_homology_decision", "final_homology_reason", "final_homology_program",
        "final_homology_sseqid", "final_homology_stitle", "final_homology_pident",
        "final_homology_alignment_length", "final_homology_evalue",
        "final_homology_bitscore", "final_homology_query_coverage",
    )
    final_fields = original_fields + [field for field in added_fields if field not in original_fields]

    enriched_rows: list[dict[str, str]] = []
    for original in rows:
        identifier = original["catalog_id"]
        classification = original["classification"]
        decision, reason = decisions[identifier]
        hit = best_hits.get(identifier)
        enriched = dict(original)
        enriched.update({
            "final_homology_decision": decision,
            "final_homology_reason": reason,
            "final_homology_program": "blastp" if classification == "functional" else "blastx",
            "final_homology_sseqid": hit["sseqid"] if hit else "",
            "final_homology_stitle": hit["stitle"] if hit else "",
            "final_homology_pident": hit["pident"] if hit else "",
            "final_homology_alignment_length": hit["length"] if hit else "",
            "final_homology_evalue": hit["evalue"] if hit else "",
            "final_homology_bitscore": hit["bitscore"] if hit else "",
            "final_homology_query_coverage": (
                f"{(abs(int(hit['qend']) - int(hit['qstart'])) + 1) / int(hit['qlen']):.6f}"
                if hit else ""
            ),
        })
        if decision == "retain":
            enriched["evidence"] = original["evidence"] + "; final best hit OR-Receptor"
        enriched_rows.append(enriched)

    for key, selected in (("catalog", retained_ids), ("excluded_catalog", excluded_ids)):
        with output_paths[key].open("wt", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=final_fields, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            for row in enriched_rows:
                if row["catalog_id"] in selected:
                    writer.writerow(row)

    with output_paths["audit"].open("wt", encoding="utf-8", newline="") as handle:
        audit_fields = (
            "catalog_id", "classification", "decision", "reason", "program",
            "sseqid", "stitle", "pident", "alignment_length", "evalue",
            "bitscore", "query_coverage",
        )
        writer = csv.DictWriter(handle, fieldnames=audit_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in enriched_rows:
            writer.writerow({
                "catalog_id": row["catalog_id"],
                "classification": row["classification"],
                "decision": row["final_homology_decision"],
                "reason": row["final_homology_reason"],
                "program": row["final_homology_program"],
                "sseqid": row["final_homology_sseqid"],
                "stitle": row["final_homology_stitle"],
                "pident": row["final_homology_pident"],
                "alignment_length": row["final_homology_alignment_length"],
                "evalue": row["final_homology_evalue"],
                "bitscore": row["final_homology_bitscore"],
                "query_coverage": row["final_homology_query_coverage"],
            })

    _write_filtered_gff3(gff3, output_paths["gff3"], retained_ids)

    fasta_handles = {
        key: output_paths[key].open("wt", encoding="utf-8", newline="\n")
        for key in (
            "functional_nucleotide", "functional_protein", "pseudogene", "truncated",
            "ambiguous", "excluded_nucleotide", "excluded_protein",
        )
    }
    try:
        for row in rows:
            identifier = row["catalog_id"]
            classification = row["classification"]
            header, sequence = nucleotide_records[classification][identifier]
            if identifier in retained_ids:
                target = "functional_nucleotide" if classification == "functional" else classification
                write_fasta_record(fasta_handles[target], header, sequence)
                if classification == "functional":
                    protein_header, protein = protein_records[identifier]
                    write_fasta_record(fasta_handles["functional_protein"], protein_header, protein)
            else:
                rejected_header = f"{header} final_homology={decisions[identifier][1]}"
                write_fasta_record(fasta_handles["excluded_nucleotide"], rejected_header, sequence)
                if classification == "functional":
                    protein_header, protein = protein_records[identifier]
                    write_fasta_record(
                        fasta_handles["excluded_protein"],
                        f"{protein_header} final_homology={decisions[identifier][1]}",
                        protein,
                    )
    finally:
        for handle in fasta_handles.values():
            handle.close()

    retained_counts = {
        classification: len(retained_ids & class_ids[classification])
        for classification in CATALOG_CLASSES
    }
    reason_counts = Counter(reason for _, reason in decisions.values())
    if sum(retained_counts.values()) + len(excluded_ids) != len(rows):
        raise ValueError("internal error: final homology partition is incomplete")

    summary: dict[str, Any] = {
        "schema_version": "1.0",
        "status": "pass",
        "stage": "NC-compatible final OR homology filter",
        "reportable_final_catalog": True,
        "species": species,
        "method": (
            "NC OR_Finder.sh final filter: blastp for functional proteins and blastx "
            "for nonfunctional nucleotide models; retain best hits annotated OR-Receptor "
            "while excluding growth-factor-receptor"
        ),
        "parameters": {"blast_evalue": str(blast_evalue), "max_target_seqs": 1},
        "counts": {
            "input_catalog_loci": len(rows),
            "input_functional": len(class_ids["functional"]),
            "input_pseudogene": len(class_ids["pseudogene"]),
            "input_truncated": len(class_ids["truncated"]),
            "input_ambiguous": len(class_ids["ambiguous"]),
            "blastp_rows": blast_rows["functional"],
            "blastx_rows": sum(blast_rows[c] for c in CATALOG_CLASSES if c != "functional"),
            "queries_with_multiple_hsps": sum(multiple_hsps.values()),
            "queries_with_best_hit": len(best_hits),
            "retained_catalog_loci": len(retained_ids),
            "excluded_non_or_loci": len(excluded_ids),
            "retained_functional": retained_counts["functional"],
            "retained_pseudogene": retained_counts["pseudogene"],
            "retained_truncated": retained_counts["truncated"],
            "retained_ambiguous": retained_counts["ambiguous"],
            **dict(sorted(reason_counts.items())),
        },
        "checks": {
            "all_catalog_loci_audited": True,
            "all_blast_result_ids_are_catalog_ids": True,
            "catalog_partition_is_complete": True,
            "retained_functional_or_present": True,
        },
        "policies": {
            "retain": "best hit contains OR-Receptor and excludes growth-factor-receptor",
            "exclude": "no hit, non-OR best hit, or growth-factor-receptor best hit",
            "rejected_records": "retain in explicit TSV and FASTA audit outputs",
        },
        "outputs": {key: path.name for key, path in output_paths.items()},
    }
    summary_path = Path(f"{output_prefix}.NC_final_OR_catalog.summary.json")
    summary["outputs"]["summary"] = summary_path.name
    summary_path.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return summary
