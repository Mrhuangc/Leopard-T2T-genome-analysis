"""Traceable classification of translated ORFs by their best GPCR database hit."""

from __future__ import annotations

from pathlib import Path
from typing import Any
import csv
import json

from or_mining.fasta import read_fasta, write_fasta_record


BLASTP_COLUMNS = (
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qlen",
    "slen",
    "stitle",
)


def _load_fasta_by_id(path: Path) -> tuple[list[str], dict[str, tuple[str, str]]]:
    order: list[str] = []
    records: dict[str, tuple[str, str]] = {}
    for header, sequence in read_fasta(path):
        identifier = header.split()[0]
        if identifier in records:
            raise ValueError(f"duplicate FASTA identifier {identifier!r} in {path}")
        order.append(identifier)
        records[identifier] = (header, sequence)
    return order, records


def classify_best_gpcr_hits(
    blastp_tsv: Path,
    protein_fasta: Path,
    nucleotide_fasta: Path,
    orf_table: Path,
    output_prefix: Path,
) -> dict[str, Any]:
    """Classify every ORF and emit retained OR sequences plus an audit table."""

    order, proteins = _load_fasta_by_id(protein_fasta)
    _, nucleotides = _load_fasta_by_id(nucleotide_fasta)
    if set(proteins) != set(nucleotides):
        raise ValueError("protein and nucleotide ORF identifiers differ")

    with orf_table.open("rt", encoding="utf-8", newline="") as handle:
        coordinate_rows = {row["orf_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    if set(proteins) != set(coordinate_rows):
        raise ValueError("ORF coordinate table identifiers differ from FASTA identifiers")

    best_hits: dict[str, dict[str, str]] = {}
    hit_row_counts: dict[str, int] = {}
    blastp_rows = 0
    with blastp_tsv.open("rt", encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, start=1):
            blastp_rows += 1
            fields = line.rstrip("\r\n").split("\t", maxsplit=len(BLASTP_COLUMNS) - 1)
            if len(fields) != len(BLASTP_COLUMNS):
                raise ValueError(
                    f"{blastp_tsv}: line {line_number}: expected 15 columns, found {len(fields)}"
                )
            row = dict(zip(BLASTP_COLUMNS, fields, strict=True))
            query = row["qseqid"]
            if query not in proteins:
                raise ValueError(f"BLASTP returned unknown query ID {query!r}")
            hit_row_counts[query] = hit_row_counts.get(query, 0) + 1
            current = best_hits.get(query)
            candidate_rank = (
                float(row["evalue"]),
                -float(row["bitscore"]),
                -int(row["length"]),
                row["sseqid"],
            )
            if current is None:
                best_hits[query] = row
            else:
                current_rank = (
                    float(current["evalue"]),
                    -float(current["bitscore"]),
                    -int(current["length"]),
                    current["sseqid"],
                )
                if candidate_rank < current_rank:
                    best_hits[query] = row

    classification_path = Path(f"{output_prefix}.gpcr_classification.tsv")
    retained_protein_path = Path(f"{output_prefix}.OR_candidates.prot.fa")
    retained_nucleotide_path = Path(f"{output_prefix}.OR_candidates.nt.fa")
    retained_coordinates_path = Path(f"{output_prefix}.OR_candidates.tsv")
    retained_ids_path = Path(f"{output_prefix}.OR_candidates.ids.txt")
    summary_path = Path(f"{output_prefix}.gpcr_classification.summary.json")
    retained_ids: list[str] = []
    reason_counts = {"retained_or": 0, "no_hit": 0, "non_or_best_hit": 0, "growth_factor": 0}

    with classification_path.open("wt", encoding="utf-8", newline="\n") as output:
        output.write(
            "orf_id\tdecision\treason\tsseqid\tstitle\tpident\talignment_length\t"
            "evalue\tbitscore\tquery_coverage\n"
        )
        for orf_id in order:
            hit = best_hits.get(orf_id)
            if hit is None:
                decision = "reject"
                reason = "no_hit"
                values = (".", ".", ".", ".", ".", ".", ".")
            else:
                annotation = f"{hit['sseqid']} {hit['stitle']}".lower()
                if "growth-factor-receptor" in annotation:
                    decision = "reject"
                    reason = "growth_factor"
                elif "or-receptor" in annotation:
                    decision = "retain"
                    reason = "retained_or"
                    retained_ids.append(orf_id)
                else:
                    decision = "reject"
                    reason = "non_or_best_hit"
                query_coverage = int(hit["length"]) / int(hit["qlen"])
                values = (
                    hit["sseqid"],
                    hit["stitle"],
                    hit["pident"],
                    hit["length"],
                    hit["evalue"],
                    hit["bitscore"],
                    f"{query_coverage:.6f}",
                )
            reason_counts[reason] += 1
            output.write(f"{orf_id}\t{decision}\t{reason}\t" + "\t".join(values) + "\n")

    with (
        retained_protein_path.open("wt", encoding="utf-8", newline="\n") as protein_out,
        retained_nucleotide_path.open("wt", encoding="utf-8", newline="\n") as nucleotide_out,
        retained_ids_path.open("wt", encoding="utf-8", newline="\n") as ids_out,
    ):
        for orf_id in retained_ids:
            write_fasta_record(protein_out, proteins[orf_id][0], proteins[orf_id][1])
            write_fasta_record(nucleotide_out, nucleotides[orf_id][0], nucleotides[orf_id][1])
            ids_out.write(f"{orf_id}\n")

    with orf_table.open("rt", encoding="utf-8", newline="") as source, retained_coordinates_path.open(
        "wt", encoding="utf-8", newline="\n"
    ) as output:
        reader = csv.DictReader(source, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("ORF coordinate table has no header")
        writer = csv.DictWriter(output, fieldnames=reader.fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        retained_set = set(retained_ids)
        for row in reader:
            if row["orf_id"] in retained_set:
                writer.writerow(row)

    summary: dict[str, Any] = {
        "schema_version": "1.0",
        "status": "pass",
        "rule": "best hit contains OR-Receptor and excludes growth-factor-receptor",
        "counts": {
            "input_orfs": len(order),
            "blastp_rows": blastp_rows,
            "blastp_best_hits": len(best_hits),
            "queries_with_multiple_hsps": sum(count > 1 for count in hit_row_counts.values()),
            **reason_counts,
        },
    }
    with summary_path.open("wt", encoding="utf-8", newline="\n") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
    return summary
