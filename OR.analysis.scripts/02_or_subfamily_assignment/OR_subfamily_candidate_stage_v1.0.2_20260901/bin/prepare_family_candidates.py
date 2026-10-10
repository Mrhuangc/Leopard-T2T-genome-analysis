#!/usr/bin/env python3
"""Prepare audited family-level datasets for felid OR subfamily phylogenies."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

CANONICAL_AA = set("ACDEFGHIKLMNPQRSTVWY")
DIAMOND_FIELDS = [
    "qseqid", "sseqid", "pident", "length", "qlen", "slen",
    "qcovhsp", "scovhsp", "evalue", "bitscore",
]


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    current: str | None = None
    chunks: list[str] = []
    for raw in path.read_text(encoding="ascii").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current is not None:
                records[current] = "".join(chunks)
            current = line[1:].split()[0]
            if not current or current in records:
                raise ValueError(f"empty or duplicate FASTA ID in {path}: {current}")
            chunks = []
        elif current is None:
            raise ValueError(f"sequence before first header in {path}")
        else:
            chunks.append(line.upper())
    if current is not None:
        records[current] = "".join(chunks)
    if not records:
        raise ValueError(f"empty FASTA: {path}")
    for identifier, sequence in records.items():
        if not sequence:
            raise ValueError(f"empty sequence: {identifier}")
        invalid = sorted(set(sequence) - CANONICAL_AA)
        if invalid:
            raise ValueError(f"non-canonical residues in {identifier}: {''.join(invalid)}")
    return records


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({field: row.get(field, "") for field in fields} for row in rows)


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for identifier, sequence in records:
            handle.write(f">{identifier}\n")
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start:start + 60] + "\n")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def expected_broad_clade(or_family: str) -> str:
    number = int(or_family[2:])
    return "alpha" if 51 <= number <= 56 else "gamma"


def family_sort_key(family: str) -> int:
    return int(family[2:])


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--query-fasta", type=Path, required=True)
    parser.add_argument("--reference-fasta", type=Path, required=True)
    parser.add_argument("--reference-metadata", type=Path, required=True)
    parser.add_argument("--broad-clade-tsv", type=Path, required=True)
    parser.add_argument("--diamond-tsv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--min-query-cover", type=float, default=60.0)
    parser.add_argument("--min-subject-cover", type=float, default=60.0)
    parser.add_argument("--max-evalue", type=float, default=1e-10)
    parser.add_argument("--candidate-score-ratio", type=float, default=0.90)
    parser.add_argument("--max-candidate-families", type=int, default=3)
    args = parser.parse_args()

    outdir = args.output_dir.resolve()
    fasta_dir = outdir / "family_fastas"
    fasta_dir.mkdir(parents=True, exist_ok=True)
    queries = read_fasta(args.query_fasta.resolve())
    references = read_fasta(args.reference_fasta.resolve())

    metadata_rows = read_tsv(args.reference_metadata.resolve())
    metadata = {row["reference_id"]: row for row in metadata_rows}
    if len(metadata) != len(metadata_rows):
        raise ValueError("duplicate reference_id in reference metadata")
    if set(metadata) != set(references):
        raise ValueError("reference FASTA IDs do not equal reference metadata IDs")

    broad_rows = read_tsv(args.broad_clade_tsv.resolve())
    broad = {row["representative_id"]: row["joint_clade"].lower() for row in broad_rows}
    if len(broad) != len(broad_rows):
        raise ValueError("duplicate representative_id in broad-clade table")
    missing_broad = sorted(set(queries) - set(broad))
    extra_broad = sorted(set(broad) - set(queries))

    hits: dict[str, list[dict[str, object]]] = defaultdict(list)
    unknown_queries: set[str] = set()
    unknown_subjects: set[str] = set()
    malformed_rows = 0
    with args.diamond_tsv.resolve().open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for values in reader:
            if not values:
                continue
            if len(values) != len(DIAMOND_FIELDS):
                malformed_rows += 1
                continue
            row = dict(zip(DIAMOND_FIELDS, values))
            query_id, subject_id = row["qseqid"], row["sseqid"]
            if query_id not in queries:
                unknown_queries.add(query_id)
                continue
            if subject_id not in metadata:
                unknown_subjects.add(subject_id)
                continue
            qcov = float(row["qcovhsp"])
            scov = float(row["scovhsp"])
            evalue = float(row["evalue"])
            if qcov < args.min_query_cover or scov < args.min_subject_cover or evalue > args.max_evalue:
                continue
            hits[query_id].append(
                {
                    **row,
                    "qcovhsp": qcov,
                    "scovhsp": scov,
                    "evalue": evalue,
                    "bitscore": float(row["bitscore"]),
                }
            )

    if malformed_rows or unknown_queries or unknown_subjects:
        raise ValueError(
            "invalid DIAMOND table: "
            f"malformed={malformed_rows}, unknown_queries={len(unknown_queries)}, "
            f"unknown_subjects={len(unknown_subjects)}"
        )

    assignment_rows: list[dict[str, object]] = []
    unresolved_rows: list[dict[str, object]] = []
    family_queries: dict[str, set[str]] = defaultdict(set)
    status_counts: dict[str, int] = defaultdict(int)

    for query_id in sorted(queries):
        query_hits = sorted(hits.get(query_id, []), key=lambda row: (-float(row["bitscore"]), str(row["sseqid"])))
        observed_broad = broad.get(query_id, "missing")
        base: dict[str, object] = {
            "representative_id": query_id,
            "joint_broad_clade": observed_broad,
            "best_reference_id": "",
            "best_reference_symbol": "",
            "best_family": "",
            "best_subfamily": "",
            "expected_broad_clade": "",
            "best_pident": "",
            "best_query_cover": "",
            "best_subject_cover": "",
            "best_evalue": "",
            "best_bitscore": "",
            "candidate_families": "",
            "candidate_subfamilies": "",
            "second_family_score_ratio": "",
            "broad_clade_concordance": "not_applicable",
            "candidate_status": "",
        }
        if not query_hits:
            base["candidate_status"] = "no_qualified_reference_hit"
            status_counts[str(base["candidate_status"])] += 1
            assignment_rows.append(base)
            unresolved_rows.append(base)
            continue

        best = query_hits[0]
        best_score = float(best["bitscore"])
        family_best: dict[str, float] = {}
        subfamily_best: dict[str, float] = {}
        for hit in query_hits:
            ref = metadata[str(hit["sseqid"])]
            family = ref["or_family"]
            subfamily = ref["or_subfamily"]
            family_best[family] = max(family_best.get(family, 0.0), float(hit["bitscore"]))
            subfamily_best[subfamily] = max(subfamily_best.get(subfamily, 0.0), float(hit["bitscore"]))
        candidate_families = sorted(
            (family for family, score in family_best.items() if score / best_score >= args.candidate_score_ratio),
            key=lambda family: (-family_best[family], family_sort_key(family)),
        )[: args.max_candidate_families]
        candidate_subfamilies = sorted(
            (subfamily for subfamily, score in subfamily_best.items() if score / best_score >= args.candidate_score_ratio),
            key=lambda subfamily: (-subfamily_best[subfamily], subfamily),
        )
        best_ref = metadata[str(best["sseqid"])]
        best_family = best_ref["or_family"]
        expected = expected_broad_clade(best_family)
        concordance = "agree" if observed_broad == expected else "disagree"
        if observed_broad not in {"alpha", "gamma"}:
            status = "not_assignable_with_human_reference"
        elif concordance == "disagree":
            status = "broad_clade_conflict"
        elif len(candidate_families) == 1:
            status = "single_family_candidate"
        else:
            status = "multiple_family_candidates"
        ranked_families = sorted(family_best, key=lambda family: (-family_best[family], family_sort_key(family)))
        second_ratio = family_best[ranked_families[1]] / best_score if len(ranked_families) > 1 else 0.0
        base.update(
            {
                "best_reference_id": best["sseqid"],
                "best_reference_symbol": best_ref["approved_symbol"],
                "best_family": best_family,
                "best_subfamily": best_ref["or_subfamily"],
                "expected_broad_clade": expected,
                "best_pident": best["pident"],
                "best_query_cover": f"{float(best['qcovhsp']):.2f}",
                "best_subject_cover": f"{float(best['scovhsp']):.2f}",
                "best_evalue": best["evalue"],
                "best_bitscore": f"{best_score:.2f}",
                "candidate_families": ";".join(candidate_families),
                "candidate_subfamilies": ";".join(candidate_subfamilies),
                "second_family_score_ratio": f"{second_ratio:.4f}",
                "broad_clade_concordance": concordance,
                "candidate_status": status,
            }
        )
        assignment_rows.append(base)
        status_counts[status] += 1
        if status in {"no_qualified_reference_hit", "not_assignable_with_human_reference", "broad_clade_conflict"}:
            unresolved_rows.append(base)
        if status != "not_assignable_with_human_reference":
            for family in candidate_families:
                family_queries[family].add(query_id)

    assignment_fields = [
        "representative_id", "joint_broad_clade", "best_reference_id",
        "best_reference_symbol", "best_family", "best_subfamily",
        "expected_broad_clade", "best_pident", "best_query_cover",
        "best_subject_cover", "best_evalue", "best_bitscore",
        "candidate_families", "candidate_subfamilies",
        "second_family_score_ratio", "broad_clade_concordance", "candidate_status",
    ]
    assignment_path = outdir / "candidate_assignments.tsv"
    unresolved_path = outdir / "unresolved_queries.tsv"
    write_tsv(assignment_path, assignment_rows, assignment_fields)
    write_tsv(unresolved_path, unresolved_rows, assignment_fields)

    reference_by_family: dict[str, list[str]] = defaultdict(list)
    for ref_id, row in metadata.items():
        reference_by_family[row["or_family"]].append(ref_id)

    manifest_rows: list[dict[str, object]] = []
    for family in sorted(family_queries, key=family_sort_key):
        ref_ids = sorted(reference_by_family[family])
        query_ids = sorted(family_queries[family])
        output = fasta_dir / f"{family}.combined.prot.fa"
        combined = [(ref_id, references[ref_id]) for ref_id in ref_ids]
        combined.extend((query_id, queries[query_id]) for query_id in query_ids)
        write_fasta(output, combined)
        manifest_rows.append(
            {
                "or_family": family,
                "expected_broad_clade": expected_broad_clade(family),
                "human_reference_sequences": len(ref_ids),
                "felid_candidate_sequences": len(query_ids),
                "combined_sequences": len(combined),
                "combined_fasta": output.relative_to(outdir).as_posix(),
                "combined_fasta_sha256": sha256(output),
            }
        )
    manifest_path = outdir / "family_manifest.tsv"
    write_tsv(
        manifest_path,
        manifest_rows,
        [
            "or_family", "expected_broad_clade", "human_reference_sequences",
            "felid_candidate_sequences", "combined_sequences", "combined_fasta",
            "combined_fasta_sha256",
        ],
    )

    checks = {
        "query_assignment_count_matches": len(assignment_rows) == len(queries),
        "query_assignment_ids_match": {str(row["representative_id"]) for row in assignment_rows} == set(queries),
        "broad_clade_ids_complete": not missing_broad,
        "broad_clade_has_no_extra_ids": not extra_broad,
        "reference_metadata_fasta_ids_match": set(metadata) == set(references),
        "family_fastas_created": bool(manifest_rows),
    }
    summary = {
        "schema_version": "1.0",
        "status": "pass" if all(checks.values()) else "fail",
        "method": "DIAMOND candidate retrieval followed by family-level phylogenetic preparation",
        "parameters": {
            "minimum_query_cover_percent": args.min_query_cover,
            "minimum_subject_cover_percent": args.min_subject_cover,
            "maximum_evalue": args.max_evalue,
            "candidate_bitscore_ratio": args.candidate_score_ratio,
            "maximum_candidate_families": args.max_candidate_families,
        },
        "counts": {
            "query_proteins": len(queries),
            "human_reference_proteins": len(references),
            "family_datasets": len(manifest_rows),
            "unresolved_or_flagged_queries": len(unresolved_rows),
            "missing_broad_clade_ids": len(missing_broad),
            "extra_broad_clade_ids": len(extra_broad),
            "candidate_status": dict(sorted(status_counts.items())),
        },
        "checks": checks,
        "interpretation": "Candidate labels are not final subfamily assignments; family phylogenies are required.",
    }
    marker = outdir / "PREPARATION_COMPLETE.json"
    marker.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(
        f"PREPARATION_{summary['status'].upper()} queries={len(queries)} "
        f"families={len(manifest_rows)} flagged={len(unresolved_rows)} output={outdir}"
    )
    return 0 if summary["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
