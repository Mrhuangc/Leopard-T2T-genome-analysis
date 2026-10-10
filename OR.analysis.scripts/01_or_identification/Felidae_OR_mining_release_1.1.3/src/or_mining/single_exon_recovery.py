"""Recover NC-compatible single-exon OR pseudogenes and genomic fragments."""

from __future__ import annotations

from bisect import bisect_left
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any
import csv
import json
import re
import subprocess

from or_mining.blast import BlastHit, _open_text, parse_blast_row
from or_mining.fasta import read_fasta, write_fasta_record


@dataclass(frozen=True, slots=True)
class Seed:
    hit: BlastHit
    expanded_start_0: int
    expanded_end_0: int


CANDIDATE_FIELDS = (
    "candidate_id", "sseqid", "strand", "query_id", "true_start_1", "true_end_1",
    "expanded_start_1", "expanded_end_1", "slen", "alignment_aa", "pident",
    "evalue", "bitscore", "faidx_region", "source_hit_count",
)
TARGETED_FIELDS = (
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart",
    "qend", "sstart", "send", "evalue", "sframe", "sseq",
)
FINAL_BLAST_FIELDS = (
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart",
    "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "stitle",
)
CLASS_PREFIXES = {
    "functional": "ORC", "pseudogene": "ORP", "truncated": "ORT", "ambiguous": "ORA",
}


def _fasta_map(path: Path, allow_empty: bool = True) -> dict[str, tuple[str, str]]:
    records: dict[str, tuple[str, str]] = {}
    for header, sequence in read_fasta(path):
        identifier = header.split()[0]
        if identifier in records:
            raise ValueError(f"duplicate FASTA ID {identifier!r} in {path}")
        if not sequence:
            raise ValueError(f"empty FASTA sequence for {identifier!r} in {path}")
        records[identifier] = (header, sequence.upper())
    if not allow_empty and not records:
        raise ValueError(f"FASTA is empty: {path}")
    return records


def _catalog_exclusion_index(
    catalog: Path, clearance: int
) -> dict[str, tuple[list[int], list[int]]]:
    intervals: dict[str, list[tuple[int, int]]] = {}
    with catalog.open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sseqid", "start_1", "end_1"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"catalog lacks required fields: {sorted(required)}")
        for row in reader:
            start_0 = max(0, int(row["start_1"]) - 1 - clearance)
            end_0 = int(row["end_1"]) + clearance
            intervals.setdefault(row["sseqid"], []).append((start_0, end_0))
    index: dict[str, tuple[list[int], list[int]]] = {}
    for subject, values in intervals.items():
        values.sort()
        starts = [value[0] for value in values]
        prefix_max: list[int] = []
        current = 0
        for _, end_0 in values:
            current = max(current, end_0)
            prefix_max.append(current)
        index[subject] = (starts, prefix_max)
    return index


def _overlaps_index(
    hit: BlastHit, index: dict[str, tuple[list[int], list[int]]]
) -> bool:
    values = index.get(hit.sseqid)
    if values is None:
        return False
    starts, prefix_max = values
    count = bisect_left(starts, hit.genome_end_0)
    return count > 0 and prefix_max[count - 1] > hit.genome_start_0


def _query_sequences(reference_or: Path, species_or: Path) -> dict[str, str]:
    combined: dict[str, str] = {}
    for path in (reference_or, species_or):
        for identifier, (_, sequence) in _fasta_map(path).items():
            previous = combined.get(identifier)
            if previous is not None and previous != sequence:
                raise ValueError(f"query ID {identifier!r} has conflicting protein sequences")
            combined[identifier] = sequence
    return combined


def build_single_exon_candidates(
    merged_hits: Path,
    catalog: Path,
    reference_or: Path,
    species_or: Path,
    species: str,
    output_prefix: Path,
    clearance: int = 100,
    flank: int = 100,
) -> dict[str, Any]:
    if clearance < 0 or flank < 0:
        raise ValueError("clearance and flank must be non-negative")
    exclusion = _catalog_exclusion_index(catalog, clearance)
    query_sequences = _query_sequences(reference_or, species_or)
    residual: list[Seed] = []
    input_hits = 0
    excluded_hits = 0
    missing_queries: set[str] = set()
    with _open_text(merged_hits) as merged_handle:
        hit_stream = (
            parse_blast_row(line) for line in merged_handle if line.strip()
        )
        for hit in hit_stream:
            input_hits += 1
            if _overlaps_index(hit, exclusion):
                excluded_hits += 1
                continue
            if hit.qseqid not in query_sequences:
                missing_queries.add(hit.qseqid)
                continue
            residual.append(
                Seed(
                    hit,
                    max(0, hit.genome_start_0 - flank),
                    min(hit.slen, hit.genome_end_0 + flank),
                )
            )
    if missing_queries:
        raise ValueError(
            "residual TBLASTN queries absent from reference/species OR FASTA: "
            + ", ".join(sorted(missing_queries)[:10])
        )

    residual.sort(key=lambda seed: (
        seed.hit.sseqid, seed.hit.strand, seed.expanded_start_0,
        seed.expanded_end_0, seed.hit.qseqid,
    ))
    groups: list[list[Seed]] = []
    current: list[Seed] = []
    current_end = -1
    for seed in residual:
        if current and (
            seed.hit.sseqid != current[0].hit.sseqid
            or seed.hit.strand != current[0].hit.strand
            or seed.expanded_start_0 > current_end
        ):
            groups.append(current)
            current = []
            current_end = -1
        current.append(seed)
        current_end = max(current_end, seed.expanded_end_0)
    if current:
        groups.append(current)

    candidates: list[dict[str, str]] = []
    selected_queries: list[tuple[str, str]] = []
    for index, group in enumerate(groups, start=1):
        best = min(group, key=lambda seed: (
            -seed.hit.length, -seed.hit.bitscore, seed.hit.evalue,
            seed.hit.qseqid, seed.hit.genome_start_0,
        ))
        hit = best.hit
        candidate_id = f"{species}_ORSE{index:06d}"
        faidx_region = f"{hit.sseqid}:{best.expanded_start_0 + 1}-{best.expanded_end_0}"
        candidates.append({
            "candidate_id": candidate_id,
            "sseqid": hit.sseqid,
            "strand": hit.strand,
            "query_id": hit.qseqid,
            "true_start_1": str(hit.genome_start_0 + 1),
            "true_end_1": str(hit.genome_end_0),
            "expanded_start_1": str(best.expanded_start_0 + 1),
            "expanded_end_1": str(best.expanded_end_0),
            "slen": str(hit.slen),
            "alignment_aa": str(hit.length),
            "pident": f"{hit.pident:.6g}",
            "evalue": f"{hit.evalue:.6g}",
            "bitscore": f"{hit.bitscore:.6g}",
            "faidx_region": faidx_region,
            "source_hit_count": str(len(group)),
        })
        selected_queries.append((candidate_id, query_sequences[hit.qseqid]))

    candidate_path = Path(f"{output_prefix}.single_exon_residual_candidates.tsv")
    faidx_path = Path(f"{output_prefix}.single_exon_residual_candidates.faidx.txt")
    query_path = Path(f"{output_prefix}.single_exon_residual_queries.prot.fa")
    summary_path = Path(f"{output_prefix}.single_exon_residual_candidates.summary.json")
    with candidate_path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CANDIDATE_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(candidates)
    faidx_path.write_text(
        "".join(row["faidx_region"] + "\n" for row in candidates), encoding="utf-8"
    )
    with query_path.open("wt", encoding="utf-8", newline="\n") as handle:
        for identifier, sequence in selected_queries:
            write_fasta_record(handle, identifier, sequence)
    summary: dict[str, Any] = {
        "schema_version": "1.0", "status": "pass", "species": species,
        "method": "NC R_script_numero2_2.R-compatible residual hit reduction",
        "parameters": {"confirmed_locus_clearance_bp": clearance, "candidate_flank_bp": flank},
        "counts": {
            "input_hits": input_hits,
            "hits_excluded_by_confirmed_catalog": excluded_hits,
            "residual_hits": len(residual),
            "residual_overlap_components": len(groups),
            "single_exon_residual_candidates": len(candidates),
        },
        "checks": {
            "all_residual_queries_resolved": True,
            "all_candidates_outside_confirmed_catalog_clearance": True,
        },
    }
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return summary


def normalize_extracted_targets(
    raw_fasta: Path, candidates_tsv: Path, output_fasta: Path, summary_path: Path
) -> dict[str, Any]:
    raw = _fasta_map(raw_fasta)
    rows = list(csv.DictReader(candidates_tsv.open("rt", encoding="utf-8", newline=""), delimiter="\t"))
    expected_regions = {row["faidx_region"] for row in rows}
    if set(raw) != expected_regions:
        raise ValueError("samtools-extracted regions and candidate manifest differ")
    total_bases = 0
    with output_fasta.open("wt", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            sequence = raw[row["faidx_region"]][1]
            expected = int(row["expanded_end_1"]) - int(row["expanded_start_1"]) + 1
            if len(sequence) != expected:
                raise ValueError(f"{row['candidate_id']}: extracted target length mismatch")
            total_bases += len(sequence)
            write_fasta_record(handle, row["candidate_id"], sequence)
    summary = {
        "schema_version": "1.0", "status": "pass",
        "counts": {"candidate_targets": len(rows), "extracted_bases": total_bases},
    }
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return summary


def _run_pair(
    candidate_id: str,
    query: str,
    target: str,
    evalue: str,
    timeout_seconds: int,
) -> tuple[str, str]:
    outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue sframe sseq"
    with TemporaryDirectory(prefix=f"{candidate_id}_") as temp_name:
        temp = Path(temp_name)
        query_path = temp / "query.fa"
        target_path = temp / "target.fa"
        query_path.write_text(f">{candidate_id}\n{query}\n", encoding="utf-8")
        target_path.write_text(f">{candidate_id}\n{target}\n", encoding="utf-8")
        base = [
            "tblastn", "-query", str(query_path), "-subject", str(target_path),
            "-evalue", str(evalue), "-outfmt", outfmt, "-num_threads", "1",
        ]
        standard = subprocess.run(
            base, check=True, capture_output=True, text=True, timeout=timeout_seconds
        ).stdout
        gapless = subprocess.run(
            base[:1] + ["-gapopen", "32767", "-gapextend", "32767"] + base[1:],
            check=True, capture_output=True, text=True, timeout=timeout_seconds,
        ).stdout
    return standard, gapless


def run_targeted_tblastn(
    candidates_tsv: Path,
    query_fasta: Path,
    target_fasta: Path,
    standard_output: Path,
    gapless_output: Path,
    summary_path: Path,
    gate_path: Path | None,
    evalue: str,
    workers: int,
    timeout_seconds: int,
    limit: int | None = None,
) -> dict[str, Any]:
    if workers < 1 or timeout_seconds < 1:
        raise ValueError("workers and timeout_seconds must be positive")
    rows = list(csv.DictReader(candidates_tsv.open("rt", encoding="utf-8", newline=""), delimiter="\t"))
    if limit is not None:
        rows = rows[:limit]
    queries = _fasta_map(query_fasta)
    targets = _fasta_map(target_fasta)
    ids = [row["candidate_id"] for row in rows]
    if not set(ids).issubset(queries) or not set(ids).issubset(targets):
        raise ValueError("selected candidates are missing from query or target FASTA")
    with ThreadPoolExecutor(max_workers=workers) as executor:
        results = list(executor.map(
            lambda identifier: _run_pair(
                identifier, queries[identifier][1], targets[identifier][1],
                evalue, timeout_seconds,
            ),
            ids,
        ))
    standard_text = "".join(result[0] for result in results)
    gapless_text = "".join(result[1] for result in results)
    standard_output.write_text(standard_text, encoding="utf-8")
    gapless_output.write_text(gapless_text, encoding="utf-8")
    standard_rows = sum(bool(line) for line in standard_text.splitlines())
    gapless_rows = sum(bool(line) for line in gapless_text.splitlines())
    status = "pass" if not ids or (standard_rows > 0 and gapless_rows > 0) else "fail"
    summary: dict[str, Any] = {
        "schema_version": "1.0", "status": status,
        "test_type": "pilot" if limit is not None else "formal_search",
        "parameters": {
            "evalue": str(evalue), "workers": workers,
            "per_pair_timeout_seconds": timeout_seconds,
        },
        "counts": {
            "selected_candidates": len(ids), "completed_candidates": len(results),
            "standard_alignment_rows": standard_rows, "gapless_alignment_rows": gapless_rows,
        },
    }
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if status != "pass":
        raise ValueError("targeted TBLASTN pilot/search produced no usable alignments")
    if gate_path is not None:
        gate_path.write_text("status=pass\n", encoding="utf-8")
    return summary


def _parse_targeted(path: Path, expected: set[str]) -> dict[str, list[dict[str, str]]]:
    rows: dict[str, list[dict[str, str]]] = {identifier: [] for identifier in expected}
    with path.open("rt", encoding="utf-8", newline="") as handle:
        for line_number, raw in enumerate(handle, start=1):
            if not raw.strip():
                continue
            fields = raw.rstrip("\r\n").split("\t", maxsplit=len(TARGETED_FIELDS) - 1)
            if len(fields) != len(TARGETED_FIELDS):
                raise ValueError(f"{path}: malformed targeted TBLASTN row {line_number}")
            row = dict(zip(TARGETED_FIELDS, fields, strict=True))
            if row["qseqid"] not in expected or row["sseqid"] != row["qseqid"]:
                raise ValueError(f"{path}: unexpected query-target pair at line {line_number}")
            rows[row["qseqid"]].append(row)
    return rows


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def classify_recovered_candidates(
    candidates_tsv: Path,
    target_fasta: Path,
    standard_tsv: Path,
    gapless_tsv: Path,
    output_prefix: Path,
) -> dict[str, Any]:
    rows = list(csv.DictReader(candidates_tsv.open("rt", encoding="utf-8", newline=""), delimiter="\t"))
    expected = {row["candidate_id"] for row in rows}
    targets = _fasta_map(target_fasta)
    if set(targets) != expected:
        raise ValueError("target FASTA and recovery candidate IDs differ")
    standard = _parse_targeted(standard_tsv, expected)
    gapless = _parse_targeted(gapless_tsv, expected)
    calls_path = Path(f"{output_prefix}.single_exon_recovery_calls.tsv")
    combined_path = Path(f"{output_prefix}.single_exon_recovery_candidates.nucl.fa")
    fasta_paths = {
        classification: Path(f"{output_prefix}.single_exon_{classification}.nucl.fa")
        for classification in ("pseudogene", "truncated", "ambiguous")
    }
    rejected_path = Path(f"{output_prefix}.single_exon_no_targeted_hit.nucl.fa")
    call_fields = (
        *CANDIDATE_FIELDS, "classification", "frameshift", "stop_codon_count",
        "distinct_subject_frames", "edge_state", "ambiguous_bases", "max_n_run",
        "targeted_hsp_count", "evidence",
    )
    output_calls: list[dict[str, str]] = []
    class_counts: Counter[str] = Counter()
    handles = {
        key: path.open("wt", encoding="utf-8", newline="\n")
        for key, path in {**fasta_paths, "rejected": rejected_path, "combined": combined_path}.items()
    }
    try:
        for row in rows:
            identifier = row["candidate_id"]
            expanded = targets[identifier][1]
            start_offset = int(row["true_start_1"]) - int(row["expanded_start_1"])
            end_offset = int(row["true_end_1"]) - int(row["expanded_start_1"]) + 1
            sequence = expanded[start_offset:end_offset]
            if row["strand"] == "-":
                sequence = _reverse_complement(sequence)
            standard_rows = standard[identifier]
            gapless_rows = gapless[identifier]
            frames = {entry["sframe"] for entry in gapless_rows}
            frameshift = len(frames) > 1
            stops = sum(entry["sseq"].count("*") for entry in standard_rows)
            n_runs = [len(value) for value in re.findall(r"N+", expanded)]
            max_n_run = max(n_runs, default=0)
            ambiguous = "N" in sequence
            edge = (
                int(row["true_start_1"]) <= 100
                or int(row["slen"]) - int(row["true_end_1"]) <= 100
                or max_n_run >= 50
            )
            if not standard_rows or not gapless_rows:
                classification = "rejected_no_targeted_hit"
                evidence = "residual seed lacked reproducible targeted TBLASTN alignment"
            elif ambiguous:
                classification = "ambiguous"
                evidence = "single-exon residual OR region contains ambiguous bases"
            elif frameshift or stops:
                classification = "pseudogene"
                defects = []
                if frameshift:
                    defects.append("multiple subject frames")
                if stops:
                    defects.append(f"{stops} translated stop(s)")
                evidence = "single-exon residual OR with LoF: " + ", ".join(defects)
            else:
                classification = "truncated"
                evidence = (
                    "single-exon residual OR fragment without detected LoF"
                    + ("; assembly-edge evidence" if edge else "")
                )
            output = dict(row)
            output.update({
                "classification": classification,
                "frameshift": str(frameshift).lower(),
                "stop_codon_count": str(stops),
                "distinct_subject_frames": ";".join(sorted(frames)),
                "edge_state": str(edge).lower(),
                "ambiguous_bases": str(ambiguous).lower(),
                "max_n_run": str(max_n_run),
                "targeted_hsp_count": str(len(standard_rows)),
                "evidence": evidence,
            })
            output_calls.append(output)
            class_counts[classification] += 1
            header = f"{identifier} {row['sseqid']}:{row['true_start_1']}-{row['true_end_1']} strand={row['strand']} class={classification}"
            if classification == "rejected_no_targeted_hit":
                write_fasta_record(handles["rejected"], header, sequence)
            else:
                write_fasta_record(handles[classification], header, sequence)
                write_fasta_record(handles["combined"], header, sequence)
    finally:
        for handle in handles.values():
            handle.close()
    with calls_path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=call_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_calls)
    summary = {
        "schema_version": "1.0", "status": "pass",
        "method": "NC single-exon residual TBLASTN frame/stop/edge classification",
        "counts": {"input_candidates": len(rows), **dict(sorted(class_counts.items()))},
        "checks": {"all_candidates_classified": sum(class_counts.values()) == len(rows)},
    }
    summary_path = Path(f"{output_prefix}.single_exon_recovery_classification.summary.json")
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return summary


def _best_final_hits(path: Path, expected: set[str]) -> dict[str, dict[str, str]]:
    best: dict[str, dict[str, str]] = {}
    with path.open("rt", encoding="utf-8", newline="") as handle:
        for line_number, raw in enumerate(handle, start=1):
            if not raw.strip():
                continue
            fields = raw.rstrip("\r\n").split("\t", maxsplit=len(FINAL_BLAST_FIELDS) - 1)
            if len(fields) != len(FINAL_BLAST_FIELDS):
                raise ValueError(f"{path}: malformed BLASTX row {line_number}")
            row = dict(zip(FINAL_BLAST_FIELDS, fields, strict=True))
            query = row["qseqid"]
            if query not in expected:
                raise ValueError(f"{path}: unknown recovered query {query!r}")
            rank = (float(row["evalue"]), -float(row["bitscore"]), -int(row["length"]), row["sseqid"])
            previous = best.get(query)
            if previous is None:
                best[query] = row
            else:
                previous_rank = (
                    float(previous["evalue"]), -float(previous["bitscore"]),
                    -int(previous["length"]), previous["sseqid"],
                )
                if rank < previous_rank:
                    best[query] = row
    return best


def _is_or_hit(hit: dict[str, str] | None) -> tuple[bool, str]:
    if hit is None:
        return False, "no_hit"
    annotation = f"{hit['sseqid']} {hit['stitle']}".lower()
    if "growth-factor-receptor" in annotation:
        return False, "growth_factor_receptor"
    if "or-receptor" in annotation:
        return True, "or_receptor_best_hit"
    return False, "non_or_best_hit"


def _overlap_components(rows: list[dict[str, str]]) -> list[list[dict[str, str]]]:
    ordered = sorted(rows, key=lambda row: (
        row["sseqid"], int(row["true_start_1"]), int(row["true_end_1"]), row["candidate_id"],
    ))
    components: list[list[dict[str, str]]] = []
    current: list[dict[str, str]] = []
    current_end = -1
    for row in ordered:
        start, end = int(row["true_start_1"]), int(row["true_end_1"])
        if current and (row["sseqid"] != current[0]["sseqid"] or start > current_end):
            components.append(current)
            current, current_end = [], -1
        current.append(row)
        current_end = max(current_end, end)
    if current:
        components.append(current)
    return components


def build_complete_catalog(
    species: str,
    existing_catalog: Path,
    existing_gff3: Path,
    existing_functional_nucleotide: Path,
    existing_functional_protein: Path,
    existing_pseudogene: Path,
    existing_truncated: Path,
    existing_ambiguous: Path,
    recovery_calls: Path,
    recovered_pseudogene: Path,
    recovered_truncated: Path,
    recovered_ambiguous: Path,
    recovery_blastx: Path,
    output_prefix: Path,
) -> dict[str, Any]:
    with existing_catalog.open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("existing catalog has no header")
        fields = list(reader.fieldnames)
        existing_rows = list(reader)
    calls = list(csv.DictReader(recovery_calls.open("rt", encoding="utf-8", newline=""), delimiter="\t"))
    eligible = [row for row in calls if row["classification"] in {"pseudogene", "truncated", "ambiguous"}]
    expected = {row["candidate_id"] for row in eligible}
    best_hits = _best_final_hits(recovery_blastx, expected)
    retained = []
    audit_rows: list[dict[str, str]] = []
    for row in eligible:
        hit = best_hits.get(row["candidate_id"])
        keep, reason = _is_or_hit(hit)
        audit_rows.append({
            "candidate_id": row["candidate_id"], "classification": row["classification"],
            "decision": "retain" if keep else "exclude", "reason": reason,
            "sseqid": hit["sseqid"] if hit else "", "stitle": hit["stitle"] if hit else "",
            "evalue": hit["evalue"] if hit else "", "bitscore": hit["bitscore"] if hit else "",
            "redundancy_decision": "pending" if keep else "not_applicable",
        })
        if keep:
            retained.append(row)
    audit_by_id = {row["candidate_id"]: row for row in audit_rows}

    selected: list[dict[str, str]] = []
    removed_redundant = 0
    for component in _overlap_components(retained):
        winner = min(component, key=lambda row: (
            -(int(row["true_end_1"]) - int(row["true_start_1"]) + 1),
            -float(best_hits[row["candidate_id"]]["bitscore"]), row["candidate_id"],
        ))
        selected.append(winner)
        for row in component:
            audit_by_id[row["candidate_id"]]["redundancy_decision"] = (
                "selected" if row is winner else f"removed_overlap_keep_{winner['candidate_id']}"
            )
            removed_redundant += row is not winner

    existing_intervals: dict[str, list[tuple[int, int, str]]] = {}
    for row in existing_rows:
        existing_intervals.setdefault(row["sseqid"], []).append(
            (int(row["start_1"]), int(row["end_1"]), row["catalog_id"])
        )
    for row in selected:
        for start, end, identifier in existing_intervals.get(row["sseqid"], []):
            if start <= int(row["true_end_1"]) and int(row["true_start_1"]) <= end:
                raise ValueError(
                    f"recovered candidate {row['candidate_id']} overlaps existing locus {identifier}"
                )

    existing_fastas = {
        "functional_nucleotide": _fasta_map(existing_functional_nucleotide),
        "functional_protein": _fasta_map(existing_functional_protein),
        "pseudogene": _fasta_map(existing_pseudogene),
        "truncated": _fasta_map(existing_truncated),
        "ambiguous": _fasta_map(existing_ambiguous),
    }
    existing_class_ids = {
        classification: {
            row["catalog_id"] for row in existing_rows if row["classification"] == classification
        }
        for classification in ("functional", "pseudogene", "truncated", "ambiguous")
    }
    if set(existing_fastas["functional_nucleotide"]) != existing_class_ids["functional"]:
        raise ValueError("existing functional nucleotide FASTA and catalog IDs differ")
    if set(existing_fastas["functional_protein"]) != existing_class_ids["functional"]:
        raise ValueError("existing functional protein FASTA and catalog IDs differ")
    for classification in ("pseudogene", "truncated", "ambiguous"):
        if set(existing_fastas[classification]) != existing_class_ids[classification]:
            raise ValueError(f"existing {classification} FASTA and catalog IDs differ")
    recovered_fastas = {
        "pseudogene": _fasta_map(recovered_pseudogene),
        "truncated": _fasta_map(recovered_truncated),
        "ambiguous": _fasta_map(recovered_ambiguous),
    }
    for classification in ("pseudogene", "truncated", "ambiguous"):
        expected_recovered = {
            row["candidate_id"] for row in eligible if row["classification"] == classification
        }
        if set(recovered_fastas[classification]) != expected_recovered:
            raise ValueError(f"recovered {classification} FASTA and calls table IDs differ")
    counters = Counter(row["classification"] for row in existing_rows)
    next_numbers: dict[str, int] = {}
    for classification, prefix in CLASS_PREFIXES.items():
        pattern = re.compile(rf"^{re.escape(species)}_{prefix}(\d+)$")
        values = [int(match.group(1)) for row in existing_rows if (match := pattern.match(row["catalog_id"]))]
        next_numbers[classification] = max(values, default=0)

    new_rows: list[dict[str, str]] = []
    new_id_by_candidate: dict[str, str] = {}
    for row in sorted(selected, key=lambda value: (
        value["sseqid"], int(value["true_start_1"]), int(value["true_end_1"]), value["candidate_id"],
    )):
        classification = row["classification"]
        next_numbers[classification] += 1
        catalog_id = f"{species}_{CLASS_PREFIXES[classification]}{next_numbers[classification]:06d}"
        new_id_by_candidate[row["candidate_id"]] = catalog_id
        hit = best_hits[row["candidate_id"]]
        new_row = {field: "" for field in fields}
        updates = {
            "catalog_id": catalog_id,
            "source_id": row["candidate_id"],
            "source_stage": "single_exon_residual_recovery",
            "classification": classification,
            "evidence": row["evidence"] + "; final best hit OR-Receptor",
            "query_id": row["query_id"],
            "sseqid": row["sseqid"],
            "start_1": row["true_start_1"],
            "end_1": row["true_end_1"],
            "strand": row["strand"],
            "genomic_span_bp": str(int(row["true_end_1"]) - int(row["true_start_1"]) + 1),
            "spliced_cds_bp": str(int(row["true_end_1"]) - int(row["true_start_1"]) + 1),
            "protein_aa": "", "exon_count": "1", "requires_7tm_gate": "false",
            "pre7tm_classification": classification,
            "interpro_7tm_status": "not_applicable",
            "final_homology_decision": "retain",
            "final_homology_reason": "or_receptor_best_hit",
            "final_homology_program": "blastx",
            "final_homology_sseqid": hit["sseqid"],
            "final_homology_stitle": hit["stitle"],
            "final_homology_pident": hit["pident"],
            "final_homology_alignment_length": hit["length"],
            "final_homology_evalue": hit["evalue"],
            "final_homology_bitscore": hit["bitscore"],
            "final_homology_query_coverage": f"{(abs(int(hit['qend']) - int(hit['qstart'])) + 1) / int(hit['qlen']):.6f}",
        }
        for key, value in updates.items():
            if key in new_row:
                new_row[key] = value
        new_rows.append(new_row)
        counters[classification] += 1

    output_paths = {
        "catalog": Path(f"{output_prefix}.complete_OR_catalog.tsv"),
        "gff3": Path(f"{output_prefix}.complete_OR_catalog.gff3"),
        "audit": Path(f"{output_prefix}.single_exon_recovery_audit.tsv"),
        "functional_nucleotide": Path(f"{output_prefix}.complete_functional_OR.nucl.fa"),
        "functional_protein": Path(f"{output_prefix}.complete_functional_OR.prot.fa"),
        "pseudogene": Path(f"{output_prefix}.complete_pseudogene_OR.nucl.fa"),
        "truncated": Path(f"{output_prefix}.complete_truncated_OR.nucl.fa"),
        "ambiguous": Path(f"{output_prefix}.complete_ambiguous_OR.nucl.fa"),
        "rejected": Path(f"{output_prefix}.single_exon_nonOR_rejected.nucl.fa"),
    }
    with output_paths["catalog"].open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(existing_rows)
        writer.writerows(new_rows)
    audit_fields = (
        "candidate_id", "classification", "decision", "reason", "sseqid", "stitle",
        "evalue", "bitscore", "redundancy_decision", "assigned_catalog_id",
    )
    with output_paths["audit"].open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=audit_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in audit_rows:
            writer.writerow({**row, "assigned_catalog_id": new_id_by_candidate.get(row["candidate_id"], "")})

    with existing_gff3.open("rt", encoding="utf-8") as source, output_paths["gff3"].open(
        "wt", encoding="utf-8", newline="\n"
    ) as output:
        existing_text = source.read()
        existing_gff_genes = sum(
            len(line.split("\t")) == 9 and line.split("\t")[2] == "gene"
            for line in existing_text.splitlines() if line and not line.startswith("#")
        )
        if existing_gff_genes != len(existing_rows):
            raise ValueError("existing GFF3 gene count and catalog row count differ")
        output.write(existing_text)
        if existing_text and not existing_text.endswith("\n"):
            output.write("\n")
        for row in selected:
            catalog_id = new_id_by_candidate[row["candidate_id"]]
            attrs = (
                f"ID={catalog_id};source_id={row['candidate_id']};"
                f"classification={row['classification']};source_stage=single_exon_residual_recovery;"
                "final_homology=pass"
            )
            output.write(
                f"{row['sseqid']}\tfelidae_or\tgene\t{row['true_start_1']}\t{row['true_end_1']}\t.\t"
                f"{row['strand']}\t.\t{attrs}\n"
            )
            output.write(
                f"{row['sseqid']}\tfelidae_or\texon\t{row['true_start_1']}\t{row['true_end_1']}\t.\t"
                f"{row['strand']}\t.\tID={catalog_id}.exon1;Parent={catalog_id}\n"
            )

    output_handles = {
        key: output_paths[key].open("wt", encoding="utf-8", newline="\n")
        for key in ("functional_nucleotide", "functional_protein", "pseudogene", "truncated", "ambiguous", "rejected")
    }
    try:
        for _, (header, sequence) in existing_fastas["functional_nucleotide"].items():
            write_fasta_record(output_handles["functional_nucleotide"], header, sequence)
        for _, (header, sequence) in existing_fastas["functional_protein"].items():
            write_fasta_record(output_handles["functional_protein"], header, sequence)
        for classification in ("pseudogene", "truncated", "ambiguous"):
            for _, (header, sequence) in existing_fastas[classification].items():
                write_fasta_record(output_handles[classification], header, sequence)
        for row in selected:
            classification = row["classification"]
            candidate_id = row["candidate_id"]
            sequence = recovered_fastas[classification][candidate_id][1]
            catalog_id = new_id_by_candidate[candidate_id]
            header = (
                f"{catalog_id} source_id={candidate_id} {row['sseqid']}:"
                f"{row['true_start_1']}-{row['true_end_1']} strand={row['strand']} "
                f"class={classification} exons=1"
            )
            write_fasta_record(output_handles[classification], header, sequence)
        rejected_ids = {
            row["candidate_id"] for row in eligible if row["candidate_id"] not in new_id_by_candidate
        }
        for row in eligible:
            if row["candidate_id"] in rejected_ids:
                sequence = recovered_fastas[row["classification"]][row["candidate_id"]][1]
                write_fasta_record(output_handles["rejected"], row["candidate_id"], sequence)
    finally:
        for handle in output_handles.values():
            handle.close()

    final_loci = len(existing_rows) + len(new_rows)
    summary = {
        "schema_version": "1.0", "status": "pass", "species": species,
        "stage": "complete NC-compatible OR catalog including single-exon nonfunctional recovery",
        "reportable_final_catalog": True,
        "counts": {
            "existing_catalog_loci": len(existing_rows),
            "recovery_candidates_with_targeted_alignment": len(eligible),
            "recovery_candidates_with_or_best_hit": len(retained),
            "recovery_redundant_overlaps_removed": removed_redundant,
            "new_catalog_loci": len(new_rows),
            "complete_catalog_loci": final_loci,
            "functional": counters["functional"],
            "pseudogene": counters["pseudogene"],
            "truncated": counters["truncated"],
            "ambiguous": counters["ambiguous"],
            "single_exon_non_or_rejected": len(eligible) - len(retained),
            "single_exon_total_rejected_or_redundant": len(eligible) - len(new_rows),
        },
        "checks": {
            "existing_catalog_preserved": True,
            "new_loci_do_not_overlap_existing_catalog": True,
            "new_overlap_components_resolved": True,
            "catalog_partition_is_complete": sum(counters.values()) == final_loci,
        },
        "outputs": {key: path.name for key, path in output_paths.items()},
    }
    summary_path = Path(f"{output_prefix}.complete_OR_catalog.summary.json")
    summary["outputs"]["summary"] = summary_path.name
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return summary
