"""Build TBLASTN-seeded query sets for targeted multi-exon OR modelling."""

from __future__ import annotations

import csv
import json
from bisect import bisect_right
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from or_mining.blast import _open_text, parse_blast_row
from or_mining.fasta import read_fasta, write_fasta_record


@dataclass(frozen=True, slots=True)
class TargetRegion:
    region_id: str
    sseqid: str
    start_0: int
    end_0: int


@dataclass(slots=True)
class QueryEvidence:
    coverage_mask: int = 0
    hsp_count: int = 0
    max_bitscore: float = 0.0
    best_evalue: float = float("inf")

    def add(self, qstart: int, qend: int, bitscore: float, evalue: float) -> None:
        width = qend - qstart + 1
        self.coverage_mask |= ((1 << width) - 1) << (qstart - 1)
        self.hsp_count += 1
        self.max_bitscore = max(self.max_bitscore, bitscore)
        self.best_evalue = min(self.best_evalue, evalue)

    @property
    def covered_residues(self) -> int:
        return self.coverage_mask.bit_count()


def _read_unique_fasta(paths: list[Path]) -> dict[str, tuple[str, str]]:
    records: dict[str, tuple[str, str]] = {}
    for path in paths:
        for header, sequence in read_fasta(path):
            identifier = header.split()[0]
            if identifier in records:
                raise ValueError(f"duplicate FASTA identifier: {identifier}")
            if not sequence:
                raise ValueError(f"empty FASTA sequence: {identifier}")
            records[identifier] = (header, sequence)
    if not records:
        raise ValueError("combined Exonerate query FASTA is empty")
    return records


def _read_regions(path: Path) -> list[TargetRegion]:
    with path.open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"region_id", "sseqid", "expanded_start_1", "expanded_end_1"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path}: missing required columns {sorted(required)}")
        regions: list[TargetRegion] = []
        for line_number, row in enumerate(reader, start=2):
            try:
                start_1 = int(row["expanded_start_1"])
                end_1 = int(row["expanded_end_1"])
            except (TypeError, ValueError) as exc:
                raise ValueError(f"{path}: line {line_number}: invalid coordinates") from exc
            if start_1 < 1 or end_1 < start_1:
                raise ValueError(f"{path}: line {line_number}: invalid interval")
            regions.append(TargetRegion(row["region_id"], row["sseqid"], start_1 - 1, end_1))
    if not regions:
        raise ValueError("multi-exon region table is empty")
    if len({region.region_id for region in regions}) != len(regions):
        raise ValueError("duplicate region_id in multi-exon region table")
    return regions


def prepare_targeted_pairs(
    reference_or: Path,
    species_or: Path,
    target_fasta: Path,
    regions_tsv: Path,
    merged_hits: Path,
    output_dir: Path,
    manifest_path: Path,
    summary_path: Path,
    top_queries_per_region: int = 3,
) -> dict[str, Any]:
    """Select top supported queries per region and write one bounded task per region."""

    if top_queries_per_region < 1:
        raise ValueError("top_queries_per_region must be at least 1")

    query_records = _read_unique_fasta([reference_or, species_or])
    target_records = _read_unique_fasta([target_fasta])
    regions = _read_regions(regions_tsv)
    if set(target_records) != {region.region_id for region in regions}:
        raise ValueError("target FASTA identifiers and multi-exon region IDs do not match")

    regions_by_subject: dict[str, list[TargetRegion]] = {}
    for region in regions:
        regions_by_subject.setdefault(region.sseqid, []).append(region)
    starts_by_subject: dict[str, list[int]] = {}
    for subject, subject_regions in regions_by_subject.items():
        subject_regions.sort(key=lambda region: (region.start_0, region.end_0, region.region_id))
        starts_by_subject[subject] = [region.start_0 for region in subject_regions]

    evidence: dict[str, dict[str, QueryEvidence]] = {region.region_id: {} for region in regions}
    input_hits = 0
    overlapping_hits = 0
    missing_query_hit_rows = 0
    missing_query_ids: set[str] = set()
    with _open_text(merged_hits) as handle:
        for line_number, line in enumerate(handle, start=1):
            try:
                hit = parse_blast_row(line)
            except ValueError as exc:
                raise ValueError(f"{merged_hits}: line {line_number}: {exc}") from exc
            input_hits += 1
            subject_regions = regions_by_subject.get(hit.sseqid)
            if not subject_regions:
                continue
            starts = starts_by_subject[hit.sseqid]
            index = max(0, bisect_right(starts, hit.genome_start_0) - 1)
            while index < len(subject_regions):
                region = subject_regions[index]
                if region.start_0 >= hit.genome_end_0:
                    break
                if hit.genome_start_0 < region.end_0 and region.start_0 < hit.genome_end_0:
                    overlapping_hits += 1
                    if hit.qseqid not in query_records:
                        missing_query_hit_rows += 1
                        missing_query_ids.add(hit.qseqid)
                    else:
                        query_evidence = evidence[region.region_id].setdefault(
                            hit.qseqid, QueryEvidence()
                        )
                        query_evidence.add(hit.qstart, hit.qend, hit.bitscore, hit.evalue)
                index += 1

    ranked_regions: list[tuple[TargetRegion, list[tuple[str, QueryEvidence]]]] = []
    regions_without_usable_queries: list[str] = []
    for region in regions:
        ranked = sorted(
            evidence[region.region_id].items(),
            key=lambda item: (
                -item[1].covered_residues,
                -item[1].max_bitscore,
                item[1].best_evalue,
                -item[1].hsp_count,
                item[0],
            ),
        )[:top_queries_per_region]
        if not ranked:
            regions_without_usable_queries.append(region.region_id)
        ranked_regions.append((region, ranked))
    if regions_without_usable_queries:
        preview = ", ".join(regions_without_usable_queries[:10])
        raise ValueError(f"candidate regions without usable query evidence: {preview}")

    ranked_regions.sort(
        key=lambda item: (
            -item[1][0][1].covered_residues,
            -item[1][0][1].max_bitscore,
            item[1][0][1].best_evalue,
            item[0].region_id,
        )
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    selected_query_records = 0
    manifest_fields = (
        "pair_id",
        "region_id",
        "sseqid",
        "expanded_start_1",
        "expanded_end_1",
        "selected_query_ids",
        "best_query_id",
        "best_query_covered_residues",
        "best_query_max_bitscore",
        "best_query_evalue",
        "best_query_hsp_count",
    )
    with manifest_path.open("wt", encoding="utf-8", newline="") as manifest_handle:
        writer = csv.DictWriter(manifest_handle, fieldnames=manifest_fields, delimiter="\t")
        writer.writeheader()
        for index, (region, ranked) in enumerate(ranked_regions, start=1):
            pair_id = f"pair_{index:06d}"
            query_path = output_dir / f"{pair_id}.query.prot.fa"
            target_path = output_dir / f"{pair_id}.target.nucl.fa"
            with query_path.open("wt", encoding="utf-8", newline="\n") as query_handle:
                for query_id, _ in ranked:
                    header, sequence = query_records[query_id]
                    write_fasta_record(query_handle, header, sequence)
                    selected_query_records += 1
            target_header, target_sequence = target_records[region.region_id]
            with target_path.open("wt", encoding="utf-8", newline="\n") as target_handle:
                write_fasta_record(target_handle, target_header, target_sequence)
            best_query_id, best = ranked[0]
            writer.writerow(
                {
                    "pair_id": pair_id,
                    "region_id": region.region_id,
                    "sseqid": region.sseqid,
                    "expanded_start_1": region.start_0 + 1,
                    "expanded_end_1": region.end_0,
                    "selected_query_ids": ",".join(query_id for query_id, _ in ranked),
                    "best_query_id": best_query_id,
                    "best_query_covered_residues": best.covered_residues,
                    "best_query_max_bitscore": f"{best.max_bitscore:.6g}",
                    "best_query_evalue": f"{best.best_evalue:.6g}",
                    "best_query_hsp_count": best.hsp_count,
                }
            )

    all_by_all_comparisons = len(query_records) * len(regions)
    summary: dict[str, Any] = {
        "schema_version": "1.0",
        "status": "pass",
        "method": "TBLASTN-seeded top-query selection per multi-exon candidate region",
        "ranking": "covered query residues, maximum bitscore, evalue, HSP count",
        "parameters": {"top_queries_per_region": top_queries_per_region},
        "counts": {
            "reference_and_species_query_records": len(query_records),
            "target_regions": len(regions),
            "input_tblastn_hits": input_hits,
            "overlapping_tblastn_hits": overlapping_hits,
            "missing_query_hit_rows": missing_query_hit_rows,
            "missing_query_ids": len(missing_query_ids),
            "regions_with_usable_queries": len(ranked_regions),
            "pair_tasks": len(ranked_regions),
            "selected_query_records": selected_query_records,
            "all_by_all_comparisons_avoided": all_by_all_comparisons
            - selected_query_records,
        },
        "missing_query_id_preview": sorted(missing_query_ids)[:20],
    }
    with summary_path.open("wt", encoding="utf-8", newline="\n") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
    return summary
