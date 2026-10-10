"""Classify strictly refined multi-exon OR models with explicit coordinates."""

from __future__ import annotations

import csv
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

from or_mining.exonerate_model_selection import (
    ExonerateModel,
    Region,
    _read_regions,
    parse_exonerate_models,
)
from or_mining.fasta import read_fasta, write_fasta_record


DNA_COMPLEMENT = str.maketrans(
    "ACGTRYKMSWBDHVNacgtrykmswbdhvn",
    "TGCAYRMKSWVHDBNtgcayrmkswvhdbn",
)
STOP_CODONS = {"TAA", "TAG", "TGA"}
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


@dataclass(frozen=True, slots=True)
class StrictCandidate:
    candidate_id: str
    region_id: str
    query_id: str
    sseqid: str
    strand: str
    refinement_start_1: int
    refinement_end_1: int
    target_sequence: str
    contig_length: int


@dataclass(frozen=True, slots=True)
class IntactOrf:
    start_0: int
    end_0: int
    nucleotide: str
    protein: str
    genomic_positions: tuple[int, ...]
    aligned_overlap_bp: int
    exon_intervals: tuple[tuple[int, int], ...]


def reverse_complement(sequence: str) -> str:
    return sequence.translate(DNA_COMPLEMENT)[::-1]


def translate(sequence: str) -> str:
    sequence = sequence.upper()
    return "".join(
        CODON_TABLE.get(sequence[offset : offset + 3], "X")
        for offset in range(0, len(sequence) - 2, 3)
    )


def _fasta_map(paths: Iterable[Path]) -> dict[str, str]:
    records: dict[str, str] = {}
    for path in paths:
        for header, sequence in read_fasta(path):
            identifier = header.split()[0]
            if identifier in records:
                raise ValueError(f"duplicate FASTA identifier: {identifier}")
            if not sequence:
                raise ValueError(f"empty FASTA sequence: {identifier}")
            records[identifier] = sequence.upper()
    return records


def _load_strict_candidates(
    manifest_path: Path,
    source_regions_path: Path,
    source_fasta_path: Path,
) -> tuple[dict[str, StrictCandidate], dict[str, Region], dict[str, set[str]]]:
    source_regions = _read_regions(source_regions_path)
    source_sequences = _fasta_map([source_fasta_path])
    if set(source_sequences) != set(source_regions):
        raise ValueError("multi-exon region FASTA and table IDs differ")

    candidates: dict[str, StrictCandidate] = {}
    parser_regions: dict[str, Region] = {}
    allowed_queries: dict[str, set[str]] = {}
    required = {
        "candidate_id",
        "region_id",
        "query_id",
        "sseqid",
        "strand",
        "refinement_start_1",
        "refinement_end_1",
        "target_bases",
    }
    with manifest_path.open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"strict manifest is missing columns: {sorted(required)}")
        for line_number, row in enumerate(reader, start=2):
            candidate_id = row["candidate_id"]
            region = source_regions.get(row["region_id"])
            if not candidate_id or candidate_id in candidates or region is None:
                raise ValueError(f"strict manifest line {line_number}: invalid candidate/region")
            if row["sseqid"] != region.sseqid or row["strand"] not in {"+", "-"}:
                raise ValueError(f"strict manifest line {line_number}: scaffold/strand mismatch")
            start_1 = int(row["refinement_start_1"])
            end_1 = int(row["refinement_end_1"])
            if start_1 < region.expanded_start_1 or end_1 > region.expanded_end_1:
                raise ValueError(f"strict manifest line {line_number}: target outside source region")
            local_start_0 = start_1 - region.expanded_start_1
            local_end_0 = end_1 - region.expanded_start_1 + 1
            target = source_sequences[region.region_id][local_start_0:local_end_0]
            if len(target) != int(row["target_bases"]):
                raise ValueError(f"strict manifest line {line_number}: target length mismatch")
            candidate = StrictCandidate(
                candidate_id=candidate_id,
                region_id=region.region_id,
                query_id=row["query_id"],
                sseqid=region.sseqid,
                strand=row["strand"],
                refinement_start_1=start_1,
                refinement_end_1=end_1,
                target_sequence=target,
                contig_length=region.slen,
            )
            candidates[candidate_id] = candidate
            parser_regions[candidate_id] = Region(
                candidate_id, region.sseqid, start_1, end_1, region.slen
            )
            allowed_queries[candidate_id] = {row["query_id"]}
    if not candidates:
        raise ValueError("strict manifest is empty")
    return candidates, parser_regions, allowed_queries


def _transcript_path(
    model: ExonerateModel,
    candidate: StrictCandidate,
    include_outer_flanks: bool,
) -> tuple[str, tuple[int, ...], set[int]]:
    exons = sorted(model.exons, key=lambda feature: feature.start_1)
    target = candidate.target_sequence
    local_positions: list[int] = []
    if model.strand == "+":
        if include_outer_flanks:
            local_positions.extend(range(1, exons[0].start_1))
        for exon in exons:
            local_positions.extend(range(exon.start_1, exon.end_1 + 1))
        if include_outer_flanks:
            local_positions.extend(range(exons[-1].end_1 + 1, len(target) + 1))
        sequence = "".join(target[position - 1] for position in local_positions)
    else:
        if include_outer_flanks:
            local_positions.extend(range(len(target), exons[-1].end_1, -1))
        for exon in reversed(exons):
            local_positions.extend(range(exon.end_1, exon.start_1 - 1, -1))
        if include_outer_flanks:
            local_positions.extend(range(exons[0].start_1 - 1, 0, -1))
        forward = "".join(target[position - 1] for position in local_positions)
        sequence = forward.translate(DNA_COMPLEMENT)

    genomic_positions = tuple(
        candidate.refinement_start_1 + position - 1 for position in local_positions
    )
    aligned_positions = {
        candidate.refinement_start_1 + position - 1
        for exon in exons
        for position in range(exon.start_1, exon.end_1 + 1)
    }
    return sequence.upper(), genomic_positions, aligned_positions


def _position_intervals(positions: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    if not positions:
        return ()
    groups: list[list[int]] = [[positions[0]]]
    for position in positions[1:]:
        if abs(position - groups[-1][-1]) == 1:
            groups[-1].append(position)
        else:
            groups.append([position])
    return tuple((min(group), max(group)) for group in groups)


def _find_intact_orf(
    sequence: str,
    genomic_positions: tuple[int, ...],
    aligned_positions: set[int],
    query_length_aa: int,
    minimum_orf_nt: int,
    minimum_model_overlap: float,
    maximum_exons: int,
) -> IntactOrf | None:
    candidates: list[IntactOrf] = []
    required_overlap = math.ceil(len(aligned_positions) * minimum_model_overlap)
    for frame in range(3):
        stops = [
            offset
            for offset in range(frame, len(sequence) - 2, 3)
            if sequence[offset : offset + 3] in STOP_CODONS
        ]
        for start_0 in range(frame, len(sequence) - 2, 3):
            if sequence[start_0 : start_0 + 3] != "ATG":
                continue
            stop_0 = next((offset for offset in stops if offset > start_0), None)
            if stop_0 is None:
                continue
            end_0 = stop_0 + 3
            nucleotide = sequence[start_0:end_0]
            if len(nucleotide) < minimum_orf_nt or set(nucleotide) - set("ACGT"):
                continue
            mapped = genomic_positions[start_0:end_0]
            overlap = sum(position in aligned_positions for position in mapped)
            intervals = _position_intervals(mapped)
            if overlap < required_overlap or not (2 <= len(intervals) <= maximum_exons):
                continue
            protein_with_stop = translate(nucleotide)
            if not protein_with_stop.endswith("*") or "*" in protein_with_stop[:-1]:
                continue
            candidates.append(
                IntactOrf(
                    start_0,
                    end_0,
                    nucleotide,
                    protein_with_stop[:-1],
                    mapped,
                    overlap,
                    intervals,
                )
            )
    if not candidates:
        return None
    return max(
        candidates,
        key=lambda orf: (
            orf.aligned_overlap_bp,
            -abs(len(orf.protein) - query_length_aa),
            len(orf.nucleotide),
            -orf.start_0,
        ),
    )


def _longest_n_run(sequence: str) -> int:
    return max((len(match.group(0)) for match in re.finditer(r"N+", sequence.upper())), default=0)


def _gff_escape(value: object) -> str:
    return str(value).replace("%", "%25").replace(";", "%3B").replace("=", "%3D")


CALL_FIELDS = (
    "candidate_id", "query_id", "sseqid", "strand", "classification", "evidence",
    "genomic_start_1", "genomic_end_1", "query_start_0", "query_end_0",
    "query_length_aa", "query_coverage", "score", "identity_percent",
    "similarity_percent", "model_exons", "model_introns", "model_spliced_bp",
    "frameshifts", "frameshift_bases", "premature_stops", "edge_evidence",
    "longest_N_run", "intact_orf", "orf_cds_bp", "orf_protein_aa", "orf_exons",
)


def classify_multiexon_models(
    strict_raw: Path,
    strict_manifest: Path,
    source_regions_tsv: Path,
    source_regions_fasta: Path,
    reference_or: Path,
    species_or: Path,
    output_prefix: Path,
    minimum_orf_nt: int = 750,
    minimum_model_overlap: float = 0.80,
    maximum_exons: int = 2,
    edge_distance_bp: int = 100,
    n_window_bp: int = 200,
    edge_n_run_bp: int = 50,
) -> dict[str, Any]:
    if minimum_orf_nt < 3 or minimum_orf_nt % 3:
        raise ValueError("minimum_orf_nt must be a positive multiple of three")
    if not 0 < minimum_model_overlap <= 1:
        raise ValueError("minimum_model_overlap must be in (0, 1]")
    if maximum_exons < 2:
        raise ValueError("maximum_exons must be at least two")

    candidates, parser_regions, allowed_queries = _load_strict_candidates(
        strict_manifest, source_regions_tsv, source_regions_fasta
    )
    query_sequences = _fasta_map([reference_or, species_or])
    models = parse_exonerate_models(strict_raw, parser_regions, allowed_queries)
    by_candidate = {model.target_id: model for model in models}
    if len(by_candidate) != len(models) or set(by_candidate) != set(candidates):
        raise ValueError("strict Exonerate models are not one-to-one with the manifest")

    calls_path = output_prefix.with_suffix(".multiexon_calls.tsv")
    gff_path = output_prefix.with_suffix(".multiexon_calls.gff3")
    summary_path = output_prefix.with_suffix(".multiexon_classification.summary.json")
    fasta_paths = {
        "functional": output_prefix.with_suffix(".multiexon_functional.nucl.fa"),
        "functional_protein": output_prefix.with_suffix(".multiexon_functional.prot.fa"),
        "pseudogene": output_prefix.with_suffix(".multiexon_pseudogene.nucl.fa"),
        "truncated": output_prefix.with_suffix(".multiexon_truncated.nucl.fa"),
        "ambiguous": output_prefix.with_suffix(".multiexon_ambiguous.nucl.fa"),
    }
    fasta_handles = {
        key: path.open("wt", encoding="utf-8", newline="\n")
        for key, path in fasta_paths.items()
    }
    counts = {"functional": 0, "pseudogene": 0, "truncated": 0, "ambiguous": 0}
    total_exons = 0
    try:
        with calls_path.open("wt", encoding="utf-8", newline="") as calls_handle, gff_path.open(
            "wt", encoding="utf-8", newline="\n"
        ) as gff_handle:
            writer = csv.DictWriter(calls_handle, fieldnames=CALL_FIELDS, delimiter="\t")
            writer.writeheader()
            gff_handle.write("##gff-version 3\n")
            gff_handle.write("##coordinate-system 1-based inclusive\n")

            for candidate_id in sorted(candidates):
                candidate = candidates[candidate_id]
                model = by_candidate[candidate_id]
                query = query_sequences.get(model.query_id)
                if query is None:
                    raise ValueError(f"missing strict query sequence: {model.query_id}")
                if model.strand != candidate.strand:
                    raise ValueError(f"{candidate_id}: strict and screened strands disagree")

                transcript, mapping, aligned_positions = _transcript_path(
                    model, candidate, include_outer_flanks=True
                )
                model_sequence, _, _ = _transcript_path(
                    model, candidate, include_outer_flanks=False
                )
                query_span = model.query_span_aa
                query_coverage = query_span / len(query)
                orf = _find_intact_orf(
                    transcript,
                    mapping,
                    aligned_positions,
                    len(query),
                    minimum_orf_nt,
                    minimum_model_overlap,
                    maximum_exons,
                )
                translated_tcs = translate(model.target_coding_sequence)
                stop_limit = max(
                    0,
                    math.floor(len(query) * 0.95)
                    - min(model.query_start_0, model.query_end_0),
                )
                premature_stops = sum(
                    residue == "*" and position <= stop_limit
                    for position, residue in enumerate(translated_tcs, start=1)
                )
                genomic_model_start = model.genomic_start_1
                genomic_model_end = model.genomic_end_1
                local_neighborhood_start = max(0, model.gene.start_1 - n_window_bp - 1)
                local_neighborhood_end = min(
                    len(candidate.target_sequence), model.gene.end_1 + n_window_bp
                )
                longest_n = _longest_n_run(
                    candidate.target_sequence[local_neighborhood_start:local_neighborhood_end]
                )
                edge_evidence = (
                    genomic_model_start <= edge_distance_bp
                    or candidate.contig_length - genomic_model_end < edge_distance_bp
                    or longest_n >= edge_n_run_bp
                )
                frame_disruption = model.frameshifts > 0 or len(model_sequence) % 3 != 0

                if orf is not None:
                    classification = "functional"
                    evidence = "intact_start_stop_orf"
                    output_sequence = orf.nucleotide
                    output_intervals = orf.exon_intervals
                    output_start = min(orf.genomic_positions)
                    output_end = max(orf.genomic_positions)
                elif frame_disruption or premature_stops:
                    classification = "pseudogene"
                    evidence_parts = []
                    if frame_disruption:
                        evidence_parts.append("frameshift_or_nontriplet_cds")
                    if premature_stops:
                        evidence_parts.append("premature_stop")
                    evidence = ",".join(evidence_parts)
                    output_sequence = model_sequence
                    output_intervals = tuple(
                        (candidate.refinement_start_1 + exon.start_1 - 1,
                         candidate.refinement_start_1 + exon.end_1 - 1)
                        for exon in sorted(model.exons, key=lambda feature: feature.start_1)
                    )
                    output_start, output_end = genomic_model_start, genomic_model_end
                elif edge_evidence or query_coverage < 0.95:
                    classification = "truncated"
                    evidence = "assembly_edge_or_partial_query" if edge_evidence else "partial_query"
                    output_sequence = model_sequence
                    output_intervals = tuple(
                        (candidate.refinement_start_1 + exon.start_1 - 1,
                         candidate.refinement_start_1 + exon.end_1 - 1)
                        for exon in sorted(model.exons, key=lambda feature: feature.start_1)
                    )
                    output_start, output_end = genomic_model_start, genomic_model_end
                else:
                    classification = "ambiguous"
                    evidence = "no_intact_orf_without_explicit_disruption"
                    output_sequence = model_sequence
                    output_intervals = tuple(
                        (candidate.refinement_start_1 + exon.start_1 - 1,
                         candidate.refinement_start_1 + exon.end_1 - 1)
                        for exon in sorted(model.exons, key=lambda feature: feature.start_1)
                    )
                    output_start, output_end = genomic_model_start, genomic_model_end

                counts[classification] += 1
                total_exons += len(output_intervals)
                header = (
                    f"{candidate_id} {candidate.sseqid}:{output_start}-{output_end} "
                    f"strand={model.strand} class={classification} exons={len(output_intervals)}"
                )
                write_fasta_record(fasta_handles[classification], header, output_sequence)
                if orf is not None:
                    write_fasta_record(fasta_handles["functional_protein"], header, orf.protein)

                writer.writerow(
                    {
                        "candidate_id": candidate_id,
                        "query_id": model.query_id,
                        "sseqid": candidate.sseqid,
                        "strand": model.strand,
                        "classification": classification,
                        "evidence": evidence,
                        "genomic_start_1": output_start,
                        "genomic_end_1": output_end,
                        "query_start_0": model.query_start_0,
                        "query_end_0": model.query_end_0,
                        "query_length_aa": len(query),
                        "query_coverage": f"{query_coverage:.6f}",
                        "score": f"{model.score:.6g}",
                        "identity_percent": "" if model.identity is None else f"{model.identity:.6g}",
                        "similarity_percent": "" if model.similarity is None else f"{model.similarity:.6g}",
                        "model_exons": model.exon_count,
                        "model_introns": model.intron_count,
                        "model_spliced_bp": len(model_sequence),
                        "frameshifts": model.frameshifts,
                        "frameshift_bases": model.frameshift_bases,
                        "premature_stops": premature_stops,
                        "edge_evidence": str(edge_evidence).lower(),
                        "longest_N_run": longest_n,
                        "intact_orf": str(orf is not None).lower(),
                        "orf_cds_bp": 0 if orf is None else len(orf.nucleotide),
                        "orf_protein_aa": 0 if orf is None else len(orf.protein),
                        "orf_exons": 0 if orf is None else len(orf.exon_intervals),
                    }
                )

                gene_attributes = (
                    f"ID={_gff_escape(candidate_id)};classification={classification};"
                    f"query={_gff_escape(model.query_id)};evidence={_gff_escape(evidence)}"
                )
                gff_handle.write(
                    f"{candidate.sseqid}\tfelidae_or\tgene\t{output_start}\t{output_end}\t"
                    f"{model.score:.6g}\t{model.strand}\t.\t{gene_attributes}\n"
                )
                transcript_order = (
                    sorted(output_intervals)
                    if model.strand == "+"
                    else sorted(output_intervals, reverse=True)
                )
                phase_by_interval: dict[tuple[int, int], str] = {}
                coding_before = 0
                for interval in transcript_order:
                    phase_by_interval[interval] = (
                        str((3 - coding_before % 3) % 3) if orf is not None else "."
                    )
                    coding_before += interval[1] - interval[0] + 1
                for exon_number, interval in enumerate(transcript_order, start=1):
                    start_1, end_1 = interval
                    exon_id = f"{candidate_id}.exon{exon_number}"
                    gff_handle.write(
                        f"{candidate.sseqid}\tfelidae_or\texon\t{start_1}\t{end_1}\t.\t"
                        f"{model.strand}\t.\tID={exon_id};Parent={candidate_id}\n"
                    )
                    gff_handle.write(
                        f"{candidate.sseqid}\tfelidae_or\tCDS\t{start_1}\t{end_1}\t.\t"
                        f"{model.strand}\t{phase_by_interval[interval]}\tParent={candidate_id}\n"
                    )
    finally:
        for handle in fasta_handles.values():
            handle.close()

    summary: dict[str, Any] = {
        "schema_version": "1.0",
        "status": "pass",
        "stage": "strict multi-exon OR sequence reconstruction and classification",
        "method": "NC-compatible intact-ORF and defect classification with explicit coordinate audit",
        "coordinates": {
            "strict_exonerate_gff": "target-local 1-based inclusive",
            "output_tsv_and_gff3": "genome-global 1-based inclusive",
            "query_ranges": "0-based half-open",
        },
        "parameters": {
            "minimum_orf_nt": minimum_orf_nt,
            "minimum_model_overlap": minimum_model_overlap,
            "maximum_exons": maximum_exons,
            "edge_distance_bp": edge_distance_bp,
            "edge_n_run_bp": edge_n_run_bp,
            "edge_n_window_bp": n_window_bp,
            "functional_rule": (
                "ATG-to-stop ORF, no internal stop/ambiguity, >= minimum length, "
                "spanning 2..maximum_exons"
            ),
            "pseudogene_rule": (
                "frameshift/non-triplet spliced model or premature stop before "
                "95% query position"
            ),
            "truncated_rule": (
                "no intact ORF/explicit disruption plus assembly-edge evidence "
                "or <95% query coverage"
            ),
        },
        "counts": {
            "input_candidates": len(models),
            "classified_candidates": sum(counts.values()),
            "functional": counts["functional"],
            "pseudogene": counts["pseudogene"],
            "truncated": counts["truncated"],
            "ambiguous": counts["ambiguous"],
            "output_exons": total_exons,
        },
        "outputs": {
            "calls": calls_path.name,
            "gff3": gff_path.name,
            **{key: path.name for key, path in fasta_paths.items()},
        },
    }
    if summary["counts"]["classified_candidates"] != len(models):
        raise ValueError("multi-exon classification did not partition every model")
    summary_path.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return summary
