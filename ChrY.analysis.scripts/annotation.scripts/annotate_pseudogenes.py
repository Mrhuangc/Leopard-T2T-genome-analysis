#!/usr/bin/env python3
"""One-file FCA/PPA ChrY pseudogene annotation and RAB9A/Y review.

All project-specific query, alignment parsing, ORF damage, classification,
calibration, and quality-gate logic is included below. Requires only the
external Python package mappy and the miniprot executable for a full run.
The CLI also reads FASTA, GFF3, BED, and JSON directly with the standard library.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
from collections import Counter, defaultdict
from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path
from urllib.parse import quote, unquote

# ===== models: self-contained historical logic =====
@dataclass(frozen=True, order=True)
class Interval:
    start0: int
    end0: int

    def __post_init__(self) -> None:
        if self.start0 < 0:
            raise ValueError(f"start0 must be >= 0, got {self.start0}")
        if self.end0 <= self.start0:
            raise ValueError(
                f"end0 must be greater than start0, got {self.start0}-{self.end0}"
            )

    @property
    def length(self) -> int:
        return self.end0 - self.start0


@dataclass(frozen=True)
class GffFeature:
    seqid: str
    source: str
    feature_type: str
    start0: int
    end0: int
    score: str
    strand: str
    phase: str
    attributes: dict[str, str] = field(compare=False)
    line_number: int = 0

    def __post_init__(self) -> None:
        Interval(self.start0, self.end0)
        if self.strand not in {"+", "-", ".", "?"}:
            raise ValueError(
                f"invalid strand {self.strand!r} at source line {self.line_number}"
            )

    @property
    def interval(self) -> Interval:
        return Interval(self.start0, self.end0)

    @property
    def feature_id(self) -> str | None:
        return self.attributes.get("ID")

    @property
    def parents(self) -> tuple[str, ...]:
        raw = self.attributes.get("Parent", "")
        return tuple(parent for parent in raw.split(",") if parent)


@dataclass(frozen=True)
class Region:
    chrom: str
    start0: int
    end0: int
    region_id: str
    region_class: str
    confidence: str
    source_region_id: str | None = None
    line_number: int = 0

    def __post_init__(self) -> None:
        Interval(self.start0, self.end0)

    @property
    def interval(self) -> Interval:
        return Interval(self.start0, self.end0)

    @property
    def length(self) -> int:
        return self.end0 - self.start0


@dataclass(frozen=True)
class TranslationStatus:
    protein: str
    frame_ok: bool
    complete_start: bool
    complete_stop: bool
    internal_stop_count: int
    ambiguous_aa_count: int

    @property
    def intact(self) -> bool:
        return (
            self.frame_ok
            and self.complete_start
            and self.complete_stop
            and self.internal_stop_count == 0
            and self.ambiguous_aa_count == 0
        )


@dataclass(frozen=True)
class CdsModel:
    gene_id: str
    transcript_id: str
    family: str
    seqid: str
    start0: int
    end0: int
    strand: str
    exons: tuple[Interval, ...]
    phases: tuple[str, ...]
    cds_sequence: str
    translation: TranslationStatus
    species: str = "unknown"
    annotation_class: str = "unknown"
    attributes: dict[str, str] = field(default_factory=dict, compare=False)

    @property
    def interval(self) -> Interval:
        return Interval(self.start0, self.end0)


@dataclass(frozen=True)
class Query:
    query_id: str
    family: str
    species: str
    cds_sequence: str
    protein_sequence: str
    source_gene_id: str


@dataclass(frozen=True)
class HomologyHit:
    query_id: str
    family: str
    target: str
    start0: int
    end0: int
    strand: str
    query_length: int
    query_start0: int
    query_end0: int
    matches: int
    alignment_length: int
    mapq: int
    is_primary: bool
    cigar: tuple[tuple[int, int], ...] = ()
    cs: str = ""
    frameshift_count: int = 0
    stop_codon_count: int = 0
    evidence_source: str = "minimap2_cdna"
    reported_identity: float | None = None

    def __post_init__(self) -> None:
        Interval(self.start0, self.end0)
        if self.query_length <= 0:
            raise ValueError("query_length must be positive")
        if not 0 <= self.query_start0 <= self.query_end0 <= self.query_length:
            raise ValueError("invalid query coordinates")
        if self.alignment_length <= 0:
            raise ValueError("alignment_length must be positive")
        if self.strand not in {"+", "-"}:
            raise ValueError(f"invalid homology-hit strand {self.strand!r}")

    @property
    def interval(self) -> Interval:
        return Interval(self.start0, self.end0)

    @property
    def query_coverage(self) -> float:
        return (self.query_end0 - self.query_start0) / self.query_length

    @property
    def identity(self) -> float:
        if self.reported_identity is not None:
            return self.reported_identity
        return self.matches / self.alignment_length


@dataclass(frozen=True)
class HomologyLocus:
    family: str
    target: str
    start0: int
    end0: int
    hits: tuple[HomologyHit, ...]

    @property
    def interval(self) -> Interval:
        return Interval(self.start0, self.end0)

    @property
    def best_hit(self) -> HomologyHit:
        return max(
            self.hits,
            key=lambda hit: (
                hit.query_coverage * hit.identity,
                hit.query_coverage,
                hit.identity,
                hit.matches,
            ),
        )


@dataclass(frozen=True)
class DamageEvidence:
    projected_cds: str
    projected_exons: tuple[Interval, ...]
    translation: TranslationStatus
    disabling_lesions: tuple[str, ...]

    @property
    def intact(self) -> bool:
        return self.translation.intact and not self.disabling_lesions


@dataclass(frozen=True)
class PseudogeneCall:
    family: str
    chrom: str
    start0: int
    end0: int
    strand: str
    classification: str
    evidence_grade: str
    best_hit: HomologyHit
    damage: DamageEvidence
    coding_overlap_fraction: float
    repeat_overlap_fraction: float
    reason: str
    region_id: str = "unassigned"
    region_class: str = "unassigned"

    @property
    def interval(self) -> Interval:
        return Interval(self.start0, self.end0)


@dataclass(frozen=True)
class PseudogeneCallSet:
    species: str
    calls: tuple[PseudogeneCall, ...]

    @property
    def reliable(self) -> tuple[PseudogeneCall, ...]:
        return tuple(call for call in self.calls if call.classification == "reliable_pseudogene")


@dataclass(frozen=True)
class CalibrationMetrics:
    reference_explicit_count: int
    recovered_explicit_count: int
    sensitivity: float
    matched_reference_ids: tuple[str, ...] = ()


@dataclass(frozen=True)
class GateReport:
    gate_pass: bool
    checks: dict[str, bool]
    messages: tuple[str, ...]


# ===== intervals: self-contained historical logic =====
def merge_intervals(intervals: Iterable[tuple[int, int] | Interval]) -> list[tuple[int, int]]:
    normalized: list[tuple[int, int]] = []
    for interval in intervals:
        if isinstance(interval, Interval):
            start0, end0 = interval.start0, interval.end0
        else:
            start0, end0 = interval
            Interval(start0, end0)
        normalized.append((start0, end0))
    if not normalized:
        return []
    normalized.sort()
    merged = [normalized[0]]
    for start0, end0 in normalized[1:]:
        previous_start, previous_end = merged[-1]
        if start0 <= previous_end:
            merged[-1] = (previous_start, max(previous_end, end0))
        else:
            merged.append((start0, end0))
    return merged


def covered_bp(intervals: Iterable[tuple[int, int] | Interval]) -> int:
    return sum(end0 - start0 for start0, end0 in merge_intervals(intervals))


def interval_overlap(a: Interval, b: Interval) -> int:
    return max(0, min(a.end0, b.end0) - max(a.start0, b.start0))


def reciprocal_overlap(a: Interval, b: Interval) -> tuple[float, float]:
    overlap = interval_overlap(a, b)
    return overlap / a.length, overlap / b.length


def clip_interval(start0: int, end0: int, chrom_length: int) -> Interval:
    if chrom_length <= 0:
        raise ValueError("chrom_length must be positive")
    clipped_start = max(0, start0)
    clipped_end = min(chrom_length, end0)
    return Interval(clipped_start, clipped_end)


# ===== io: self-contained historical logic =====
def read_fasta(path: Path | str) -> dict[str, str]:
    path = Path(path)
    records: dict[str, list[str]] = {}
    current: str | None = None
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                if not current:
                    raise ValueError(f"empty FASTA identifier at line {line_number}")
                if current in records:
                    raise ValueError(f"duplicate FASTA identifier {current!r} at line {line_number}")
                records[current] = []
            else:
                if current is None:
                    raise ValueError(f"sequence before FASTA header at line {line_number}")
                records[current].append(line.upper())
    return {name: "".join(parts) for name, parts in records.items()}


def parse_attributes(text: str) -> dict[str, str]:
    if text in {"", "."}:
        return {}
    attributes: dict[str, str] = {}
    for field in text.strip().strip(";").split(";"):
        if not field:
            continue
        if "=" in field:
            key, value = field.split("=", 1)
        elif " " in field:
            key, value = field.split(" ", 1)
            value = value.strip('"')
        else:
            key, value = field, ""
        attributes[unquote(key)] = unquote(value)
    return attributes


def read_gff3(path: Path | str) -> list[GffFeature]:
    path = Path(path)
    features: list[GffFeature] = []
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(
                    f"{path}: line {line_number} has {len(fields)} columns, expected 9"
                )
            seqid, source, feature_type, start, end, score, strand, phase, attrs = fields
            try:
                start1 = int(start)
                end1 = int(end)
            except ValueError as error:
                raise ValueError(f"{path}: invalid coordinate at line {line_number}") from error
            if start1 < 1 or end1 < start1:
                raise ValueError(f"{path}: invalid interval at line {line_number}")
            features.append(
                GffFeature(
                    seqid=seqid,
                    source=source,
                    feature_type=feature_type,
                    start0=start1 - 1,
                    end0=end1,
                    score=score,
                    strand=strand,
                    phase=phase,
                    attributes=parse_attributes(attrs),
                    line_number=line_number,
                )
            )
    return features


def validate_gff_parentage(features: list[GffFeature]) -> None:
    ids = {feature.feature_id for feature in features if feature.feature_id}
    for feature in features:
        for parent in feature.parents:
            if parent not in ids:
                raise ValueError(
                    f"line {feature.line_number}: Parent {parent!r} is missing from GFF3"
                )


def read_bed_regions(path: Path | str) -> list[Region]:
    path = Path(path)
    regions: list[Region] = []
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if fields[0].lower() == "chrom" and fields[1:3] == ["start0", "end0"]:
                continue
            if len(fields) < 6:
                raise ValueError(
                    f"{path}: line {line_number} has {len(fields)} columns, expected >=6"
                )
            chrom, start, end, region_id, region_class, confidence = fields[:6]
            source_region_id = fields[6] if len(fields) > 6 else None
            try:
                start0, end0 = int(start), int(end)
            except ValueError as error:
                raise ValueError(f"{path}: invalid BED coordinate at line {line_number}") from error
            regions.append(
                Region(
                    chrom=chrom,
                    start0=start0,
                    end0=end0,
                    region_id=region_id,
                    region_class=region_class,
                    confidence=confidence,
                    source_region_id=source_region_id,
                    line_number=line_number,
                )
            )
    return regions


def validate_partition(regions: list[Region], chrom_length: int) -> None:
    if not regions:
        raise ValueError("partition is empty")
    ordered = sorted(regions, key=lambda region: (region.chrom, region.start0))
    chromosomes = {region.chrom for region in ordered}
    if len(chromosomes) != 1:
        raise ValueError(f"partition has multiple chromosomes: {sorted(chromosomes)}")
    if ordered[0].start0 != 0:
        raise ValueError(f"partition starts at {ordered[0].start0}, expected 0")
    previous_end = 0
    for region in ordered:
        if region.start0 != previous_end:
            raise ValueError(
                f"partition gap or overlap before {region.region_id}: "
                f"expected {previous_end}, got {region.start0}"
            )
        previous_end = region.end0
    if previous_end != chrom_length:
        raise ValueError(
            f"partition end {previous_end} does not match chromosome length {chrom_length}"
        )


def load_config(path: Path | str) -> dict:
    path = Path(path).resolve()
    with path.open("r", encoding="utf-8") as handle:
        config = json.load(handle)
    config["_config_path"] = str(path)
    config["_config_dir"] = str(path.parent)
    return config


def resolve_config_path(config: dict, value: str | None) -> Path | None:
    if value is None:
        return None
    return (Path(config["_config_dir"]) / value).resolve()


# ===== sequences: self-contained historical logic =====
_COMPLEMENT = str.maketrans(
    "ACGTRYSWKMBDHVNacgtryswkmbdhvn",
    "TGCAYRSWMKVHDBNtgcayrswmkvhdbn",
)
# Standard genetic code, NCBI codon order T/C/A/G on each position.
_BASES = "TCAG"
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_CODON_TABLE = dict(zip((a + b + c for a in _BASES for b in _BASES for c in _BASES), _AA))
_FORWARD_TABLE = {codon: amino for codon, amino in _CODON_TABLE.items() if amino != "*"}
_STOP_CODONS = {codon for codon, amino in _CODON_TABLE.items() if amino == "*"}


def reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


def translate_cds(sequence: str) -> TranslationStatus:
    sequence = sequence.upper()
    frame_ok = len(sequence) % 3 == 0
    codons = [sequence[index : index + 3] for index in range(0, len(sequence) - 2, 3)]
    amino_acids: list[str] = []
    for codon in codons:
        if codon in _STOP_CODONS:
            amino_acids.append("*")
        else:
            amino_acids.append(_FORWARD_TABLE.get(codon, "X"))
    protein = "".join(amino_acids)
    complete_start = bool(codons) and codons[0] == "ATG"
    complete_stop = frame_ok and bool(codons) and codons[-1] in _STOP_CODONS
    internal_stop_count = protein[:-1].count("*") if protein else 0
    return TranslationStatus(
        protein=protein,
        frame_ok=frame_ok,
        complete_start=complete_start,
        complete_stop=complete_stop,
        internal_stop_count=internal_stop_count,
        ambiguous_aa_count=protein.count("X"),
    )


def normalize_family(attributes: dict[str, str], gene_id: str) -> str:
    family = (
        attributes.get("gene_symbol")
        or attributes.get("Name")
        or attributes.get("gene")
        or gene_id
    )
    family = family.strip()
    if family.upper() == "CCDC71LY":
        return "CCDC71L"
    return family


def extract_cds_models(
    genome: str,
    features: list[GffFeature],
    species: str = "unknown",
) -> dict[str, CdsModel]:
    genes = {
        feature.feature_id: feature
        for feature in features
        if feature.feature_type in {"gene", "pseudogene"} and feature.feature_id
    }
    transcripts = {
        feature.feature_id: feature
        for feature in features
        if feature.feature_type in {"mRNA", "transcript"} and feature.feature_id
    }
    transcripts_by_gene: dict[str, list[str]] = defaultdict(list)
    for transcript_id, transcript in transcripts.items():
        for parent in transcript.parents:
            transcripts_by_gene[parent].append(transcript_id)
    cds_by_transcript: dict[str, list[GffFeature]] = defaultdict(list)
    for feature in features:
        if feature.feature_type != "CDS":
            continue
        for parent in feature.parents:
            cds_by_transcript[parent].append(feature)

    models: dict[str, CdsModel] = {}
    for gene_id, gene in genes.items():
        candidates = transcripts_by_gene.get(gene_id, [])
        if not candidates and gene_id in cds_by_transcript:
            candidates = [gene_id]
        candidates = [candidate for candidate in candidates if cds_by_transcript.get(candidate)]
        if not candidates:
            continue
        transcript_id = max(
            candidates,
            key=lambda candidate: sum(
                feature.end0 - feature.start0 for feature in cds_by_transcript[candidate]
            ),
        )
        cds_features = sorted(
            cds_by_transcript[transcript_id],
            key=lambda feature: feature.start0,
            reverse=gene.strand == "-",
        )
        sequence_parts: list[str] = []
        exons: list[Interval] = []
        phases: list[str] = []
        for cds in cds_features:
            if cds.end0 > len(genome):
                raise ValueError(
                    f"{gene_id}: CDS ends at {cds.end0}, beyond chromosome length {len(genome)}"
                )
            part = genome[cds.start0 : cds.end0]
            if gene.strand == "-":
                part = reverse_complement(part)
            sequence_parts.append(part)
            exons.append(cds.interval)
            phases.append(cds.phase)
        cds_sequence = "".join(sequence_parts)
        models[gene_id] = CdsModel(
            gene_id=gene_id,
            transcript_id=transcript_id,
            family=normalize_family(gene.attributes, gene_id),
            seqid=gene.seqid,
            start0=gene.start0,
            end0=gene.end0,
            strand=gene.strand,
            exons=tuple(exons),
            phases=tuple(phases),
            cds_sequence=cds_sequence,
            translation=translate_cds(cds_sequence),
            species=species,
            annotation_class=gene.attributes.get("annotation_class", "unknown"),
            attributes=dict(gene.attributes),
        )
    return models


def select_nonredundant_queries(
    models: Iterable[CdsModel], max_per_family: int = 4
) -> list[Query]:
    by_family: dict[str, list[CdsModel]] = defaultdict(list)
    for model in models:
        if not model.cds_sequence:
            continue
        by_family[model.family].append(model)
    queries: list[Query] = []
    for family in sorted(by_family):
        models_for_family = sorted(
            by_family[family],
            key=lambda model: (
                model.translation.intact,
                len(model.cds_sequence),
                model.species,
                model.gene_id,
            ),
            reverse=True,
        )
        seen_sequences: set[str] = set()
        retained: list[CdsModel] = []
        best_by_species: dict[str, CdsModel] = {}
        for model in models_for_family:
            best_by_species.setdefault(model.species, model)
        for model in sorted(best_by_species.values(), key=lambda item: item.species):
            if model.cds_sequence in seen_sequences:
                continue
            seen_sequences.add(model.cds_sequence)
            retained.append(model)
        for model in models_for_family:
            if len(retained) >= max_per_family:
                break
            if model.cds_sequence in seen_sequences:
                continue
            seen_sequences.add(model.cds_sequence)
            retained.append(model)
        for index, model in enumerate(retained, start=1):
            queries.append(
                Query(
                    query_id=f"{family}|{model.species}|q{index}|{model.gene_id}",
                    family=family,
                    species=model.species,
                    cds_sequence=model.cds_sequence,
                    protein_sequence=model.translation.protein.rstrip("*"),
                    source_gene_id=model.gene_id,
                )
            )
    return queries


def derive_pseudogene_queries(
    genome: str,
    features: list[GffFeature],
    *,
    species: str = "PPA_PG",
    min_family_copies: int = 2,
    max_per_family: int = 2,
) -> list[Query]:
    genes = {
        feature.feature_id: feature
        for feature in features
        if feature.feature_type == "pseudogene" and feature.feature_id
    }
    transcript_to_gene: dict[str, str] = {}
    for feature in features:
        if feature.feature_type != "pseudogenic_transcript" or not feature.feature_id:
            continue
        if feature.parents:
            transcript_to_gene[feature.feature_id] = feature.parents[0]
    exons_by_gene: dict[str, list[GffFeature]] = defaultdict(list)
    for feature in features:
        if feature.feature_type != "pseudogenic_exon":
            continue
        for parent in feature.parents:
            gene_id = transcript_to_gene.get(parent)
            if gene_id:
                exons_by_gene[gene_id].append(feature)

    sequences_by_family: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for gene_id, gene in genes.items():
        exons = sorted(
            exons_by_gene.get(gene_id, []),
            key=lambda exon: exon.start0,
            reverse=gene.strand == "-",
        )
        if not exons:
            continue
        parts = []
        for exon in exons:
            part = genome[exon.start0 : exon.end0]
            if gene.strand == "-":
                part = reverse_complement(part)
            parts.append(part)
        sequence = "".join(parts)
        family = normalize_family(gene.attributes, gene_id)
        sequences_by_family[family].append((gene_id, sequence))

    queries: list[Query] = []
    for family, family_records in sorted(sequences_by_family.items()):
        if len(family_records) < min_family_copies:
            continue
        length_bins: dict[int, list[tuple[str, str]]] = defaultdict(list)
        for record in family_records:
            length_bins[(len(record[1]) // 50) * 50].append(record)
        selected_bins = sorted(
            length_bins.items(),
            key=lambda item: (len(item[1]), max(len(sequence) for _, sequence in item[1])),
            reverse=True,
        )[:max_per_family]
        for group_index, (_, records) in enumerate(selected_bins, start=1):
            modal_length = Counter(len(sequence) for _, sequence in records).most_common(1)[0][0]
            modal_records = [record for record in records if len(record[1]) == modal_length]
            consensus = "".join(
                Counter(sequence[position] for _, sequence in modal_records).most_common(1)[0][0]
                for position in range(modal_length)
            )
            frame_candidates: list[tuple[tuple[int, int, int, int], int, str, TranslationStatus]] = []
            for frame in range(3):
                candidate = consensus[frame:]
                candidate = candidate[: len(candidate) // 3 * 3]
                if not candidate:
                    continue
                status = translate_cds(candidate)
                score = (
                    -status.internal_stop_count,
                    int(status.complete_start),
                    -status.ambiguous_aa_count,
                    len(status.protein),
                )
                frame_candidates.append((score, frame, candidate, status))
            if not frame_candidates:
                continue
            _, frame, candidate, status = max(frame_candidates)
            representative_id = modal_records[0][0]
            query_id = (
                f"{family}|{species}|pg_group{group_index}_len{modal_length}_frame{frame}"
            )
            queries.append(
                Query(
                    query_id=query_id,
                    family=family,
                    species=species,
                    cds_sequence=candidate,
                    protein_sequence=status.protein.rstrip("*").replace("*", "X"),
                    source_gene_id=representative_id,
                )
            )
    return queries


def _write_wrapped_fasta(handle, identifier: str, sequence: str, width: int = 60) -> None:
    handle.write(f">{identifier}\n")
    for index in range(0, len(sequence), width):
        handle.write(sequence[index : index + width] + "\n")


def write_queries(queries: list[Query], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    cds_path = output_dir / "query_cds.fa"
    protein_path = output_dir / "query_proteins.fa"
    manifest_path = output_dir / "query_manifest.tsv"
    with (
        cds_path.open("w", encoding="utf-8") as cds_handle,
        protein_path.open("w", encoding="utf-8") as protein_handle,
        manifest_path.open("w", encoding="utf-8", newline="") as manifest_handle,
    ):
        writer = csv.writer(manifest_handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "query_id",
                "family",
                "species",
                "source_gene_id",
                "cds_length",
                "protein_length",
            ]
        )
        for query in sorted(queries, key=lambda item: (item.family, item.species, item.query_id)):
            _write_wrapped_fasta(cds_handle, query.query_id, query.cds_sequence)
            _write_wrapped_fasta(protein_handle, query.query_id, query.protein_sequence)
            writer.writerow(
                [
                    query.query_id,
                    query.family,
                    query.species,
                    query.source_gene_id,
                    len(query.cds_sequence),
                    len(query.protein_sequence),
                ]
            )


# ===== homology: self-contained historical logic =====
def build_miniprot_command(
    *,
    executable: Path,
    target_fasta: Path,
    query_proteins: Path,
    threads: int = 4,
) -> list[str]:
    return [
        str(executable),
        "-t",
        str(threads),
        "-N",
        "500",
        "--outn",
        "500",
        "--outs",
        "0.25",
        "--outc",
        "0.20",
        "--gff",
        "--aln",
        str(target_fasta),
        str(query_proteins),
    ]


def run_miniprot(
    *,
    executable: Path,
    target_fasta: Path,
    query_proteins: Path,
    output_gff: Path,
    log_path: Path,
    threads: int = 4,
) -> None:
    command = build_miniprot_command(
        executable=executable,
        target_fasta=target_fasta,
        query_proteins=query_proteins,
        threads=threads,
    )
    output_gff.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with (
        output_gff.open("w", encoding="utf-8") as output_handle,
        log_path.open("w", encoding="utf-8") as log_handle,
    ):
        subprocess.run(
            command,
            check=True,
            stdout=output_handle,
            stderr=log_handle,
            text=True,
        )


def map_cdna_queries(
    target_fasta: Path,
    queries: list[Query],
    preset: str = "splice:hq",
    best_n: int = 500,
) -> list[HomologyHit]:
    import mappy as mp  # external minimap2 Python binding; loaded only for RAB9A/Y review

    aligner = mp.Aligner(str(target_fasta), preset=preset, best_n=best_n)
    if not aligner:
        raise RuntimeError(f"failed to build minimap2 index for {target_fasta}")
    hits: list[HomologyHit] = []
    for query in queries:
        for alignment in aligner.map(query.cds_sequence, cs=True):
            hits.append(
                HomologyHit(
                    query_id=query.query_id,
                    family=query.family,
                    target=alignment.ctg,
                    start0=alignment.r_st,
                    end0=alignment.r_en,
                    strand="+" if alignment.strand == 1 else "-",
                    query_length=len(query.cds_sequence),
                    query_start0=alignment.q_st,
                    query_end0=alignment.q_en,
                    matches=alignment.mlen,
                    alignment_length=alignment.blen,
                    mapq=alignment.mapq,
                    is_primary=bool(alignment.is_primary),
                    cigar=tuple(alignment.cigar or ()),
                    cs=alignment.cs or "",
                )
            )
    return hits


def write_hits_tsv(hits: list[HomologyHit], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "query_id",
        "family",
        "target",
        "start0",
        "end0",
        "strand",
        "query_length",
        "query_start0",
        "query_end0",
        "query_coverage",
        "matches",
        "alignment_length",
        "identity",
        "mapq",
        "is_primary",
        "cigar",
        "cs",
        "frameshift_count",
        "stop_codon_count",
        "evidence_source",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for hit in hits:
            writer.writerow(
                {
                    "query_id": hit.query_id,
                    "family": hit.family,
                    "target": hit.target,
                    "start0": hit.start0,
                    "end0": hit.end0,
                    "strand": hit.strand,
                    "query_length": hit.query_length,
                    "query_start0": hit.query_start0,
                    "query_end0": hit.query_end0,
                    "query_coverage": f"{hit.query_coverage:.6f}",
                    "matches": hit.matches,
                    "alignment_length": hit.alignment_length,
                    "identity": f"{hit.identity:.6f}",
                    "mapq": hit.mapq,
                    "is_primary": str(hit.is_primary).lower(),
                    "cigar": ",".join(f"{length}:{operation}" for length, operation in hit.cigar),
                    "cs": hit.cs,
                    "frameshift_count": hit.frameshift_count,
                    "stop_codon_count": hit.stop_codon_count,
                    "evidence_source": hit.evidence_source,
                }
            )


def _cigar_from_exons(exons: list[Interval]) -> tuple[tuple[int, int], ...]:
    if not exons:
        return ()
    ordered = sorted(exons, key=lambda exon: exon.start0)
    cigar: list[tuple[int, int]] = []
    previous_end: int | None = None
    for exon in ordered:
        if previous_end is not None:
            gap = exon.start0 - previous_end
            if gap < 0:
                raise ValueError("miniprot CDS intervals overlap within one projection")
            if gap:
                cigar.append((gap, 3))
        cigar.append((exon.length, 0))
        previous_end = exon.end0
    return tuple(cigar)


def parse_miniprot_gff(
    path: Path, query_protein_lengths: dict[str, int]
) -> list[HomologyHit]:
    features = read_gff3(path)
    cds_by_parent: dict[str, list[Interval]] = defaultdict(list)
    for feature in features:
        if feature.feature_type != "CDS":
            continue
        for parent in feature.parents:
            cds_by_parent[parent].append(feature.interval)
    hits: list[HomologyHit] = []
    for feature in features:
        if feature.feature_type != "mRNA":
            continue
        projection_id = feature.attributes.get("ID", "")
        target_fields = feature.attributes.get("Target", "").split()
        if len(target_fields) < 3:
            raise ValueError(
                f"{path}: miniprot mRNA {projection_id!r} has invalid Target attribute"
            )
        query_id = target_fields[0]
        if query_id not in query_protein_lengths:
            raise ValueError(f"{path}: no protein length for miniprot query {query_id!r}")
        query_start_aa = int(target_fields[1])
        query_end_aa = int(target_fields[2])
        if query_start_aa < 1 or query_end_aa < query_start_aa:
            raise ValueError(f"{path}: invalid miniprot Target range for {query_id}")
        aligned_nt = (query_end_aa - query_start_aa + 1) * 3
        identity = float(feature.attributes.get("Identity", "0"))
        family = query_id.split("|", 1)[0]
        exons = cds_by_parent.get(projection_id, [feature.interval])
        hits.append(
            HomologyHit(
                query_id=query_id,
                family=family,
                target=feature.seqid,
                start0=feature.start0,
                end0=feature.end0,
                strand=feature.strand,
                query_length=query_protein_lengths[query_id] * 3,
                query_start0=(query_start_aa - 1) * 3,
                query_end0=query_end_aa * 3,
                matches=round(identity * aligned_nt),
                alignment_length=aligned_nt,
                mapq=0,
                is_primary=feature.attributes.get("Rank", "1") == "1",
                cigar=_cigar_from_exons(exons),
                cs="",
                frameshift_count=int(feature.attributes.get("Frameshift", "0")),
                stop_codon_count=int(feature.attributes.get("StopCodon", "0")),
                evidence_source="miniprot_protein",
                reported_identity=identity,
            )
        )
    return hits


def cluster_family_hits(
    hits: list[HomologyHit], max_gap: int = 2000
) -> list[HomologyLocus]:
    grouped: dict[tuple[str, str], list[HomologyHit]] = defaultdict(list)
    for hit in hits:
        grouped[(hit.family, hit.target)].append(hit)
    loci: list[HomologyLocus] = []
    for (family, target), family_hits in sorted(grouped.items()):
        ordered = sorted(family_hits, key=lambda hit: (hit.start0, hit.end0))
        current: list[HomologyHit] = []
        current_start = 0
        current_end = 0
        for hit in ordered:
            if not current or hit.start0 > current_end + max_gap:
                if current:
                    loci.append(
                        HomologyLocus(
                            family=family,
                            target=target,
                            start0=current_start,
                            end0=current_end,
                            hits=tuple(current),
                        )
                    )
                current = [hit]
                current_start = hit.start0
                current_end = hit.end0
            else:
                current.append(hit)
                current_end = max(current_end, hit.end0)
        if current:
            loci.append(
                HomologyLocus(
                    family=family,
                    target=target,
                    start0=current_start,
                    end0=current_end,
                    hits=tuple(current),
                )
            )
    return sorted(loci, key=lambda locus: (locus.target, locus.start0, locus.family))


def project_target_exons(hit: HomologyHit) -> list[Interval]:
    if not hit.cigar:
        return [hit.interval]
    reference_position = hit.start0
    exon_start = reference_position
    exons: list[Interval] = []
    for length, operation in hit.cigar:
        if operation in {0, 2, 7, 8}:  # M, D, =, X consume target reference
            reference_position += length
        elif operation == 3:  # N: intron on target reference
            if reference_position > exon_start:
                exons.append(Interval(exon_start, reference_position))
            reference_position += length
            exon_start = reference_position
        elif operation in {1, 4, 5, 6}:  # I/S/H/P do not consume target
            continue
        else:
            raise ValueError(f"unsupported CIGAR operation {operation}")
    if reference_position > exon_start:
        exons.append(Interval(exon_start, reference_position))
    if not exons:
        exons = [hit.interval]
    return exons


def assess_projected_orf(hit: HomologyHit, target_sequence: str) -> DamageEvidence:
    exons = project_target_exons(hit)
    projected = "".join(target_sequence[exon.start0 : exon.end0] for exon in exons)
    if hit.strand == "-":
        projected = reverse_complement(projected)
        exons = list(reversed(exons))
    translation = translate_cds(projected)
    lesions: list[str] = []
    if not translation.frame_ok:
        lesions.append("frameshift_length")
    if translation.internal_stop_count:
        lesions.append("internal_stop")
    if not translation.complete_start:
        lesions.append("missing_start")
    if not translation.complete_stop:
        lesions.append("missing_stop")
    if translation.ambiguous_aa_count:
        lesions.append("ambiguous_sequence")
    if hit.frameshift_count:
        lesions.append("frameshift_alignment")
    if hit.stop_codon_count:
        lesions.append("inframe_stop_alignment")
    return DamageEvidence(
        projected_cds=projected,
        projected_exons=tuple(exons),
        translation=translation,
        disabling_lesions=tuple(lesions),
    )


# ===== pseudogenes: self-contained historical logic =====
DEFAULT_THRESHOLDS = {
    "min_query_coverage": 0.50,
    "min_identity": 0.70,
    "intact_min_query_coverage": 0.90,
    "intact_min_identity": 0.80,
    "accepted_coding_overlap": 0.80,
    "annotation_conflict_overlap": 0.20,
    "ambiguous_reciprocal_overlap": 0.80,
    "ambiguous_score_fraction": 0.95,
    "strong_truncation_max_query_coverage": 0.75,
    "max_repeat_fraction_for_truncation": 0.80,
}


def _as_interval(record: object) -> Interval:
    if isinstance(record, Interval):
        return record
    if isinstance(record, GffFeature):
        return record.interval
    if isinstance(record, tuple) and len(record) >= 2:
        return Interval(int(record[0]), int(record[1]))
    if hasattr(record, "interval"):
        return record.interval
    raise TypeError(f"cannot obtain interval from {type(record).__name__}")


def _overlap_fraction(interval: Interval, records: list[object]) -> float:
    overlaps: list[tuple[int, int]] = []
    for record in records:
        other = _as_interval(record)
        start0 = max(interval.start0, other.start0)
        end0 = min(interval.end0, other.end0)
        if end0 > start0:
            overlaps.append((start0, end0))
    return covered_bp(overlaps) / interval.length if overlaps else 0.0


def _coding_overlap_fraction(locus: HomologyLocus, coding_models: list[CdsModel]) -> float:
    if not coding_models:
        return 0.0
    return max(
        (
            max(reciprocal_overlap(locus.interval, model.interval))
            for model in coding_models
            if model.seqid == locus.target and interval_overlap(locus.interval, model.interval)
        ),
        default=0.0,
    )


def _assign_region(locus: HomologyLocus, regions: list[Region]) -> tuple[str, str]:
    overlaps = [
        (interval_overlap(locus.interval, region.interval), region)
        for region in regions
        if region.chrom == locus.target
    ]
    overlap, region = max(overlaps, default=(0, None), key=lambda item: item[0])
    if overlap == 0 or region is None:
        return "unassigned", "unassigned"
    return region.region_id, region.region_class


def _independent_query_sources(locus: HomologyLocus) -> set[str]:
    sources: set[str] = set()
    for hit in locus.hits:
        parts = hit.query_id.split("|")
        sources.add(parts[1] if len(parts) > 2 else hit.query_id)
    return sources


def _ambiguous_locus_indices(
    loci: list[HomologyLocus], thresholds: dict[str, float]
) -> set[int]:
    ambiguous: set[int] = set()
    for left_index, left in enumerate(loci):
        left_score = left.best_hit.query_coverage * left.best_hit.identity
        for right_index in range(left_index + 1, len(loci)):
            right = loci[right_index]
            if left.target != right.target or left.family == right.family:
                continue
            left_overlap, right_overlap = reciprocal_overlap(left.interval, right.interval)
            if min(left_overlap, right_overlap) < thresholds["ambiguous_reciprocal_overlap"]:
                continue
            right_score = right.best_hit.query_coverage * right.best_hit.identity
            if min(left_score, right_score) >= thresholds["ambiguous_score_fraction"] * max(
                left_score, right_score
            ):
                ambiguous.update({left_index, right_index})
    return ambiguous


def classify_pseudogene_loci(
    loci: list[HomologyLocus],
    coding_models: list[CdsModel],
    repeats: list[object],
    target_sequence: str,
    *,
    species: str,
    regions: list[Region] | None = None,
    thresholds: dict[str, float] | None = None,
) -> PseudogeneCallSet:
    effective = dict(DEFAULT_THRESHOLDS)
    if thresholds:
        effective.update(thresholds)
    regions = regions or []
    ambiguous_indices = _ambiguous_locus_indices(loci, effective)
    calls: list[PseudogeneCall] = []
    for index, locus in enumerate(loci):
        hit = locus.best_hit
        damage = assess_projected_orf(hit, target_sequence)
        coding_overlap = _coding_overlap_fraction(locus, coding_models)
        repeat_overlap = _overlap_fraction(locus.interval, repeats)
        region_id, region_class = _assign_region(locus, regions)
        classification = "candidate_fragment"
        grade = "low"
        reason = "homology_without_decisive_disabling_lesion"

        if index in ambiguous_indices:
            classification = "ambiguous_repetitive_hit"
            reason = "equivalent_overlapping_hits_from_multiple_families"
        elif coding_overlap >= effective["accepted_coding_overlap"]:
            classification = "accepted_coding_locus"
            grade = "high"
            reason = "overlaps_accepted_coding_model"
        elif coding_overlap >= effective["annotation_conflict_overlap"]:
            classification = "annotation_overlap_conflict"
            grade = "moderate"
            reason = "substantial_overlap_with_accepted_coding_model"
        elif (
            hit.query_coverage < effective["min_query_coverage"]
            or hit.identity < effective["min_identity"]
        ):
            classification = "subthreshold_homology"
            reason = "below_query_coverage_or_identity_threshold"
        elif (
            hit.query_coverage >= effective["intact_min_query_coverage"]
            and hit.identity >= effective["intact_min_identity"]
            and damage.intact
        ):
            classification = "unannotated_intact_candidate"
            grade = "high"
            reason = "full_length_intact_ORF_without_accepted_model"
        elif any(
            lesion
            in {
                "internal_stop",
                "frameshift_alignment",
                "inframe_stop_alignment",
            }
            for lesion in damage.disabling_lesions
        ):
            classification = "reliable_pseudogene"
            grade = "high"
            reason = "explicit_disabling_lesion"
        elif "frameshift_length" in damage.disabling_lesions:
            classification = "candidate_fragment"
            grade = "moderate"
            reason = "length_or_terminal_defect_without_localized_lesion"
        elif (
            hit.evidence_source == "miniprot_protein"
            and hit.query_coverage >= 0.90
            and hit.identity >= effective["min_identity"]
            and "missing_stop" in damage.disabling_lesions
        ):
            classification = "candidate_fragment"
            grade = "moderate"
            reason = "full_protein_projection_missing_terminal_codon"
        elif (
            hit.query_coverage <= effective["strong_truncation_max_query_coverage"]
            and len(_independent_query_sources(locus)) >= 2
            and repeat_overlap < effective["max_repeat_fraction_for_truncation"]
        ):
            classification = "reliable_pseudogene"
            grade = "moderate"
            reason = "concordant_strong_truncation"

        calls.append(
            PseudogeneCall(
                family=locus.family,
                chrom=locus.target,
                start0=hit.start0,
                end0=hit.end0,
                strand=hit.strand,
                classification=classification,
                evidence_grade=grade,
                best_hit=hit,
                damage=damage,
                coding_overlap_fraction=coding_overlap,
                repeat_overlap_fraction=repeat_overlap,
                reason=reason,
                region_id=region_id,
                region_class=region_class,
            )
        )
    return PseudogeneCallSet(
        species=species,
        calls=tuple(sorted(calls, key=lambda call: (call.chrom, call.start0, call.family))),
    )


def _escape(value: object) -> str:
    return quote(str(value), safe="._|:+-,")


def _format_attributes(items: list[tuple[str, object]]) -> str:
    return ";".join(f"{key}={_escape(value)}" for key, value in items)


def write_pseudogene_gff3(
    callset: PseudogeneCallSet, path: Path, *, chrom_length: int
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    family_counts: dict[str, int] = {}
    with path.open("w", encoding="utf-8") as handle:
        handle.write("##gff-version 3\n")
        handle.write(f"##sequence-region ChrY 1 {chrom_length}\n")
        handle.write(
            f"# Reliable {callset.species} ChrY pseudogenes reconstructed by "
            "standardized CDS homology projection.\n"
        )
        for call in callset.reliable:
            family_counts[call.family] = family_counts.get(call.family, 0) + 1
            serial = family_counts[call.family]
            call_id = f"{callset.species}.ChrY.PG.{call.family}.{serial:03d}"
            transcript_id = f"{call_id}.pt1"
            gene_attrs = _format_attributes(
                [
                    ("ID", call_id),
                    ("Name", call.family),
                    ("gene_symbol", call.family),
                    ("gene_biotype", "pseudogene"),
                    ("annotation_class", "reliable_pseudogene"),
                    ("evidence_grade", call.evidence_grade),
                    ("best_query", call.best_hit.query_id),
                    ("query_coverage", f"{call.best_hit.query_coverage:.4f}"),
                    ("identity", f"{call.best_hit.identity:.4f}"),
                    ("disabling_lesions", ",".join(call.damage.disabling_lesions)),
                    ("reason", call.reason),
                    ("sequence_class", call.region_class),
                    ("region_id", call.region_id),
                    ("coding_overlap_fraction", f"{call.coding_overlap_fraction:.4f}"),
                    ("repeat_overlap_fraction", f"{call.repeat_overlap_fraction:.4f}"),
                ]
            )
            handle.write(
                "\t".join(
                    [
                        call.chrom,
                        "felid_chry_homology",
                        "pseudogene",
                        str(call.start0 + 1),
                        str(call.end0),
                        ".",
                        call.strand,
                        ".",
                        gene_attrs,
                    ]
                )
                + "\n"
            )
            handle.write(
                "\t".join(
                    [
                        call.chrom,
                        "felid_chry_homology",
                        "pseudogenic_transcript",
                        str(call.start0 + 1),
                        str(call.end0),
                        ".",
                        call.strand,
                        ".",
                        _format_attributes(
                            [
                                ("ID", transcript_id),
                                ("Parent", call_id),
                                ("Name", call.family),
                                ("Target", call.best_hit.query_id),
                            ]
                        ),
                    ]
                )
                + "\n"
            )
            for exon_index, exon in enumerate(
                sorted(call.damage.projected_exons, key=lambda item: item.start0), start=1
            ):
                handle.write(
                    "\t".join(
                        [
                            call.chrom,
                            "felid_chry_homology",
                            "pseudogenic_exon",
                            str(exon.start0 + 1),
                            str(exon.end0),
                            ".",
                            call.strand,
                            ".",
                            _format_attributes(
                                [
                                    ("ID", f"{transcript_id}.pexon{exon_index}"),
                                    ("Parent", transcript_id),
                                ]
                            ),
                        ]
                    )
                    + "\n"
                )


def write_calls_tsv(callset: PseudogeneCallSet, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "species",
        "family",
        "chrom",
        "start0",
        "end0",
        "strand",
        "classification",
        "evidence_grade",
        "reason",
        "best_query",
        "query_coverage",
        "identity",
        "mapq",
        "disabling_lesions",
        "projected_cds_length",
        "projected_protein",
        "coding_overlap_fraction",
        "repeat_overlap_fraction",
        "region_id",
        "region_class",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for call in callset.calls:
            writer.writerow(
                {
                    "species": callset.species,
                    "family": call.family,
                    "chrom": call.chrom,
                    "start0": call.start0,
                    "end0": call.end0,
                    "strand": call.strand,
                    "classification": call.classification,
                    "evidence_grade": call.evidence_grade,
                    "reason": call.reason,
                    "best_query": call.best_hit.query_id,
                    "query_coverage": f"{call.best_hit.query_coverage:.6f}",
                    "identity": f"{call.best_hit.identity:.6f}",
                    "mapq": call.best_hit.mapq,
                    "disabling_lesions": ",".join(call.damage.disabling_lesions),
                    "projected_cds_length": len(call.damage.projected_cds),
                    "projected_protein": call.damage.translation.protein,
                    "coding_overlap_fraction": f"{call.coding_overlap_fraction:.6f}",
                    "repeat_overlap_fraction": f"{call.repeat_overlap_fraction:.6f}",
                    "region_id": call.region_id,
                    "region_class": call.region_class,
                }
            )


def _is_explicit_damage_status(status: str) -> bool:
    normalized = status.lower()
    return any(
        token in normalized
        for token in ("internal_stop", "localized_frameshift", "frameshift_event", "truncation_breakpoint")
    )


def compare_to_reference(
    callset: PseudogeneCallSet,
    reference_pseudogenes: list[GffFeature],
    *,
    explicit_damage_only: bool = True,
) -> CalibrationMetrics:
    references = [
        feature
        for feature in reference_pseudogenes
        if feature.feature_type == "pseudogene"
        and (
            not explicit_damage_only
            or _is_explicit_damage_status(feature.attributes.get("orf_status", ""))
        )
    ]
    matched: list[str] = []
    reliable = callset.reliable
    for reference in references:
        reference_family = normalize_family(
            reference.attributes, reference.feature_id or "unknown"
        )
        for call in reliable:
            if call.chrom != reference.seqid or call.family != reference_family:
                continue
            left, right = reciprocal_overlap(call.best_hit.interval, reference.interval)
            if min(left, right) >= 0.50:
                matched.append(reference.feature_id or f"line_{reference.line_number}")
                break
    total = len(references)
    recovered = len(matched)
    sensitivity = recovered / total if total else 1.0
    return CalibrationMetrics(
        reference_explicit_count=total,
        recovered_explicit_count=recovered,
        sensitivity=sensitivity,
        matched_reference_ids=tuple(matched),
    )


def validate_pseudogene_callset(
    callset: PseudogeneCallSet,
    *,
    chrom_length: int,
    calibration: CalibrationMetrics | None,
    search_manifest: dict[str, object],
) -> GateReport:
    checks: dict[str, bool] = {}
    checks["search_completed"] = bool(search_manifest.get("search_completed"))
    checks["query_set_nonempty"] = int(search_manifest.get("query_count", 0)) > 0
    checks["target_length_matches"] = int(search_manifest.get("target_length", -1)) == chrom_length
    if calibration is not None:
        checks["ppa_explicit_damage_sensitivity"] = calibration.sensitivity >= 0.75
    checks["coordinates_valid"] = all(
        0 <= call.start0 < call.end0 <= chrom_length for call in callset.calls
    )
    checks["no_substantial_coding_overlap"] = all(
        call.coding_overlap_fraction < DEFAULT_THRESHOLDS["annotation_conflict_overlap"]
        for call in callset.reliable
    )
    checks["evidence_complete"] = all(
        bool(call.best_hit.query_id)
        and bool(call.reason)
        and (
            call.reason == "concordant_strong_truncation"
            or any(
                lesion
                in {
                    "internal_stop",
                    "frameshift_alignment",
                    "inframe_stop_alignment",
                }
                for lesion in call.damage.disabling_lesions
            )
        )
        for call in callset.reliable
    )
    messages = tuple(name for name, passed in checks.items() if not passed)
    return GateReport(gate_pass=all(checks.values()), checks=checks, messages=messages)


def write_gate_report(report: GateReport, path: Path, **metadata: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "gate_pass": report.gate_pass,
        "checks": report.checks,
        "failed_checks": list(report.messages),
        **metadata,
    }
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def classify_rab9ay_evidence(
    hit: HomologyHit, damage: DamageEvidence, synteny_coherent: bool
) -> str:
    if (
        hit.query_coverage >= 0.90
        and hit.identity >= 0.85
        and damage.intact
        and synteny_coherent
    ):
        return "missed_intact_coding_gene"
    if (
        hit.query_coverage >= 0.50
        and hit.identity >= 0.70
        and (
            any(
                lesion
                in {
                    "frameshift_length",
                    "internal_stop",
                    "frameshift_alignment",
                    "inframe_stop_alignment",
                }
                for lesion in damage.disabling_lesions
            )
        )
    ):
        return "RAB9AY_pseudogene"
    if hit.query_coverage >= 0.50 and hit.identity >= 0.70:
        return "low_confidence_fragment"
    if hit.identity >= 0.70:
        return "subthreshold_RAB9A_like_fragment"
    return "no_reliable_homolog_detected"


# ===== pseudogene_workflow: self-contained historical logic =====
def _is_effectively_empty_gff(path: Path) -> bool:
    with path.open("r", encoding="utf-8") as handle:
        return not any(line.strip() and not line.startswith("#") for line in handle)


def optional_reference_pseudogenes(config: dict, species: str) -> Path | None:
    """Optional historical reference, never a prerequisite annotation result.

    Accept the old chry_pseudogene_gff field for explicit historical replays.
    An explicitly configured but nonexistent file is an input error.
    """
    values = config["species"][species]
    value = values.get("reference_pseudogene_gff", values.get("chry_pseudogene_gff"))
    if value is None or (isinstance(value, str) and not value.strip()):
        return None
    if not isinstance(value, str):
        raise ValueError(f"{species}: reference_pseudogene_gff must be a path or null")
    path = resolve_config_path(config, value)
    if path is None or not path.is_file():
        raise FileNotFoundError(f"{species} optional reference_pseudogene_gff: {path}")
    return path


def validate_inputs(config_path: Path, output: Path | None = None) -> list[dict[str, object]]:
    config = load_config(config_path)
    rows: list[dict[str, object]] = []
    required_keys = {
        "chry_fasta",
        "partition_bed",
        "chry_high_confidence_gff",
        "chry_provisional_gff",
        "genome_gff",
        "genome_cds_fasta",
        "genome_pep_fasta",
        "repeat_gff",
    }
    for species, values in config["species"].items():
        missing_keys = required_keys - values.keys()
        if missing_keys:
            raise ValueError(f"{species}: missing config keys {sorted(missing_keys)}")
        resolved: dict[str, Path] = {}
        for key in sorted(required_keys):
            path = resolve_config_path(config, values[key])
            assert path is not None
            if not path.is_file():
                raise FileNotFoundError(f"{species} {key}: {path}")
            resolved[key] = path

        fasta = read_fasta(resolved["chry_fasta"])
        if len(fasta) != 1:
            raise ValueError(f"{species}: expected one ChrY FASTA record, found {len(fasta)}")
        chrom, sequence = next(iter(fasta.items()))
        regions = read_bed_regions(resolved["partition_bed"])
        validate_partition(regions, len(sequence))
        high_confidence = read_gff3(resolved["chry_high_confidence_gff"])
        provisional = read_gff3(resolved["chry_provisional_gff"])
        reference_path = optional_reference_pseudogenes(config, species)
        pseudogenes = read_gff3(reference_path) if reference_path is not None else []
        rows.append(
            {
                "species": species,
                "chrom": chrom,
                "chrom_length": len(sequence),
                "n_count": sequence.count("N"),
                "region_count": len(regions),
                "high_confidence_genes": sum(
                    feature.feature_type == "gene" for feature in high_confidence
                ),
                "provisional_genes": sum(
                    feature.feature_type == "gene" for feature in provisional
                ),
                "reference_pseudogene_supplied": reference_path is not None,
                "reference_pseudogene_records": sum(
                    feature.feature_type in {"pseudogene", "gene"} for feature in pseudogenes
                ),
                "status": "valid",
            }
        )
    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        columns = list(rows[0])
        with output.open("w", encoding="utf-8") as handle:
            handle.write("\t".join(columns) + "\n")
            for row in rows:
                handle.write("\t".join(str(row[column]) for column in columns) + "\n")
    return rows


def build_queries(config_path: Path, output_dir: Path) -> list[Query]:
    config = load_config(config_path)
    models = []
    y_sequences: dict[str, str] = {}
    for species, values in config["species"].items():
        fasta_path = resolve_config_path(config, values["chry_fasta"])
        high_path = resolve_config_path(config, values["chry_high_confidence_gff"])
        provisional_path = resolve_config_path(config, values["chry_provisional_gff"])
        assert fasta_path is not None and high_path is not None and provisional_path is not None
        genome = next(iter(read_fasta(fasta_path).values()))
        y_sequences[species] = genome
        features = read_gff3(high_path) + read_gff3(provisional_path)
        models.extend(extract_cds_models(genome, features, species=species).values())
    queries = select_nonredundant_queries(models, max_per_family=4)

    ppa_pseudogene_path = optional_reference_pseudogenes(config, "PPA")
    if ppa_pseudogene_path is not None:
        queries.extend(
            derive_pseudogene_queries(
                y_sequences["PPA"],
                read_gff3(ppa_pseudogene_path),
                species="PPA_PG",
                min_family_copies=2,
                max_per_family=2,
            )
        )

    ppa = config["species"]["PPA"]
    ppa_cds_path = resolve_config_path(config, ppa["genome_cds_fasta"])
    ppa_pep_path = resolve_config_path(config, ppa["genome_pep_fasta"])
    assert ppa_cds_path is not None and ppa_pep_path is not None
    ppa_cds = read_fasta(ppa_cds_path)
    ppa_pep = read_fasta(ppa_pep_path)
    if "RAB9A" in ppa_cds:
        queries.append(
            Query(
                query_id="RAB9AY|PPA_X|RAB9A",
                family="RAB9AY",
                species="PPA_X",
                cds_sequence=ppa_cds["RAB9A"],
                protein_sequence=ppa_pep.get("RAB9A", "").rstrip("*"),
                source_gene_id="RAB9A",
            )
        )
    for species, values in config["species"].items():
        genome_gff_path = resolve_config_path(config, values["genome_gff"])
        cds_path = resolve_config_path(config, values["genome_cds_fasta"])
        pep_path = resolve_config_path(config, values["genome_pep_fasta"])
        assert genome_gff_path is not None and cds_path is not None and pep_path is not None
        genome_features = read_gff3(genome_gff_path)
        cds_records = read_fasta(cds_path)
        pep_records = read_fasta(pep_path)
        for feature in genome_features:
            if feature.feature_type != "gene" or feature.seqid == "ChrY":
                continue
            gene_id = feature.feature_id or ""
            name = feature.attributes.get("Name", gene_id)
            if not name.startswith("CCDC71L"):
                continue
            if gene_id not in cds_records or gene_id not in pep_records:
                continue
            queries.append(
                Query(
                    query_id=f"CCDC71L|{species}_A2|{gene_id}",
                    family="CCDC71L",
                    species=f"{species}_A2",
                    cds_sequence=cds_records[gene_id],
                    protein_sequence=pep_records[gene_id].rstrip("*"),
                    source_gene_id=gene_id,
                )
            )
    unique: dict[tuple[str, str], Query] = {}
    for query in queries:
        if query.family in {"UNRESOLVED", "XG", "ERVK-25"}:
            continue
        unique.setdefault((query.family, query.cds_sequence), query)
    final_queries = list(unique.values())
    write_queries(final_queries, output_dir)
    return final_queries


def _load_queries(query_dir: Path) -> list[Query]:
    cds = read_fasta(query_dir / "query_cds.fa")
    proteins = read_fasta(query_dir / "query_proteins.fa")
    queries: list[Query] = []
    with (query_dir / "query_manifest.tsv").open("r", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            query_id = row["query_id"]
            queries.append(
                Query(
                    query_id=query_id,
                    family=row["family"],
                    species=row["species"],
                    cds_sequence=cds[query_id],
                    protein_sequence=proteins[query_id],
                    source_gene_id=row["source_gene_id"],
                )
            )
    return queries


def _species_analysis_data(config: dict, species: str) -> dict[str, object]:
    values = config["species"][species]
    fasta_path = resolve_config_path(config, values["chry_fasta"])
    high_path = resolve_config_path(config, values["chry_high_confidence_gff"])
    provisional_path = resolve_config_path(config, values["chry_provisional_gff"])
    repeat_path = resolve_config_path(config, values["repeat_gff"])
    partition_path = resolve_config_path(config, values["partition_bed"])
    assert all(
        path is not None
        for path in (fasta_path, high_path, provisional_path, repeat_path, partition_path)
    )
    fasta_records = read_fasta(fasta_path)
    chrom, sequence = next(iter(fasta_records.items()))
    coding_features = read_gff3(high_path) + read_gff3(provisional_path)
    return {
        "chrom": chrom,
        "sequence": sequence,
        "fasta_path": fasta_path,
        "coding_features": coding_features,
        "coding_models": list(
            extract_cds_models(sequence, coding_features, species=species).values()
        ),
        "repeats": read_gff3(repeat_path),
        "regions": read_bed_regions(partition_path),
    }


def _nearest_gene_context(
    models: list[CdsModel], start0: int, end0: int, n_each_side: int = 2
) -> dict[str, str]:
    def label(model: CdsModel) -> str:
        return f"{model.family}|{model.gene_id}|{model.start0}-{model.end0}|{model.strand}"

    left = sorted(
        (model for model in models if model.end0 <= start0),
        key=lambda model: start0 - model.end0,
    )[:n_each_side]
    overlap = sorted(
        (
            model
            for model in models
            if min(model.end0, end0) > max(model.start0, start0)
        ),
        key=lambda model: model.start0,
    )
    right = sorted(
        (model for model in models if model.start0 >= end0),
        key=lambda model: model.start0 - end0,
    )[:n_each_side]
    return {
        "left_genes": ";".join(label(model) for model in reversed(left)),
        "overlapping_genes": ";".join(label(model) for model in overlap),
        "right_genes": ";".join(label(model) for model in right),
    }


def _write_rab9ay_review(
    *,
    output_dir: Path,
    queries: list[Query],
    fca_data: dict[str, object],
    ppa_data: dict[str, object],
    fca_protein_hits,
) -> dict[str, object]:
    fca_sequence = fca_data["sequence"]
    assert isinstance(fca_sequence, str)
    fca_fasta = fca_data["fasta_path"]
    assert isinstance(fca_fasta, Path)
    rab_queries = [query for query in queries if query.family == "RAB9AY"]
    nucleotide_hits = map_cdna_queries(
        fca_fasta, rab_queries, preset="splice:hq", best_n=50
    )
    write_hits_tsv(nucleotide_hits, output_dir / "RAB9AY.nucleotide_hits.tsv")
    protein_hits = [hit for hit in fca_protein_hits if hit.family == "RAB9AY"]
    write_hits_tsv(protein_hits, output_dir / "RAB9AY.protein_hits.tsv")
    all_evidence = protein_hits or nucleotide_hits
    best_hit = max(
        all_evidence,
        key=lambda hit: (hit.query_coverage * hit.identity, hit.alignment_length),
        default=None,
    )

    ppa_models = ppa_data["coding_models"]
    assert isinstance(ppa_models, list)
    ppa_rab_models = [model for model in ppa_models if model.family == "RAB9AY"]
    ppa_sequence = ppa_data["sequence"]
    assert isinstance(ppa_sequence, str)
    synteny_hits = []
    synteny_query_start0 = 0
    synteny_query_end0 = 0
    if ppa_rab_models:
        ppa_rab = ppa_rab_models[0]
        synteny_query_start0 = max(0, ppa_rab.start0 - 100_000)
        synteny_query_end0 = min(len(ppa_sequence), ppa_rab.end0 + 100_000)
        synteny_sequence = ppa_sequence[synteny_query_start0:synteny_query_end0]
        synteny_query = Query(
            query_id="RAB9AY_SYNTENY|PPA|plusminus100kb",
            family="RAB9AY_SYNTENY",
            species="PPA",
            cds_sequence=synteny_sequence,
            protein_sequence="",
            source_gene_id=ppa_rab.gene_id,
        )
        synteny_hits = [
            hit
            for hit in map_cdna_queries(
                fca_fasta, [synteny_query], preset="asm20", best_n=50
            )
            if hit.alignment_length >= 5_000 and hit.identity >= 0.60
        ]
    write_hits_tsv(synteny_hits, output_dir / "RAB9AY.synteny_hits.tsv")
    query_covered = sum(
        end0 - start0
        for start0, end0 in merge_intervals(
            [(hit.query_start0, hit.query_end0) for hit in synteny_hits]
        )
    )
    synteny_query_length = max(1, synteny_query_end0 - synteny_query_start0)
    synteny_query_coverage = query_covered / synteny_query_length
    synteny_coherent = bool(
        best_hit
        and synteny_query_coverage >= 0.50
        and any(interval_overlap(best_hit.interval, hit.interval) > 0 for hit in synteny_hits)
    )

    if best_hit is not None:
        damage = assess_projected_orf(best_hit, fca_sequence)
        classification = classify_rab9ay_evidence(best_hit, damage, synteny_coherent)
        fca_context = _nearest_gene_context(
            fca_data["coding_models"], best_hit.start0, best_hit.end0
        )
        region_id = "unassigned"
        region_class = "unassigned"
        overlapping_regions = [
            (interval_overlap(best_hit.interval, region.interval), region)
            for region in fca_data["regions"]
        ]
        overlap, region = max(overlapping_regions, key=lambda item: item[0])
        if overlap:
            region_id, region_class = region.region_id, region.region_class
    else:
        damage = None
        classification = "no_reliable_homolog_detected"
        fca_context = {"left_genes": "", "overlapping_genes": "", "right_genes": ""}
        region_id = region_class = "unassigned"

    if ppa_rab_models:
        ppa_context = _nearest_gene_context(
            ppa_models, ppa_rab_models[0].start0, ppa_rab_models[0].end0
        )
    else:
        ppa_context = {"left_genes": "", "overlapping_genes": "", "right_genes": ""}

    summary = {
        "classification": classification,
        "fca_chrom": best_hit.target if best_hit else "",
        "fca_start0": best_hit.start0 if best_hit else "",
        "fca_end0": best_hit.end0 if best_hit else "",
        "fca_strand": best_hit.strand if best_hit else "",
        "best_query": best_hit.query_id if best_hit else "",
        "protein_query_coverage": f"{best_hit.query_coverage:.6f}" if best_hit else "",
        "protein_identity": f"{best_hit.identity:.6f}" if best_hit else "",
        "projected_cds_length": len(damage.projected_cds) if damage else 0,
        "disabling_lesions": ",".join(damage.disabling_lesions) if damage else "",
        "region_id": region_id,
        "region_class": region_class,
        "synteny_query_start0_PPA": synteny_query_start0,
        "synteny_query_end0_PPA": synteny_query_end0,
        "synteny_query_coverage": f"{synteny_query_coverage:.6f}",
        "synteny_coherent": str(synteny_coherent).lower(),
        "FCA_left_genes": fca_context["left_genes"],
        "FCA_overlapping_genes": fca_context["overlapping_genes"],
        "FCA_right_genes": fca_context["right_genes"],
        "PPA_left_genes": ppa_context["left_genes"],
        "PPA_overlapping_genes": ppa_context["overlapping_genes"],
        "PPA_right_genes": ppa_context["right_genes"],
        "interpretation": (
            "High-identity syntenic RAB9A/Y-like remnant below the pre-specified 50% coverage threshold; "
            "not counted as a reliable pseudogene."
            if classification == "subthreshold_RAB9A_like_fragment"
            else classification
        ),
    }
    review_path = output_dir / "RAB9AY.review.tsv"
    with review_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary), delimiter="\t")
        writer.writeheader()
        writer.writerow(summary)

    with (output_dir / "RAB9AY.projected_sequences.fa").open(
        "w", encoding="utf-8"
    ) as handle:
        if damage and best_hit:
            handle.write(
                f">FCA_RAB9A_like_fragment|ChrY:{best_hit.start0 + 1}-{best_hit.end0}|{best_hit.strand}\n"
            )
            handle.write(damage.projected_cds + "\n")
            local_start = max(0, best_hit.start0 - 500)
            local_end = min(len(fca_sequence), best_hit.end0 + 500)
            handle.write(f">FCA_RAB9AY_local_context|ChrY:{local_start + 1}-{local_end}\n")
            handle.write(fca_sequence[local_start:local_end] + "\n")
        for query in rab_queries:
            handle.write(f">{query.query_id}\n{query.cds_sequence}\n")

    with (output_dir / "RAB9AY.hits.gff3").open("w", encoding="utf-8") as handle:
        handle.write("##gff-version 3\n")
        handle.write(f"##sequence-region ChrY 1 {len(fca_sequence)}\n")
        for index, hit in enumerate(protein_hits + nucleotide_hits, start=1):
            attributes = (
                f"ID=RAB9AY.hit.{index};Name=RAB9AY;Target={hit.query_id};"
                f"query_coverage={hit.query_coverage:.4f};identity={hit.identity:.4f};"
                f"evidence_source={hit.evidence_source}"
            )
            handle.write(
                "\t".join(
                    [
                        hit.target,
                        hit.evidence_source,
                        "match",
                        str(hit.start0 + 1),
                        str(hit.end0),
                        ".",
                        hit.strand,
                        ".",
                        attributes,
                    ]
                )
                + "\n"
            )
    return summary


def run_pseudogene_analysis(
    *,
    config_path: Path,
    output_dir: Path,
    query_dir: Path,
    miniprot_executable: Path,
    threads: int = 4,
    reuse_alignments: bool = False,
) -> dict[str, object]:
    config = load_config(config_path)
    if not (query_dir / "query_manifest.tsv").is_file():
        build_queries(config_path, query_dir)
    queries = _load_queries(query_dir)
    protein_lengths = {query.query_id: len(query.protein_sequence) for query in queries}
    output_dir.mkdir(parents=True, exist_ok=True)
    species_data = {
        species: _species_analysis_data(config, species) for species in ("PPA", "FCA")
    }
    callsets = {}
    protein_hits_by_species = {}
    for species in ("PPA", "FCA"):
        miniprot_gff = output_dir / f"{species}.miniprot.gff"
        log_path = output_dir / f"{species}.miniprot.log"
        if not (reuse_alignments and miniprot_gff.is_file() and miniprot_gff.stat().st_size):
            run_miniprot(
                executable=miniprot_executable,
                target_fasta=species_data[species]["fasta_path"],
                query_proteins=query_dir / "query_proteins.fa",
                output_gff=miniprot_gff,
                log_path=log_path,
                threads=threads,
            )
        hits = parse_miniprot_gff(miniprot_gff, protein_lengths)
        hits = [
            hit
            for hit in hits
            if hit.query_coverage >= 0.20 and hit.identity >= 0.65
        ]
        protein_hits_by_species[species] = hits
        write_hits_tsv(hits, output_dir / f"{species}.miniprot_hits.tsv")
        loci = cluster_family_hits(hits, max_gap=50)
        data = species_data[species]
        callset = classify_pseudogene_loci(
            loci,
            data["coding_models"],
            data["repeats"],
            data["sequence"],
            species=species,
            regions=data["regions"],
        )
        callsets[species] = callset
        write_calls_tsv(callset, output_dir / f"{species}.pseudogene_all_calls.tsv")

    ppa_reference_path = optional_reference_pseudogenes(config, "PPA")
    reference_metrics = (compare_to_reference(
        callsets["PPA"], read_gff3(ppa_reference_path), explicit_damage_only=True
    ) if ppa_reference_path is not None else None)
    # No usable reference means calibration was not assessed, never 100% sensitivity.
    calibration = (reference_metrics if reference_metrics is not None
                   and reference_metrics.reference_explicit_count > 0 else None)
    calibration_info = {
        "status": ("evaluated" if calibration is not None else
                   "not_assessed_no_reference" if ppa_reference_path is None else
                   "not_assessed_no_explicit_damage_reference"),
        "reference_count": calibration.reference_explicit_count if calibration else 0,
        "recovered_count": calibration.recovered_explicit_count if calibration else None,
        "sensitivity": calibration.sensitivity if calibration else None,
        "reference_definition": "PPA internal_stop/localized disabling-event subset",
    }
    with (output_dir / "PPA.calibration.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "reference_definition",
                "reference_count",
                "recovered_count",
                "sensitivity",
                "matched_reference_ids",
                "status",
            ]
        )
        writer.writerow(
            [
                "PPA internal_stop or localized disabling-event pseudogenes",
                calibration_info["reference_count"],
                calibration_info["recovered_count"] if calibration else "NA",
                f"{calibration.sensitivity:.6f}" if calibration else "NA",
                ",".join(calibration.matched_reference_ids) if calibration else "",
                calibration_info["status"],
            ]
        )

    search_manifests: dict[str, dict[str, object]] = {}
    species_gates = {}
    for species in ("PPA", "FCA"):
        sequence = species_data[species]["sequence"]
        assert isinstance(sequence, str)
        search_manifest = {
            "search_completed": True,
            "query_count": len(queries),
            "target_length": len(sequence),
            "protein_search": "miniprot",
            "nucleotide_cross_check": "mappy/minimap2",
            "miniprot_command_profile": (
                "-N 500 --outn 500 --outs 0.25 --outc 0.20 --gff --aln"
            ),
        }
        search_manifests[species] = search_manifest
        species_gates[species] = validate_pseudogene_callset(
            callsets[species],
            chrom_length=len(sequence),
            calibration=calibration,
            search_manifest=search_manifest,
        )
        reliable_gff = output_dir / f"{species}.ChrY.reliable_pseudogenes.gff3"
        write_pseudogene_gff3(
            callsets[species], reliable_gff, chrom_length=len(sequence)
        )
        reliable_features = read_gff3(reliable_gff)
        validate_gff_parentage(reliable_features)
        candidate_calls = PseudogeneCallSet(
            species=species,
            calls=tuple(
                call
                for call in callsets[species].calls
                if call.classification
                not in {"reliable_pseudogene", "accepted_coding_locus"}
            ),
        )
        write_calls_tsv(
            candidate_calls,
            output_dir / f"{species}.ChrY.pseudogene_candidates.tsv",
        )

    combined_checks = {
        f"{species}.{name}": passed
        for species in ("PPA", "FCA")
        for name, passed in species_gates[species].checks.items()
    }
    gate = GateReport(
        gate_pass=all(combined_checks.values()),
        checks=combined_checks,
        messages=tuple(name for name, passed in combined_checks.items() if not passed),
    )
    rab9ay = _write_rab9ay_review(
        output_dir=output_dir,
        queries=queries,
        fca_data=species_data["FCA"],
        ppa_data=species_data["PPA"],
        fca_protein_hits=protein_hits_by_species["FCA"],
    )
    write_gate_report(
        gate,
        output_dir / "quality_gate.json",
        shared_query_count=len(queries),
        shared_query_manifest=str(query_dir / "query_manifest.tsv"),
        preliminary_hit_filter={"query_coverage": 0.20, "identity": 0.65},
        standardized_calling_thresholds=DEFAULT_THRESHOLDS,
        species_quality_gates={
            species: {
                "gate_pass": species_gates[species].gate_pass,
                "checks": species_gates[species].checks,
                "failed_checks": list(species_gates[species].messages),
                "search_manifest": search_manifests[species],
            }
            for species in ("PPA", "FCA")
        },
        calibration=calibration_info,
        fca_reliable_pseudogene_count=len(callsets["FCA"].reliable),
        fca_reliable_family_counts=dict(
            Counter(call.family for call in callsets["FCA"].reliable)
        ),
        ppa_reliable_pseudogene_count=len(callsets["PPA"].reliable),
        ppa_reliable_family_counts=dict(
            Counter(call.family for call in callsets["PPA"].reliable)
        ),
        rab9ay_classification=rab9ay["classification"],
    )
    summary_rows = []
    for species in ("PPA", "FCA"):
        counts = Counter(call.classification for call in callsets[species].calls)
        families = Counter(call.family for call in callsets[species].reliable)
        summary_rows.append(
            {
                "species": species,
                "total_loci": len(callsets[species].calls),
                "accepted_coding": counts["accepted_coding_locus"],
                "reliable_pseudogene": counts["reliable_pseudogene"],
                "candidate_fragment": counts["candidate_fragment"],
                "subthreshold": counts["subthreshold_homology"],
                "unannotated_intact_candidate": counts["unannotated_intact_candidate"],
                "reliable_family_counts": ",".join(
                    f"{family}:{count}" for family, count in sorted(families.items())
                ),
            }
        )
    with (output_dir / "pseudogene_summary.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(summary_rows[0]), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(summary_rows)
    result = {
        "gate_pass": gate.gate_pass,
        "failed_checks": list(gate.messages),
        "ppa_calibration_status": calibration_info["status"],
        "ppa_calibration_sensitivity": calibration_info["sensitivity"],
        "fca_reliable_pseudogene_count": len(callsets["FCA"].reliable),
        "fca_reliable_family_counts": dict(
            Counter(call.family for call in callsets["FCA"].reliable)
        ),
        "ppa_reliable_pseudogene_count": len(callsets["PPA"].reliable),
        "ppa_reliable_family_counts": dict(
            Counter(call.family for call in callsets["PPA"].reliable)
        ),
        "rab9ay_classification": rab9ay["classification"],
    }
    return result

# ===== command-line execution =====
def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path, help="FCA/PPA input JSON")
    parser.add_argument("--out", required=True, type=Path, help="Output root")
    parser.add_argument("--miniprot", default="miniprot", type=Path)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--reuse-alignments", action="store_true",
                        help="Reuse nonempty {FCA,PPA}.miniprot.gff in 02_pseudogenes")
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args(argv)
    if args.threads < 1:
        parser.error("--threads must be positive")

    inventory = args.out / "00_input_inventory.tsv"
    rows = validate_inputs(args.config, inventory)
    if args.validate_only:
        print(json.dumps({"input_inventory": str(inventory), "species": rows},
                         ensure_ascii=False, indent=2))
        return 0

    query_dir = args.out / "01_queries"
    queries = build_queries(args.config, query_dir)
    if not queries:
        raise ValueError("No eligible protein/CDS queries were built from the supplied GFF3s")
    result = run_pseudogene_analysis(
        config_path=args.config,
        output_dir=args.out / "02_pseudogenes",
        query_dir=query_dir,
        miniprot_executable=args.miniprot,
        threads=args.threads,
        reuse_alignments=args.reuse_alignments,
    )
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if result["gate_pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
