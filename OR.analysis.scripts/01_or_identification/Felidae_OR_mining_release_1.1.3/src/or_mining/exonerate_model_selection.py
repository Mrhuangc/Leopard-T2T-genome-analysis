"""Parse targeted Exonerate output and select NC-compatible multi-exon loci.

This stage is deliberately a *screen*, not final gene classification.  It
reproduces the published OR_Finder ordering: retain spliced alignments, reject
probable chimeras using reduced TBLASTN segments, prefer longer query coverage
in descending tiers, and keep the shortest non-overlapping genomic model in
each tier.  Selected loci are expanded for a later strict best-query
refinement.
"""

from __future__ import annotations

import csv
import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

from or_mining.blast import _open_text, parse_blast_row
from or_mining.fasta import read_fasta


DNA_LINE = re.compile(r"^[ACGTRYKMSWBDHVNacgtrykmswbdhvn]+$")
ATTRIBUTE_NUMBER = re.compile(
    r"(?:^|;\s*)(?P<key>identity|similarity|frameshifts)\s+(?P<value>[0-9.]+)",
    re.IGNORECASE,
)
DEFAULT_QUERY_TIERS = (300, 270, 240, 210, 180, 150)


@dataclass(frozen=True, slots=True)
class Region:
    region_id: str
    sseqid: str
    expanded_start_1: int
    expanded_end_1: int
    slen: int

    @property
    def width(self) -> int:
        return self.expanded_end_1 - self.expanded_start_1 + 1


@dataclass(frozen=True, slots=True)
class GffFeature:
    kind: str
    start_1: int
    end_1: int
    strand: str
    attributes: str

    @property
    def width(self) -> int:
        return self.end_1 - self.start_1 + 1


@dataclass(slots=True)
class ExonerateModel:
    block_index: int
    query_id: str
    query_start_0: int
    query_end_0: int
    target_id: str
    target_start_0: int
    target_end_0: int
    strand: str
    score: float
    gene: GffFeature
    exons: tuple[GffFeature, ...]
    introns: tuple[GffFeature, ...]
    target_coding_sequence: str
    identity: float | None
    similarity: float | None
    frameshifts: int
    frameshift_bases: int
    region: Region
    tblastn_segments: int | None = None
    filter_status: str = "unclassified"
    selection_tier_aa: int | None = None
    candidate_id: str | None = None

    @property
    def query_span_aa(self) -> int:
        return abs(self.query_end_0 - self.query_start_0)

    @property
    def exon_count(self) -> int:
        return len(self.exons)

    @property
    def intron_count(self) -> int:
        return len(self.introns)

    @property
    def local_start_1(self) -> int:
        return self.gene.start_1

    @property
    def local_end_1(self) -> int:
        return self.gene.end_1

    @property
    def genomic_start_1(self) -> int:
        return self.region.expanded_start_1 + self.local_start_1 - 1

    @property
    def genomic_end_1(self) -> int:
        return self.region.expanded_start_1 + self.local_end_1 - 1

    @property
    def genomic_span_bp(self) -> int:
        return self.genomic_end_1 - self.genomic_start_1 + 1


def _read_regions(path: Path) -> dict[str, Region]:
    regions: dict[str, Region] = {}
    with path.open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "region_id",
            "sseqid",
            "expanded_start_1",
            "expanded_end_1",
            "slen",
        }
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path}: missing required columns {sorted(required)}")
        for line_number, row in enumerate(reader, start=2):
            try:
                start_1 = int(row["expanded_start_1"])
                end_1 = int(row["expanded_end_1"])
                slen = int(row["slen"])
            except (TypeError, ValueError) as exc:
                raise ValueError(f"{path}: line {line_number}: invalid coordinates") from exc
            region_id = row["region_id"]
            if not region_id or start_1 < 1 or end_1 < start_1 or end_1 > slen:
                raise ValueError(f"{path}: line {line_number}: invalid region")
            if region_id in regions:
                raise ValueError(f"{path}: duplicate region_id {region_id!r}")
            regions[region_id] = Region(region_id, row["sseqid"], start_1, end_1, slen)
    if not regions:
        raise ValueError("multi-exon region table is empty")
    return regions


def _read_manifest(path: Path) -> dict[str, set[str]]:
    selected: dict[str, set[str]] = {}
    with path.open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"region_id", "selected_query_ids"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path}: missing required columns {sorted(required)}")
        for line_number, row in enumerate(reader, start=2):
            region_id = row["region_id"]
            query_ids = {value for value in row["selected_query_ids"].split(",") if value}
            if not region_id or not query_ids:
                raise ValueError(f"{path}: line {line_number}: empty region/query set")
            if region_id in selected:
                raise ValueError(f"{path}: duplicate region_id {region_id!r}")
            selected[region_id] = query_ids
    return selected


def _validate_target_fasta(path: Path, regions: dict[str, Region]) -> None:
    seen: set[str] = set()
    for header, sequence in read_fasta(path):
        identifier = header.split()[0]
        region = regions.get(identifier)
        if region is None:
            raise ValueError(f"{path}: target FASTA has unknown region {identifier!r}")
        if identifier in seen:
            raise ValueError(f"{path}: duplicate FASTA identifier {identifier!r}")
        if len(sequence) != region.width:
            raise ValueError(
                f"{path}: {identifier}: expected {region.width} bp, got {len(sequence)}"
            )
        seen.add(identifier)
    if seen != set(regions):
        missing = sorted(set(regions) - seen)[:10]
        raise ValueError(f"{path}: missing target regions: {missing}")


def _attribute_numbers(attributes: str) -> dict[str, float]:
    values: dict[str, float] = {}
    for match in ATTRIBUTE_NUMBER.finditer(attributes):
        values[match.group("key").lower()] = float(match.group("value"))
    return values


def _parse_gff_line(line: str, target_id: str, strand: str) -> GffFeature | None:
    fields = line.rstrip("\r\n").split("\t")
    if len(fields) != 9 or fields[0].startswith("#"):
        return None
    if fields[0] != target_id:
        raise ValueError(f"GFF target {fields[0]!r} disagrees with vulgar target {target_id!r}")
    try:
        start_1 = int(fields[3])
        end_1 = int(fields[4])
    except ValueError as exc:
        raise ValueError(f"invalid GFF coordinates: {line.rstrip()}") from exc
    if start_1 < 1 or end_1 < start_1:
        raise ValueError(f"invalid GFF interval: {line.rstrip()}")
    if fields[6] != strand:
        raise ValueError(f"GFF strand {fields[6]!r} disagrees with vulgar strand {strand!r}")
    return GffFeature(fields[2].lower(), start_1, end_1, strand, fields[8])


def _parse_block(
    lines: list[str],
    block_index: int,
    regions: dict[str, Region],
    manifest: dict[str, set[str]],
) -> ExonerateModel:
    vulgar_lines = [line for line in lines if line.startswith("vulgar:")]
    if len(vulgar_lines) != 1:
        raise ValueError(f"alignment block {block_index}: expected one vulgar record")
    fields = vulgar_lines[0].split()
    if len(fields) < 10:
        raise ValueError(f"alignment block {block_index}: malformed vulgar record")
    try:
        query_id = fields[1]
        query_start_0 = int(fields[2])
        query_end_0 = int(fields[3])
        target_id = fields[5]
        target_start_0 = int(fields[6])
        target_end_0 = int(fields[7])
        strand = fields[8]
        score = float(fields[9])
    except ValueError as exc:
        raise ValueError(f"alignment block {block_index}: malformed vulgar values") from exc
    if strand not in {"+", "-"}:
        raise ValueError(f"alignment block {block_index}: invalid target strand {strand!r}")
    operation_fields = fields[10:]
    if len(operation_fields) % 3:
        raise ValueError(f"alignment block {block_index}: malformed vulgar operations")
    frameshift_bases = 0
    for offset in range(0, len(operation_fields), 3):
        operation, _query_length, target_length = operation_fields[offset : offset + 3]
        try:
            target_operation_length = int(target_length)
        except ValueError as exc:
            raise ValueError(
                f"alignment block {block_index}: malformed vulgar operation length"
            ) from exc
        if operation == "F":
            frameshift_bases += target_operation_length
    region = regions.get(target_id)
    if region is None:
        raise ValueError(f"alignment block {block_index}: unknown target {target_id!r}")
    allowed_queries = manifest.get(target_id)
    if allowed_queries is None or query_id not in allowed_queries:
        raise ValueError(
            f"alignment block {block_index}: query {query_id!r} was not selected for {target_id!r}"
        )

    in_gff = False
    after_gff = False
    features: list[GffFeature] = []
    coding_parts: list[str] = []
    for line in lines:
        if line.startswith("# --- START OF GFF DUMP ---"):
            in_gff = True
            after_gff = False
            continue
        if line.startswith("# --- END OF GFF DUMP ---"):
            in_gff = False
            after_gff = True
            continue
        if in_gff:
            feature = _parse_gff_line(line, target_id, strand)
            if feature is not None:
                features.append(feature)
        elif after_gff:
            sequence_line = line.strip()
            if sequence_line and DNA_LINE.fullmatch(sequence_line):
                coding_parts.append(sequence_line.upper())

    genes = [feature for feature in features if feature.kind == "gene"]
    exons = tuple(feature for feature in features if feature.kind == "exon")
    introns = tuple(feature for feature in features if feature.kind == "intron")
    cds = tuple(feature for feature in features if feature.kind == "cds")
    if len(genes) != 1 or not exons or len(cds) != len(exons):
        raise ValueError(
            f"alignment block {block_index}: expected one gene and paired CDS/exon features"
        )
    gene = genes[0]
    if gene.end_1 > region.width or any(feature.end_1 > region.width for feature in features):
        raise ValueError(f"alignment block {block_index}: GFF coordinates exceed target region")
    if introns and len(exons) != len(introns) + 1:
        raise ValueError(f"alignment block {block_index}: inconsistent exon/intron structure")
    if not introns and len(exons) != 1:
        raise ValueError(f"alignment block {block_index}: multiple exons without an intron")
    coding_sequence = "".join(coding_parts)
    expected_coding_bases = sum(feature.width for feature in cds)
    expected_tcs_bases = expected_coding_bases - frameshift_bases
    if not coding_sequence or len(coding_sequence) != expected_tcs_bases:
        raise ValueError(
            f"alignment block {block_index}: %tcs length {len(coding_sequence)} "
            f"does not match frameshift-adjusted CDS length {expected_tcs_bases}"
        )

    gene_values = _attribute_numbers(gene.attributes)
    frameshifts = sum(
        int(_attribute_numbers(exon.attributes).get("frameshifts", 0)) for exon in exons
    )
    return ExonerateModel(
        block_index=block_index,
        query_id=query_id,
        query_start_0=query_start_0,
        query_end_0=query_end_0,
        target_id=target_id,
        target_start_0=target_start_0,
        target_end_0=target_end_0,
        strand=strand,
        score=score,
        gene=gene,
        exons=exons,
        introns=introns,
        target_coding_sequence=coding_sequence,
        identity=gene_values.get("identity"),
        similarity=gene_values.get("similarity"),
        frameshifts=frameshifts,
        frameshift_bases=frameshift_bases,
        region=region,
    )


def parse_exonerate_models(
    path: Path,
    regions: dict[str, Region],
    manifest: dict[str, set[str]],
) -> list[ExonerateModel]:
    models: list[ExonerateModel] = []
    current: list[str] | None = None
    with path.open("rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("C4 Alignment:"):
                if current is not None:
                    models.append(_parse_block(current, len(models) + 1, regions, manifest))
                current = [line]
            elif current is not None:
                current.append(line)
    if current is not None:
        models.append(_parse_block(current, len(models) + 1, regions, manifest))
    if not models:
        raise ValueError("Exonerate result contains no alignment blocks")
    return models


def _reduce_intervals(intervals: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    reduced: list[tuple[int, int]] = []
    for start_0, end_0 in sorted(intervals):
        if not reduced or start_0 > reduced[-1][1]:
            reduced.append((start_0, end_0))
        else:
            reduced[-1] = (reduced[-1][0], max(reduced[-1][1], end_0))
    return reduced


def add_tblastn_segment_counts(
    models: list[ExonerateModel],
    merged_hits: Path,
    clearance_bp: int = 50,
    minimum_segment_bp: int = 50,
) -> None:
    if clearance_bp < 0 or minimum_segment_bp < 1:
        raise ValueError("TBLASTN clearance and minimum segment width are invalid")
    intervals: dict[int, list[tuple[int, int]]] = {model.block_index: [] for model in models}
    models_by_subject: dict[str, list[ExonerateModel]] = {}
    for model in models:
        models_by_subject.setdefault(model.region.sseqid, []).append(model)
    with _open_text(merged_hits) as handle:
        for line_number, line in enumerate(handle, start=1):
            try:
                hit = parse_blast_row(line)
            except ValueError as exc:
                raise ValueError(f"{merged_hits}: line {line_number}: {exc}") from exc
            for model in models_by_subject.get(hit.sseqid, []):
                if hit.strand != model.strand:
                    continue
                model_start_0 = model.genomic_start_1 - 1
                model_end_0 = model.genomic_end_1
                if (
                    hit.genome_start_0 >= max(0, model_start_0 - clearance_bp)
                    and hit.genome_end_0 <= min(model.region.slen, model_end_0 + clearance_bp)
                ):
                    intervals[model.block_index].append(
                        (hit.genome_start_0, hit.genome_end_0)
                    )
    for model in models:
        reduced = _reduce_intervals(intervals[model.block_index])
        model.tblastn_segments = sum(
            end_0 - start_0 > minimum_segment_bp for start_0, end_0 in reduced
        )


def _overlap(left: ExonerateModel, right: ExonerateModel) -> bool:
    return (
        left.region.sseqid == right.region.sseqid
        and left.strand == right.strand
        and left.genomic_start_1 <= right.genomic_end_1
        and right.genomic_start_1 <= left.genomic_end_1
    )


def _overlap_components(models: list[ExonerateModel]) -> list[list[ExonerateModel]]:
    ordered = sorted(
        models,
        key=lambda model: (
            model.region.sseqid,
            model.strand,
            model.genomic_start_1,
            model.genomic_end_1,
            model.block_index,
        ),
    )
    components: list[list[ExonerateModel]] = []
    current: list[ExonerateModel] = []
    current_end = -1
    current_key: tuple[str, str] | None = None
    for model in ordered:
        key = (model.region.sseqid, model.strand)
        if not current or key != current_key or model.genomic_start_1 > current_end:
            if current:
                components.append(current)
            current = [model]
            current_key = key
            current_end = model.genomic_end_1
        else:
            current.append(model)
            current_end = max(current_end, model.genomic_end_1)
    if current:
        components.append(current)
    return components


def select_models(
    models: list[ExonerateModel],
    maximum_exons: int = 2,
    minimum_intron_bp: int = 50,
    maximum_intron_bp: int = 10_000,
    query_tiers: tuple[int, ...] = DEFAULT_QUERY_TIERS,
) -> list[ExonerateModel]:
    if maximum_exons < 2:
        raise ValueError("maximum_exons must be at least 2")
    if minimum_intron_bp < 1 or maximum_intron_bp < minimum_intron_bp:
        raise ValueError("invalid intron bounds")
    if not query_tiers or tuple(sorted(query_tiers, reverse=True)) != query_tiers:
        raise ValueError("query tiers must be non-empty and descending")

    eligible: list[ExonerateModel] = []
    for model in models:
        if model.intron_count == 0:
            model.filter_status = "rejected_single_exon"
        elif model.query_span_aa <= query_tiers[-1]:
            model.filter_status = "rejected_short_query_alignment"
        elif model.exon_count > maximum_exons:
            model.filter_status = "rejected_exon_limit"
        elif any(
            intron.width < minimum_intron_bp or intron.width > maximum_intron_bp
            for intron in model.introns
        ):
            model.filter_status = "rejected_intron_length"
        elif model.tblastn_segments is None:
            raise ValueError("TBLASTN segment counts were not assigned")
        elif model.tblastn_segments == 0:
            model.filter_status = "rejected_no_tblastn_support"
        elif model.tblastn_segments > model.exon_count:
            model.filter_status = "rejected_probable_chimera"
        else:
            model.filter_status = "eligible"
            eligible.append(model)

    selected: list[ExonerateModel] = []
    for tier in query_tiers:
        tier_models = [
            model
            for model in eligible
            if model.query_span_aa > tier
            and not any(_overlap(model, accepted) for accepted in selected)
        ]
        for component in _overlap_components(tier_models):
            best = min(
                component,
                key=lambda model: (
                    model.genomic_span_bp,
                    -model.query_span_aa,
                    -model.score,
                    model.frameshifts,
                    -(model.identity or 0.0),
                    model.query_id,
                    model.block_index,
                ),
            )
            if not any(_overlap(best, accepted) for accepted in selected):
                best.filter_status = "selected_for_refinement"
                best.selection_tier_aa = tier
                selected.append(best)

    selected_ids = {model.block_index for model in selected}
    for model in eligible:
        if model.block_index not in selected_ids:
            model.filter_status = "rejected_overlapping_better_model"
    selected.sort(
        key=lambda model: (
            model.region.sseqid,
            model.genomic_start_1,
            model.genomic_end_1,
            model.strand,
        )
    )
    for index, model in enumerate(selected, start=1):
        model.candidate_id = f"ORMEC{index:06d}"
    return selected


def _format_optional(value: float | int | None) -> str:
    if value is None:
        return "."
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def _model_row(model: ExonerateModel, refinement_flank_bp: int) -> dict[str, str]:
    refinement_start_1 = max(1, model.genomic_start_1 - refinement_flank_bp)
    refinement_end_1 = min(model.region.slen, model.genomic_end_1 + refinement_flank_bp)
    return {
        "candidate_id": model.candidate_id or ".",
        "raw_block_index": str(model.block_index),
        "filter_status": model.filter_status,
        "selection_tier_aa": _format_optional(model.selection_tier_aa),
        "region_id": model.target_id,
        "query_id": model.query_id,
        "sseqid": model.region.sseqid,
        "strand": model.strand,
        "local_start_1": str(model.local_start_1),
        "local_end_1": str(model.local_end_1),
        "genomic_start_1": str(model.genomic_start_1),
        "genomic_end_1": str(model.genomic_end_1),
        "refinement_start_1": str(refinement_start_1),
        "refinement_end_1": str(refinement_end_1),
        "query_start_0": str(model.query_start_0),
        "query_end_0": str(model.query_end_0),
        "query_span_aa": str(model.query_span_aa),
        "exon_count": str(model.exon_count),
        "intron_count": str(model.intron_count),
        "genomic_span_bp": str(model.genomic_span_bp),
        "target_coding_bases": str(len(model.target_coding_sequence)),
        "exonerate_score": f"{model.score:.6g}",
        "identity_percent": _format_optional(model.identity),
        "similarity_percent": _format_optional(model.similarity),
        "frameshifts": str(model.frameshifts),
        "frameshift_bases": str(model.frameshift_bases),
        "tblastn_segments": _format_optional(model.tblastn_segments),
    }


MODEL_FIELDS = (
    "candidate_id",
    "raw_block_index",
    "filter_status",
    "selection_tier_aa",
    "region_id",
    "query_id",
    "sseqid",
    "strand",
    "local_start_1",
    "local_end_1",
    "genomic_start_1",
    "genomic_end_1",
    "refinement_start_1",
    "refinement_end_1",
    "query_start_0",
    "query_end_0",
    "query_span_aa",
    "exon_count",
    "intron_count",
    "genomic_span_bp",
    "target_coding_bases",
    "exonerate_score",
    "identity_percent",
    "similarity_percent",
    "frameshifts",
    "frameshift_bases",
    "tblastn_segments",
)


def _write_table(path: Path, models: list[ExonerateModel], refinement_flank_bp: int) -> None:
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MODEL_FIELDS, delimiter="\t")
        writer.writeheader()
        for model in models:
            writer.writerow(_model_row(model, refinement_flank_bp))


def _write_selected_gff3(path: Path, models: list[ExonerateModel]) -> None:
    with path.open("wt", encoding="utf-8", newline="\n") as handle:
        handle.write("##gff-version 3\n")
        for model in models:
            assert model.candidate_id is not None
            attributes = (
                f"ID={model.candidate_id};query={model.query_id};"
                f"source_region={model.target_id};exonerate_score={model.score:.6g};"
                f"query_span_aa={model.query_span_aa};frameshifts={model.frameshifts}"
            )
            handle.write(
                f"{model.region.sseqid}\tFelidae_OR_mining\tgene\t"
                f"{model.genomic_start_1}\t{model.genomic_end_1}\t{model.score:.6g}\t"
                f"{model.strand}\t.\t{attributes}\n"
            )
            ordered_exons = sorted(
                model.exons,
                key=lambda feature: feature.start_1,
                reverse=model.strand == "-",
            )
            for exon_number, exon in enumerate(ordered_exons, start=1):
                start_1 = model.region.expanded_start_1 + exon.start_1 - 1
                end_1 = model.region.expanded_start_1 + exon.end_1 - 1
                exon_id = f"{model.candidate_id}.exon{exon_number}"
                handle.write(
                    f"{model.region.sseqid}\tFelidae_OR_mining\texon\t{start_1}\t{end_1}\t.\t"
                    f"{model.strand}\t.\tID={exon_id};Parent={model.candidate_id}\n"
                )


def select_exonerate_candidates(
    raw_exonerate: Path,
    regions_tsv: Path,
    target_fasta: Path,
    manifest_tsv: Path,
    merged_hits: Path,
    output_prefix: Path,
    maximum_exons: int = 2,
    minimum_intron_bp: int = 50,
    maximum_intron_bp: int = 10_000,
    refinement_flank_bp: int = 600,
    blast_clearance_bp: int = 50,
    minimum_blast_segment_bp: int = 50,
) -> dict[str, Any]:
    if refinement_flank_bp < 0:
        raise ValueError("refinement_flank_bp must be non-negative")
    regions = _read_regions(regions_tsv)
    manifest = _read_manifest(manifest_tsv)
    if set(manifest) != set(regions):
        raise ValueError("targeted-pair manifest and region table IDs differ")
    _validate_target_fasta(target_fasta, regions)
    models = parse_exonerate_models(raw_exonerate, regions, manifest)
    support_models = [
        model
        for model in models
        if model.intron_count > 0
        and model.query_span_aa > DEFAULT_QUERY_TIERS[-1]
        and model.exon_count <= maximum_exons
        and all(
            minimum_intron_bp <= intron.width <= maximum_intron_bp
            for intron in model.introns
        )
    ]
    add_tblastn_segment_counts(
        support_models,
        merged_hits,
        clearance_bp=blast_clearance_bp,
        minimum_segment_bp=minimum_blast_segment_bp,
    )
    selected = select_models(
        models,
        maximum_exons=maximum_exons,
        minimum_intron_bp=minimum_intron_bp,
        maximum_intron_bp=maximum_intron_bp,
    )
    if not selected:
        raise ValueError("no multi-exon candidates passed the NC-compatible screen")

    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    all_models_path = Path(f"{output_prefix}.multiexon_models.all.tsv")
    candidates_path = Path(f"{output_prefix}.multiexon_refinement_candidates.tsv")
    gff3_path = Path(f"{output_prefix}.multiexon_screened_models.gff3")
    summary_path = Path(f"{output_prefix}.multiexon_model_selection.summary.json")
    _write_table(all_models_path, models, refinement_flank_bp)
    _write_table(candidates_path, selected, refinement_flank_bp)
    _write_selected_gff3(gff3_path, selected)

    status_counts: dict[str, int] = {}
    for model in models:
        status_counts[model.filter_status] = status_counts.get(model.filter_status, 0) + 1
    summary: dict[str, Any] = {
        "schema_version": "1.0",
        "status": "pass",
        "stage": "multi-exon candidate screening before strict local refinement",
        "method": "NC R_script_exonerate_parsing_new.R-compatible tiered non-overlap selection",
        "coordinates": {
            "exonerate_gff": "target-local 1-based inclusive",
            "output": "genome-global 1-based inclusive",
        },
        "parameters": {
            "query_alignment_tiers_aa_strictly_greater_than": list(DEFAULT_QUERY_TIERS),
            "maximum_exons": maximum_exons,
            "minimum_intron_bp": minimum_intron_bp,
            "maximum_intron_bp": maximum_intron_bp,
            "refinement_flank_bp": refinement_flank_bp,
            "tblastn_model_clearance_bp": blast_clearance_bp,
            "minimum_reduced_tblastn_segment_bp_strictly_greater_than": (
                minimum_blast_segment_bp
            ),
            "chimera_rule": "reduced TBLASTN segments must be <= exon count",
        },
        "counts": {
            "input_alignment_blocks": len(models),
            "input_target_regions": len(regions),
            "models_with_introns": sum(model.intron_count > 0 for model in models),
            "selected_refinement_candidates": len(selected),
            **dict(sorted(status_counts.items())),
        },
        "outputs": {
            "all_models": all_models_path.name,
            "refinement_candidates": candidates_path.name,
            "screened_gff3": gff3_path.name,
        },
    }
    with summary_path.open("wt", encoding="utf-8", newline="\n") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
    return summary


def main() -> int:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("raw_exonerate", type=Path)
    parser.add_argument("regions_tsv", type=Path)
    parser.add_argument("target_fasta", type=Path)
    parser.add_argument("manifest_tsv", type=Path)
    parser.add_argument("merged_hits", type=Path)
    parser.add_argument("--output-prefix", required=True, type=Path)
    parser.add_argument("--maximum-exons", type=int, default=2)
    parser.add_argument("--minimum-intron-bp", type=int, default=50)
    parser.add_argument("--maximum-intron-bp", type=int, default=10_000)
    parser.add_argument("--refinement-flank-bp", type=int, default=600)
    parser.add_argument("--blast-clearance-bp", type=int, default=50)
    parser.add_argument("--minimum-blast-segment-bp", type=int, default=50)
    args = parser.parse_args()
    select_exonerate_candidates(
        args.raw_exonerate,
        args.regions_tsv,
        args.target_fasta,
        args.manifest_tsv,
        args.merged_hits,
        args.output_prefix,
        maximum_exons=args.maximum_exons,
        minimum_intron_bp=args.minimum_intron_bp,
        maximum_intron_bp=args.maximum_intron_bp,
        refinement_flank_bp=args.refinement_flank_bp,
        blast_clearance_bp=args.blast_clearance_bp,
        minimum_blast_segment_bp=args.minimum_blast_segment_bp,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
