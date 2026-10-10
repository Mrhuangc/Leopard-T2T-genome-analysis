"""Build an audited pre-7TM OR locus catalog.

The published NC workflow applies an InterProScan IPR017452 gate after the
single- and multi-exon searches.  This module therefore consolidates loci but
deliberately does not call the functional-candidate count a final functional
gene count.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any
from urllib.parse import unquote
import csv
import json
import re

from or_mining.fasta import read_fasta, write_fasta_record
from or_mining.multiexon_classification import translate


CLASSIFICATIONS = ("functional_candidate", "pseudogene", "truncated", "ambiguous")


@dataclass(frozen=True, slots=True)
class Feature:
    feature_type: str
    start_1: int
    end_1: int
    phase: str


@dataclass(frozen=True, slots=True)
class Locus:
    source_id: str
    source_stage: str
    classification: str
    evidence: str
    query_id: str
    sseqid: str
    start_1: int
    end_1: int
    strand: str
    nucleotide: str
    protein: str | None
    features: tuple[Feature, ...]

    @property
    def exon_count(self) -> int:
        return sum(feature.feature_type == "exon" for feature in self.features)


def _fasta_map(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    for header, sequence in read_fasta(path):
        identifier = header.split()[0]
        if identifier in records:
            raise ValueError(f"duplicate FASTA identifier {identifier!r} in {path}")
        sequence = sequence.upper()
        if not sequence:
            raise ValueError(f"empty FASTA sequence for {identifier!r} in {path}")
        records[identifier] = sequence
    return records


def _read_fai(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) < 2:
            raise ValueError(f"{path}: malformed FAI line {line_number}")
        identifier, raw_length = fields[:2]
        if identifier in lengths:
            raise ValueError(f"{path}: duplicate contig {identifier!r}")
        length = int(raw_length)
        if length < 1:
            raise ValueError(f"{path}: non-positive length for {identifier!r}")
        lengths[identifier] = length
    if not lengths:
        raise ValueError(f"genome FAI is empty: {path}")
    return lengths


def _validate_interval(locus_id: str, contig: str, start_1: int, end_1: int,
                       contig_lengths: dict[str, int]) -> None:
    if contig not in contig_lengths:
        raise ValueError(f"{locus_id}: contig {contig!r} is absent from the genome FAI")
    if not 1 <= start_1 <= end_1 <= contig_lengths[contig]:
        raise ValueError(
            f"{locus_id}: interval {contig}:{start_1}-{end_1} is outside 1-{contig_lengths[contig]}"
        )


def _validate_translation(identifier: str, nucleotide: str, protein: str) -> None:
    if len(nucleotide) % 3:
        raise ValueError(f"{identifier}: functional CDS length is not divisible by three")
    translated = translate(nucleotide)
    if translated.endswith("*"):
        translated = translated[:-1]
    if "*" in translated:
        raise ValueError(f"{identifier}: functional CDS contains an internal stop")
    if translated != protein:
        raise ValueError(f"{identifier}: nucleotide and protein sequences disagree")


def _load_single_exon(
    nucleotide_fasta: Path,
    protein_fasta: Path,
    coordinates_tsv: Path,
    contig_lengths: dict[str, int],
) -> list[Locus]:
    nucleotides = _fasta_map(nucleotide_fasta)
    proteins = _fasta_map(protein_fasta)
    if set(nucleotides) != set(proteins):
        raise ValueError("single-exon nucleotide and protein FASTA identifiers differ")
    loci: list[Locus] = []
    observed: set[str] = set()
    with coordinates_tsv.open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"orf_id", "sseqid", "start_1", "end_1", "strand"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"single-exon coordinate table lacks fields: {sorted(required)}")
        for row in reader:
            identifier = row["orf_id"]
            if identifier in observed:
                raise ValueError(f"duplicate single-exon coordinate ID {identifier!r}")
            observed.add(identifier)
            if identifier not in nucleotides:
                raise ValueError(f"single-exon coordinate ID {identifier!r} is absent from FASTA")
            start_1, end_1 = int(row["start_1"]), int(row["end_1"])
            strand = row["strand"]
            if strand not in {"+", "-"}:
                raise ValueError(f"{identifier}: invalid strand {strand!r}")
            _validate_interval(identifier, row["sseqid"], start_1, end_1, contig_lengths)
            nucleotide, protein = nucleotides[identifier], proteins[identifier]
            if end_1 - start_1 + 1 != len(nucleotide):
                raise ValueError(f"{identifier}: genomic span and nucleotide length differ")
            _validate_translation(identifier, nucleotide, protein)
            loci.append(
                Locus(
                    identifier, "tree_filtered_single_exon", "functional_candidate",
                    "GPCR best-hit plus NC-compatible OR-clade membership", "",
                    row["sseqid"], start_1, end_1, strand, nucleotide, protein,
                    (Feature("exon", start_1, end_1, "."),
                     Feature("CDS", start_1, end_1, "0")),
                )
            )
    if observed != set(nucleotides):
        missing = sorted(set(nucleotides) - observed)
        raise ValueError(f"single-exon coordinate table omitted IDs: {missing[:5]}")
    return loci


def _parse_attributes(raw: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for item in raw.split(";"):
        if not item:
            continue
        if "=" not in item:
            raise ValueError(f"malformed GFF3 attribute {item!r}")
        key, value = item.split("=", 1)
        attributes[key] = unquote(value)
    return attributes


def _load_gff3(path: Path) -> tuple[dict[str, tuple[str, int, int, str, str]],
                                    dict[str, list[Feature]]]:
    genes: dict[str, tuple[str, int, int, str, str]] = {}
    children: dict[str, list[Feature]] = {}
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 9:
            raise ValueError(f"{path}: expected 9 GFF3 fields at line {line_number}")
        contig, _, feature_type, raw_start, raw_end, _, strand, phase, raw_attributes = fields
        start_1, end_1 = int(raw_start), int(raw_end)
        attributes = _parse_attributes(raw_attributes)
        if feature_type == "gene":
            identifier = attributes.get("ID")
            if not identifier or identifier in genes:
                raise ValueError(f"{path}: missing/duplicate gene ID at line {line_number}")
            genes[identifier] = (
                contig, start_1, end_1, strand, attributes.get("classification", "")
            )
        elif feature_type in {"exon", "CDS"}:
            parent = attributes.get("Parent")
            if not parent:
                raise ValueError(f"{path}: {feature_type} lacks Parent at line {line_number}")
            children.setdefault(parent, []).append(Feature(feature_type, start_1, end_1, phase))
    return genes, children


def _load_multiexon(
    calls_tsv: Path,
    gff3: Path,
    functional_nucleotide: Path,
    functional_protein: Path,
    pseudogene_nucleotide: Path,
    truncated_nucleotide: Path,
    ambiguous_nucleotide: Path,
    contig_lengths: dict[str, int],
) -> list[Locus]:
    class_fastas = {
        "functional": _fasta_map(functional_nucleotide),
        "pseudogene": _fasta_map(pseudogene_nucleotide),
        "truncated": _fasta_map(truncated_nucleotide),
        "ambiguous": _fasta_map(ambiguous_nucleotide),
    }
    proteins = _fasta_map(functional_protein)
    genes, children = _load_gff3(gff3)
    loci: list[Locus] = []
    observed: set[str] = set()
    with calls_tsv.open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "candidate_id", "query_id", "sseqid", "strand", "classification", "evidence",
            "genomic_start_1", "genomic_end_1",
        }
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"multi-exon calls table lacks fields: {sorted(required)}")
        for row in reader:
            identifier = row["candidate_id"]
            classification = row["classification"]
            if classification not in class_fastas:
                raise ValueError(f"{identifier}: invalid multi-exon class {classification!r}")
            if identifier in observed:
                raise ValueError(f"duplicate multi-exon call ID {identifier!r}")
            observed.add(identifier)
            if identifier not in class_fastas[classification]:
                raise ValueError(f"{identifier}: absent from its {classification} FASTA")
            start_1, end_1 = int(row["genomic_start_1"]), int(row["genomic_end_1"])
            strand = row["strand"]
            if strand not in {"+", "-"}:
                raise ValueError(f"{identifier}: invalid strand {strand!r}")
            _validate_interval(identifier, row["sseqid"], start_1, end_1, contig_lengths)
            gene = genes.get(identifier)
            if gene != (row["sseqid"], start_1, end_1, strand, classification):
                raise ValueError(f"{identifier}: calls table and GFF3 gene disagree")
            features = tuple(children.get(identifier, []))
            exons = [feature for feature in features if feature.feature_type == "exon"]
            cds = [feature for feature in features if feature.feature_type == "CDS"]
            if len(exons) < 2 or len(exons) != len(cds):
                raise ValueError(f"{identifier}: expected matching multi-exon exon/CDS features")
            exon_intervals = sorted((feature.start_1, feature.end_1) for feature in exons)
            cds_intervals = sorted((feature.start_1, feature.end_1) for feature in cds)
            if exon_intervals != cds_intervals:
                raise ValueError(f"{identifier}: exon and CDS intervals differ")
            if any(not (start_1 <= left <= right <= end_1) for left, right in exon_intervals):
                raise ValueError(f"{identifier}: child feature lies outside the gene")
            nucleotide = class_fastas[classification][identifier]
            spliced_length = sum(right - left + 1 for left, right in exon_intervals)
            if len(nucleotide) != spliced_length:
                raise ValueError(f"{identifier}: spliced exon length and nucleotide FASTA differ")
            protein = proteins.get(identifier)
            if classification == "functional":
                if protein is None:
                    raise ValueError(f"{identifier}: functional multi-exon protein is absent")
                _validate_translation(identifier, nucleotide, protein)
                output_class = "functional_candidate"
            else:
                if protein is not None:
                    raise ValueError(f"{identifier}: nonfunctional model appears in protein FASTA")
                output_class = classification
            loci.append(
                Locus(
                    identifier, "strict_multiexon", output_class, row["evidence"], row["query_id"],
                    row["sseqid"], start_1, end_1, strand, nucleotide, protein, features,
                )
            )
    expected = set().union(*(set(records) for records in class_fastas.values()))
    if observed != expected or observed != set(genes) or observed != set(children):
        raise ValueError("multi-exon calls, FASTA and GFF3 identifier sets differ")
    if set(proteins) != set(class_fastas["functional"]):
        raise ValueError("multi-exon functional nucleotide and protein identifier sets differ")
    return loci


def _natural_key(value: str) -> tuple[object, ...]:
    return tuple(
        int(part) if part.isdigit() else part.lower()
        for part in re.split(r"([0-9]+)", value)
    )


def _locus_key(locus: Locus) -> tuple[object, ...]:
    return (*_natural_key(locus.sseqid), locus.start_1, locus.end_1, locus.strand, locus.source_id)


def _overlap(left: Locus, right: Locus) -> int:
    if left.sseqid != right.sseqid:
        return 0
    return max(0, min(left.end_1, right.end_1) - max(left.start_1, right.start_1) + 1)


def _audit_overlaps(loci: list[Locus]) -> list[tuple[Locus, Locus, int]]:
    overlaps: list[tuple[Locus, Locus, int]] = []
    ordered = sorted(loci, key=_locus_key)
    for index, left in enumerate(ordered):
        for right in ordered[index + 1:]:
            if right.sseqid != left.sseqid:
                if _natural_key(right.sseqid) > _natural_key(left.sseqid):
                    break
                continue
            if right.start_1 > left.end_1:
                break
            width = _overlap(left, right)
            if width:
                overlaps.append((left, right, width))
    return overlaps


def _resolve_single_exon_redundancy(
    loci: list[Locus],
    overlaps: list[tuple[Locus, Locus, int]],
) -> tuple[list[Locus], dict[str, str], int]:
    """Reproduce the NC longest-model rule only for overlapping single-exon ORFs.

    The published ``Remove_redundancy.R`` reduces overlapping intervals and
    retains the longest member of each component.  We constrain that rule to
    tree-filtered single-exon ORFs.  Overlaps involving multi-exon calls remain
    fatal because they indicate a failure of the upstream exclusion/model
    selection contracts.
    """

    single_stage = "tree_filtered_single_exon"
    parent: dict[str, str] = {
        locus.source_id: locus.source_id
        for locus in loci
        if locus.source_stage == single_stage
    }

    def find(identifier: str) -> str:
        while parent[identifier] != identifier:
            parent[identifier] = parent[parent[identifier]]
            identifier = parent[identifier]
        return identifier

    def union(left_id: str, right_id: str) -> None:
        left_root, right_root = find(left_id), find(right_id)
        if left_root != right_root:
            parent[right_root] = left_root

    for left, right, _ in overlaps:
        if left.source_stage == right.source_stage == single_stage:
            union(left.source_id, right.source_id)

    by_id = {locus.source_id: locus for locus in loci}
    components: dict[str, list[Locus]] = {}
    for identifier in parent:
        components.setdefault(find(identifier), []).append(by_id[identifier])

    removed_to_selected: dict[str, str] = {}
    redundancy_components = 0
    for component in components.values():
        if len(component) == 1:
            continue
        redundancy_components += 1
        selected = min(
            component,
            key=lambda locus: (
                -(locus.end_1 - locus.start_1 + 1),
                -(len(locus.protein) if locus.protein is not None else 0),
                locus.source_id,
            ),
        )
        for locus in component:
            if locus.source_id != selected.source_id:
                removed_to_selected[locus.source_id] = selected.source_id

    retained = [locus for locus in loci if locus.source_id not in removed_to_selected]
    return retained, removed_to_selected, redundancy_components


CATALOG_FIELDS = (
    "catalog_id", "source_id", "source_stage", "classification", "evidence", "query_id",
    "sseqid", "start_1", "end_1", "strand", "genomic_span_bp", "spliced_cds_bp",
    "protein_aa", "exon_count", "requires_7tm_gate",
)


def build_preliminary_catalog(
    species: str,
    genome_fai: Path,
    single_nucleotide: Path,
    single_protein: Path,
    single_coordinates: Path,
    multiexon_calls: Path,
    multiexon_gff3: Path,
    multiexon_functional_nucleotide: Path,
    multiexon_functional_protein: Path,
    multiexon_pseudogene_nucleotide: Path,
    multiexon_truncated_nucleotide: Path,
    multiexon_ambiguous_nucleotide: Path,
    output_prefix: Path,
) -> dict[str, Any]:
    if not species or any(character.isspace() for character in species):
        raise ValueError("species must be a non-empty identifier without whitespace")
    contig_lengths = _read_fai(genome_fai)
    single = _load_single_exon(
        single_nucleotide, single_protein, single_coordinates, contig_lengths
    )
    multiexon = _load_multiexon(
        multiexon_calls, multiexon_gff3, multiexon_functional_nucleotide,
        multiexon_functional_protein, multiexon_pseudogene_nucleotide,
        multiexon_truncated_nucleotide, multiexon_ambiguous_nucleotide, contig_lengths,
    )
    loci = single + multiexon
    identifiers = [locus.source_id for locus in loci]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("single- and multi-exon source identifiers overlap")

    overlaps = _audit_overlaps(loci)
    invalid_overlaps = [
        overlap for overlap in overlaps
        if not (
            overlap[0].source_stage
            == overlap[1].source_stage
            == "tree_filtered_single_exon"
        )
    ]
    resolved_loci, removed_to_selected, redundancy_components = (
        _resolve_single_exon_redundancy(loci, overlaps)
    )
    overlap_path = Path(f"{output_prefix}.pre7tm_overlap_audit.tsv")
    with overlap_path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow((
            "left_id", "right_id", "sseqid", "left_start_1", "left_end_1",
            "right_start_1", "right_end_1", "overlap_bp", "left_strand",
            "right_strand", "left_stage", "right_stage", "action", "selected_id",
            "removed_ids",
        ))
        for left, right, width in overlaps:
            if left.source_stage == right.source_stage == "tree_filtered_single_exon":
                selected_candidates = {
                    removed_to_selected.get(left.source_id, left.source_id),
                    removed_to_selected.get(right.source_id, right.source_id),
                }
                if len(selected_candidates) != 1:
                    raise ValueError("internal error resolving a single-exon overlap component")
                selected_id = selected_candidates.pop()
                removed_ids = sorted(
                    identifier for identifier, winner in removed_to_selected.items()
                    if winner == selected_id
                )
                action = "collapse_single_exon_component_keep_longest"
            else:
                selected_id, removed_ids, action = "", [], "fatal_cross_stage_or_multiexon_conflict"
            writer.writerow((
                left.source_id, right.source_id, left.sseqid,
                left.start_1, left.end_1, right.start_1, right.end_1, width,
                left.strand, right.strand, left.source_stage, right.source_stage,
                action, selected_id, ";".join(removed_ids),
            ))
    if invalid_overlaps:
        first = invalid_overlaps[0]
        raise ValueError(
            "pre-7TM catalog contains an overlap outside the published single-exon "
            f"redundancy case (first: {first[0].source_id}/{first[1].source_id}, "
            f"{first[2]} bp); see the overlap audit"
        )

    remaining_overlaps = _audit_overlaps(resolved_loci)
    if remaining_overlaps:
        raise ValueError("internal error: overlaps remain after single-exon redundancy resolution")

    ordered = sorted(resolved_loci, key=_locus_key)
    counters = {classification: 0 for classification in CLASSIFICATIONS}
    prefixes = {
        "functional_candidate": "ORC", "pseudogene": "ORP",
        "truncated": "ORT", "ambiguous": "ORA",
    }
    catalog_ids: dict[str, str] = {}
    for locus in ordered:
        counters[locus.classification] += 1
        catalog_ids[locus.source_id] = (
            f"{species}_{prefixes[locus.classification]}{counters[locus.classification]:06d}"
        )

    catalog_path = Path(f"{output_prefix}.pre7tm_OR_catalog.tsv")
    id_map_path = Path(f"{output_prefix}.pre7tm_OR_id_map.tsv")
    gff_path = Path(f"{output_prefix}.pre7tm_OR_catalog.gff3")
    fasta_paths = {
        "functional_nucleotide": Path(f"{output_prefix}.pre7tm_functional_candidates.nucl.fa"),
        "functional_protein": Path(f"{output_prefix}.pre7tm_functional_candidates.prot.fa"),
        "pseudogene": Path(f"{output_prefix}.pre7tm_pseudogene.nucl.fa"),
        "truncated": Path(f"{output_prefix}.pre7tm_truncated.nucl.fa"),
        "ambiguous": Path(f"{output_prefix}.pre7tm_ambiguous.nucl.fa"),
    }
    fasta_handles = {key: path.open("wt", encoding="utf-8", newline="\n")
                     for key, path in fasta_paths.items()}
    try:
        with catalog_path.open("wt", encoding="utf-8", newline="") as catalog_handle, \
             id_map_path.open("wt", encoding="utf-8", newline="") as map_handle, \
             gff_path.open("wt", encoding="utf-8", newline="\n") as gff_handle:
            writer = csv.DictWriter(catalog_handle, fieldnames=CATALOG_FIELDS, delimiter="\t",
                                    lineterminator="\n")
            writer.writeheader()
            map_handle.write("catalog_id\tsource_id\tsource_stage\n")
            gff_handle.write("##gff-version 3\n##coordinate-system 1-based inclusive\n")
            for locus in ordered:
                catalog_id = catalog_ids[locus.source_id]
                map_handle.write(f"{catalog_id}\t{locus.source_id}\t{locus.source_stage}\n")
                writer.writerow({
                    "catalog_id": catalog_id, "source_id": locus.source_id,
                    "source_stage": locus.source_stage, "classification": locus.classification,
                    "evidence": locus.evidence, "query_id": locus.query_id,
                    "sseqid": locus.sseqid, "start_1": locus.start_1, "end_1": locus.end_1,
                    "strand": locus.strand, "genomic_span_bp": locus.end_1 - locus.start_1 + 1,
                    "spliced_cds_bp": len(locus.nucleotide),
                    "protein_aa": "" if locus.protein is None else len(locus.protein),
                    "exon_count": locus.exon_count,
                    "requires_7tm_gate": str(
                        locus.classification == "functional_candidate"
                    ).lower(),
                })
                attributes = (
                    f"ID={catalog_id};source_id={locus.source_id};"
                    f"classification={locus.classification};source_stage={locus.source_stage}"
                )
                gff_handle.write(
                    f"{locus.sseqid}\tfelidae_or\tgene\t{locus.start_1}\t{locus.end_1}\t.\t"
                    f"{locus.strand}\t.\t{attributes}\n"
                )
                child_order = sorted(
                    locus.features,
                    key=lambda feature: (
                        feature.start_1 if locus.strand == "+" else -feature.end_1,
                        0 if feature.feature_type == "exon" else 1,
                    ),
                )
                exon_number = 0
                for feature in child_order:
                    if feature.feature_type == "exon":
                        exon_number += 1
                        child_attributes = f"ID={catalog_id}.exon{exon_number};Parent={catalog_id}"
                    else:
                        child_attributes = f"Parent={catalog_id}"
                    gff_handle.write(
                        f"{locus.sseqid}\tfelidae_or\t{feature.feature_type}\t"
                        f"{feature.start_1}\t{feature.end_1}\t.\t{locus.strand}\t"
                        f"{feature.phase}\t{child_attributes}\n"
                    )
                header = (
                    f"{catalog_id} source_id={locus.source_id} {locus.sseqid}:"
                    f"{locus.start_1}-{locus.end_1} strand={locus.strand} "
                    f"class={locus.classification} exons={locus.exon_count}"
                )
                if locus.classification == "functional_candidate":
                    write_fasta_record(
                        fasta_handles["functional_nucleotide"], header, locus.nucleotide
                    )
                    assert locus.protein is not None
                    write_fasta_record(fasta_handles["functional_protein"], header, locus.protein)
                else:
                    write_fasta_record(
                        fasta_handles[locus.classification], header, locus.nucleotide
                    )
    finally:
        for handle in fasta_handles.values():
            handle.close()

    summary_path = Path(f"{output_prefix}.pre7tm_OR_catalog.summary.json")
    summary: dict[str, Any] = {
        "schema_version": "1.1",
        "status": "pass",
        "stage": "pre-7TM OR locus catalog consolidation",
        "reportable_final_catalog": False,
        "next_required_gate": "InterProScan IPR017452 on every functional candidate protein",
        "method": (
            "coordinate-aware single/multi-exon union; NC Remove_redundancy.R-compatible "
            "longest-model selection within overlapping tree-filtered single-exon components"
        ),
        "coordinates": {
            "input_and_output": "genome 1-based inclusive",
            "gff3": "1-based inclusive",
        },
        "counts": {
            "input_single_exon_functional_candidates": len(single),
            "input_multiexon_models": len(multiexon),
            "functional_candidates_requiring_7tm": counters["functional_candidate"],
            "pseudogene": counters["pseudogene"],
            "truncated": counters["truncated"],
            "ambiguous": counters["ambiguous"],
            "catalog_loci": len(ordered),
            "input_overlapping_locus_pairs": len(overlaps),
            "single_exon_redundancy_components": redundancy_components,
            "redundant_single_exon_models_removed": len(removed_to_selected),
            "invalid_cross_stage_or_multiexon_overlaps": len(invalid_overlaps),
            "remaining_overlapping_locus_pairs": len(remaining_overlaps),
            "sequence_only_deduplications": 0,
        },
        "policies": {
            "locus_identity": (
                "genomic interval and strand; identical sequences at distinct loci are retained"
            ),
            "overlap": (
                "overlapping tree-filtered single-exon components retain the longest genomic "
                "model; any overlap involving a multi-exon call is fatal"
            ),
            "functional_label": "candidate until IPR017452 gate passes",
        },
        "outputs": {
            "catalog": catalog_path.name, "gff3": gff_path.name,
            "id_map": id_map_path.name, "overlap_audit": overlap_path.name,
            **{key: path.name for key, path in fasta_paths.items()},
        },
    }
    summary_path.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return summary
