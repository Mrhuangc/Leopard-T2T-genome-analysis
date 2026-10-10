#!/usr/bin/env python3
"""Portable ChrY ampliconic-unit analysis and reporting pipeline."""

import argparse
from collections import Counter
import colorsys
import csv
from dataclasses import dataclass
import gzip
import math
from pathlib import Path
import re
import shutil
import statistics
import subprocess
import sys
import tempfile


PLOT_COLORS = {
    "forward_unit": "#4C78A8",
    "reverse_unit": "#E45756",
    "coding_tspy": "#2CA25F",
    "coding_ccdc71l": "#377EB8",
    "coding_cul4by": "#E69F00",
    "coding_gene": "#17A2A4",
    "pseudogene": "#9467BD",
    "exon": "#B2182B",
    "LINE/L1": "#E9A3C9",
    "RC/Helitron": "#756BB1",
    "SINE": "#80B1D3",
    "LTR": "#FDB462",
    "DNA transposon": "#8DD3C7",
    "satellite": "#FB8072",
    "simple repeat": "#B8E186",
    "low complexity": "#F6C85F",
    "other/unknown repeat": "#BDBDBD",
    "simple_repeat": "#B8E186",
    "low_complexity": "#F6C85F",
    "helitron_body": "#756BB1",
    "stainedglass": "#F28E2B",
    "repeat_rich_gap": "#8C2D2D",
    "unassigned_flank": "#BDBDBD",
}


GENE_COLOR_OVERRIDES = {
    "TSPY": "#2CA25F",
    "CCDC71L": "#377EB8",
    "CUL4BY": "#E69F00",
}


GENE_COLOR_PALETTE = [
    "#D55E00",
    "#CC79A7",
    "#009E73",
    "#56B4E9",
    "#B07AA1",
    "#E15759",
    "#59A14F",
    "#EDC948",
    "#4E79A7",
    "#F28E2B",
    "#76B7B2",
    "#FF9DA7",
    "#9C755F",
    "#BAB0AC",
    "#8C564B",
    "#BCBD22",
    "#17BECF",
]


ALL_UNITS_MAX_PNG_PIXELS = 30_000_000


REPEAT_PLOT_ORDER = [
    "RC/Helitron",
    "LINE/L1",
    "SINE",
    "LTR",
    "DNA transposon",
    "satellite",
    "simple repeat",
    "low complexity",
    "other/unknown repeat",
]


@dataclass(frozen=True)
class Region:
    chrom: str
    start1: int
    end1: int

    def __post_init__(self):
        if not self.chrom or self.start1 < 1 or self.end1 < self.start1:
            raise ValueError("invalid region")


@dataclass
class Feature:
    chrom: str
    source: str
    feature_type: str
    start0: int
    end0: int
    strand: str
    phase: str
    attributes: dict[str, str]
    score: str = "."
    raw_attributes: str | None = None


@dataclass
class RepeatHit:
    hit_id: str
    start0: int
    end0: int
    strand: str
    consensus_start1: int
    consensus_end1: int
    kimura: float | None


@dataclass(frozen=True)
class StainedGlassHit:
    query_chrom: str
    query_start0: int
    query_end0: int
    reference_chrom: str
    reference_start0: int
    reference_end0: int
    identity: float
    strand: str


@dataclass
class GeneModel:
    gene_id: str
    symbol: str
    strand: str
    cds_blocks: list[tuple[int, int, str]]
    is_pseudogene: bool


@dataclass
class Unit:
    unit_id: str
    chrom: str
    start0: int
    end0: int
    strand: str
    arm: str
    status: str


@dataclass(frozen=True)
class Hairpin:
    start0: int
    end0: int
    stem_left: str
    loop: str
    stem_right: str
    distance_to_terminal_end: int


@dataclass(frozen=True)
class TerminalPair:
    five_prime_offset: int
    five_prime_motif: str
    three_prime_offset: int
    three_prime_motif: str
    hairpin: Hairpin | None
    score: float
    five_prime_context: str = ""
    three_prime_context: str = ""


@dataclass(frozen=True)
class ORF:
    strand: str
    frame: int
    start0: int
    end0: int
    length_aa: int


@dataclass(frozen=True)
class HLEEvidence:
    terminal_support: float
    hairpin_support: float
    longest_orf_aa: int
    rep_domain: bool
    helicase_domain: bool
    independent_loci: int

    def __post_init__(self):
        if not 0.0 <= self.terminal_support <= 1.0 or not 0.0 <= self.hairpin_support <= 1.0:
            raise ValueError("HLE support fractions must be between zero and one")
        if self.longest_orf_aa < 0 or self.independent_loci < 0:
            raise ValueError("HLE counts cannot be negative")


@dataclass(frozen=True)
class ProjectedFeature:
    unit_id: str
    feature_type: str
    start0: int
    end0: int
    strand: str
    attributes: dict[str, str]
    source: str = "."
    phase: str = "."
    genomic_start0: int | None = None
    genomic_end0: int | None = None
    clipped: bool = False


_REGION_RE = re.compile(r"^(?P<chrom>[^:\s]+):(?P<start>[0-9]+)-(?P<end>[0-9]+)(?:\s|$)")


def parse_region_header(text: str) -> Region:
    """Parse a ``chrom:start-end`` region from a FASTA title or header."""
    match = _REGION_RE.match(text.strip().lstrip(">"))
    if match is None:
        raise ValueError("FASTA title must contain coordinates as chrom:start-end")
    try:
        return Region(match.group("chrom"), int(match.group("start")), int(match.group("end")))
    except ValueError as exc:
        raise ValueError("invalid FASTA region coordinates") from exc


def read_single_fasta(path: Path) -> tuple[str, str, Region | None]:
    """Read exactly one FASTA record, requiring coordinates in its title."""
    lines = Path(path).read_text().splitlines()
    headers = [i for i, line in enumerate(lines) if line.startswith(">")]
    if len(headers) != 1:
        raise ValueError("expected a single FASTA record")
    if any(line.strip() for line in lines[: headers[0]]):
        raise ValueError("FASTA contains nonblank preamble before record")
    title = lines[headers[0]][1:].strip()
    sequence = "".join(line.strip() for line in lines[headers[0] + 1 :] if line.strip())
    if not sequence:
        raise ValueError("FASTA sequence is empty")
    if any(base not in "ACGTNacgtn" for base in sequence):
        raise ValueError("FASTA sequence contains invalid bases")
    region = parse_region_header(title)
    closed_length = region.end1 - region.start1 + 1
    if len(sequence) == closed_length - 1:
        # Some regional FASTA exporters retain a BED half-open title.
        region = Region(region.chrom, region.start1 + 1, region.end1)
    elif len(sequence) != closed_length:
        raise ValueError("FASTA length does not match title coordinates")
    return title, sequence, region


def _parse_attributes(raw: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    if raw == "." or raw == "":
        return attrs
    for item in raw.split(";"):
        if not item:
            continue
        if "=" not in item:
            raise ValueError("malformed GFF attributes")
        key, value = item.split("=", 1)
        if not key:
            raise ValueError("malformed GFF attributes")
        attrs[key] = value
    return attrs


def read_gff(path: Path) -> list[Feature]:
    """Read GFF3 rows and convert 1-based closed coordinates to half-open."""
    features: list[Feature] = []
    for line_number, line in enumerate(Path(path).read_text().splitlines(), 1):
        if not line or line.startswith("#"):
            continue
        columns = line.split("\t")
        if len(columns) != 9:
            raise ValueError(f"malformed GFF row {line_number}: expected 9 columns")
        chrom, source, feature_type, start, end, score, strand, phase, raw_attrs = columns
        try:
            start1, end1 = int(start), int(end)
        except ValueError as exc:
            raise ValueError(f"malformed GFF row {line_number}: invalid coordinates") from exc
        try:
            score_valid = score == "." or math.isfinite(float(score))
        except ValueError:
            score_valid = False
        if (
            not chrom
            or not source
            or not feature_type
            or start1 < 1
            or end1 < start1
            or not score_valid
            or strand not in "+-."
            or phase not in "012."
        ):
            raise ValueError(f"malformed GFF row {line_number}")
        features.append(
            Feature(
                chrom,
                source,
                feature_type,
                start1 - 1,
                end1,
                strand,
                phase,
                _parse_attributes(raw_attrs),
                score,
                raw_attrs,
            )
        )
    return features


def stainedglass_hit_overlaps_region(hit: StainedGlassHit, region: Region) -> bool:
    """Return whether both half-open StainedGlass windows overlap the target."""
    region_start0 = region.start1 - 1
    return (
        hit.query_chrom == region.chrom
        and hit.reference_chrom == region.chrom
        and hit.query_start0 < region.end1
        and hit.query_end0 > region_start0
        and hit.reference_start0 < region.end1
        and hit.reference_end0 > region_start0
    )


def read_stainedglass_bed(path: Path, region: Region | None = None) -> list[StainedGlassHit]:
    """Read StainedGlass eight-column window alignments, optionally region-filtered."""
    opener = gzip.open if str(path).endswith(".gz") else open
    hits: list[StainedGlassHit] = []
    with opener(path, "rt") as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            columns = line.split("\t")
            if len(columns) < 8:
                raise ValueError(f"malformed StainedGlass BED row {line_number}: expected 8 columns")
            query_chrom, query_start, query_end, reference_chrom, reference_start, reference_end, identity, strand = columns[:8]
            try:
                query_start0, query_end0 = int(query_start), int(query_end)
                reference_start0, reference_end0 = int(reference_start), int(reference_end)
                identity_value = float(identity)
            except ValueError as exc:
                raise ValueError(f"malformed StainedGlass BED row {line_number}") from exc
            if (
                not query_chrom
                or not reference_chrom
                or query_start0 < 0
                or reference_start0 < 0
                or query_end0 <= query_start0
                or reference_end0 <= reference_start0
                or not math.isfinite(identity_value)
                or strand not in "+-"
            ):
                raise ValueError(f"malformed StainedGlass BED row {line_number}")
            hit = StainedGlassHit(
                query_chrom,
                query_start0,
                query_end0,
                reference_chrom,
                reference_start0,
                reference_end0,
                identity_value,
                strand,
            )
            if region is not None and not stainedglass_hit_overlaps_region(hit, region):
                continue
            hits.append(hit)
    return hits


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def extract_genomic(seq: str, fasta_region: Region, start1: int, end1: int, strand: str = "+") -> str:
    """Extract a genomic 1-based closed interval from a regional FASTA."""
    if strand not in "+-":
        raise ValueError("strand must be '+' or '-'")
    if start1 < fasta_region.start1 or end1 > fasta_region.end1 or end1 < start1:
        raise ValueError("coordinates outside FASTA region")
    expected_length = fasta_region.end1 - fasta_region.start1 + 1
    if len(seq) != expected_length:
        raise ValueError("FASTA sequence length does not match region")
    fragment = seq[start1 - fasta_region.start1 : end1 - fasta_region.start1 + 1]
    return fragment if strand == "+" else revcomp(fragment)


_CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
_STOP_CODONS = {"TAA", "TAG", "TGA"}


def _feature_is_pseudogene(feature: Feature) -> bool:
    if "pseudogene" in feature.feature_type.lower():
        return True
    if any(key.lower() == "pseudo" and value.lower() in {"true", "1", "yes"} for key, value in feature.attributes.items()):
        return True
    return any("pseudogene" in value.lower() for value in feature.attributes.values())


def collect_gene_models(features: list[Feature]) -> list[GeneModel]:
    """Group CDS features into genomic gene models using their GFF parent links."""
    by_id = {feature.attributes["ID"]: feature for feature in features if "ID" in feature.attributes}

    def ancestor(feature_id: str) -> list[Feature]:
        lineage: list[Feature] = []
        seen: set[str] = set()
        while feature_id in by_id and feature_id not in seen:
            seen.add(feature_id)
            feature = by_id[feature_id]
            lineage.append(feature)
            parents = feature.attributes.get("Parent", "").split(",")
            if not parents or not parents[0]:
                break
            feature_id = parents[0]
        return lineage

    grouped: dict[str, dict[str, object]] = {}
    for feature in features:
        if feature.feature_type.lower() != "cds":
            continue
        parent_ids = [parent for parent in feature.attributes.get("Parent", feature.attributes.get("ID", "")).split(",") if parent]
        for parent_id in parent_ids:
            lineage = ancestor(parent_id)
            root = lineage[-1] if lineage else None
            gene_id = root.attributes.get("ID", parent_id) if root else parent_id
            group = grouped.setdefault(
                gene_id,
                {"blocks": [], "lineage": lineage, "strand": feature.strand},
            )
            group["blocks"].append((feature.start0, feature.end0, feature.phase))
            if group["strand"] != feature.strand:
                raise ValueError(f"inconsistent CDS strands for {gene_id}")

    models: list[GeneModel] = []
    for gene_id, group in grouped.items():
        lineage = group["lineage"]
        metadata = lineage[-1] if lineage else None
        attrs = metadata.attributes if metadata else {}
        symbol = next((attrs[key] for key in ("Name", "gene_name", "gene", "locus_tag") if attrs.get(key)), gene_id)
        is_pseudogene = any(_feature_is_pseudogene(feature) for feature in lineage)
        models.append(GeneModel(gene_id, symbol, group["strand"], sorted(group["blocks"]), is_pseudogene))
    return sorted(models, key=lambda model: (model.cds_blocks[0][0], model.gene_id))


def extract_cds(model: GeneModel, sequence: str, region: Region) -> str:
    """Extract genomic CDS blocks, then orient the assembled transcript once."""
    if model.strand not in "+-":
        raise ValueError("gene model strand must be '+' or '-'")
    genomic_cds = "".join(
        extract_genomic(sequence, region, start0 + 1, end0)
        for start0, end0, _phase in sorted(model.cds_blocks)
    )
    return genomic_cds if model.strand == "+" else revcomp(genomic_cds)


def translate_dna(seq: str) -> str:
    """Translate complete codons with the standard genetic code; unknown codons are X."""
    sequence = seq.upper()
    return "".join(_CODON_TABLE.get(sequence[offset : offset + 3], "X") for offset in range(0, len(sequence) - 2, 3))


def summarize_cds_qc(cds_records: list[tuple[str, str]]) -> dict[str, object]:
    """Summarize coding-coordinate compatibility for complete coding models."""
    genes_tested = len(cds_records)
    atg_starts = 0
    terminal_stops = 0
    internal_stop_models = 0
    internal_stop_codons = 0
    compatible_terminal_models = 0
    for _gene_id, raw_cds in cds_records:
        cds = raw_cds.upper()
        in_frame = len(cds) % 3 == 0
        starts_at_atg = cds.startswith("ATG")
        terminal_stop = in_frame and len(cds) >= 3 and cds[-3:] in _STOP_CODONS
        internal_stops = sum(cds[offset : offset + 3] in _STOP_CODONS for offset in range(0, max(len(cds) - 3, 0), 3))
        atg_starts += starts_at_atg
        terminal_stops += terminal_stop
        internal_stop_models += internal_stops > 0
        internal_stop_codons += internal_stops
        compatible_terminal_models += starts_at_atg and terminal_stop

    compatible_terminal_fraction = compatible_terminal_models / genes_tested if genes_tested else 0.0
    reasons: list[str] = []
    if genes_tested == 0:
        status = "downgrade"
        reasons.append("no_complete_coding_models")
    else:
        status = "pass"
        if internal_stop_models >= 2:
            status = "fail"
            reasons.append("systematic_internal_stops")
        if compatible_terminal_fraction < 0.8:
            status = "fail"
            reasons.append("incompatible_terminal_structure")
        if not reasons:
            reasons.append("compatible_terminal_structure")
    return {
        "genes_tested": genes_tested,
        "atg_starts": atg_starts,
        "terminal_stops": terminal_stops,
        "internal_stop_models": internal_stop_models,
        "internal_stop_codons": internal_stop_codons,
        "compatible_terminal_models": compatible_terminal_models,
        "compatible_terminal_fraction": compatible_terminal_fraction,
        "coordinate_status": status,
        "coordinate_reason": ";".join(reasons),
        "downgrade": status == "downgrade",
    }


def run_input_qc(sequence: str, region: Region, repeats: list[Feature], genes: list[GeneModel]) -> dict[str, object]:
    """Run deterministic coordinate QC while excluding pseudogene models from coding checks."""
    coding_genes = [gene for gene in genes if not gene.is_pseudogene]
    qc = summarize_cds_qc([(gene.gene_id, extract_cds(gene, sequence, region)) for gene in coding_genes])
    qc.update(
        {
            "region": f"{region.chrom}:{region.start1}-{region.end1}",
            "sequence_length": len(sequence),
            "repeat_features": len(repeats),
            "gene_models": len(genes),
            "coding_models": len(coding_genes),
            "pseudogene_models": len(genes) - len(coding_genes),
        }
    )
    return qc


def summarize_pseudogene_qc(features: list[Feature], fasta_region: Region) -> dict[str, object]:
    """Count pseudogene annotations and verify that their coordinates are extractable."""
    pseudogene_features = [
        feature
        for feature in features
        if "pseudogen" in feature.feature_type.lower() or _feature_is_pseudogene(feature)
    ]
    top_level = [feature for feature in pseudogene_features if feature.feature_type.lower() == "pseudogene"]
    fasta_start0 = fasta_region.start1 - 1
    coordinates_valid = all(
        feature.chrom == fasta_region.chrom
        and fasta_start0 <= feature.start0 < feature.end0 <= fasta_region.end1
        for feature in pseudogene_features
    )
    return {
        "pseudogene_models": len(top_level),
        "pseudogene_features": len(pseudogene_features),
        "pseudogene_coordinate_status": "pass" if coordinates_valid else "fail",
    }


def _consensus_bounds(hit: RepeatHit) -> tuple[int, int]:
    """Return the closed, orientation-neutral consensus interval for a hit."""
    return min(hit.consensus_start1, hit.consensus_end1), max(hit.consensus_start1, hit.consensus_end1)


def _consensus_span(hit: RepeatHit) -> int:
    start1, end1 = _consensus_bounds(hit)
    return end1 - start1 + 1


def select_repeat_hits(features: list[Feature], family_name: str) -> list[RepeatHit]:
    """Select valid, exactly named repeat-family annotations from a GFF feature list."""
    hits: list[RepeatHit] = []
    for index, feature in enumerate(features, 1):
        attrs = feature.attributes
        if (
            attrs.get("Name") != family_name
            or feature.strand not in "+-"
            or feature.end0 <= feature.start0
        ):
            continue
        try:
            consensus_start1 = int(attrs["tstart"])
            consensus_end1 = int(attrs["tend"])
        except (KeyError, TypeError, ValueError):
            continue
        if consensus_start1 < 1 or consensus_end1 < 1:
            continue
        kimura: float | None
        try:
            parsed_kimura = float(attrs["kimura80"])
            kimura = parsed_kimura if math.isfinite(parsed_kimura) else None
        except (KeyError, TypeError, ValueError):
            kimura = None
        hit_id = attrs.get("ID") or f"repeat_{index}"
        hits.append(
            RepeatHit(
                hit_id,
                feature.start0,
                feature.end0,
                feature.strand,
                consensus_start1,
                consensus_end1,
                kimura,
            )
        )
    return sorted(hits, key=lambda hit: (hit.start0, hit.end0, hit.hit_id))


def select_focal_family_features(
    features: list[Feature], family_name: str
) -> list[Feature]:
    """Select raw exact-name annotations for a non-filtering plot overlay.

    Unlike family-only anchor discovery, an overlay does not require consensus
    coordinates or a resolved strand, and adjacent GFF records stay separate.
    """
    return sorted(
        (
            feature
            for feature in features
            if feature.attributes.get("Name") == family_name
            and feature.end0 > feature.start0
        ),
        key=lambda feature: (
            feature.start0,
            feature.end0,
            feature.attributes.get("ID", ""),
        ),
    )


def _can_merge_fragments(
    left: RepeatHit,
    right: RepeatHit,
    max_genomic_gap: int,
    max_consensus_gap: int,
) -> bool:
    """Check adjacency in genomic order and in the strand-aware consensus order."""
    if left.strand != right.strand or right.start0 < left.end0:
        return False
    genomic_gap = right.start0 - left.end0
    if genomic_gap > max_genomic_gap:
        return False
    left_low, left_high = _consensus_bounds(left)
    right_low, right_high = _consensus_bounds(right)
    consensus_gap = (
        right_low - left_high - 1
        if left.strand == "+"
        else left_low - right_high - 1
    )
    return 0 <= consensus_gap <= max_consensus_gap


def _merge_hits(left: RepeatHit, right: RepeatHit) -> RepeatHit:
    """Combine adjacent pieces into a single family hit."""
    left_span = _consensus_span(left)
    right_span = _consensus_span(right)
    values = [(left.kimura, left_span), (right.kimura, right_span)]
    known_weight = sum(span for value, span in values if value is not None)
    kimura = (
        sum(value * span for value, span in values if value is not None) / known_weight
        if known_weight
        else None
    )
    consensus_start1 = min(_consensus_bounds(left)[0], _consensus_bounds(right)[0])
    consensus_end1 = max(_consensus_bounds(left)[1], _consensus_bounds(right)[1])
    return RepeatHit(
        f"{left.hit_id},{right.hit_id}",
        left.start0,
        right.end0,
        left.strand,
        consensus_start1,
        consensus_end1,
        kimura,
    )


def merge_family_fragments(
    hits: list[RepeatHit], max_genomic_gap: int = 100, max_consensus_gap: int = 100
) -> list[RepeatHit]:
    """Merge adjacent pieces only when genomic and consensus directions agree.

    Opposite-strand overlaps are intentionally kept separate because they can mark an
    inversion junction rather than a fragmented repeat annotation.
    """
    if max_genomic_gap < 0 or max_consensus_gap < 0:
        raise ValueError("fragment gap limits must be non-negative")
    merged: list[RepeatHit] = []
    for hit in sorted(hits, key=lambda item: (item.start0, item.end0, item.hit_id)):
        if hit.strand not in "+-" or hit.end0 <= hit.start0 or _consensus_bounds(hit)[0] < 1:
            continue
        if merged and _can_merge_fragments(merged[-1], hit, max_genomic_gap, max_consensus_gap):
            merged[-1] = _merge_hits(merged[-1], hit)
        else:
            merged.append(hit)
    return merged


def _family_span(hits: list[RepeatHit]) -> int:
    if not hits:
        return 0
    return max(_consensus_bounds(hit)[1] for hit in hits) - min(_consensus_bounds(hit)[0] for hit in hits) + 1


def _primary_anchors(hits: list[RepeatHit]) -> list[RepeatHit]:
    """Keep complete or near-complete hits as homologous anchor starts."""
    span = _family_span(hits)
    if not span:
        return []
    minimum_span = 0.9 * span
    anchors = [
        hit
        for hit in hits
        if (
            hit.strand in "+-"
            and hit.end0 > hit.start0
            and _consensus_bounds(hit)[0] >= 1
            and _consensus_span(hit) >= minimum_span
        )
    ]
    return sorted(anchors, key=lambda hit: (hit.start0, hit.end0, hit.hit_id))


def _period_summary(starts: list[int]) -> dict[str, int | float | None]:
    intervals = [right - left for left, right in zip(starts, starts[1:]) if right > left]
    if not intervals:
        return {"median": None, "mad": None, "min": None, "max": None, "n_anchors": len(starts), "n_intervals": 0}
    median = statistics.median(intervals)
    return {
        "median": median,
        "mad": statistics.median([abs(interval - median) for interval in intervals]),
        "min": min(intervals),
        "max": max(intervals),
        "n_anchors": len(starts),
        "n_intervals": len(intervals),
    }


def estimate_period(hits: list[RepeatHit]) -> dict[str, object]:
    """Estimate an array period from same-strand, near-full family anchors.

    With fewer than two usable anchors on a strand its statistics are ``None`` and
    it contributes no interval to the overall estimate.  This avoids treating the
    long interval across a central inversion as a tandem-repeat period.
    """
    anchors = _primary_anchors(hits)
    by_strand = {
        strand: _period_summary([hit.start0 for hit in anchors if hit.strand == strand])
        for strand in "+-"
    }
    all_intervals: list[int] = []
    for strand in "+-":
        starts = [hit.start0 for hit in anchors if hit.strand == strand]
        all_intervals.extend(right - left for left, right in zip(starts, starts[1:]) if right > left)
    if all_intervals:
        median = statistics.median(all_intervals)
        overall = {
            "median": median,
            "mad": statistics.median([abs(interval - median) for interval in all_intervals]),
            "min": min(all_intervals),
            "max": max(all_intervals),
            "n_anchors": len(anchors),
            "n_intervals": len(all_intervals),
        }
    else:
        overall = {
            "median": None,
            "mad": None,
            "min": None,
            "max": None,
            "n_anchors": len(anchors),
            "n_intervals": 0,
        }
    return {
        **overall,
        "consensus_span": _family_span(hits),
        "by_strand": by_strand,
        "plus_median": by_strand["+"]["median"],
        "minus_median": by_strand["-"]["median"],
    }


def estimate_stainedglass_period(
    hits: list[StainedGlassHit],
    region: Region,
    min_identity: float = 95.0,
    period_override: int | None = None,
) -> dict[str, object]:
    """Estimate the smallest coherent repeat period from StainedGlass lag harmonics."""
    regional = [
        hit
        for hit in hits
        if stainedglass_hit_overlaps_region(hit, region)
    ]
    spans = [hit.query_end0 - hit.query_start0 for hit in regional]
    if not spans:
        raise ValueError("no StainedGlass alignments overlap the target region")
    window_size = int(statistics.median(spans))
    if window_size < 1:
        raise ValueError("StainedGlass window size is invalid")

    pair_identity: dict[tuple[int, int], float] = {}
    for hit in regional:
        if hit.identity < min_identity or hit.query_start0 == hit.reference_start0:
            continue
        left, right = sorted((hit.query_start0, hit.reference_start0))
        pair_identity[(left, right)] = max(pair_identity.get((left, right), -math.inf), hit.identity)
    lag_counts = Counter(
        max(1, int(round((right - left) / window_size)))
        for left, right in pair_identity
        if right > left
    )
    if not lag_counts and period_override is None:
        raise ValueError(f"no off-diagonal StainedGlass alignments reach {min_identity:g}% identity")
    n_windows = max(1, math.ceil((region.end1 - region.start1 + 1) / window_size))
    strengths = {
        lag: min(1.0, count / max(n_windows - lag, 1))
        for lag, count in lag_counts.items()
    }
    lag_rows = [
        {
            "lag_windows": lag,
            "lag_bp": lag * window_size,
            "high_identity_pairs": lag_counts[lag],
            "normalized_support": strengths[lag],
        }
        for lag in sorted(lag_counts)
    ]

    maximum_lag = max(lag_counts, default=max(2, n_windows // 2))
    candidate_rows: list[dict[str, object]] = []
    for candidate in range(2, maximum_lag // 2 + 1):
        harmonics = list(range(1, min(12, maximum_lag // candidate) + 1))
        if len(harmonics) < 2:
            continue
        harmonic_support = [
            max(strengths.get(candidate * multiple + delta, 0.0) for delta in (-1, 0, 1))
            for multiple in harmonics
        ]
        weights = [1.0 / multiple for multiple in harmonics]
        raw_score = sum(value * weight for value, weight in zip(harmonic_support, weights)) / sum(weights)
        observed = sum(value >= 0.02 for value in harmonic_support)
        coherence = min(1.0, observed / 3.0)
        exact_fundamental = strengths.get(candidate, 0.0)
        score = raw_score * coherence * (0.25 + 0.75 * exact_fundamental)
        candidate_rows.append(
            {
                "candidate_period_bp": candidate * window_size,
                "score": score,
                "observed_harmonics": observed,
                "fundamental_support": exact_fundamental,
                "harmonic_support": ",".join(f"{value:.6f}" for value in harmonic_support),
            }
        )
    candidate_rows.sort(key=lambda row: (-float(row["score"]), int(row["candidate_period_bp"])))
    if period_override is not None:
        if period_override < 1:
            raise ValueError("period override must be positive")
        period_bp = int(period_override)
        method = "stainedglass_manual_period"
    elif candidate_rows and float(candidate_rows[0]["score"]) > 0:
        best_score = float(candidate_rows[0]["score"])
        eligible = [
            row
            for row in candidate_rows
            if float(row["score"]) >= 0.50 * best_score and int(row["observed_harmonics"]) >= 2
        ]
        chosen = min(eligible or candidate_rows[:1], key=lambda row: int(row["candidate_period_bp"]))
        period_bp = int(chosen["candidate_period_bp"])
        method = "stainedglass_lag_harmonics"
    else:
        raise ValueError("StainedGlass lag spectrum does not support a repeat period; use --period-bp")
    return {
        "median": period_bp,
        "mad": window_size,
        "min": period_bp,
        "max": period_bp,
        "n_anchors": len(pair_identity),
        "n_intervals": len(lag_counts),
        "consensus_span": None,
        "by_strand": {"+": {}, "-": {}},
        "plus_median": None,
        "minus_median": None,
        "window_size": window_size,
        "min_identity": min_identity,
        "method": method,
        "lag_rows": lag_rows,
        "candidate_scores": candidate_rows,
    }


def locate_inversion_junction(hits: list[RepeatHit]) -> int | None:
    """Locate one well-supported genomic ``+`` to ``-`` arm transition."""
    anchors = _primary_anchors(hits)
    if len([hit for hit in anchors if hit.strand == "+"]) < 2 or len([hit for hit in anchors if hit.strand == "-"]) < 2:
        return None
    strands = [hit.strand for hit in anchors]
    transitions = [index for index, pair in enumerate(zip(strands, strands[1:])) if pair[0] != pair[1]]
    if len(transitions) != 1:
        return None
    transition = transitions[0]
    if strands[: transition + 1] != ["+"] * (transition + 1) or strands[transition + 1 :] != ["-"] * (len(strands) - transition - 1):
        return None
    left, right = anchors[transition], anchors[transition + 1]
    return (left.end0 + right.start0) // 2


def _period_value(period: int | float | dict[str, object]) -> float:
    value = period.get("median") if isinstance(period, dict) else period
    if not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
        raise ValueError("period must be a positive finite number or estimate with a median")
    return float(value)


def build_units(hits: list[RepeatHit], period: int | float | dict[str, object], region: Region) -> list[Unit]:
    """Build region-clipped units from successive near-full homologous anchors.

    A fractional period is rounded half-up only for the extrapolated boundary
    widths; observed anchor-to-anchor boundaries remain their exact coordinates.
    """
    anchor_hits = _primary_anchors(hits)
    if not anchor_hits:
        return []
    period_bp = math.floor(_period_value(period) + 0.5)
    if period_bp < 1:
        raise ValueError("period must round to at least one base pair")
    region_start0, region_end0 = region.start1 - 1, region.end1
    junction = locate_inversion_junction(anchor_hits)

    boundaries: list[tuple[int, int, RepeatHit, str]] = []
    first, last = anchor_hits[0], anchor_hits[-1]
    left_start0, left_end0 = first.start0 - period_bp, first.start0
    while left_end0 > region_start0:
        boundaries.append((left_start0, left_end0, first, "complete"))
        left_start0, left_end0 = left_start0 - period_bp, left_start0
    for left, right in zip(anchor_hits, anchor_hits[1:]):
        if right.start0 <= left.start0:
            continue
        status = "junction_partial" if junction is not None and left.strand == "+" and right.strand == "-" else "complete"
        boundaries.append((left.start0, right.start0, left, status))
    right_start0, right_end0 = last.start0, last.start0 + period_bp
    while right_start0 < region_end0:
        boundaries.append((right_start0, right_end0, last, "complete"))
        right_start0, right_end0 = right_end0, right_end0 + period_bp

    units: list[Unit] = []
    for start0, end0, anchor, status in sorted(boundaries, key=lambda boundary: (boundary[0], boundary[1])):
        clipped_start0, clipped_end0 = max(start0, region_start0), min(end0, region_end0)
        if clipped_start0 >= clipped_end0:
            continue
        if status == "complete" and start0 < region_start0:
            status = "partial_left"
        elif status == "complete" and end0 > region_end0:
            status = "partial_right"
        arm = "unassigned"
        if junction is not None:
            if status == "junction_partial":
                arm = "junction"
            elif end0 <= junction:
                arm = "left"
            elif start0 >= junction:
                arm = "right"
        units.append(Unit(f"unit_{len(units) + 1:03d}", region.chrom, clipped_start0, clipped_end0, anchor.strand, arm, status))
    return units


def build_stainedglass_units(period: int | float | dict[str, object], region: Region) -> list[Unit]:
    """Tile a user-defined repeat-rich region into period-sized units."""
    period_bp = math.floor(_period_value(period) + 0.5)
    region_start0 = region.start1 - 1
    units: list[Unit] = []
    for start0 in range(region_start0, region.end1, period_bp):
        end0 = min(start0 + period_bp, region.end1)
        status = "complete" if end0 - start0 == period_bp else "partial_right"
        units.append(
            Unit(
                f"unit_{len(units) + 1:03d}",
                region.chrom,
                start0,
                end0,
                "+",
                "forward",
                status,
            )
        )
    return units


def clip_units_to_region(units: list[Unit], region: Region) -> list[Unit]:
    """Clip refined units to the exact target interval without filling flanks."""
    region_start0, region_end0 = region.start1 - 1, region.end1
    clipped: list[Unit] = []
    for unit in sorted(units, key=lambda item: (item.start0, item.end0, item.unit_id)):
        if unit.chrom != region.chrom:
            raise ValueError("unit chromosome does not match target region")
        start0 = max(region_start0, unit.start0)
        end0 = min(region_end0, unit.end0)
        if end0 <= start0:
            continue
        partial_left = unit.start0 < region_start0 or unit.status in {"partial_left", "partial_both"}
        partial_right = unit.end0 > region_end0 or unit.status in {"partial_right", "partial_both"}
        if partial_left and partial_right:
            status = "partial_both"
        elif partial_left:
            status = "partial_left"
        elif partial_right:
            status = "partial_right"
        else:
            status = unit.status
        clipped.append(
            Unit(
                unit.unit_id,
                unit.chrom,
                start0,
                end0,
                unit.strand,
                unit.arm,
                status,
            )
        )
    return clipped


def summarize_refined_unit_lengths(units: list[Unit]) -> dict[str, int | float | str]:
    """Summarize complete units, falling back safely when every unit is partial."""
    if not units:
        raise ValueError("cannot summarize an empty unit set")
    complete_lengths = [
        unit.end0 - unit.start0 for unit in units if unit.status == "complete"
    ]
    basis = "complete_units" if complete_lengths else "all_units"
    lengths = complete_lengths or [unit.end0 - unit.start0 for unit in units]
    return {
        "basis": basis,
        "count": len(lengths),
        "minimum_bp": min(lengths),
        "maximum_bp": max(lengths),
        "median_bp": statistics.median(lengths),
    }


def reconcile_boundary_rows_with_units(
    rows: list[dict[str, object]], units: list[Unit], region: Region
) -> list[dict[str, object]]:
    """Make reported final boundaries agree with target-clipped unit endpoints."""
    units_by_id = {unit.unit_id: unit for unit in units}
    region_start0, region_end0 = region.start1 - 1, region.end1
    reconciled: list[dict[str, object]] = []
    for source in rows:
        row = dict(source)
        unit = units_by_id.get(str(row.get("unit_id", "")))
        role = str(row.get("boundary_role", ""))
        current = int(row["refined_boundary0"])
        if unit is not None and role == "start":
            final = unit.start0
        elif unit is not None and role == "end":
            final = unit.end0
        elif unit is not None and role == "shared":
            final = unit.start0
        else:
            final = min(max(current, region_start0), region_end0)
        if final != current:
            row["refined_boundary0"] = final
            row["adjustment_bp"] = final - int(row["coarse_boundary0"])
            suffix = "target_clip" if final in {region_start0, region_end0} else "unit_reconcile"
            method = str(row.get("method", ""))
            if suffix not in method.split("+"):
                row["method"] = f"{method}+{suffix}" if method else suffix
            if suffix == "target_clip":
                row["confidence"] = "low"
        reconciled.append(row)
    return reconciled


def discover_stainedglass_unit_blocks(
    hits: list[StainedGlassHit],
    region: Region,
    *,
    expected_period_bp: int | float | dict[str, object],
    minimum_chain_windows: int = 6,
    high_identity: float = 95.0,
    medium_identity: float = 80.0,
) -> tuple[list[Unit], list[dict[str, object]], dict[str, object]]:
    """Detect discrete homologous blocks without tiling the whole target region.

    Long, constant-diagonal StainedGlass chains vote for block endpoints. Windows
    with unusually many distant partners are retained as gap evidence but cannot
    create a unit unless a multi-window chain also spans the interval.
    """
    if not hits or minimum_chain_windows < 2 or medium_identity > high_identity:
        raise ValueError("invalid StainedGlass block-detection inputs")
    spans = [hit.query_end0 - hit.query_start0 for hit in hits if hit.query_end0 > hit.query_start0]
    if not spans:
        raise ValueError("StainedGlass hits do not contain positive-width windows")
    window_size = int(statistics.median(spans))
    period_bp = _period_value(expected_period_bp)
    period_windows = max(2, int(round(period_bp / window_size)))
    chain_minimum = min(
        minimum_chain_windows,
        max(2, int(math.ceil(period_windows * 0.45))),
    )
    region_start0, region_end0 = region.start1 - 1, region.end1

    pair_identity: dict[tuple[int, int, str], float] = {}
    for hit in hits:
        if hit.identity < medium_identity or not stainedglass_hit_overlaps_region(hit, region):
            continue
        query_window = hit.query_start0 // window_size
        reference_window = hit.reference_start0 // window_size
        if query_window == reference_window:
            continue
        key = (query_window, reference_window, hit.strand)
        pair_identity[key] = max(pair_identity.get(key, -math.inf), hit.identity)
    if not pair_identity:
        raise ValueError("no off-diagonal StainedGlass chains overlap the target region")

    high_degree: Counter[int] = Counter()
    for (query_window, _reference_window, _strand), identity in pair_identity.items():
        if identity >= high_identity:
            high_degree[query_window] += 1
    active_degrees = sorted(value for value in high_degree.values() if value > 0)
    median_degree = statistics.median(active_degrees) if active_degrees else 0.0
    hub_threshold = max(5, int(math.ceil(median_degree * 1.75)))
    hub_windows = {window for window, degree in high_degree.items() if degree >= hub_threshold}

    signature_diagonals: dict[int, dict[str, set[int]]] = {}
    for (query_window, reference_window, strand), identity in pair_identity.items():
        if identity < max(medium_identity, high_identity - 5.0):
            continue
        intercept = reference_window - query_window if strand == "+" else reference_window + query_window
        signature_diagonals.setdefault(query_window, {"+": set(), "-": set()})[strand].add(intercept)

    diagonals: dict[tuple[str, int], list[tuple[int, int, float]]] = {}
    for (query_window, reference_window, strand), identity in pair_identity.items():
        if query_window in hub_windows or reference_window in hub_windows:
            continue
        intercept = reference_window - query_window if strand == "+" else reference_window + query_window
        diagonals.setdefault((strand, intercept), []).append((query_window, reference_window, identity))

    endpoint_votes: Counter[int] = Counter()
    chain_coverage: Counter[int] = Counter()
    coherent_chains = 0
    maximum_step = 2
    for points in diagonals.values():
        ordered = sorted(points, key=lambda point: (point[0], point[1]))
        runs: list[list[tuple[int, int, float]]] = []
        current: list[tuple[int, int, float]] = []
        for point in ordered:
            if current and point[0] - current[-1][0] > maximum_step:
                runs.append(current)
                current = []
            current.append(point)
        if current:
            runs.append(current)
        for run in runs:
            span_windows = run[-1][0] - run[0][0] + 1
            high_windows = sum(identity >= high_identity for _query, _reference, identity in run)
            nonhub_windows = sum(query not in hub_windows for query, _reference, _identity in run)
            if (
                span_windows < chain_minimum
                or high_windows < max(2, chain_minimum // 2)
                or nonhub_windows < 1
            ):
                continue
            coherent_chains += 1
            reference_windows = [reference for _query, reference, _identity in run]
            weight = min(span_windows, period_windows)
            for endpoint in (
                run[0][0],
                run[-1][0] + 1,
                min(reference_windows),
                max(reference_windows) + 1,
            ):
                endpoint_votes[endpoint] += weight
            for window in range(run[0][0], run[-1][0] + 1):
                chain_coverage[window] += 1
            for window in range(min(reference_windows), max(reference_windows) + 1):
                chain_coverage[window] += 1

    if not endpoint_votes:
        raise ValueError("no coherent multi-window StainedGlass chains were detected")

    region_start_window = region_start0 // window_size
    region_end_window = math.ceil(region_end0 / window_size)
    active_windows = sorted(chain_coverage)
    if not active_windows:
        raise ValueError("coherent chains did not cover any target windows")
    array_start = max(region_start_window, active_windows[0])
    array_end = min(region_end_window, active_windows[-1] + 1)

    def jaccard(left: set[int], right: set[int]) -> float:
        union = left | right
        return len(left & right) / len(union) if union else 1.0

    smoothed_votes = {
        window: sum(endpoint_votes[index] for index in range(window - 2, window + 3))
        for window in range(array_start, array_end + 1)
    }
    maximum_vote = max(smoothed_votes.values(), default=1) or 1
    boundary_scores: dict[int, float] = {}
    for window in range(array_start, array_end + 1):
        previous = signature_diagonals.get(window - 1, {"+": set(), "-": set()})
        current = signature_diagonals.get(window, {"+": set(), "-": set()})
        discontinuity = 1.0 - (
            jaccard(previous["+"], current["+"]) + jaccard(previous["-"], current["-"])
        ) / 2.0
        boundary_scores[window] = 2.0 * smoothed_votes[window] / maximum_vote + discontinuity

    minimum_hub_run = max(2, int(math.ceil(period_windows * 0.18)))
    hub_runs: list[tuple[int, int]] = []
    current_hubs: list[int] = []
    for window in range(array_start, array_end):
        if window in hub_windows:
            current_hubs.append(window)
        elif current_hubs:
            if len(current_hubs) >= minimum_hub_run:
                hub_runs.append((current_hubs[0], current_hubs[-1] + 1))
            current_hubs = []
    if current_hubs and len(current_hubs) >= minimum_hub_run:
        hub_runs.append((current_hubs[0], current_hubs[-1] + 1))

    # A gap without a high-copy TE plateau still remains a gap when coherent
    # chains stop completely. Use this fallback only when no TE-hub run exists,
    # avoiding proliferation of small gaps inside repeat-rich units.
    gap_cores: list[tuple[int, int, str]] = [
        (start, end, "repeat_rich_gap") for start, end in hub_runs
    ]
    if not gap_cores:
        inactive_runs: list[int] = []
        minimum_inactive = max(2, int(math.ceil(period_windows * 0.20)))
        for window in range(array_start, array_end):
            if chain_coverage[window] == 0:
                inactive_runs.append(window)
            elif inactive_runs:
                if len(inactive_runs) >= minimum_inactive:
                    gap_cores.append((inactive_runs[0], inactive_runs[-1] + 1, "inter_unit_gap"))
                inactive_runs = []
        if inactive_runs and len(inactive_runs) >= minimum_inactive:
            gap_cores.append((inactive_runs[0], inactive_runs[-1] + 1, "inter_unit_gap"))

    gap_intervals: list[tuple[int, int, str]] = []
    search_left = max(1, int(math.ceil(period_windows * 0.25)))
    search_right = max(1, int(math.ceil(period_windows * 0.10)))
    for core_start, core_end, interval_type in gap_cores:
        if interval_type == "inter_unit_gap":
            left, right = core_start, core_end
        else:
            left_candidates = range(max(array_start, core_start - search_left), core_start + 1)
            left = max(left_candidates, key=lambda window: (boundary_scores.get(window, 0.0), window))
            right_low = max(left + 1, core_end - search_left)
            right_high = min(array_end, core_end + search_right)
            right_candidates = range(right_low, right_high + 1)
            right = max(right_candidates, key=lambda window: (boundary_scores.get(window, 0.0), -window))
        if gap_intervals and left <= gap_intervals[-1][1]:
            previous_start, previous_end, previous_type = gap_intervals.pop()
            gap_intervals.append((previous_start, max(previous_end, right), previous_type))
        else:
            gap_intervals.append((left, right, interval_type))

    minimum_unit_windows = max(2, int(math.ceil(period_windows * 0.60)))
    maximum_unit_windows = max(minimum_unit_windows, int(math.floor(period_windows * 1.40)))

    def partition_segment(start: int, end: int) -> list[tuple[int, int]]:
        length = end - start
        if length < minimum_unit_windows:
            return []
        minimum_count = max(1, int(math.ceil(length / maximum_unit_windows)))
        maximum_count = max(1, int(math.floor(length / minimum_unit_windows)))
        rounded_count = max(1, int(math.floor(length / period_windows + 0.5)))
        unit_count = min(max(rounded_count, minimum_count), maximum_count)
        states: dict[int, tuple[float, list[int]]] = {start: (0.0, [start])}
        for unit_index in range(1, unit_count + 1):
            next_states: dict[int, tuple[float, list[int]]] = {}
            remaining = unit_count - unit_index
            for previous, (score, path) in states.items():
                low = max(
                    previous + minimum_unit_windows,
                    end - remaining * maximum_unit_windows,
                )
                high = min(
                    previous + maximum_unit_windows,
                    end - remaining * minimum_unit_windows,
                )
                for boundary in range(low, high + 1):
                    if unit_index == unit_count and boundary != end:
                        continue
                    unit_length = boundary - previous
                    length_score = 1.5 * math.exp(
                        -((unit_length - period_windows) / max(1.0, 0.28 * period_windows)) ** 2
                    )
                    endpoint_score = 0.0 if boundary == end else boundary_scores.get(boundary, 0.0)
                    candidate = (score + length_score + endpoint_score, path + [boundary])
                    if boundary not in next_states or candidate[0] > next_states[boundary][0]:
                        next_states[boundary] = candidate
            states = next_states
        if end not in states:
            return []
        path = states[end][1]
        return list(zip(path, path[1:]))

    segments: list[tuple[int, int]] = []
    segment_start = array_start
    for gap_start, gap_end, _interval_type in gap_intervals:
        if segment_start < gap_start:
            segments.append((segment_start, gap_start))
        segment_start = gap_end
    if segment_start < array_end:
        segments.append((segment_start, array_end))

    units: list[Unit] = []
    gap_rows: list[dict[str, object]] = []

    def add_gap(start0: int, end0: int, interval_type: str) -> None:
        start0 = max(region_start0, start0)
        end0 = min(region_end0, end0)
        if end0 <= start0:
            return
        first_window = start0 // window_size
        last_window = math.ceil(end0 / window_size)
        hubs = sum(window in hub_windows for window in range(first_window, last_window))
        if interval_type == "inter_unit_gap" and hubs:
            interval_type = "repeat_rich_gap"
        gap_rows.append(
            {
                "gap_id": f"gap_{len(gap_rows) + 1:03d}",
                "chrom": region.chrom,
                "start0": start0,
                "end0": end0,
                "length_bp": end0 - start0,
                "interval_type": interval_type,
                "hub_windows": hubs,
            }
        )

    add_gap(region_start0, array_start * window_size, "left_flank")
    selected_boundaries: set[int] = set()
    for segment_start, segment_end in segments:
        for left_window, right_window in partition_segment(segment_start, segment_end):
            raw_left, raw_right = left_window * window_size, right_window * window_size
            left = max(region_start0, raw_left)
            right = min(region_end0, raw_right)
            if right <= left:
                continue
            if raw_left < region_start0 and raw_right > region_end0:
                status = "partial_both"
            elif raw_left < region_start0:
                status = "partial_left"
            elif raw_right > region_end0:
                status = "partial_right"
            else:
                status = "complete"
            units.append(
                Unit(
                    f"unit_{len(units) + 1:03d}",
                    region.chrom,
                    left,
                    right,
                    "+",
                    "unassigned",
                    status,
                )
            )
            selected_boundaries.update((left, right))
    for gap_start, gap_end, interval_type in gap_intervals:
        start0 = max(region_start0, gap_start * window_size)
        end0 = min(region_end0, gap_end * window_size)
        add_gap(start0, end0, interval_type)
        if end0 > start0:
            selected_boundaries.update((start0, end0))
    add_gap(array_end * window_size, region_end0, "right_flank")
    units.sort(key=lambda unit: (unit.start0, unit.end0))
    for index, unit in enumerate(units, 1):
        unit.unit_id = f"unit_{index:03d}"
    gap_rows.sort(key=lambda row: (int(row["start0"]), int(row["end0"])))
    for index, row in enumerate(gap_rows, 1):
        row["gap_id"] = f"gap_{index:03d}"
    return units, gap_rows, {
        "window_size": window_size,
        "period_windows": period_windows,
        "coherent_chains": coherent_chains,
        "repeat_hub_windows": len(hub_windows),
        "hub_degree_threshold": hub_threshold,
        "candidate_boundaries": sorted(selected_boundaries),
        "repeat_hub_runs": [(start * window_size, end * window_size) for start, end in hub_runs],
        "whole_region_tiled": False,
    }


def _refine_discrete_stainedglass_units(
    sequence: str,
    fasta_region: Region,
    coarse_units: list[Unit],
    repeats: list[Feature],
    coding: list[Feature],
    pseudogenes: list[Feature],
    *,
    k: int,
    sample_step: int,
    minimum_peak_support: int,
) -> tuple[list[Unit], list[dict[str, object]]]:
    """Refine true gap edges independently while retaining contiguous subarrays."""
    groups: list[list[Unit]] = []
    for unit in coarse_units:
        if not groups or groups[-1][-1].end0 < unit.start0:
            groups.append([])
        groups[-1].append(unit)

    refined_units: list[Unit] = []
    boundary_rows: list[dict[str, object]] = []
    for group in groups:
        global_offset = len(refined_units)
        group_units, group_rows = refine_stainedglass_units(
            sequence,
            fasta_region,
            group,
            repeats,
            coding,
            pseudogenes,
            k=k,
            sample_step=sample_step,
            minimum_peak_support=minimum_peak_support,
        )
        refined_units.extend(group_units)
        for local_boundary, row in enumerate(group_rows, 1):
            annotated = dict(row)
            annotated["unit_id"] = f"unit_{global_offset + local_boundary + 1:03d}"
            annotated["boundary_role"] = "shared"
            boundary_rows.append(annotated)
    for index, unit in enumerate(refined_units, 1):
        unit.unit_id = f"unit_{index:03d}"

    region_start0 = fasta_region.start1 - 1
    coarse_period = int(statistics.median(unit.end0 - unit.start0 for unit in coarse_units))
    # StainedGlass already constrains a gap edge to one 10-kb-scale window.
    # Keep sequence/annotation refinement local so a recurrent TE inside the gap
    # cannot pull the endpoint across a substantial fraction of the unit.
    search_radius = max(k, int(round(coarse_period * 0.15)))
    support_span = max(k, min(2000, int(round(coarse_period * 0.05))))

    kmer_occurrences: dict[str, list[int]] = {}
    for relative_start0 in range(0, len(sequence) - k + 1, sample_step):
        kmer = sequence[relative_start0 : relative_start0 + k].upper()
        if not set(kmer) <= {"A", "C", "G", "T"}:
            continue
        canonical = min(kmer, revcomp(kmer))
        kmer_occurrences.setdefault(canonical, []).append(region_start0 + relative_start0)
    maximum_occurrences = max(80, len(coarse_units) * 8)
    position_support: dict[int, int] = {}
    for positions in kmer_occurrences.values():
        if len(positions) > maximum_occurrences:
            continue
        for position0 in positions:
            position_support[position0] = max(position_support.get(position0, 0), len(positions))

    feature_edges = {
        position0
        for feature in [
            *repeats,
            *(item for item in coding if item.feature_type.lower() == "gene"),
            *(item for item in pseudogenes if item.feature_type.lower() == "pseudogene"),
        ]
        for position0 in (feature.start0, feature.end0)
    }

    def average_copy_support(start0: int, end0: int) -> float:
        positions = range(start0, max(start0, end0), sample_step)
        values = [
            min(1.0, position_support.get(position0, 0) / minimum_peak_support)
            for position0 in positions
        ]
        return sum(values) / len(values) if values else 0.0

    def refine_transition(coarse0: int, low0: int, high0: int, role: str) -> tuple[int, float, bool]:
        candidates = set(range(low0, high0 + 1, sample_step))
        candidates.update(edge for edge in feature_edges if low0 <= edge <= high0)
        best: tuple[float, int, bool] | None = None
        for candidate0 in sorted(candidates):
            left_support = average_copy_support(max(region_start0, candidate0 - support_span), candidate0)
            right_support = average_copy_support(candidate0, min(fasta_region.end1, candidate0 + support_span))
            transition = left_support - right_support if role == "end" else right_support - left_support
            on_annotation_edge = candidate0 in feature_edges
            score = (
                transition
                + (0.75 if on_annotation_edge else 0.0)
                - 0.30 * abs(candidate0 - coarse0) / max(1, search_radius)
            )
            item = (score, candidate0, on_annotation_edge)
            if best is None or item[0] > best[0] or (
                item[0] == best[0] and abs(item[1] - coarse0) < abs(best[1] - coarse0)
            ):
                best = item
        if best is None:
            return coarse0, 0.0, False
        return best[1], best[0], best[2]

    for index in range(len(refined_units) - 1):
        left_source = coarse_units[index]
        right_source = coarse_units[index + 1]
        if left_source.end0 >= right_source.start0:
            continue
        left = refined_units[index]
        right = refined_units[index + 1]
        refined_end0, end_support, end_annotation = refine_transition(
            left_source.end0,
            max(left.start0 + k, left_source.end0 - search_radius),
            min(right_source.start0 - k, left_source.end0 + search_radius),
            "end",
        )
        refined_start0, start_support, start_annotation = refine_transition(
            right_source.start0,
            max(refined_end0 + k, right_source.start0 - search_radius),
            min(right.end0 - k, right_source.start0 + search_radius),
            "start",
        )
        refined_units[index] = Unit(
            left.unit_id, left.chrom, left.start0, refined_end0, left.strand, left.arm, left.status
        )
        refined_units[index + 1] = Unit(
            right.unit_id, right.chrom, refined_start0, right.end0, right.strand, right.arm, right.status
        )
        for unit, role, coarse0, refined0, support, annotation in (
            (left, "end", left_source.end0, refined_end0, end_support, end_annotation),
            (right, "start", right_source.start0, refined_start0, start_support, start_annotation),
        ):
            boundary_rows.append(
                {
                    "boundary_id": f"{unit.unit_id}_{role}",
                    "unit_id": unit.unit_id,
                    "boundary_role": role,
                    "coarse_boundary0": coarse0,
                    "sequence_target0": refined0,
                    "refined_boundary0": refined0,
                    "adjustment_bp": refined0 - coarse0,
                    "method": "copy_support_transition" + ("+annotation_edge" if annotation else ""),
                    "sequence_support": support,
                    "forward_lag_bp": "",
                    "forward_support": "",
                    "reverse_center0": "",
                    "reverse_support": "",
                    "annotation_signature": "feature_edge" if annotation else "",
                    "annotation_recurrence": "",
                    "annotations_crossed": sum(
                        feature.start0 < refined0 < feature.end0
                        for feature in [*repeats, *coding, *pseudogenes]
                    ),
                    "confidence": "high" if annotation and support > 0 else "medium" if support > 0 else "low",
                }
            )

    first_source, first = coarse_units[0], refined_units[0]
    if first_source.start0 > region_start0:
        refined_start0, support, annotation = refine_transition(
            first_source.start0,
            max(region_start0, first_source.start0 - search_radius),
            min(first.end0 - k, first_source.start0 + search_radius),
            "start",
        )
        refined_units[0] = Unit(
            first.unit_id, first.chrom, refined_start0, first.end0, first.strand, first.arm, first.status
        )
        boundary_rows.append(
            {
                "boundary_id": f"{first.unit_id}_start",
                "unit_id": first.unit_id,
                "boundary_role": "start",
                "coarse_boundary0": first_source.start0,
                "sequence_target0": refined_start0,
                "refined_boundary0": refined_start0,
                "adjustment_bp": refined_start0 - first_source.start0,
                "method": "copy_support_transition" + ("+annotation_edge" if annotation else ""),
                "sequence_support": support,
                "forward_lag_bp": "",
                "forward_support": "",
                "reverse_center0": "",
                "reverse_support": "",
                "annotation_signature": "feature_edge" if annotation else "",
                "annotation_recurrence": "",
                "annotations_crossed": 0,
                "confidence": "high" if annotation and support > 0 else "medium" if support > 0 else "low",
            }
        )

    last_source, last = coarse_units[-1], refined_units[-1]
    if last_source.end0 < fasta_region.end1:
        refined_end0, support, annotation = refine_transition(
            last_source.end0,
            max(last.start0 + k, last_source.end0 - search_radius),
            min(fasta_region.end1, last_source.end0 + search_radius),
            "end",
        )
        refined_units[-1] = Unit(
            last.unit_id, last.chrom, last.start0, refined_end0, last.strand, last.arm, last.status
        )
        boundary_rows.append(
            {
                "boundary_id": f"{last.unit_id}_end",
                "unit_id": last.unit_id,
                "boundary_role": "end",
                "coarse_boundary0": last_source.end0,
                "sequence_target0": refined_end0,
                "refined_boundary0": refined_end0,
                "adjustment_bp": refined_end0 - last_source.end0,
                "method": "copy_support_transition" + ("+annotation_edge" if annotation else ""),
                "sequence_support": support,
                "forward_lag_bp": "",
                "forward_support": "",
                "reverse_center0": "",
                "reverse_support": "",
                "annotation_signature": "feature_edge" if annotation else "",
                "annotation_recurrence": "",
                "annotations_crossed": 0,
                "confidence": "high" if annotation and support > 0 else "medium" if support > 0 else "low",
            }
        )
    boundary_rows.sort(key=lambda row: (int(row["refined_boundary0"]), str(row["boundary_id"])))
    for boundary_index, row in enumerate(boundary_rows, 1):
        row["boundary_id"] = f"boundary_{boundary_index:03d}"
    return refined_units, boundary_rows


def refine_stainedglass_units(
    sequence: str,
    fasta_region: Region,
    coarse_units: list[Unit],
    repeats: list[Feature],
    coding: list[Feature],
    pseudogenes: list[Feature],
    *,
    k: int = 21,
    sample_step: int = 20,
    minimum_peak_support: int = 10,
) -> tuple[list[Unit], list[dict[str, object]]]:
    """Refine coarse StainedGlass boundaries with exact k-mer diagonals and annotation edges.

    StainedGlass defines only the coarse phase. Forward exact-kmer diagonals estimate
    local copy length, reverse-complement anti-diagonals estimate inversion centers,
    and recurrent repeat/gene edges provide base-level snapping and low-confidence
    fallback positions.
    """
    if len(coarse_units) < 1 or k < 3 or sample_step < 1 or minimum_peak_support < 1:
        raise ValueError("invalid boundary-refinement inputs")
    if any(left.end0 < right.start0 for left, right in zip(coarse_units, coarse_units[1:])):
        return _refine_discrete_stainedglass_units(
            sequence,
            fasta_region,
            coarse_units,
            repeats,
            coding,
            pseudogenes,
            k=k,
            sample_step=sample_step,
            minimum_peak_support=minimum_peak_support,
        )
    coarse_period = int(statistics.median(unit.end0 - unit.start0 for unit in coarse_units[:-1] or coarse_units))
    if coarse_period < k:
        raise ValueError("coarse unit period is shorter than the refinement k-mer")
    window_half = max(k, int(round(coarse_period * 0.65)))
    search_radius = max(sample_step, int(round(coarse_period * 0.15)))
    lag_low = int(round(coarse_period * 0.75))
    lag_high = int(round(coarse_period * 1.25))
    region_start0 = fasta_region.start1 - 1

    kmer_index: dict[str, list[tuple[int, bool]]] = {}
    for relative_start0 in range(0, len(sequence) - k + 1, sample_step):
        kmer = sequence[relative_start0 : relative_start0 + k].upper()
        if not set(kmer) <= {"A", "C", "G", "T"}:
            continue
        reverse = revcomp(kmer)
        if kmer == reverse:
            continue
        canonical = min(kmer, reverse)
        kmer_index.setdefault(canonical, []).append((region_start0 + relative_start0, kmer <= reverse))

    features = [
        *repeats,
        *(feature for feature in coding if feature.feature_type.lower() == "gene"),
        *(feature for feature in pseudogenes if feature.feature_type.lower() == "pseudogene"),
    ]
    coarse_boundaries = [unit.start0 for unit in coarse_units[1:]]
    signature_coverage: dict[tuple[str, str, str], set[int]] = {}
    signature_distances: dict[tuple[str, str, str], dict[int, int]] = {}
    edges_by_boundary: dict[
        int, list[tuple[int, tuple[str, str, str], tuple[str, str, str]]]
    ] = {}
    for boundary_index, boundary0 in enumerate(coarse_boundaries, 1):
        nearby: list[
            tuple[int, tuple[str, str, str], tuple[str, str, str]]
        ] = []
        phase_unit = coarse_units[boundary_index]
        for feature in features:
            symbol = feature.attributes.get("Name") or feature.attributes.get("gene_symbol") or feature.feature_type
            for side, position0 in (("start", feature.start0), ("end", feature.end0)):
                if abs(position0 - boundary0) <= search_radius:
                    raw_signature = (feature.feature_type, symbol, side)
                    oriented_side = (
                        side
                        if phase_unit.strand == "+"
                        else "end" if side == "start" else "start"
                    )
                    phase_signature = (
                        feature.feature_type,
                        symbol,
                        f"oriented_{oriented_side}",
                    )
                    nearby.append((position0, raw_signature, phase_signature))
                    signature_coverage.setdefault(phase_signature, set()).add(boundary_index)
                    distances = signature_distances.setdefault(phase_signature, {})
                    distances[boundary_index] = min(
                        distances.get(boundary_index, math.inf),
                        abs(position0 - boundary0),
                    )
        edges_by_boundary[boundary_index] = nearby

    # A tandem array has one biological phase, but several recurrent feature
    # edges can lie near every coarse grid point.  Without a cross-boundary
    # constraint, independently choosing the nearest edge can alternate between
    # the start and end of the same feature and shift every downstream unit.
    # Lock weak-evidence boundaries to one broadly recurrent, generic feature
    # signature.  Exact sequence support can still override this fallback.
    minimum_phase_anchor_coverage = max(
        3, int(math.ceil(len(coarse_boundaries) * 0.50))
    )

    def phase_anchor_reliability(signature: tuple[str, str, str]) -> int:
        feature_type = signature[0].lower()
        if feature_type in {"gene", "pseudogene"} or feature_type.endswith("_gene"):
            return 2
        if feature_type in {"simple_repeat", "low_complexity", "satellite"}:
            return 0
        return 1

    phase_anchor_signature = max(
        (
            signature
            for signature, covered_boundaries in signature_coverage.items()
            if len(covered_boundaries) >= minimum_phase_anchor_coverage
        ),
        key=lambda signature: (
            len(signature_coverage[signature]),
            phase_anchor_reliability(signature),
            -statistics.median(signature_distances[signature].values()),
            -signature_distances[signature].get(1, math.inf),
            signature,
        ),
        default=None,
    )

    refined_boundaries = [coarse_units[0].start0]
    boundary_rows: list[dict[str, object]] = []
    max_occurrences = 80
    for boundary_index, coarse_boundary0 in enumerate(coarse_boundaries, 1):
        left_low = max(fasta_region.start1 - 1, coarse_boundary0 - window_half)
        right_high = min(fasta_region.end1, coarse_boundary0 + window_half)
        forward_lags: Counter[int] = Counter()
        reverse_centers: Counter[int] = Counter()
        for occurrences in kmer_index.values():
            if len(occurrences) > max_occurrences:
                continue
            left = [item for item in occurrences if left_low <= item[0] < coarse_boundary0]
            right = [item for item in occurrences if coarse_boundary0 <= item[0] < right_high]
            for left_start0, left_orientation in left:
                for right_start0, right_orientation in right:
                    if left_orientation == right_orientation:
                        lag = right_start0 - left_start0
                        if lag_low <= lag <= lag_high:
                            forward_lags[lag] += 1
                    else:
                        center0 = int(round((left_start0 + right_start0 + k) / 2))
                        if abs(center0 - coarse_boundary0) <= search_radius:
                            reverse_centers[center0] += 1

        forward_lag, forward_support = min(
            forward_lags.items(),
            key=lambda item: (-item[1], abs(item[0] - coarse_period), item[0]),
            default=(coarse_period, 0),
        )
        reverse_center0, reverse_support = min(
            reverse_centers.items(),
            key=lambda item: (-item[1], abs(item[0] - coarse_boundary0), item[0]),
            default=(coarse_boundary0, 0),
        )
        if reverse_support >= minimum_peak_support and reverse_support > forward_support:
            sequence_target0 = reverse_center0
            method = "reverse_kmer_center"
            sequence_support = reverse_support
        elif forward_support >= minimum_peak_support:
            propagated = refined_boundaries[-1] + forward_lag
            coarse_based = coarse_boundary0 - coarse_period + forward_lag
            sequence_target0 = (
                propagated
                if abs(propagated - coarse_boundary0) <= search_radius
                else coarse_based
            )
            method = "forward_kmer_lag"
            sequence_support = forward_support
        else:
            sequence_target0 = coarse_boundary0
            method = "coarse_only"
            sequence_support = max(forward_support, reverse_support)

        strong_sequence = sequence_support >= minimum_peak_support
        snap_radius = max(sample_step, min(2000 if strong_sequence else 10000, search_radius))
        candidates: list[
            tuple[
                float,
                int,
                tuple[str, str, str] | None,
                tuple[str, str, str] | None,
                int,
            ]
        ] = []
        candidates.append((0.0, sequence_target0, None, None, 0))
        for position0, raw_signature, phase_signature in edges_by_boundary.get(boundary_index, []):
            distance = abs(position0 - sequence_target0)
            if distance > snap_radius:
                continue
            recurrence = len(signature_coverage.get(phase_signature, set()))
            crossing = sum(feature.start0 < position0 < feature.end0 for feature in features)
            distance_weight = 5.0 if strong_sequence else 1.0
            score = distance_weight * distance / snap_radius + 0.08 * crossing - 2.0 * recurrence / max(1, len(coarse_boundaries))
            if phase_anchor_signature is not None and phase_signature != phase_anchor_signature:
                score += 1.0 if strong_sequence else 3.0
            candidates.append(
                (score, position0, raw_signature, phase_signature, recurrence)
            )
        (
            _score,
            refined_boundary0,
            signature,
            selected_phase_signature,
            annotation_recurrence,
        ) = min(
            candidates,
            key=lambda item: (item[0], abs(item[1] - sequence_target0), item[1]),
        )
        lower = refined_boundaries[-1] + max(k, coarse_period // 2)
        upper = min(fasta_region.end1 - k, coarse_boundary0 + search_radius)
        if lower >= upper or not lower <= refined_boundary0 <= upper:
            refined_boundary0 = min(max(sequence_target0, lower), max(lower, upper))
            signature = None
            selected_phase_signature = None
            annotation_recurrence = 0
        if signature is not None:
            edge_method = (
                "phase_anchor_edge"
                if selected_phase_signature == phase_anchor_signature
                else "repeat_edge"
            )
            method = (
                method + f"+{edge_method}"
                if method != "coarse_only"
                else edge_method
            )
        confidence = (
            "high"
            if sequence_support >= 25
            else "medium"
            if strong_sequence or annotation_recurrence >= max(3, len(coarse_boundaries) // 2)
            else "low"
        )
        annotations_crossed = sum(feature.start0 < refined_boundary0 < feature.end0 for feature in features)
        refined_boundaries.append(refined_boundary0)
        boundary_rows.append(
            {
                "boundary_id": f"boundary_{boundary_index:03d}",
                "coarse_boundary0": coarse_boundary0,
                "sequence_target0": sequence_target0,
                "refined_boundary0": refined_boundary0,
                "adjustment_bp": refined_boundary0 - coarse_boundary0,
                "method": method,
                "sequence_support": sequence_support,
                "forward_lag_bp": forward_lag if forward_support else "",
                "forward_support": forward_support,
                "reverse_center0": reverse_center0 if reverse_support else "",
                "reverse_support": reverse_support,
                "annotation_signature": "|".join(signature) if signature else "",
                "phase_anchor_signature": "|".join(phase_anchor_signature) if phase_anchor_signature else "",
                "annotation_recurrence": annotation_recurrence,
                "annotations_crossed": annotations_crossed,
                "confidence": confidence,
            }
        )
    refined_boundaries.append(coarse_units[-1].end0)

    refined_units: list[Unit] = []
    for index, (start0, end0) in enumerate(zip(refined_boundaries, refined_boundaries[1:])):
        if end0 <= start0:
            raise RuntimeError("boundary refinement produced non-increasing units")
        source = coarse_units[min(index, len(coarse_units) - 1)]
        refined_units.append(
            Unit(
                f"unit_{index + 1:03d}",
                source.chrom,
                start0,
                end0,
                source.strand,
                source.arm,
                source.status,
            )
        )
    return refined_units, boundary_rows


def infer_stainedglass_unit_orientations(
    hits: list[StainedGlassHit], units: list[Unit], min_identity: float = 95.0
) -> tuple[dict[str, str], dict[str, object]]:
    """Propagate same/opposite StainedGlass strand relations across a unit graph."""
    if not units:
        return {}, {"disconnected_units": 0, "reference_unit_id": "", "unit_links": {}}

    def unit_index(start0: int, end0: int) -> int | None:
        midpoint = (start0 + end0 - 1) // 2
        for index, unit in enumerate(units):
            if unit.start0 <= midpoint < unit.end0:
                return index
        return None

    relations: dict[tuple[int, int], Counter[str]] = {}
    for hit in hits:
        if hit.identity < min_identity:
            continue
        left = unit_index(hit.query_start0, hit.query_end0)
        right = unit_index(hit.reference_start0, hit.reference_end0)
        if left is None or right is None or left == right:
            continue
        pair = tuple(sorted((left, right)))
        relations.setdefault(pair, Counter())[hit.strand] += 1

    graph: dict[int, list[tuple[int, int, int]]] = {index: [] for index in range(len(units))}
    constraints: list[tuple[int, int, int, int]] = []
    unit_links = Counter()
    for (left, right), counts in relations.items():
        relation = 0 if counts["+"] >= counts["-"] else 1
        support = abs(counts["+"] - counts["-"])
        if support == 0:
            continue
        graph[left].append((right, relation, support))
        graph[right].append((left, relation, support))
        constraints.append((left, right, relation, support))
        unit_links[left] += support
        unit_links[right] += support
    reference = min(range(len(units)), key=lambda index: (-unit_links[index], units[index].unit_id))

    parity: dict[int, int] = {}
    primary_component: set[int] = set()
    roots = [reference] + [index for index in range(len(units)) if index != reference]
    for root in roots:
        if root in parity:
            continue
        parity[root] = 0
        queue = [root]
        component: set[int] = {root}
        while queue:
            current = queue.pop(0)
            for neighbor, relation, _support in sorted(graph[current], key=lambda edge: (-edge[2], edge[0])):
                expected = parity[current] ^ relation
                if neighbor in parity:
                    continue
                parity[neighbor] = expected
                component.add(neighbor)
                queue.append(neighbor)
        if root == reference:
            primary_component = component

    # A spanning tree can lock in weak root edges. Optimize the full signed graph
    # while holding the arbitrary reference orientation fixed.
    for _iteration in range(max(1, len(units) * 4)):
        best_delta = 0
        best_index: int | None = None
        for index in range(len(units)):
            if index == reference:
                continue
            delta = 0
            for left, right, relation, support in constraints:
                if index not in (left, right):
                    continue
                satisfied = (parity[left] ^ parity[right]) == relation
                delta += -support if satisfied else support
            if delta > best_delta:
                best_delta, best_index = delta, index
        if best_index is None:
            break
        parity[best_index] ^= 1

    conflicts = sum((parity[left] ^ parity[right]) != relation for left, right, relation, _support in constraints)
    total_weight = sum(support for _left, _right, _relation, support in constraints)
    satisfied_weight = sum(
        support
        for left, right, relation, support in constraints
        if (parity[left] ^ parity[right]) == relation
    )
    orientations = {units[index].unit_id: "+" if parity[index] == 0 else "-" for index in range(len(units))}
    return orientations, {
        "reference_unit_id": units[reference].unit_id,
        "disconnected_units": len(units) - len(primary_component),
        "orientation_conflicts": conflicts,
        "orientation_edges": len(constraints),
        "orientation_weight_agreement": satisfied_weight / total_weight if total_weight else 0.0,
        "unit_links": {units[index].unit_id: unit_links[index] for index in range(len(units))},
    }


def oriented_unit_sequence(unit: Unit, sequence: str, region: Region) -> str:
    """Extract a unit in the biological direction of its repeat anchor."""
    if unit.chrom != region.chrom:
        raise ValueError("unit chromosome does not match FASTA region")
    if unit.strand not in "+-":
        raise ValueError("unit strand must be '+' or '-'")
    return extract_genomic(sequence, region, unit.start0 + 1, unit.end0, unit.strand)


def unit_phase_offset(unit: Unit, hits: list[RepeatHit]) -> int:
    """Return the oriented start of the dominant, same-strand family body."""
    candidates: list[tuple[int, int, int, int]] = []
    for hit in hits:
        genomic_start0 = max(unit.start0, hit.start0)
        genomic_end0 = min(unit.end0, hit.end0)
        if genomic_start0 >= genomic_end0:
            continue
        if unit.strand == "+":
            projected_start0 = genomic_start0 - unit.start0
        else:
            projected_start0 = unit.end0 - genomic_end0
        candidates.append(
            (
                int(hit.strand == unit.strand),
                genomic_end0 - genomic_start0,
                _consensus_span(hit),
                projected_start0,
            )
        )
    if not candidates:
        return 0
    best = max(candidates, key=lambda item: (item[0], item[1], item[2], -item[3]))
    return best[3]


def circular_rephase(sequence: str, offset: int) -> str:
    """Rotate one complete tandem unit so its dominant family body begins at zero."""
    if not sequence:
        return sequence
    if not 0 <= offset < len(sequence):
        raise ValueError("phase offset is outside the unit sequence")
    return sequence[offset:] + sequence[:offset]


def _directional_kmers(sequence: str, k: int) -> set[str]:
    sequence = sequence.upper()
    if k < 1 or len(sequence) < k:
        return set()
    return {
        kmer
        for start0 in range(len(sequence) - k + 1)
        if set(kmer := sequence[start0 : start0 + k]) <= {"A", "C", "G", "T"}
    }


def oriented_kmer_support(reference: str, query: str, k: int = 21) -> dict[str, int | float]:
    """Compare orientation-sensitive k-mer sets in the observed and reverse directions."""
    if k < 1:
        raise ValueError("k must be positive")
    reference_kmers = _directional_kmers(reference, k)

    def jaccard(candidate: str) -> float:
        candidate_kmers = _directional_kmers(candidate, k)
        union = reference_kmers | candidate_kmers
        return len(reference_kmers & candidate_kmers) / len(union) if union else 0.0

    return {
        "fraction": jaccard(query),
        "reverse_fraction": jaccard(revcomp(query)),
        "reference_kmers": len(reference_kmers),
        "k": k,
    }


def anchor_support(reference: str, query: str, k: int = 21, step: int = 1000, radius: int = 250) -> dict[str, int | float]:
    """Measure orientation-sensitive support at regularly spaced k-mer anchors.

    Anchors containing ambiguity codes are not informative.  Candidate matches
    may move within ``radius``, but a nonzero displacement is accepted only
    when at least two anchors and 80% of all tested anchors support the same
    shift.  This tolerates boundary offsets without mistaking reverse-order
    k-mer matches for oriented support.
    """
    if k < 1 or step < 1 or radius < 0:
        raise ValueError("k and step must be positive and radius non-negative")
    reference = reference.upper()
    query = query.upper()
    tested = 0
    offset_counts: dict[int, int] = {}
    if len(reference) >= k:
        for start0 in range(0, len(reference) - k + 1, step):
            anchor = reference[start0 : start0 + k]
            if set(anchor) <= {"A", "C", "G", "T"}:
                tested += 1
                low = max(0, start0 - radius)
                high = min(len(query) - k, start0 + radius)
                offsets_for_anchor: set[int] = set()
                for query_start0 in range(low, high + 1):
                    if query[query_start0 : query_start0 + k] == anchor:
                        offsets_for_anchor.add(query_start0 - start0)
                for offset in offsets_for_anchor:
                    offset_counts[offset] = offset_counts.get(offset, 0) + 1
    selected_offset = 0
    if offset_counts:
        best_offset = min(offset_counts, key=lambda offset: (-offset_counts[offset], abs(offset), offset))
        if best_offset == 0 or (tested >= 2 and offset_counts[best_offset] / tested >= 0.8):
            selected_offset = best_offset
    matched = offset_counts.get(selected_offset, 0)
    return {
        "matched": matched,
        "tested": tested,
        "fraction": matched / tested if tested else 0.0,
        "k": k,
        "step": step,
        "radius": radius,
    }


def _canonical_kmers(sequence: str, k: int) -> set[str]:
    sequence = sequence.upper()
    if len(sequence) < k:
        return set()
    kmers: set[str] = set()
    for start0 in range(len(sequence) - k + 1):
        kmer = sequence[start0 : start0 + k]
        if set(kmer) <= {"A", "C", "G", "T"}:
            kmers.add(min(kmer, revcomp(kmer)))
    return kmers


def _jaccard(left: set[str], right: set[str]) -> float:
    if not left and not right:
        return 1.0
    return len(left & right) / len(left | right)


def choose_medoid(records: dict[str, str], k: int = 21) -> str:
    """Return the deterministic k-mer medoid record name, never a consensus."""
    if k < 1:
        raise ValueError("k must be positive")
    if not records:
        raise ValueError("medoid selection requires at least one record")
    names = sorted(records)
    kmer_sets = {name: _canonical_kmers(records[name], k) for name in names}
    median_length = statistics.median(len(records[name]) for name in names)
    ranked: list[tuple[float, float, str]] = []
    for name in names:
        similarities = [_jaccard(kmer_sets[name], kmer_sets[other]) for other in names if other != name]
        score = statistics.median(similarities) if similarities else 1.0
        ranked.append((-score, abs(len(records[name]) - median_length), name))
    return min(ranked)[2]


def _record_items(records: dict[str, str]) -> list[tuple[str, str]]:
    items = list(records.items())
    for name, sequence in items:
        if not name or any(character.isspace() for character in name) or ">" in name:
            raise ValueError("FASTA record names must be nonempty and contain no whitespace or '>'")
        if any(character.isspace() for character in sequence):
            raise ValueError("FASTA sequences must not contain whitespace")
    return items


def write_fasta(records: dict[str, str], path: Path) -> None:
    """Write mapping-ordered FASTA records, wrapping sequence lines at 80 bp."""
    lines: list[str] = []
    for name, sequence in _record_items(records):
        lines.append(f">{name}")
        lines.extend(sequence[start0 : start0 + 80] for start0 in range(0, len(sequence), 80))
    Path(path).write_text("\n".join(lines) + ("\n" if lines else ""))


def _read_fasta_text(text: str) -> dict[str, str]:
    records: dict[str, str] = {}
    name: str | None = None
    chunks: list[str] = []
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records[name] = "".join(chunks)
            name = line[1:].split()[0] if line[1:].strip() else ""
            if not name or name in records:
                raise RuntimeError("MAFFT produced invalid or duplicate FASTA record names")
            chunks = []
        elif name is None:
            raise RuntimeError("MAFFT produced malformed FASTA output")
        else:
            chunks.append(line.upper())
    if name is not None:
        records[name] = "".join(chunks)
    return records


def run_mafft_consensus(records: dict[str, str], executable: str | Path) -> tuple[str, list[str]]:
    """Run an explicitly requested MAFFT and build a gap-aware majority consensus."""
    items = _record_items(records)
    if not items:
        raise ValueError("MAFFT consensus requires at least one record")
    if any(not sequence for _name, sequence in items):
        raise ValueError("MAFFT consensus does not accept empty sequences")
    with tempfile.TemporaryDirectory(prefix="chry-mafft-") as directory:
        input_path = Path(directory) / "records.fa"
        write_fasta(records, input_path)
        try:
            result = subprocess.run(
                [str(executable), "--quiet", str(input_path)],
                text=True,
                capture_output=True,
                check=False,
            )
        except OSError as exc:
            raise RuntimeError(f"MAFFT failed to start: {exc}") from exc
    if result.returncode != 0:
        detail = result.stderr.strip() or f"exit status {result.returncode}"
        raise RuntimeError(f"MAFFT failed: {detail}")
    output_records = _read_fasta_text(result.stdout)
    input_names = [name for name, _sequence in items]
    if set(output_records) != set(input_names):
        raise RuntimeError("MAFFT output record names do not match input records")
    aligned = [output_records[name] for name in input_names]
    lengths = {len(sequence) for sequence in aligned}
    if len(lengths) != 1:
        raise RuntimeError("MAFFT produced unequal alignment lengths")
    if not lengths or next(iter(lengths)) == 0:
        raise RuntimeError("MAFFT produced an empty alignment")

    consensus: list[str] = []
    for column in zip(*aligned):
        bases = [base for base in column if base != "-"]
        if not bases:
            continue
        counts = {base: bases.count(base) for base in set(bases)}
        maximum = max(counts.values())
        tied = {base for base, count in counts.items() if count == maximum}
        reference_base = next((base for base in column if base != "-"), None)
        consensus.append(reference_base if reference_base in tied else min(tied))
    return "".join(consensus), aligned


def find_hairpins(
    seq: str,
    end0: int,
    stem_min: int = 4,
    loop_min: int = 3,
    loop_max: int = 12,
) -> list[Hairpin]:
    """Find short inverted repeats in the 100 bp preceding a terminal motif."""
    if not 0 <= end0 <= len(seq):
        raise ValueError("hairpin terminal coordinate is outside the sequence")
    if stem_min < 1 or loop_min < 0 or loop_max < loop_min:
        raise ValueError("invalid hairpin stem or loop bounds")
    sequence = seq.upper()
    window_start0 = max(0, end0 - 100)
    hairpins: list[Hairpin] = []
    for stem_length in range(stem_min, 7):
        for left_start0 in range(window_start0, end0):
            left_end0 = left_start0 + stem_length
            if left_end0 > end0:
                break
            stem_left = sequence[left_start0:left_end0]
            if not stem_left or not set(stem_left) <= {"A", "C", "G", "T"}:
                continue
            expected_right = revcomp(stem_left)
            for loop_length in range(loop_min, loop_max + 1):
                right_start0 = left_end0 + loop_length
                right_end0 = right_start0 + stem_length
                if right_end0 > end0:
                    break
                if sequence[right_start0:right_end0] != expected_right:
                    continue
                hairpins.append(
                    Hairpin(
                        left_start0,
                        right_end0,
                        stem_left,
                        sequence[left_end0:right_start0],
                        expected_right,
                        end0 - right_end0,
                    )
                )
    return sorted(
        hairpins,
        key=lambda hairpin: (
            hairpin.distance_to_terminal_end,
            -len(hairpin.stem_left),
            hairpin.start0,
            hairpin.loop,
        ),
    )


def find_terminal_pairs(seq: str, flank: int = 100) -> list[TerminalPair]:
    """Enumerate inward-shifted 5' TC and 3' CTRR terminal candidates."""
    if flank < 0:
        raise ValueError("terminal flank must be non-negative")
    sequence = seq.upper()
    five_limit = min(len(sequence) - 1, flank)
    five_starts = [start0 for start0 in range(max(five_limit, 0)) if sequence[start0 : start0 + 2] == "TC"]
    three_start = max(0, len(sequence) - flank - 4)
    three_starts = [
        start0
        for start0 in range(three_start, max(len(sequence) - 3, 0))
        if sequence[start0 : start0 + 2] == "CT" and set(sequence[start0 + 2 : start0 + 4]) <= {"A", "G"}
    ]
    pairs: list[TerminalPair] = []
    scale = max(flank, 1)
    for five_start0 in five_starts:
        for three_start0 in three_starts:
            if five_start0 + 2 > three_start0:
                continue
            hairpins = find_hairpins(sequence, three_start0)
            hairpin = min(
                hairpins,
                key=lambda candidate: (
                    candidate.distance_to_terminal_end,
                    -len(candidate.stem_left),
                    candidate.start0,
                ),
                default=None,
            )
            three_offset = len(sequence) - (three_start0 + 4)
            boundary_score = max(0.0, 1.0 - five_start0 / scale) + max(0.0, 1.0 - three_offset / scale)
            hairpin_score = 2.0 + 0.1 * len(hairpin.stem_left) - 0.01 * hairpin.distance_to_terminal_end if hairpin else 0.0
            context_start = max(0, five_start0 - 10)
            context_end = min(len(sequence), three_start0 + 14)
            pairs.append(
                TerminalPair(
                    five_start0,
                    "TC",
                    three_offset,
                    sequence[three_start0 : three_start0 + 4],
                    hairpin,
                    boundary_score + hairpin_score,
                    sequence[context_start : min(len(sequence), five_start0 + 12)],
                    sequence[max(0, three_start0 - 10) : context_end],
                )
            )
    return sorted(
        pairs,
        key=lambda pair: (
            pair.five_prime_offset,
            pair.three_prime_offset,
            pair.three_prime_motif,
        ),
    )


def scan_six_frame_orfs(seq: str) -> list[ORF]:
    """Return maximal stop-free codon runs on both strands.

    ORF length is structural evidence only; this function deliberately makes
    no claim about Rep/HUH or helicase domains.
    """
    sequence = seq.upper()
    orfs: list[ORF] = []
    for strand, oriented in (("+", sequence), ("-", revcomp(sequence))):
        for frame in range(3):
            codon_end0 = frame + ((len(oriented) - frame) // 3) * 3
            segment_start0 = frame
            for codon_start0 in range(frame, codon_end0, 3):
                if oriented[codon_start0 : codon_start0 + 3] not in _STOP_CODONS:
                    continue
                if codon_start0 > segment_start0:
                    start0, end0 = segment_start0, codon_start0
                    if strand == "-":
                        start0, end0 = len(sequence) - end0, len(sequence) - start0
                    orfs.append(ORF(strand, frame, start0, end0, (end0 - start0) // 3))
                segment_start0 = codon_start0 + 3
            if codon_end0 > segment_start0:
                start0, end0 = segment_start0, codon_end0
                if strand == "-":
                    start0, end0 = len(sequence) - end0, len(sequence) - start0
                orfs.append(ORF(strand, frame, start0, end0, (end0 - start0) // 3))
    return sorted(orfs, key=lambda orf: (orf.strand != "+", orf.frame, orf.start0, orf.end0))


def classify_hle_evidence(evidence: HLEEvidence) -> str:
    """Apply the conservative structural/domain evidence hierarchy."""
    strong_structure = evidence.terminal_support >= 0.5 and evidence.hairpin_support >= 0.5
    if not strong_structure:
        return "composite ampliconic repeat"
    if evidence.rep_domain and evidence.helicase_domain and evidence.independent_loci > 0:
        return "autonomous HLE"
    if evidence.rep_domain and evidence.helicase_domain:
        return "putative autonomous HLE"
    if evidence.independent_loci > 0:
        return "putative non-autonomous HLE"
    return "Helitron-derived composite repeat"


def project_feature_to_unit(feature: Feature, unit: Unit) -> ProjectedFeature | None:
    """Clip a genomic feature to a unit and project into oriented coordinates."""
    if unit.strand not in "+-" or unit.end0 <= unit.start0:
        raise ValueError("unit must have a valid interval and strand")
    if feature.end0 < feature.start0:
        raise ValueError("feature interval is invalid")
    if feature.chrom != unit.chrom:
        return None
    genomic_start0 = max(feature.start0, unit.start0)
    genomic_end0 = min(feature.end0, unit.end0)
    if genomic_start0 >= genomic_end0:
        return None
    if unit.strand == "+":
        start0 = genomic_start0 - unit.start0
        end0 = genomic_end0 - unit.start0
        strand = feature.strand
    else:
        start0 = unit.end0 - genomic_end0
        end0 = unit.end0 - genomic_start0
        strand = {"+": "-", "-": "+", ".": "."}.get(feature.strand)
        if strand is None:
            raise ValueError("feature strand must be '+', '-', or '.'")
    return ProjectedFeature(
        unit.unit_id,
        feature.feature_type,
        start0,
        end0,
        strand,
        dict(feature.attributes),
        feature.source,
        feature.phase,
        genomic_start0,
        genomic_end0,
        genomic_start0 != feature.start0 or genomic_end0 != feature.end0,
    )


def circular_rephase_intervals(start0: int, end0: int, length: int, offset: int) -> list[tuple[int, int]]:
    """Rotate a half-open interval and split it if it crosses the new unit boundary."""
    if length < 1 or not 0 <= start0 < end0 <= length or not 0 <= offset < length:
        raise ValueError("invalid circular interval or phase offset")
    if end0 - start0 == length:
        return [(0, length)]
    rotated_start0 = (start0 - offset) % length
    rotated_end0 = (end0 - offset) % length
    if rotated_start0 < rotated_end0:
        return [(rotated_start0, rotated_end0)]
    return [
        (interval_start0, interval_end0)
        for interval_start0, interval_end0 in ((0, rotated_end0), (rotated_start0, length))
        if interval_end0 > interval_start0
    ]


def rephase_projected_features(
    features: list[ProjectedFeature], units: list[Unit], phase_offsets: dict[str, int]
) -> list[ProjectedFeature]:
    """Move annotations into the same phase as the oriented unit FASTA records."""
    units_by_id = {unit.unit_id: unit for unit in units}
    rephased: list[ProjectedFeature] = []
    for feature in features:
        unit = units_by_id[feature.unit_id]
        length = unit.end0 - unit.start0
        offset = phase_offsets.get(feature.unit_id, 0)
        for start0, end0 in circular_rephase_intervals(feature.start0, feature.end0, length, offset):
            rephased.append(
                ProjectedFeature(
                    feature.unit_id,
                    feature.feature_type,
                    start0,
                    end0,
                    feature.strand,
                    dict(feature.attributes),
                    feature.source,
                    feature.phase,
                    feature.genomic_start0,
                    feature.genomic_end0,
                    feature.clipped,
                )
            )
    return rephased


def _is_tspy(feature: ProjectedFeature) -> bool:
    if feature.feature_type.lower() not in {"gene", "pseudogene"}:
        return False
    return any("TSPY" in str(value).upper() for value in feature.attributes.values())


def cluster_tspy_phases(projected: list[ProjectedFeature], tolerance: int = 750) -> dict[str, float]:
    """Infer phase centers from modes of oriented TSPY feature midpoints."""
    if tolerance < 0:
        raise ValueError("phase clustering tolerance must be non-negative")
    observations = sorted((feature.start0 + feature.end0) / 2 for feature in projected if _is_tspy(feature))
    if not observations:
        return {}
    clusters: list[list[float]] = []
    for midpoint in observations:
        eligible = [
            (abs(midpoint - statistics.median(cluster)), index)
            for index, cluster in enumerate(clusters)
            if abs(midpoint - statistics.median(cluster)) <= tolerance
        ]
        if eligible:
            clusters[min(eligible)[1]].append(midpoint)
        else:
            clusters.append([midpoint])
    centers = sorted(float(statistics.median(cluster)) for cluster in clusters)
    if len(centers) > 26:
        raise ValueError("more than 26 TSPY phase modes are not supported")
    return {chr(ord("A") + index): center for index, center in enumerate(centers)}


def _projected_gene_state(feature: ProjectedFeature) -> str:
    feature_type = feature.feature_type.lower()
    pseudo_attribute = any(
        key.lower() == "pseudo" and str(value).lower() in {"true", "1", "yes"}
        for key, value in feature.attributes.items()
    )
    pseudo_text = any("pseudogene" in str(value).lower() for value in feature.attributes.values())
    if "pseudogene" in feature_type or pseudo_attribute or pseudo_text:
        return "pseudogene"
    if feature_type in {"gene", "mrna", "transcript"} or feature_type.endswith("_gene"):
        return "coding"
    return "fragment"


def format_subtype(states: dict[str, str]) -> str:
    """Format phase-specific states without conflating coding and pseudogene calls."""
    if not states:
        return "no_TSPY"
    return "__".join(f"{phase}_{states[phase]}" for phase in sorted(states))


def classify_unit_subtype(
    unit: Unit,
    projected: list[ProjectedFeature],
    phase_centers: dict[str, float],
) -> str:
    """Assign each observed unit TSPY call to its nearest inferred phase."""
    if not phase_centers:
        return format_subtype({})
    states_by_phase: dict[str, set[str]] = {phase: set() for phase in phase_centers}
    for feature in projected:
        if feature.unit_id != unit.unit_id or not _is_tspy(feature):
            continue
        midpoint = (feature.start0 + feature.end0) / 2
        phase = min(phase_centers, key=lambda name: (abs(midpoint - phase_centers[name]), name))
        states_by_phase[phase].add(_projected_gene_state(feature))
    priority = {"coding": 0, "pseudogene": 1, "fragment": 2}
    states: dict[str, str] = {}
    for phase, observed_states in states_by_phase.items():
        if not observed_states:
            states[phase] = "missing"
        else:
            states[phase] = "+".join(sorted(observed_states, key=lambda state: (priority[state], state)))
    return format_subtype(states)


def classify_family_free_tspy_subtype(unit: Unit, projected: list[ProjectedFeature]) -> str:
    """Describe macro-unit TSPY content by coding/pseudogene counts, without phase letters."""
    coding = 0
    pseudogene = 0
    for feature in projected:
        if feature.unit_id != unit.unit_id or not _is_tspy(feature):
            continue
        if feature.genomic_start0 is not None and feature.genomic_end0 is not None:
            genomic_midpoint = (feature.genomic_start0 + feature.genomic_end0) // 2
            if not unit.start0 <= genomic_midpoint < unit.end0:
                continue
        state = _projected_gene_state(feature)
        coding += int(state == "coding")
        pseudogene += int(state == "pseudogene")
    if coding == 0 and pseudogene == 0:
        return "no_TSPY"
    return f"TSPY_coding_{coding}__pseudogene_{pseudogene}"


def classify_structural_units(
    units: list[Unit],
    hits: list[StainedGlassHit],
    *,
    min_identity: float = 95.0,
    min_reciprocal_coverage: float = 0.75,
    minimum_class_copies: int = 2,
    minimum_length_ratio: float = 0.85,
) -> dict[str, str]:
    """Assign conservative whole-unit structural classes from reciprocal coverage.

    Strand and gene annotations are deliberately excluded. A distinct class is
    reported only when the same structural state is supported by multiple copies;
    unsupported singleton/partial states remain variants of the dominant class.
    """
    if not units:
        return {}
    if not 0 < min_reciprocal_coverage <= 1 or minimum_class_copies < 1:
        raise ValueError("invalid structural-class thresholds")
    spans = [hit.query_end0 - hit.query_start0 for hit in hits if hit.query_end0 > hit.query_start0]
    window_size = int(statistics.median(spans)) if spans else 1

    def locate(start0: int, end0: int) -> int | None:
        midpoint = (start0 + end0 - 1) // 2
        return next(
            (index for index, unit in enumerate(units) if unit.start0 <= midpoint < unit.end0),
            None,
        )

    pair_support: dict[tuple[int, int], dict[int, set[int]]] = {}
    for hit in hits:
        if hit.identity < min_identity:
            continue
        query_index = locate(hit.query_start0, hit.query_end0)
        reference_index = locate(hit.reference_start0, hit.reference_end0)
        if query_index is None or reference_index is None or query_index == reference_index:
            continue
        left, right = sorted((query_index, reference_index))
        support = pair_support.setdefault((left, right), {left: set(), right: set()})
        if query_index == left:
            support[left].add(hit.query_start0)
            support[right].add(hit.reference_start0)
        else:
            support[left].add(hit.reference_start0)
            support[right].add(hit.query_start0)

    parent = list(range(len(units)))

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(left: int, right: int) -> None:
        left_root, right_root = find(left), find(right)
        if left_root != right_root:
            parent[right_root] = left_root

    for (left, right), support in pair_support.items():
        left_length = units[left].end0 - units[left].start0
        right_length = units[right].end0 - units[right].start0
        length_ratio = min(left_length, right_length) / max(left_length, right_length)
        if length_ratio < minimum_length_ratio:
            continue
        left_windows = max(1, math.ceil(left_length / window_size))
        right_windows = max(1, math.ceil(right_length / window_size))
        reciprocal_coverage = min(
            len(support[left]) / left_windows,
            len(support[right]) / right_windows,
        )
        if reciprocal_coverage >= min_reciprocal_coverage:
            union(left, right)

    components: dict[int, list[int]] = {}
    for index in range(len(units)):
        components.setdefault(find(index), []).append(index)
    reproducible = [
        members for members in components.values() if len(members) >= minimum_class_copies
    ]
    if not reproducible:
        reproducible = [max(components.values(), key=lambda members: (len(members), -members[0]))]
    reproducible.sort(key=lambda members: (-len(members), units[min(members)].start0))
    class_by_root: dict[int, str] = {}
    for class_index, members in enumerate(reproducible, 1):
        label = f"U{class_index}"
        for member in members:
            class_by_root[find(member)] = label
    dominant = "U1"
    return {
        unit.unit_id: class_by_root.get(find(index), f"{dominant}_variant")
        for index, unit in enumerate(units)
    }


def project_family_bodies_to_units(
    units: list[Unit], hits: list[RepeatHit]
) -> dict[str, list[tuple[int, int]]]:
    """Project selected family-body intervals into each unit's orientation."""
    bodies_by_unit: dict[str, list[tuple[int, int]]] = {unit.unit_id: [] for unit in units}
    for unit in units:
        for hit in hits:
            genomic_start0 = max(hit.start0, unit.start0)
            genomic_end0 = min(hit.end0, unit.end0)
            if genomic_start0 >= genomic_end0:
                continue
            if unit.strand == "+":
                start0 = genomic_start0 - unit.start0
                end0 = genomic_end0 - unit.start0
            else:
                start0 = unit.end0 - genomic_end0
                end0 = unit.end0 - genomic_start0
            bodies_by_unit[unit.unit_id].append((start0, end0))
        bodies_by_unit[unit.unit_id].sort()
    return bodies_by_unit


def normalize_tspy_to_family_bodies(
    projected: list[ProjectedFeature],
    bodies_by_unit: dict[str, list[tuple[int, int]]],
) -> list[ProjectedFeature]:
    """Express top-level TSPY annotations relative to their associated family body."""
    normalized: list[ProjectedFeature] = []
    for feature in projected:
        if not _is_tspy(feature):
            continue
        bodies = bodies_by_unit.get(feature.unit_id, [])
        if not bodies:
            continue
        midpoint = (feature.start0 + feature.end0) / 2

        def body_rank(body: tuple[int, int]) -> tuple[float, int, int]:
            start0, end0 = body
            if start0 <= midpoint < end0:
                distance = 0.0
            else:
                distance = min(abs(midpoint - start0), abs(midpoint - end0))
            return distance, start0, end0

        body_start0, _body_end0 = min(bodies, key=body_rank)
        normalized.append(
            ProjectedFeature(
                feature.unit_id,
                feature.feature_type,
                feature.start0 - body_start0,
                feature.end0 - body_start0,
                feature.strand,
                dict(feature.attributes),
                feature.source,
                feature.phase,
                feature.genomic_start0,
                feature.genomic_end0,
                feature.clipped,
            )
        )
    return normalized


def infer_body_relative_tspy_phases(
    units: list[Unit],
    projected: list[ProjectedFeature],
    hits: list[RepeatHit],
    tolerance: int = 750,
) -> tuple[dict[str, float], dict[str, str]]:
    """Infer complete-unit body-relative phases, then classify every unit."""
    bodies_by_unit = project_family_bodies_to_units(units, hits)
    normalized = normalize_tspy_to_family_bodies(projected, bodies_by_unit)
    complete_ids = {unit.unit_id for unit in units if unit.status == "complete"}
    complete_observations = [
        feature for feature in normalized if feature.unit_id in complete_ids
    ]
    coding_seeds = [
        feature
        for feature in complete_observations
        if feature.feature_type.lower() == "gene"
        and _projected_gene_state(feature) == "coding"
        and not feature.clipped
    ]
    phase_seeds = coding_seeds or [
        feature
        for feature in complete_observations
        if _projected_gene_state(feature) == "pseudogene"
    ]
    phase_centers = cluster_tspy_phases(phase_seeds, tolerance=tolerance)
    subtypes = {
        unit.unit_id: classify_unit_subtype(unit, normalized, phase_centers)
        for unit in units
    }
    return phase_centers, subtypes


def _repeat_class(feature: ProjectedFeature) -> str | None:
    class_value: str | None = None
    for key, value in feature.attributes.items():
        if key.lower() in {"class", "repeat_class", "classification"} and value:
            class_value = str(value)
            break
    raw = class_value or feature.feature_type
    lowered = raw.lower().replace("-", "_").replace(" ", "_")
    feature_lower = feature.feature_type.lower()
    if not class_value and not any(
        token in feature_lower
        for token in ("repeat", "line", "sine", "ltr", "helitron", "transpos", "satellite", "low_complexity")
    ):
        return None
    if "simple" in lowered and "repeat" in lowered:
        return "simple_repeat"
    if "low" in lowered and "complex" in lowered:
        return "low_complexity"
    for label in ("LINE", "SINE", "LTR", "Helitron", "satellite"):
        if label.lower() in lowered:
            return label
    safe = re.sub(r"[^A-Za-z0-9]+", "_", raw).strip("_")
    return safe or "other_TE"


def _is_nested_repeat(feature: ProjectedFeature) -> bool:
    for key, value in feature.attributes.items():
        if "nested" in key.lower() and str(value).lower() not in {"", "0", "false", "no", "none"}:
            return True
        if key.lower() in {"relation", "relationship"} and "nested" in str(value).lower():
            return True
    return False


def _union_length(intervals: list[tuple[int, int]]) -> int:
    if not intervals:
        return 0
    merged_length = 0
    current_start0, current_end0 = sorted(intervals)[0]
    for start0, end0 in sorted(intervals)[1:]:
        if start0 <= current_end0:
            current_end0 = max(current_end0, end0)
        else:
            merged_length += current_end0 - current_start0
            current_start0, current_end0 = start0, end0
    return merged_length + current_end0 - current_start0


def summarize_inter_unit_gaps(
    units: list[Unit],
    region: Region,
    repeats: list[Feature],
    coarse_gaps: list[dict[str, object]],
) -> list[dict[str, object]]:
    """Report unassigned flanks/gaps and their non-double-counted repeat content."""
    if not units:
        return []
    ordered = sorted(units, key=lambda unit: (unit.start0, unit.end0))
    intervals: list[tuple[int, int, str]] = []
    region_start0 = region.start1 - 1
    if region_start0 < ordered[0].start0:
        intervals.append((region_start0, ordered[0].start0, "left_flank"))
    for left, right in zip(ordered, ordered[1:]):
        if left.end0 < right.start0:
            intervals.append((left.end0, right.start0, "inter_unit_gap"))
    if ordered[-1].end0 < region.end1:
        intervals.append((ordered[-1].end0, region.end1, "right_flank"))

    rows: list[dict[str, object]] = []
    for start0, end0, default_type in intervals:
        best_coarse: dict[str, object] | None = None
        best_overlap = 0
        for coarse in coarse_gaps:
            overlap = max(
                0,
                min(end0, int(coarse["end0"])) - max(start0, int(coarse["start0"])),
            )
            if overlap > best_overlap:
                best_overlap, best_coarse = overlap, coarse
        interval_type = default_type
        hub_windows = 0
        if default_type == "inter_unit_gap" and best_coarse is not None:
            interval_type = str(best_coarse.get("interval_type", default_type))
            hub_windows = int(best_coarse.get("hub_windows", 0))

        repeat_intervals: list[tuple[int, int]] = []
        class_intervals: dict[str, list[tuple[int, int]]] = {}
        for feature in repeats:
            if feature.chrom != region.chrom:
                continue
            overlap_start = max(start0, feature.start0)
            overlap_end = min(end0, feature.end0)
            if overlap_start >= overlap_end:
                continue
            interval = (overlap_start, overlap_end)
            repeat_intervals.append(interval)
            class_intervals.setdefault(feature.feature_type, []).append(interval)
        repeat_bp = _union_length(repeat_intervals)
        class_bp = {
            repeat_class: _union_length(class_hits)
            for repeat_class, class_hits in class_intervals.items()
        }
        rows.append(
            {
                "gap_id": f"gap_{len(rows) + 1:03d}",
                "chrom": region.chrom,
                "start1": start0 + 1,
                "end1": end0,
                "length_bp": end0 - start0,
                "interval_type": interval_type,
                "repeat_bp": repeat_bp,
                "repeat_fraction": repeat_bp / (end0 - start0),
                "top_repeat_classes": ";".join(
                    f"{name}:{bp}"
                    for name, bp in sorted(class_bp.items(), key=lambda item: (-item[1], item[0]))[:5]
                ),
                "hub_windows": hub_windows,
            }
        )
    return rows


def build_repeat_inventory(
    repeats: list[Feature],
    units: list[Unit],
    gap_rows: list[dict[str, object]],
    region: Region,
) -> list[dict[str, object]]:
    """Retain every target-overlapping repeat GFF record exactly once.

    Unit-level projection can legitimately split one annotation at a boundary.
    This inventory is the lossless companion table: it keeps the original repeat
    identity and records every unit, inter-unit gap, or flank that it overlaps.
    """
    region_start0 = region.start1 - 1
    ordered_units = sorted(units, key=lambda unit: (unit.start0, unit.end0, unit.unit_id))
    rows: list[dict[str, object]] = []
    for index, feature in enumerate(repeats, 1):
        start0 = max(region_start0, feature.start0)
        end0 = min(region.end1, feature.end0)
        if feature.chrom != region.chrom or end0 <= start0:
            continue
        assignments: list[tuple[int, str]] = []
        retained_intervals: list[tuple[int, int]] = []
        for unit in ordered_units:
            overlap_start0 = max(start0, unit.start0)
            overlap_end0 = min(end0, unit.end0)
            if overlap_end0 > overlap_start0:
                assignments.append((overlap_start0, unit.unit_id))
                retained_intervals.append((overlap_start0, overlap_end0))
        for gap in gap_rows:
            gap_start0 = int(gap["start1"]) - 1
            gap_end0 = int(gap["end1"])
            overlap_start0 = max(start0, gap_start0)
            overlap_end0 = min(end0, gap_end0)
            if overlap_end0 > overlap_start0:
                assignments.append(
                    (
                        overlap_start0,
                        f"{gap['gap_id']}:{gap['interval_type']}",
                    )
                )
                retained_intervals.append((overlap_start0, overlap_end0))
        retained_bp = _union_length(retained_intervals)
        if retained_bp < end0 - start0:
            assignments.append((start0, "unassigned"))
        repeat_id = feature.attributes.get("ID") or f"repeat_{index:06d}"
        repeat_class = _repeat_class(feature) or feature.feature_type or "other_repeat"
        rows.append(
            {
                "repeat_id": repeat_id,
                "chrom": feature.chrom,
                "original_start1": feature.start0 + 1,
                "original_end1": feature.end0,
                "original_length_bp": feature.end0 - feature.start0,
                "target_start1": start0 + 1,
                "target_end1": end0,
                "target_overlap_bp": end0 - start0,
                "start1": start0 + 1,
                "end1": end0,
                "length_bp": end0 - start0,
                "source": feature.source,
                "feature_type": feature.feature_type,
                "repeat_class": repeat_class,
                "repeat_family": feature.attributes.get("Name", ""),
                "plot_category": repeat_plot_category(feature),
                "strand": feature.strand,
                "score": feature.score,
                "phase": feature.phase,
                "assigned_intervals": ";".join(
                    label for _position0, label in sorted(set(assignments))
                ),
                "retained_bp": end0 - start0,
                "attributes": _serialize_attributes(feature.attributes),
                "raw_attributes": (
                    feature.raw_attributes
                    if feature.raw_attributes is not None
                    else _serialize_attributes(feature.attributes)
                ),
            }
        )
    return rows


def summarize_repeat_retention(
    repeats: list[Feature], inventory_rows: list[dict[str, object]], region: Region
) -> dict[str, int | str]:
    """Cross-check record counts, target-overlap bases, and drawable categories."""
    region_start0 = region.start1 - 1
    expected = [
        feature
        for feature in repeats
        if feature.chrom == region.chrom
        and feature.start0 < region.end1
        and feature.end0 > region_start0
    ]
    target_overlap_bp = sum(
        min(feature.end0, region.end1) - max(feature.start0, region_start0)
        for feature in expected
    )
    inventory_overlap_bp = sum(
        int(row.get("target_overlap_bp", row.get("length_bp", 0)))
        for row in inventory_rows
    )
    retained_bp = sum(int(row.get("retained_bp", 0)) for row in inventory_rows)
    unmapped = sum(
        repeat_plot_category(feature) not in PLOT_COLORS for feature in expected
    ) + sum(
        str(row.get("plot_category", "")) not in PLOT_COLORS
        for row in inventory_rows
    )
    passed = (
        len(inventory_rows) == len(expected)
        and inventory_overlap_bp == target_overlap_bp
        and retained_bp == target_overlap_bp
        and unmapped == 0
    )
    return {
        "repeat_input_features": len(expected),
        "repeat_inventory_features": len(inventory_rows),
        "repeat_target_overlap_bp": target_overlap_bp,
        "repeat_inventory_overlap_bp": inventory_overlap_bp,
        "repeat_inventory_retained_bp": retained_bp,
        "repeat_plot_unmapped_features": unmapped,
        "repeat_retention_status": "pass" if passed else "fail",
    }


def collapse_repeat_coverage(projected: list[ProjectedFeature]) -> dict[str, int]:
    """Union repeat intervals per class/unit and report nested hits separately."""
    intervals_by_class_unit: dict[tuple[str, str], list[tuple[int, int]]] = {}
    all_intervals_by_unit: dict[str, list[tuple[int, int]]] = {}
    nested_by_class: dict[str, int] = {}
    nested_total = 0
    for feature in projected:
        repeat_class = _repeat_class(feature)
        if repeat_class is None or feature.end0 <= feature.start0:
            continue
        interval = (feature.start0, feature.end0)
        intervals_by_class_unit.setdefault((repeat_class, feature.unit_id), []).append(interval)
        all_intervals_by_unit.setdefault(feature.unit_id, []).append(interval)
        if _is_nested_repeat(feature):
            nested_total += 1
            nested_by_class[repeat_class] = nested_by_class.get(repeat_class, 0) + 1

    classes = sorted({repeat_class for repeat_class, _unit_id in intervals_by_class_unit})
    result: dict[str, int] = {}
    for repeat_class in classes:
        result[f"{repeat_class}_bp"] = sum(
            _union_length(intervals)
            for (candidate_class, _unit_id), intervals in intervals_by_class_unit.items()
            if candidate_class == repeat_class
        )
        result[f"{repeat_class}_nested_hits"] = nested_by_class.get(repeat_class, 0)
    result["total_repeat_bp"] = sum(_union_length(intervals) for intervals in all_intervals_by_unit.values())
    result["nested_hit_count"] = nested_total
    return result


def _plot_modules():
    """Import Matplotlib with a server-safe backend at first plotting use."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import patches

    return plt, patches


def _plot_file_is_complete(path: Path, kind: str) -> bool:
    """Check the container markers that are lost when a plot write is truncated."""
    path = Path(path)
    if not path.is_file() or path.stat().st_size < 12:
        return False
    with path.open("rb") as handle:
        if kind == "png":
            if handle.read(8) != b"\x89PNG\r\n\x1a\n":
                return False
            handle.seek(-12, 2)
            return handle.read(12) == b"\x00\x00\x00\x00IEND\xaeB\x60\x82"
        if kind == "pdf":
            if not handle.read(5).startswith(b"%PDF"):
                return False
            handle.seek(max(path.stat().st_size - 1024, 0))
            return b"%%EOF" in handle.read()
    raise ValueError("plot kind must be 'pdf' or 'png'")


def _save_plot_pair(fig, pdf_path: Path, png_path: Path, dpi: int) -> None:
    if dpi < 1:
        raise ValueError("plot dpi must be positive")
    pdf_path, png_path = Path(pdf_path), Path(png_path)
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    png_path.parent.mkdir(parents=True, exist_ok=True)
    pdf_tmp = pdf_path.with_name(f".{pdf_path.name}.tmp")
    png_tmp = png_path.with_name(f".{png_path.name}.tmp")
    try:
        fig.savefig(pdf_tmp, format="pdf", bbox_inches="tight")
        fig.savefig(png_tmp, format="png", dpi=dpi, bbox_inches="tight")
        if not _plot_file_is_complete(pdf_tmp, "pdf") or not _plot_file_is_complete(png_tmp, "png"):
            raise RuntimeError("plot renderer produced an incomplete PDF/PNG file")
        pdf_tmp.replace(pdf_path)
        png_tmp.replace(png_path)
    finally:
        pdf_tmp.unlink(missing_ok=True)
        png_tmp.unlink(missing_ok=True)


def _bounded_png_dpi(fig, requested_dpi: int, max_pixels: int) -> int:
    """Cap raster resolution while leaving the paired PDF fully vectorized."""
    if requested_dpi < 1 or max_pixels < 1:
        raise ValueError("PNG dpi and pixel budget must be positive")
    width, height = fig.get_size_inches()
    if width <= 0 or height <= 0:
        raise ValueError("figure dimensions must be positive")
    maximum_dpi = max(1, math.floor(math.sqrt(max_pixels / (width * height))))
    return min(requested_dpi, maximum_dpi)


def repeat_plot_category(feature: Feature | ProjectedFeature) -> str:
    """Map every known repeat annotation to one broad, drawable category.

    This function is called only for records originating from the repeat GFF.
    Unknown or newly introduced repeat types deliberately fall back to a visible
    ``other/unknown repeat`` category instead of being silently discarded.
    """
    raw = " ".join(
        [feature.feature_type, feature.source]
        + [str(value) for value in feature.attributes.values()]
    ).lower().replace("-", "_")
    if "simple" in raw and "repeat" in raw:
        return "simple repeat"
    if "low" in raw and "complex" in raw:
        return "low complexity"
    if "helitron" in raw or "rc/" in raw or "rolling_circle" in raw:
        return "RC/Helitron"
    if "sine" in raw:
        return "SINE"
    if "line" in raw or "l1" in raw:
        return "LINE/L1"
    if "ltr" in raw or any(token in raw for token in ("erv", "gypsy", "copia")):
        return "LTR"
    if "satellite" in raw:
        return "satellite"
    if (
        "dna/" in raw
        or "dna_transpos" in raw
        or "tir" in raw
        or any(token in raw for token in ("tcmar", "piggybac", "hat", "cacta", "mule"))
    ):
        return "DNA transposon"
    return "other/unknown repeat"


def _observed_repeat_plot_categories(
    features: list[Feature] | list[ProjectedFeature],
) -> list[str]:
    observed = {repeat_plot_category(feature) for feature in features}
    return [category for category in REPEAT_PLOT_ORDER if category in observed]


def _raw_repeat_label(feature: Feature | ProjectedFeature) -> str | None:
    """Backward-compatible classifier for callers that have not separated repeats."""
    return repeat_plot_category(feature) if _repeat_class(feature) is not None else None


def _feature_symbol(feature: Feature | ProjectedFeature) -> str:
    for key in ("gene_symbol", "Name", "gene_name", "gene", "locus_tag", "ID"):
        value = str(feature.attributes.get(key, "")).strip()
        if value and value != ".":
            return value
    start0 = (
        feature.genomic_start0
        if isinstance(feature, ProjectedFeature) and feature.genomic_start0 is not None
        else feature.start0
    )
    end0 = (
        feature.genomic_end0
        if isinstance(feature, ProjectedFeature) and feature.genomic_end0 is not None
        else feature.end0
    )
    return f"unnamed_gene_{start0 + 1}_{end0}"


def _is_top_level_coding_gene(feature: Feature | ProjectedFeature) -> bool:
    kind = feature.feature_type.lower()
    return kind == "gene" or (kind.endswith("_gene") and "pseudo" not in kind)


def _coding_gene_symbol(feature: Feature | ProjectedFeature) -> str:
    return str(_feature_symbol(feature)).strip().upper() or "GENE"


def _gene_color_candidates():
    """Yield a finite deterministic pool, beginning with accessible named colors."""
    yield from GENE_COLOR_PALETTE
    for lightness in (0.42, 0.55, 0.32, 0.68):
        for saturation in (0.72, 0.88, 0.56):
            for hue_degree in range(360):
                red, green, blue = colorsys.hls_to_rgb(
                    hue_degree / 360.0, lightness, saturation
                )
                yield (
                    f"#{round(red * 255):02X}"
                    f"{round(green * 255):02X}"
                    f"{round(blue * 255):02X}"
                )


def build_coding_gene_color_map(
    features: list[Feature] | list[ProjectedFeature],
) -> dict[str, str]:
    """Assign every observed coding-gene symbol a deterministic unique color."""
    symbols = sorted(
        {
            _coding_gene_symbol(feature)
            for feature in features
            if _is_top_level_coding_gene(feature)
        }
    )
    colors: dict[str, str] = {
        symbol: GENE_COLOR_OVERRIDES[symbol]
        for symbol in symbols
        if symbol in GENE_COLOR_OVERRIDES
    }
    used = set(colors.values())
    candidates = iter(_gene_color_candidates())
    for symbol in symbols:
        if symbol in colors:
            continue
        for color in candidates:
            if color not in used:
                break
        else:
            raise ValueError(
                f"too many distinct coding genes to assign unique colors: {len(symbols)}"
            )
        colors[symbol] = color
        used.add(color)
    return colors


def stainedglass_window_support(
    hits: list[StainedGlassHit], region: Region, min_identity: float = 95.0
) -> list[dict[str, int | float]]:
    """Summarize high-identity off-diagonal links per genomic window."""
    counts: Counter[tuple[int, int]] = Counter()
    identities: dict[tuple[int, int], float] = {}
    for hit in hits:
        if hit.identity < min_identity or hit.query_start0 == hit.reference_start0:
            continue
        for key in ((hit.query_start0, hit.query_end0), (hit.reference_start0, hit.reference_end0)):
            if key[0] < region.end1 and key[1] > region.start1 - 1:
                counts[key] += 1
                identities[key] = max(identities.get(key, 0.0), hit.identity)
    return [
        {"start0": start0, "end0": end0, "high_identity_links": counts[(start0, end0)], "max_identity": identities[(start0, end0)]}
        for start0, end0 in sorted(counts)
    ]


def _add_compact_legend(ax, patches, labels: list[tuple[str, str]], columns: int = 4) -> None:
    handles = [patches.Patch(facecolor=color, edgecolor="#333333", label=label) for label, color in labels]
    ax.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.18),
        ncol=columns,
        frameon=False,
        fontsize=7,
        handlelength=1.2,
        columnspacing=1.0,
    )


def plot_region_overview(data: dict[str, object], pdf_path: Path, png_path: Path, dpi: int = 600) -> None:
    """Draw the array on a true genomic axis with units and source annotations."""
    plt, patches = _plot_modules()
    region = data["region"]
    units = list(data.get("units", []))
    features = list(data.get("features", []))
    repeat_features = list(data.get("repeat_features", []))
    if not repeat_features:
        repeat_features = [feature for feature in features if _raw_repeat_label(feature) is not None]
    family_hits = list(data.get("family_hits", []))
    stainedglass_windows = list(data.get("stainedglass_windows", []))
    gaps = list(data.get("gaps", []))
    gene_colors = dict(data.get("gene_colors") or build_coding_gene_color_map(features))
    observed_repeat_categories = _observed_repeat_plot_categories(repeat_features)
    repeat_lane_height = 0.30 / max(1, len(observed_repeat_categories))
    repeat_lane_y = {
        category: 0.20 + index * repeat_lane_height
        for index, category in enumerate(observed_repeat_categories)
    }
    fig, ax = plt.subplots(figsize=(11.5, 4.2))

    for boundary0 in data.get("coarse_boundaries", []):
        ax.axvline(int(boundary0) + 0.5, color="#AFAFAF", linewidth=0.45, linestyle="--", zorder=0)

    for gap in gaps:
        start0, end0 = int(gap["start1"]) - 1, int(gap["end1"])
        interval_type = str(gap["interval_type"])
        color = (
            PLOT_COLORS["repeat_rich_gap"]
            if interval_type == "repeat_rich_gap"
            else PLOT_COLORS["unassigned_flank"]
        )
        ax.add_patch(
            patches.Rectangle(
                (start0 + 0.5, 2.92),
                end0 - start0,
                0.66,
                facecolor=color,
                edgecolor="none",
                alpha=0.38,
                zorder=0,
            )
        )

    for unit in units:
        color = PLOT_COLORS["forward_unit" if unit.strand == "+" else "reverse_unit"]
        width = unit.end0 - unit.start0
        head = min(max(width * 0.16, 1), max(width / 2, 1))
        start = unit.start0 + 0.5 if unit.strand == "+" else unit.end0 + 0.5
        dx = width if unit.strand == "+" else -width
        arrow = patches.FancyArrow(
            start,
            3.25,
            dx,
            0,
            width=0.44,
            head_width=0.72,
            head_length=head,
            length_includes_head=True,
            facecolor=color,
            edgecolor="#333333",
            linewidth=0.7,
            linestyle="--" if unit.status != "complete" else "-",
        )
        ax.add_patch(arrow)

    for hit in family_hits:
        ax.add_patch(
            patches.Rectangle(
                (hit.start0 + 0.5, 2.33),
                hit.end0 - hit.start0,
                0.14 if stainedglass_windows else 0.36,
                facecolor=PLOT_COLORS["helitron_body"],
                edgecolor="none",
                alpha=0.9,
            )
        )

    maximum_links = max((int(row["high_identity_links"]) for row in stainedglass_windows), default=1)
    for row in stainedglass_windows:
        height_scale = 0.22 if family_hits else 0.35
        height = 0.10 + height_scale * int(row["high_identity_links"]) / maximum_links
        ax.add_patch(
            patches.Rectangle(
                (int(row["start0"]) + 0.5, 2.49 if family_hits else 2.35),
                int(row["end0"]) - int(row["start0"]),
                height,
                facecolor=PLOT_COLORS["stainedglass"],
                edgecolor="none",
                alpha=0.9,
            )
        )

    # Draw long repeat bodies first, then shorter/nested annotations on top.
    # Every record from the repeat GFF maps to a visible fallback category.
    for feature in sorted(
        repeat_features,
        key=lambda item: (-(item.end0 - item.start0), item.start0, item.end0),
    ):
        if feature.chrom != region.chrom or feature.end0 <= region.start1 - 1 or feature.start0 >= region.end1:
            continue
        start0, end0 = max(feature.start0, region.start1 - 1), min(feature.end0, region.end1)
        repeat_label = repeat_plot_category(feature)
        ax.add_patch(
            patches.Rectangle(
                (start0 + 0.5, repeat_lane_y[repeat_label]),
                end0 - start0,
                repeat_lane_height,
                facecolor=PLOT_COLORS[repeat_label],
                edgecolor="#333333",
                linewidth=0.25,
            )
        )

    for feature in features:
        if feature.chrom != region.chrom or feature.end0 <= region.start1 - 1 or feature.start0 >= region.end1:
            continue
        start0, end0 = max(feature.start0, region.start1 - 1), min(feature.end0, region.end1)
        kind = feature.feature_type.lower()
        if kind in {"exon", "cds"}:
            y, height, color = 0.75, 0.28, PLOT_COLORS["exon"]
        elif "pseudogene" in kind or _feature_is_pseudogene(feature):
            y, height, color = 1.25, 0.45, PLOT_COLORS["pseudogene"]
        elif _is_top_level_coding_gene(feature):
            y, height = 1.25, 0.45
            color = gene_colors[_coding_gene_symbol(feature)]
        else:
            continue
        ax.add_patch(
            patches.Rectangle(
                (start0 + 0.5, y),
                end0 - start0,
                height,
                facecolor=color,
                edgecolor="#333333",
                linewidth=0.4,
            )
        )

    label_step = max(1, math.ceil(len(units) / 12))
    for index, unit in enumerate(units):
        if index % label_step == 0 or index == len(units) - 1:
            ax.text(
                unit.start0 + 1,
                3.72,
                f"{unit.start0 + 1:,}",
                rotation=45,
                ha="left",
                va="bottom",
                fontsize=6,
            )
    ax.set_xlim(region.start1 - 0.5, region.end1 + 0.5)
    ax.set_ylim(0, 4.35)
    evidence_label = (
        "family + StainedGlass"
        if family_hits and stainedglass_windows
        else "family body"
        if family_hits
        else "StainedGlass links"
    )
    ax.set_yticks([0.35, 0.9, 1.48, 2.53, 3.25], ["repeats", "exons", "genes", evidence_label, "units"])
    ax.set_xlabel(f"{region.chrom} genomic coordinate (1-based)")
    ax.set_title(
        "ChrY ampliconic complex: refined units (dashed = StainedGlass coarse grid)"
        if data.get("coarse_boundaries")
        else "ChrY ampliconic array overview"
    )
    ax.ticklabel_format(axis="x", style="plain", useOffset=False)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.grid(axis="x", color="#DDDDDD", linewidth=0.4, zorder=0)
    legend_labels = [
        ("forward unit", PLOT_COLORS["forward_unit"]),
        ("reverse unit", PLOT_COLORS["reverse_unit"]),
    ]
    legend_labels.extend((symbol, gene_colors[symbol]) for symbol in sorted(gene_colors))
    legend_labels.extend(
        [
            ("pseudogene", PLOT_COLORS["pseudogene"]),
            ("exon/CDS", PLOT_COLORS["exon"]),
        ]
    )
    legend_labels.extend(
        (category, PLOT_COLORS[category])
        for category in observed_repeat_categories
    )
    if family_hits:
        legend_labels.append(("focal repeat family", PLOT_COLORS["helitron_body"]))
    if stainedglass_windows:
        legend_labels.append(("StainedGlass high-identity links", PLOT_COLORS["stainedglass"]))
    if gaps:
        legend_labels.extend(
            [
                ("repeat-rich inter-unit gap", PLOT_COLORS["repeat_rich_gap"]),
                ("unassigned flank/gap", PLOT_COLORS["unassigned_flank"]),
            ]
        )
    _add_compact_legend(
        ax,
        patches,
        legend_labels,
        5,
    )
    _save_plot_pair(fig, pdf_path, png_path, dpi)
    plt.close(fig)


def _unit_gene_color_key(feature: Feature | ProjectedFeature) -> str:
    """Return the subtype-plot color key for a top-level coding gene."""
    symbol = _feature_symbol(feature).upper()
    if symbol.startswith("TSPY"):
        return "coding_tspy"
    if symbol.startswith("CCDC71L"):
        return "coding_ccdc71l"
    if symbol.startswith("CUL4BY"):
        return "coding_cul4by"
    return "coding_gene"


def _coding_gene_copy_profile(
    unit_id: str,
    projected: list[ProjectedFeature],
) -> Counter[str]:
    """Count distinct top-level coding genes for representative selection."""
    profile: Counter[str] = Counter()
    seen: set[tuple[object, ...]] = set()
    for feature in projected:
        kind = feature.feature_type.lower()
        if feature.unit_id != unit_id or not (
            kind == "gene" or (kind.endswith("_gene") and "pseudo" not in kind)
        ):
            continue
        symbol = _feature_symbol(feature).upper()
        feature_key: tuple[object, ...] = (
            feature.attributes.get("ID") or "",
            symbol,
            feature.genomic_start0 if feature.genomic_start0 is not None else feature.start0,
            feature.genomic_end0 if feature.genomic_end0 is not None else feature.end0,
        )
        if feature_key in seen:
            continue
        seen.add(feature_key)
        profile[symbol] += 1
    return profile


def select_subtype_representatives(
    units: list[Unit],
    subtypes: dict[str, str],
    projected: list[ProjectedFeature] | None = None,
) -> list[Unit]:
    """Choose complete representatives with modal gene content, then median length."""
    projected = projected or []
    grouped: dict[str, list[Unit]] = {}
    for unit in units:
        grouped.setdefault(subtypes.get(unit.unit_id, "unclassified"), []).append(unit)
    representatives: list[Unit] = []
    for subtype in dict.fromkeys(subtypes.get(unit.unit_id, "unclassified") for unit in units):
        candidates = grouped[subtype]
        complete = [unit for unit in candidates if unit.status == "complete"]
        eligible = complete or candidates
        median_length = statistics.median(unit.end0 - unit.start0 for unit in eligible)
        profiles = {
            unit.unit_id: _coding_gene_copy_profile(unit.unit_id, projected)
            for unit in eligible
        }
        symbols = sorted({symbol for profile in profiles.values() for symbol in profile})
        modal_counts: dict[str, int] = {}
        for symbol in symbols:
            frequencies = Counter(profiles[unit.unit_id][symbol] for unit in eligible)
            maximum_frequency = max(frequencies.values())
            modal_counts[symbol] = min(
                count for count, frequency in frequencies.items() if frequency == maximum_frequency
            )

        def representative_key(unit: Unit) -> tuple[float, float, str]:
            copy_distance = sum(
                abs(profiles[unit.unit_id][symbol] - modal_counts[symbol])
                for symbol in symbols
            )
            return (
                float(copy_distance),
                abs((unit.end0 - unit.start0) - median_length),
                unit.unit_id,
            )

        representatives.append(
            min(eligible, key=representative_key)
        )
    return representatives


def plot_unit_subtypes(data: dict[str, object], pdf_path: Path, png_path: Path, dpi: int = 600) -> None:
    """Draw representative, all, or single oriented units on relative coordinates."""
    plt, patches = _plot_modules()
    units = list(data.get("units", []))
    projected = list(data.get("projected", []))
    projected_repeats = list(data.get("projected_repeats", []))
    if not projected_repeats:
        projected_repeats = [feature for feature in projected if _raw_repeat_label(feature) is not None]
    subtypes = dict(data.get("subtypes", {}))
    family_hits = list(data.get("family_hits", []))
    phase_offsets = dict(data.get("phase_offsets", {}))
    plot_mode = str(data.get("unit_plot_mode", "representatives"))
    if plot_mode == "representatives":
        representatives = select_subtype_representatives(
            units,
            subtypes,
            projected,
        )
    elif plot_mode in {"all", "single"}:
        representatives = units
    else:
        raise ValueError("unit_plot_mode must be representatives, all, or single")
    if not representatives and units:
        representatives = [units[0]]
    height = max(3.2, 1.15 * max(len(representatives), 1) + 1.8)
    fig, ax = plt.subplots(figsize=(11.5, height))
    max_length = max((unit.end0 - unit.start0 for unit in representatives), default=1)
    representative_ids = {unit.unit_id for unit in representatives}
    representative_repeats = [
        feature for feature in projected_repeats if feature.unit_id in representative_ids
    ]
    observed_repeat_categories = _observed_repeat_plot_categories(
        representative_repeats
    )
    repeat_legend_source = (
        projected_repeats
        if plot_mode == "representatives"
        else representative_repeats
    )
    array_repeat_category_counts = Counter(
        repeat_plot_category(feature) for feature in repeat_legend_source
    )
    array_repeat_categories = [
        category
        for category in REPEAT_PLOT_ORDER
        if category in array_repeat_category_counts
    ]
    repeat_lane_height = 0.36 / max(1, len(observed_repeat_categories))
    repeat_lane_y_offset = {
        category: -0.18 + index * repeat_lane_height
        for index, category in enumerate(observed_repeat_categories)
    }
    gene_colors = dict(data.get("gene_colors") or build_coding_gene_color_map(projected))
    displayed_gene_symbols = sorted(
        {
            _coding_gene_symbol(feature)
            for feature in projected
            if feature.unit_id in representative_ids
            and _is_top_level_coding_gene(feature)
        }
    )

    for row, unit in enumerate(representatives):
        y = len(representatives) - row
        length = unit.end0 - unit.start0
        ax.plot([0, length], [y, y], color="#777777", linewidth=0.8, zorder=0)
        for hit in family_hits:
            start0, end0 = max(hit.start0, unit.start0), min(hit.end0, unit.end0)
            if start0 >= end0:
                continue
            if unit.strand == "+":
                relative_start, relative_end = start0 - unit.start0, end0 - unit.start0
            else:
                relative_start, relative_end = unit.end0 - end0, unit.end0 - start0
            for phased_start, phased_end in circular_rephase_intervals(
                relative_start,
                relative_end,
                length,
                int(phase_offsets.get(unit.unit_id, 0)),
            ):
                ax.add_patch(
                    patches.Rectangle(
                        (phased_start, y - 0.24),
                        phased_end - phased_start,
                        0.48,
                        facecolor="none",
                        edgecolor=PLOT_COLORS["helitron_body"],
                        linewidth=1.2,
                        fill=False,
                        zorder=4,
                    )
                )
        unit_repeats = sorted(
            (feature for feature in projected_repeats if feature.unit_id == unit.unit_id),
            key=lambda item: (-(item.end0 - item.start0), item.start0, item.end0),
        )
        for feature in unit_repeats:
            if feature.end0 <= feature.start0:
                continue
            repeat_label = repeat_plot_category(feature)
            ax.add_patch(
                patches.Rectangle(
                    (feature.start0, y + repeat_lane_y_offset[repeat_label]),
                    feature.end0 - feature.start0,
                    repeat_lane_height,
                    facecolor=PLOT_COLORS[repeat_label],
                    edgecolor="#333333",
                    linewidth=0.25,
                )
            )
        for feature in projected:
            if feature.unit_id != unit.unit_id:
                continue
            kind = feature.feature_type.lower()
            if kind in {"exon", "cds"}:
                color, vertical, feature_height = PLOT_COLORS["exon"], y - 0.14, 0.28
            elif "pseudogene" in kind or _feature_is_pseudogene(
                Feature(".", feature.source, feature.feature_type, feature.start0, feature.end0, feature.strand, feature.phase, feature.attributes)
            ):
                color, vertical, feature_height = PLOT_COLORS["pseudogene"], y - 0.28, 0.56
            elif _is_top_level_coding_gene(feature):
                color = gene_colors[_coding_gene_symbol(feature)]
                vertical, feature_height = y - 0.28, 0.56
            else:
                continue
            ax.add_patch(
                patches.Rectangle(
                    (feature.start0, vertical),
                    feature.end0 - feature.start0,
                    feature_height,
                    facecolor=color,
                    edgecolor="#333333",
                    linewidth=0.35,
                )
            )
            if kind in {"gene", "pseudogene"} or kind.endswith("_gene"):
                label_color = (
                    PLOT_COLORS["pseudogene"]
                    if "pseudogene" in kind
                    else gene_colors[_coding_gene_symbol(feature)]
                )
                ax.text(
                    (feature.start0 + feature.end0) / 2,
                    y + 0.34,
                    _feature_symbol(feature),
                    rotation=45,
                    ha="left",
                    va="bottom",
                    fontsize=6,
                    color=label_color,
                    fontweight="bold",
                )
        if plot_mode == "representatives":
            row_label = f"{subtypes.get(unit.unit_id, 'unclassified')} ({unit.unit_id})"
        elif plot_mode == "all":
            row_label = (
                f"{unit.unit_id} | {subtypes.get(unit.unit_id, 'unclassified')} | "
                f"{unit.chrom}:{unit.start0 + 1:,}-{unit.end0:,} | {unit.strand} | {unit.status}"
            )
        else:
            row_label = f"{subtypes.get(unit.unit_id, 'unclassified')} | {unit.status}"
        ax.text(-0.02 * max_length, y, row_label, ha="right", va="center", fontsize=7)

    ax.set_xlim(0, max_length)
    ax.set_ylim(0.35, max(len(representatives), 1) + 0.75)
    ax.set_yticks([])
    ax.set_xlabel("Oriented unit coordinate (bp)")
    if plot_mode == "all":
        title = "All oriented ampliconic units"
    elif plot_mode == "single" and representatives:
        unit = representatives[0]
        title = (
            f"{unit.unit_id}: {unit.chrom}:{unit.start0 + 1:,}-{unit.end0:,} "
            f"({unit.strand}; {unit.status})"
        )
    else:
        title = (
            "Representative ampliconic structural classes"
            if str(data.get("analysis_mode", "")).startswith("stainedglass")
            or not family_hits
            else "Representative ampliconic unit subtypes"
        )
    ax.set_title(title)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.grid(axis="x", color="#DDDDDD", linewidth=0.4, zorder=0)
    subtype_legend = [
        *((symbol, gene_colors[symbol]) for symbol in displayed_gene_symbols),
        ("pseudogene", PLOT_COLORS["pseudogene"]),
        ("exon/CDS", PLOT_COLORS["exon"]),
    ]
    subtype_legend.extend(
        (
            category
            if category in observed_repeat_categories
            else f"{category} (other units; n={array_repeat_category_counts[category]})",
            PLOT_COLORS[category],
        )
        for category in array_repeat_categories
    )
    if family_hits:
        subtype_legend.insert(0, ("focal repeat family", PLOT_COLORS["helitron_body"]))
    _add_compact_legend(
        ax,
        patches,
        subtype_legend,
        4,
    )
    output_dpi = (
        _bounded_png_dpi(fig, dpi, ALL_UNITS_MAX_PNG_PIXELS)
        if plot_mode == "all"
        else dpi
    )
    _save_plot_pair(fig, pdf_path, png_path, output_dpi)
    plt.close(fig)


def plot_all_units(
    data: dict[str, object],
    pdf_path: Path,
    png_path: Path,
    unit_plot_dir: Path,
    prefix: str,
    dpi: int = 600,
) -> list[Path]:
    """Write one combined all-unit figure plus one PDF/PNG pair per unit."""
    prefix = validate_output_prefix(prefix)
    units = list(data.get("units", []))
    if not units:
        raise ValueError("cannot plot all units when no units are available")
    safe_unit_ids = [
        re.sub(r"[^A-Za-z0-9._-]+", "_", unit.unit_id) for unit in units
    ]
    if len(set(safe_unit_ids)) != len(safe_unit_ids):
        raise ValueError("colliding unit plot filenames after unit-ID sanitization")
    global_gene_colors = dict(
        data.get("gene_colors")
        or build_coding_gene_color_map(list(data.get("projected", [])))
    )
    combined_data = dict(data)
    combined_data["unit_plot_mode"] = "all"
    combined_data["gene_colors"] = global_gene_colors
    plot_unit_subtypes(combined_data, pdf_path, png_path, dpi=dpi)

    unit_plot_dir = Path(unit_plot_dir)
    unit_plot_dir.mkdir(parents=True, exist_ok=True)
    generated: list[Path] = []
    for unit, safe_unit_id in zip(units, safe_unit_ids):
        unit_pdf = unit_plot_dir / f"{prefix}.{safe_unit_id}.structure.pdf"
        unit_png = unit_plot_dir / f"{prefix}.{safe_unit_id}.structure.png"
        single_data = dict(data)
        single_data["units"] = [unit]
        single_data["unit_plot_mode"] = "single"
        single_data["gene_colors"] = global_gene_colors
        plot_unit_subtypes(single_data, unit_pdf, unit_png, dpi=dpi)
        generated.extend([unit_pdf, unit_png])

    incomplete = [
        path
        for path in generated
        if not _plot_file_is_complete(path, path.suffix.lstrip("."))
    ]
    if incomplete:
        raise RuntimeError(f"individual unit plots are incomplete: {incomplete}")

    expected_names = {path.name for path in generated}
    owned_name = re.compile(
        rf"^{re.escape(prefix)}\.[A-Za-z0-9._-]+\.structure\.(?:pdf|png)$"
    )
    for existing in unit_plot_dir.iterdir():
        if (
            existing.is_file()
            and owned_name.fullmatch(existing.name)
            and existing.name not in expected_names
        ):
            existing.unlink()
    return generated


def plot_hle_boundary_evidence(data: dict[str, object], pdf_path: Path, png_path: Path, dpi: int = 600) -> None:
    """Draw original/corrected family boundaries and conservative HLE evidence."""
    plt, patches = _plot_modules()
    row = dict(data.get("hle_row", {}))
    fig, ax = plt.subplots(figsize=(11.5, 4.2))
    original_start = int(row.get("original_start1", 1))
    original_end = int(row.get("original_end1", original_start + 1))
    left_offset = int(row.get("corrected_start_offset", 0))
    right_offset = int(row.get("corrected_end_offset", 0))
    corrected_start = original_start + left_offset
    corrected_end = original_end - right_offset
    span = max(original_end - original_start + 1, 1)
    ax.plot([original_start, original_end], [2.8, 2.8], color="#777777", linewidth=7, solid_capstyle="butt")
    ax.plot([corrected_start, corrected_end], [1.9, 1.9], color=PLOT_COLORS["helitron_body"], linewidth=12, solid_capstyle="butt")
    ax.scatter([original_start, original_end], [2.8, 2.8], marker="|", s=240, color="#222222")
    ax.scatter([corrected_start, corrected_end], [1.9, 1.9], marker="|", s=280, color="#222222")
    ax.text(original_start, 3.12, f"original 5′ {original_start:,}", ha="left", fontsize=8)
    ax.text(original_end, 3.12, f"original 3′ {original_end:,}", ha="right", fontsize=8)
    ax.text(corrected_start, 2.15, f"+{left_offset} bp  {row.get('five_prime_motif') or 'motif absent'}", ha="left", fontsize=8)
    ax.text(corrected_end, 2.15, f"{row.get('three_prime_motif') or 'motif absent'}  −{right_offset} bp", ha="right", fontsize=8)
    hairpin = "–".join(
        str(row.get(key) or "?")
        for key in ("hairpin_stem_left", "hairpin_loop", "hairpin_stem_right")
    )
    ax.text((corrected_start + corrected_end) / 2, 1.45, f"candidate hairpin: {hairpin}", ha="center", fontsize=9)
    missing: list[str] = []
    if not row.get("rep_domain"):
        missing.append("Rep/HUH domain")
    if not row.get("helicase_domain"):
        missing.append("helicase domain")
    if int(row.get("independent_loci", 0)) == 0:
        missing.append("independent insertion locus")
    ax.text(
        original_start,
        0.7,
        f"Classification: {row.get('classification', 'not classified')}\nMissing autonomous evidence: {', '.join(missing) if missing else 'none'}",
        ha="left",
        va="center",
        fontsize=9,
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "#F2F2F2", "edgecolor": "#999999"},
    )
    padding = max(span * 0.04, 1)
    ax.set_xlim(original_start - padding, original_end + padding)
    ax.set_ylim(0.1, 3.65)
    ax.set_yticks([1.9, 2.8], ["corrected", "Earl Grey original"])
    ax.set_xlabel("Representative family-body coordinate (1-based)")
    ax.set_title("Helitron-like boundary evidence")
    ax.ticklabel_format(axis="x", style="plain", useOffset=False)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.grid(axis="x", color="#DDDDDD", linewidth=0.4, zorder=0)
    _save_plot_pair(fig, pdf_path, png_path, dpi)
    plt.close(fig)


def plot_stainedglass_period_evidence(
    period: dict[str, object], pdf_path: Path, png_path: Path, dpi: int = 600
) -> None:
    """Plot the regional high-identity lag spectrum and selected fundamental period."""
    plt, _patches = _plot_modules()
    rows = list(period.get("lag_rows", []))
    fig, ax = plt.subplots(figsize=(11.5, 4.2))
    x = [int(row["lag_bp"]) for row in rows]
    y = [float(row["normalized_support"]) for row in rows]
    ax.plot(x, y, color="#4C78A8", linewidth=0.8)
    ax.scatter(x, y, color="#4C78A8", s=7)
    selected = int(float(period["median"]))
    region_limit = max(x, default=selected)
    for harmonic in range(1, region_limit // selected + 1):
        ax.axvline(
            selected * harmonic,
            color=PLOT_COLORS["stainedglass"],
            linewidth=1.2 if harmonic == 1 else 0.5,
            alpha=0.9 if harmonic == 1 else 0.35,
        )
    ax.text(selected, max(y, default=1.0), f" selected {selected:,} bp", ha="left", va="top", fontsize=8)
    ax.set_xlim(0, region_limit * 1.01 if region_limit else selected * 1.01)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Window-pair genomic lag (bp)")
    ax.set_ylabel("Normalized ≥95% identity support")
    ax.set_title("StainedGlass repeat-period evidence")
    ax.ticklabel_format(axis="x", style="plain", useOffset=False)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(color="#DDDDDD", linewidth=0.4)
    _save_plot_pair(fig, pdf_path, png_path, dpi)
    plt.close(fig)


UNIT_COLUMNS = [
    "unit_id",
    "chrom",
    "start1",
    "end1",
    "strand",
    "arm",
    "status",
    "length_bp",
    "family_hit_count",
    "family_body_bp",
    "period_deviation_bp",
    "tspy_phase_A",
    "tspy_phase_B",
    "subtype",
    "anchor_support",
    "sequence_record",
]

HLE_COLUMNS = [
    "candidate_id",
    "original_start1",
    "original_end1",
    "corrected_start_offset",
    "corrected_end_offset",
    "five_prime_motif",
    "three_prime_motif",
    "five_prime_TC_support",
    "three_prime_CTRR_support",
    "hairpin_support",
    "hairpin_stem_left",
    "hairpin_loop",
    "hairpin_stem_right",
    "target_context",
    "longest_orf_aa",
    "rep_domain",
    "helicase_domain",
    "independent_loci",
    "classification",
]

SUBTYPE_COLUMNS = [
    "unit_id",
    "tspy_phase_A",
    "tspy_phase_B",
    "subtype",
    "phase_A_center_bp",
    "phase_B_center_bp",
]

REPEAT_COMPONENT_COLUMNS = [
    "unit_id",
    "feature_type",
    "repeat_class",
    "relative_start0",
    "relative_end0",
    "strand",
    "genomic_start1",
    "genomic_end1",
    "clipped",
    "nested",
    "attributes",
]

REPEAT_INVENTORY_COLUMNS = [
    "repeat_id",
    "chrom",
    "original_start1",
    "original_end1",
    "original_length_bp",
    "target_start1",
    "target_end1",
    "target_overlap_bp",
    "start1",
    "end1",
    "length_bp",
    "source",
    "feature_type",
    "repeat_class",
    "repeat_family",
    "plot_category",
    "strand",
    "score",
    "phase",
    "assigned_intervals",
    "retained_bp",
    "attributes",
    "raw_attributes",
]

STAINEDGLASS_PERIOD_COLUMNS = [
    "rank",
    "candidate_period_bp",
    "score",
    "observed_harmonics",
    "fundamental_support",
    "harmonic_support",
    "selected",
]

GENE_COMPOSITION_COLUMNS = [
    "unit_id",
    "coding_genes",
    "coding_gene_count",
    "pseudogenes",
    "pseudogene_count",
]

BOUNDARY_REFINEMENT_COLUMNS = [
    "boundary_id",
    "unit_id",
    "boundary_role",
    "coarse_boundary0",
    "sequence_target0",
    "refined_boundary0",
    "adjustment_bp",
    "method",
    "sequence_support",
    "forward_lag_bp",
    "forward_support",
    "reverse_center0",
    "reverse_support",
    "annotation_signature",
    "phase_anchor_signature",
    "annotation_recurrence",
    "annotations_crossed",
    "confidence",
]

GAP_COLUMNS = [
    "gap_id",
    "chrom",
    "start1",
    "end1",
    "length_bp",
    "interval_type",
    "repeat_bp",
    "repeat_fraction",
    "top_repeat_classes",
    "hub_windows",
]


def _raw_fasta_record(path: Path) -> tuple[str, str]:
    lines = Path(path).read_text().splitlines()
    headers = [index for index, line in enumerate(lines) if line.startswith(">")]
    if len(headers) != 1:
        raise ValueError("expected a single FASTA record")
    if any(line.strip() for line in lines[: headers[0]]):
        raise ValueError("FASTA contains nonblank preamble before record")
    title = lines[headers[0]][1:].strip()
    sequence = "".join(line.strip() for line in lines[headers[0] + 1 :] if line.strip())
    if not sequence:
        raise ValueError("FASTA sequence is empty")
    if any(base not in "ACGTNacgtn" for base in sequence):
        raise ValueError("FASTA sequence contains invalid bases")
    return title, sequence


def read_analysis_fasta(path: Path, target_region: Region, region_origin: int | None) -> tuple[str, str, Region]:
    """Read a regional FASTA, allowing only the CLI's explicit origin fallback."""
    try:
        title, sequence, fasta_region = read_single_fasta(path)
    except ValueError as exc:
        if region_origin is None or "coordinates" not in str(exc):
            raise
        raw_title, raw_sequence = _raw_fasta_record(path)
        if ":" in raw_title:
            raise
        if region_origin < 1:
            raise ValueError("region origin must be a positive 1-based coordinate") from exc
        title, sequence = raw_title, raw_sequence
        fasta_region = Region(target_region.chrom, region_origin, region_origin + len(sequence) - 1)
    if fasta_region is None:
        raise ValueError("FASTA region is unavailable")
    if (
        fasta_region.chrom != target_region.chrom
        or target_region.start1 < fasta_region.start1
        or target_region.end1 > fasta_region.end1
    ):
        raise ValueError("target region is outside the available FASTA coordinates")
    return title, sequence, fasta_region


def filter_features_to_region(features: list[Feature], region: Region) -> list[Feature]:
    """Retain source annotations on the target chromosome that overlap the target."""
    start0 = region.start1 - 1
    return [
        feature
        for feature in features
        if feature.chrom == region.chrom and feature.start0 < region.end1 and feature.end0 > start0
    ]


def _model_overlaps_region(model: GeneModel, region: Region) -> bool:
    start0 = region.start1 - 1
    return any(block_start < region.end1 and block_end > start0 for block_start, block_end, _phase in model.cds_blocks)


def _complete_coding_models(
    features: list[Feature], target_region: Region, fasta_region: Region
) -> list[GeneModel]:
    fasta_start0 = fasta_region.start1 - 1
    models = collect_gene_models([feature for feature in features if feature.chrom == target_region.chrom])
    return [
        model
        for model in models
        if (
            not model.is_pseudogene
            and model.cds_blocks
            and _model_overlaps_region(model, target_region)
            and all(fasta_start0 <= start0 < end0 <= fasta_region.end1 for start0, end0, _phase in model.cds_blocks)
        )
    ]


def _read_independent_heliano_loci(path: Path | None, region: Region, family_name: str) -> int:
    if path is None or not Path(path).is_file():
        return 0
    count = 0
    try:
        lines = Path(path).read_text().splitlines()
    except OSError:
        return 0
    region_start0 = region.start1 - 1
    for line in lines:
        if not line.strip() or line.startswith("#"):
            continue
        columns = line.split("\t")
        if len(columns) < 3:
            continue
        try:
            start0, end0 = int(columns[1]), int(columns[2])
        except ValueError:
            continue
        if end0 <= start0:
            continue
        if not any(identifier.casefold() == family_name.casefold() for identifier in columns[3:]):
            continue
        evidence_tags = {value.casefold() for value in columns[3:]}
        if not ({"distinct_flanks=true", "independent_locus=true"} & evidence_tags):
            continue
        if columns[0] != region.chrom or end0 <= region_start0 or start0 >= region.end1:
            count += 1
    return count


def build_hle_evidence(
    family_name: str,
    hits: list[RepeatHit],
    sequence: str,
    fasta_region: Region,
    heliano_bed: Path | None,
    target_region: Region | None = None,
) -> tuple[dict[str, object], dict[str, str]]:
    """Aggregate one modal paired terminal structure across near-full family copies."""
    records: dict[str, str] = {}
    record_hits: dict[str, RepeatHit] = {}
    fasta_start0 = fasta_region.start1 - 1
    evidence_hits = _primary_anchors(hits) or hits
    for index, hit in enumerate(evidence_hits, 1):
        if hit.start0 < fasta_start0 or hit.end0 > fasta_region.end1 or hit.end0 <= hit.start0:
            continue
        name = hit.hit_id or f"family_copy_{index:03d}"
        if name in records:
            name = f"{name}_{index:03d}"
        records[name] = extract_genomic(sequence, fasta_region, hit.start0 + 1, hit.end0, hit.strand)
        record_hits[name] = hit

    best_pairs: dict[str, TerminalPair] = {}
    for name, body in records.items():
        pairs = find_terminal_pairs(body)
        if pairs:
            best_pairs[name] = max(pairs, key=lambda pair: pair.score)

    copy_count = len(records)
    tc_support = 0.0
    ctrr_support = 0.0
    hairpin_support = 0.0
    representative_name = next(iter(records), None)
    representative_pair: TerminalPair | None = None
    corrected_records: dict[str, str] = {}
    if best_pairs:
        signatures = Counter(
            (
                pair.five_prime_offset,
                pair.three_prime_offset,
                pair.five_prime_motif,
                pair.three_prime_motif,
                pair.hairpin.stem_left if pair.hairpin else "",
                pair.hairpin.loop if pair.hairpin else "",
                pair.hairpin.stem_right if pair.hairpin else "",
            )
            for pair in best_pairs.values()
        )
        modal = min(signatures, key=lambda signature: (-signatures[signature], signature))
        representative_name = next(
            name
            for name, pair in best_pairs.items()
            if (
                pair.five_prime_offset,
                pair.three_prime_offset,
                pair.five_prime_motif,
                pair.three_prime_motif,
                pair.hairpin.stem_left if pair.hairpin else "",
                pair.hairpin.loop if pair.hairpin else "",
                pair.hairpin.stem_right if pair.hairpin else "",
            )
            == modal
        )
        representative_pair = best_pairs[representative_name]
        representative_structure = (
            representative_pair.five_prime_motif,
            representative_pair.three_prime_motif,
            representative_pair.hairpin.stem_left if representative_pair.hairpin else "",
            representative_pair.hairpin.loop if representative_pair.hairpin else "",
            representative_pair.hairpin.stem_right if representative_pair.hairpin else "",
        )
        modal_names = [
            name
            for name, pair in best_pairs.items()
            if (
                pair.five_prime_motif,
                pair.three_prime_motif,
                pair.hairpin.stem_left if pair.hairpin else "",
                pair.hairpin.loop if pair.hairpin else "",
                pair.hairpin.stem_right if pair.hairpin else "",
            )
            == representative_structure
        ]
        paired_support = len(modal_names) / copy_count if copy_count else 0.0
        tc_support = paired_support
        ctrr_support = paired_support
        hairpin_support = paired_support if representative_pair.hairpin else 0.0
        for name in modal_names:
            pair = best_pairs[name]
            body = records[name]
            end0 = len(body) - pair.three_prime_offset
            if pair.five_prime_offset < end0:
                corrected_records[name] = body[pair.five_prime_offset:end0]

    representative_records = corrected_records or records
    longest_orf_aa = 0
    for body in representative_records.values():
        orfs = scan_six_frame_orfs(body)
        longest_orf_aa = max(longest_orf_aa, max((orf.length_aa for orf in orfs), default=0))

    representative_hit = record_hits.get(representative_name) if representative_name else None
    independent_loci = _read_independent_heliano_loci(
        heliano_bed, target_region or fasta_region, family_name
    )
    evidence = HLEEvidence(
        terminal_support=min(tc_support, ctrr_support),
        hairpin_support=hairpin_support,
        longest_orf_aa=longest_orf_aa,
        rep_domain=False,
        helicase_domain=False,
        independent_loci=independent_loci,
    )
    hairpin = representative_pair.hairpin if representative_pair else None
    row: dict[str, object] = {
        "candidate_id": family_name,
        "original_start1": representative_hit.start0 + 1 if representative_hit else "",
        "original_end1": representative_hit.end0 if representative_hit else "",
        "corrected_start_offset": representative_pair.five_prime_offset if representative_pair else 0,
        "corrected_end_offset": representative_pair.three_prime_offset if representative_pair else 0,
        "five_prime_motif": representative_pair.five_prime_motif if representative_pair else "",
        "three_prime_motif": representative_pair.three_prime_motif if representative_pair else "",
        "five_prime_TC_support": tc_support,
        "three_prime_CTRR_support": ctrr_support,
        "hairpin_support": hairpin_support,
        "hairpin_stem_left": hairpin.stem_left if hairpin else "",
        "hairpin_loop": hairpin.loop if hairpin else "",
        "hairpin_stem_right": hairpin.stem_right if hairpin else "",
        "target_context": (
            f"5prime={representative_pair.five_prime_context};3prime={representative_pair.three_prime_context}"
            if representative_pair
            else ""
        ),
        "longest_orf_aa": longest_orf_aa,
        "rep_domain": False,
        "helicase_domain": False,
        "independent_loci": independent_loci,
        "classification": classify_hle_evidence(evidence),
    }
    return row, representative_records


def _write_tsv(path: Path, columns: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _serialize_attributes(attributes: dict[str, str]) -> str:
    return ";".join(f"{key}={attributes[key]}" for key in sorted(attributes))


def _phase_state(subtype: str, phase: str) -> str:
    prefix = f"{phase}_"
    return next((part[len(prefix) :] for part in subtype.split("__") if part.startswith(prefix)), "missing")


def _overlap_union_length(hit_intervals: list[tuple[int, int]]) -> int:
    return _union_length([(start0, end0) for start0, end0 in hit_intervals if end0 > start0])


def _resolve_executable(requested: str | Path | None) -> str | None:
    if requested is None:
        return None
    requested_text = str(requested)
    candidate = Path(requested_text)
    if candidate.is_file():
        return str(candidate)
    return shutil.which(requested_text)


def _dependency_statuses(args: argparse.Namespace, mafft_status: str) -> dict[str, str]:
    minimap_requested = getattr(args, "minimap2", None)
    if minimap_requested is None:
        minimap_status = "unused (not requested)"
    elif _resolve_executable(minimap_requested) is None:
        minimap_status = "unavailable (requested executable not found)"
    else:
        minimap_status = "unused (available; core run does not require self-alignment)"
    heliano_bed = getattr(args, "heliano_bed", None)
    if heliano_bed is not None and Path(heliano_bed).is_file():
        heliano_status = "used (readable candidate BED supplied; only explicitly validated distinct-flank loci count)"
    elif heliano_bed is not None:
        heliano_status = "unavailable (--heliano-bed is not readable)"
    else:
        heliano_status = "unavailable (no --heliano-bed supplied)"
    return {
        "MAFFT": mafft_status,
        "minimap2": minimap_status,
        "HELIANO": heliano_status,
        "TEsorter": "unused (no domain-results input; Rep/helicase evidence remains false)",
    }


def _candidate_records(records: dict[str, str], units: list[Unit]) -> dict[str, str]:
    complete = {unit.unit_id: records[unit.unit_id] for unit in units if unit.status == "complete" and unit.unit_id in records}
    return complete or dict(records)


def _summary_text(
    family_name: str,
    hle_row: dict[str, object],
    period: dict[str, object],
    units: list[Unit],
    subtypes: dict[str, str],
    dependency_statuses: dict[str, str],
    reconstruction_mode: str,
) -> str:
    arm_counts = Counter(unit.arm for unit in units)
    junction_count = sum(unit.status == "junction_partial" for unit in units)
    coding_units = sum("coding" in subtype for subtype in subtypes.values())
    pseudogene_units = sum("pseudogene" in subtype for subtype in subtypes.values())
    present: list[str] = []
    if float(hle_row["five_prime_TC_support"]) > 0:
        present.append(f"5′ TC ({float(hle_row['five_prime_TC_support']):.1%} of family-body copies)")
    if float(hle_row["three_prime_CTRR_support"]) > 0:
        present.append(f"3′ CTRR ({float(hle_row['three_prime_CTRR_support']):.1%})")
    if float(hle_row["hairpin_support"]) > 0:
        present.append(f"terminal hairpin ({float(hle_row['hairpin_support']):.1%})")
    if int(hle_row["longest_orf_aa"]) > 0:
        present.append(f"longest six-frame ORF {hle_row['longest_orf_aa']} aa")
    missing: list[str] = []
    if not hle_row["rep_domain"]:
        missing.append("Rep/HUH domain")
    if not hle_row["helicase_domain"]:
        missing.append("helicase domain")
    if int(hle_row["independent_loci"]) == 0:
        missing.append("independent insertion loci")
    missing.append("target-site duplication / distinct-flank validation")
    dependencies = "\n".join(f"- {name}: {status}" for name, status in dependency_statuses.items())
    period_value = period.get("median")
    period_label = f"{float(period_value):g} bp" if isinstance(period_value, (int, float)) else "not estimable"
    mad_value = period.get("mad")
    mad_label = f"{float(mad_value):g} bp" if isinstance(mad_value, (int, float)) else "not estimable"
    plus_label = f"{float(period['plus_median']):g} bp" if isinstance(period.get("plus_median"), (int, float)) else "not estimable"
    minus_label = f"{float(period['minus_median']):g} bp" if isinstance(period.get("minus_median"), (int, float)) else "not estimable"
    outlier_threshold = max(3.0 * float(mad_value or 0.0), 500.0)
    period_outliers = sum(
        unit.status == "complete"
        and isinstance(period_value, (int, float))
        and abs((unit.end0 - unit.start0) - float(period_value)) > outlier_threshold
        for unit in units
    )
    return (
        f"# {family_name} ChrY ampliconic-unit analysis\n\n"
        f"Classification: **{hle_row['classification']}**.\n\n"
        f"Evidence present: {', '.join(present) if present else 'none detected'}.\n\n"
        f"Evidence missing: {', '.join(missing) if missing else 'none'}.\n\n"
        f"Period and array structure: median {period_label}, MAD {mad_label}, "
        f"plus/minus medians {plus_label}/{minus_label}, {period_outliers} outlier units; {len(units)} units; "
        f"left arm {arm_counts.get('left', 0)}, right arm {arm_counts.get('right', 0)}, "
        f"junction {arm_counts.get('junction', 0)}; {junction_count} junction-partial units.\n\n"
        f"TSPY coding/pseudogene states: {coding_units} units with coding calls; "
        f"{pseudogene_units} units with pseudogene calls.\n\n"
        f"Sequence reconstruction mode: `{reconstruction_mode}`.\n\n"
        f"Optional dependency downgrades/status:\n\n{dependencies}\n\n"
        "Array copy count is not a transposition count. Tandem or palindromic copies in this array are not independent insertion events.\n"
    )


def _format_symbol_counts(counter: Counter[str]) -> str:
    return ";".join(f"{symbol}:{counter[symbol]}" for symbol in sorted(counter))


def _unit_gene_composition(units: list[Unit], projected: list[ProjectedFeature]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for unit in units:
        coding: Counter[str] = Counter()
        pseudo: Counter[str] = Counter()
        for feature in projected:
            if feature.unit_id != unit.unit_id:
                continue
            kind = feature.feature_type.lower()
            if feature.genomic_start0 is not None and feature.genomic_end0 is not None:
                genomic_midpoint = (feature.genomic_start0 + feature.genomic_end0) // 2
                if not unit.start0 <= genomic_midpoint < unit.end0:
                    continue
            if kind == "pseudogene" or _projected_gene_state(feature) == "pseudogene" and kind.endswith("_gene"):
                pseudo[_feature_symbol(feature)] += 1
            elif kind == "gene" or kind.endswith("_gene"):
                coding[_feature_symbol(feature)] += 1
        rows.append(
            {
                "unit_id": unit.unit_id,
                "coding_genes": _format_symbol_counts(coding),
                "coding_gene_count": sum(coding.values()),
                "pseudogenes": _format_symbol_counts(pseudo),
                "pseudogene_count": sum(pseudo.values()),
            }
        )
    return rows


def _family_free_summary_text(
    period: dict[str, object],
    units: list[Unit],
    structural_classes: dict[str, str],
    gene_states: dict[str, str],
    gap_rows: list[dict[str, object]],
    coding: list[Feature],
    pseudogenes: list[Feature],
    orientation_info: dict[str, object],
    boundary_rows: list[dict[str, object]],
    analysis_mode: str = "stainedglass_family_free",
    focal_family_name: str | None = None,
    focal_family_hit_count: int = 0,
) -> str:
    coding_symbols = Counter(
        _feature_symbol(feature)
        for feature in coding
        if feature.feature_type.lower() == "gene" or feature.feature_type.lower().endswith("_gene")
    )
    pseudo_symbols = Counter(
        _feature_symbol(feature) for feature in pseudogenes if feature.feature_type.lower() == "pseudogene"
    )
    complete = sum(unit.status == "complete" for unit in units)
    partial = len(units) - complete
    orientation_pattern = "".join(unit.strand for unit in units)
    coding_units = sum("coding_0" not in state and "coding" in state for state in gene_states.values())
    pseudogene_units = sum("pseudogene_0" not in state and "pseudogene" in state for state in gene_states.values())
    length_stats = summarize_refined_unit_lengths(units)
    length_scope = (
        "complete-unit"
        if length_stats["basis"] == "complete_units"
        else "available-unit (no complete units)"
    )
    confidence_counts = Counter(str(row["confidence"]) for row in boundary_rows)
    class_counts = Counter(structural_classes.values())
    inter_unit_gaps = [row for row in gap_rows if row["interval_type"] == "repeat_rich_gap"]
    flank_rows = [row for row in gap_rows if "flank" in str(row["interval_type"])]
    assigned_bp = sum(unit.end0 - unit.start0 for unit in units)
    region_bp = int(period.get("target_region_bp", assigned_bp))
    top_period = int(period["candidate_scores"][0]["candidate_period_bp"]) if period.get("candidate_scores") else int(period["median"])
    focal_family_text = (
        f"Focal repeat family: `{focal_family_name}` "
        f"({focal_family_hit_count} raw exact-name annotations), "
        "used as an overlay while StainedGlass retained control of unit discovery. "
        "All other repeat annotations were retained independently.\n\n"
        if focal_family_name
        else "No focal repeat family was supplied; all repeat annotations were retained without family filtering.\n\n"
    )
    return (
        "# ChrY repeat-rich ampliconic-complex analysis\n\n"
        f"Analysis mode: `{analysis_mode}`.\n\n"
        f"StainedGlass period: **{int(float(period['median'])):,} bp** "
        f"({period['method']}; window {int(period['window_size']):,} bp; "
        f"identity threshold {float(period['min_identity']):g}%).\n\n"
        f"The strongest higher-order lag peak is {top_period:,} bp; {int(float(period['median'])):,} bp is retained only as the smallest coherent fundamental used to seed boundary searches.\n\n"
        f"Array structure: {len(units)} discrete units ({complete} complete, {partial} boundary partial); "
        f"the requested interval is not forced to tile. Units occupy {assigned_bp:,}/{region_bp:,} bp "
        f"({assigned_bp / max(1, region_bp):.1%}); {len(inter_unit_gaps)} repeat-rich inter-unit gaps "
        f"and {len(flank_rows)} unassigned flanks are retained. "
        f"refined {length_scope} lengths {int(length_stats['minimum_bp']):,}–"
        f"{int(length_stats['maximum_bp']):,} bp "
        f"(median {int(float(length_stats['median_bp'])):,} bp). "
        f"orientation pattern `{orientation_pattern}`. Reference/root: "
        f"{orientation_info.get('reference_unit_id', 'NA')}; disconnected units: "
        f"{orientation_info.get('disconnected_units', 0)}; orientation conflicts: "
        f"{orientation_info.get('orientation_conflicts', 0)}/{orientation_info.get('orientation_edges', 0)}; "
        f"weighted agreement: {float(orientation_info.get('orientation_weight_agreement', 0.0)):.1%}.\n\n"
        f"Structural classes (orientation and TSPY state excluded): "
        f"{_format_symbol_counts(class_counts)}.\n\n"
        f"Boundary refinement confidence: high {confidence_counts.get('high', 0)}, "
        f"medium {confidence_counts.get('medium', 0)}, low {confidence_counts.get('low', 0)}.\n\n"
        f"Coding genes: {_format_symbol_counts(coding_symbols) or 'none'}.\n\n"
        f"Reliable pseudogenes: {_format_symbol_counts(pseudo_symbols) or 'none'}.\n\n"
        f"TSPY unit states: {coding_units} units with coding calls; "
        f"{pseudogene_units} units with pseudogene calls.\n\n"
        f"{focal_family_text}"
        "RC/Helitron autonomy and Helitron termini were not evaluated by the StainedGlass workflow. "
        f"StainedGlass supplied only {int(period['window_size']):,}-bp-scale coarse evidence. "
        "Sustained high-degree windows were treated as repeat hubs, not unit anchors; "
        "gap-facing ends were refined independently from exact sequence-copy support and annotation edges. "
        "Low-confidence boundaries remain provisional.\n"
    )


def _run_stainedglass_analysis(args: argparse.Namespace) -> dict[str, Path]:
    """Run StainedGlass-guided analysis, optionally with a focal-family overlay."""
    stainedglass_path = getattr(args, "stainedglass_bed", None)
    if stainedglass_path is None:
        raise ValueError("--stainedglass-bed is required for StainedGlass analysis")
    target_region = parse_region_header(str(args.target_region).replace(",", ""))
    title, sequence, fasta_region = read_analysis_fasta(
        Path(args.fasta), target_region, getattr(args, "region_origin", None)
    )
    repeat_source = read_gff(Path(args.repeat_gff))
    coding_source = read_gff(Path(args.coding_gff3))
    pseudogene_source = read_gff(Path(args.pseudogene_gff3))
    repeats = filter_features_to_region(repeat_source, target_region)
    coding = filter_features_to_region(coding_source, target_region)
    pseudogenes = filter_features_to_region(pseudogene_source, target_region)
    focal_family_name = getattr(args, "family", None)
    focal_family_hits = (
        select_focal_family_features(repeats, focal_family_name)
        if focal_family_name
        else []
    )
    if focal_family_name and not focal_family_hits:
        raise ValueError(
            f"no repeat annotations found for exact family name {focal_family_name!r}"
        )
    analysis_mode = "stainedglass_hybrid" if focal_family_name else "stainedglass_family_free"
    complete_models = _complete_coding_models(coding_source, target_region, fasta_region)
    qc = run_input_qc(sequence, fasta_region, repeats, complete_models)
    qc.update(summarize_pseudogene_qc(pseudogenes, fasta_region))
    qc["gene_models"] = int(qc["coding_models"]) + int(qc["pseudogene_models"])
    if qc["coordinate_status"] == "fail":
        raise RuntimeError(f"coding-coordinate QC failed: {qc['coordinate_reason']}")
    if qc["pseudogene_coordinate_status"] == "fail":
        raise RuntimeError("pseudogene-coordinate QC failed")

    stainedglass_hits = read_stainedglass_bed(Path(stainedglass_path), target_region)
    period = estimate_stainedglass_period(
        stainedglass_hits,
        target_region,
        period_override=getattr(args, "period_bp", None),
    )
    period["target_region_bp"] = target_region.end1 - target_region.start1 + 1
    coarse_units, coarse_gap_rows, detection_info = discover_stainedglass_unit_blocks(
        stainedglass_hits,
        target_region,
        expected_period_bp=period,
    )
    coarse_orientations, _coarse_orientation_info = infer_stainedglass_unit_orientations(
        stainedglass_hits, coarse_units
    )
    coarse_units = [
        Unit(
            unit.unit_id,
            unit.chrom,
            unit.start0,
            unit.end0,
            coarse_orientations[unit.unit_id],
            unit.arm,
            unit.status,
        )
        for unit in coarse_units
    ]
    units, boundary_rows = refine_stainedglass_units(
        sequence,
        fasta_region,
        coarse_units,
        repeats,
        coding,
        pseudogenes,
    )
    units = clip_units_to_region(units, target_region)
    if not units:
        raise ValueError("no refined StainedGlass units overlap the target region")
    boundary_rows = reconcile_boundary_rows_with_units(
        boundary_rows, units, target_region
    )
    orientations, orientation_info = infer_stainedglass_unit_orientations(stainedglass_hits, units)
    units = [
        Unit(
            unit.unit_id,
            unit.chrom,
            unit.start0,
            unit.end0,
            orientations[unit.unit_id],
            "forward" if orientations[unit.unit_id] == "+" else "reverse",
            unit.status,
        )
        for unit in units
    ]
    phase_offsets = {unit.unit_id: 0 for unit in units}
    unit_records = {
        unit.unit_id: oriented_unit_sequence(unit, sequence, fasta_region)
        for unit in units
    }
    complete_ids = {unit.unit_id for unit in units if unit.status == "complete"}
    reference_name = str(orientation_info["reference_unit_id"])
    if reference_name not in complete_ids:
        reference_name = next((unit.unit_id for unit in units if unit.status == "complete"), units[0].unit_id)
    reference_sequence = unit_records[reference_name]

    projected: list[ProjectedFeature] = []
    projected_repeats: list[ProjectedFeature] = []
    for unit in units:
        for feature in [*repeats, *coding, *pseudogenes]:
            candidate = project_feature_to_unit(feature, unit)
            if candidate is not None:
                projected.append(candidate)
                if feature in repeats:
                    projected_repeats.append(candidate)
    assigned_genes: list[ProjectedFeature] = []
    for feature in [*coding, *pseudogenes]:
        kind = feature.feature_type.lower()
        if kind not in {"gene", "pseudogene"} and not kind.endswith("_gene"):
            continue
        midpoint = (feature.start0 + feature.end0) // 2
        unit = next((candidate for candidate in units if candidate.start0 <= midpoint < candidate.end0), None)
        if unit is not None:
            candidate = project_feature_to_unit(feature, unit)
            if candidate is not None:
                assigned_genes.append(candidate)
    phase_centers: dict[str, float] = {}
    gene_states = {
        unit.unit_id: classify_family_free_tspy_subtype(unit, assigned_genes)
        for unit in units
    }
    subtypes = classify_structural_units(units, stainedglass_hits)
    projected = rephase_projected_features(projected, units, phase_offsets)
    projected_repeats = rephase_projected_features(projected_repeats, units, phase_offsets)
    gap_rows = summarize_inter_unit_gaps(units, target_region, repeats, coarse_gap_rows)

    unit_links = dict(orientation_info.get("unit_links", {}))
    maximum_links = max((int(value) for value in unit_links.values()), default=1)
    unit_rows: list[dict[str, object]] = []
    subtype_rows: list[dict[str, object]] = []
    for unit in units:
        subtype = subtypes[unit.unit_id]
        gene_state = gene_states[unit.unit_id]
        intersecting_family_hits = [
            hit
            for hit in focal_family_hits
            if hit.start0 < unit.end0 and hit.end0 > unit.start0
        ]
        family_intervals = [
            (max(hit.start0, unit.start0), min(hit.end0, unit.end0))
            for hit in intersecting_family_hits
        ]
        unit_rows.append(
            {
                "unit_id": unit.unit_id,
                "chrom": unit.chrom,
                "start1": unit.start0 + 1,
                "end1": unit.end0,
                "strand": unit.strand,
                "arm": unit.arm,
                "status": unit.status,
                "length_bp": unit.end0 - unit.start0,
                "family_hit_count": len(intersecting_family_hits),
                "family_body_bp": _overlap_union_length(family_intervals),
                "period_deviation_bp": (unit.end0 - unit.start0) - float(period["median"]),
                "tspy_phase_A": gene_state,
                "tspy_phase_B": "",
                "subtype": subtype,
                "anchor_support": int(unit_links.get(unit.unit_id, 0)) / maximum_links,
                "sequence_record": unit.unit_id,
            }
        )
        subtype_rows.append(
            {
                "unit_id": unit.unit_id,
                "tspy_phase_A": gene_state,
                "tspy_phase_B": "",
                "subtype": subtype,
                "phase_A_center_bp": phase_centers.get("A", ""),
                "phase_B_center_bp": phase_centers.get("B", ""),
            }
        )

    zero_length_repeat_projections = sum(
        feature.end0 <= feature.start0 for feature in projected_repeats
    )
    repeat_rows = [
        {
            "unit_id": feature.unit_id,
            "feature_type": feature.feature_type,
            "repeat_class": _repeat_class(feature) or feature.feature_type or "other_repeat",
            "relative_start0": feature.start0,
            "relative_end0": feature.end0,
            "strand": feature.strand,
            "genomic_start1": (feature.genomic_start0 + 1) if feature.genomic_start0 is not None else "",
            "genomic_end1": feature.genomic_end0 if feature.genomic_end0 is not None else "",
            "clipped": feature.clipped,
            "nested": _is_nested_repeat(feature),
            "attributes": _serialize_attributes(feature.attributes),
        }
        for feature in projected_repeats
        if feature.end0 > feature.start0
    ]
    gene_rows = _unit_gene_composition(units, assigned_genes)
    period_rows = [
        {
            **row,
            "rank": rank,
            "selected": int(row["candidate_period_bp"]) == int(float(period["median"])),
        }
        for rank, row in enumerate(period["candidate_scores"], 1)
    ]

    repeat_inventory_rows = build_repeat_inventory(repeats, units, gap_rows, target_region)
    repeat_retention_qc = summarize_repeat_retention(
        repeats, repeat_inventory_rows, target_region
    )
    disconnected = int(orientation_info["disconnected_units"])
    conflicts = int(orientation_info["orientation_conflicts"])
    orientation_agreement = float(orientation_info["orientation_weight_agreement"])
    unit_length_stats = summarize_refined_unit_lengths(units)
    qc.update(
        {
            "analysis_mode": analysis_mode,
            "target_region": f"{target_region.chrom}:{target_region.start1}-{target_region.end1}",
            "fasta_record": title,
            "stainedglass_hits": len(stainedglass_hits),
            "stainedglass_window_size_bp": period["window_size"],
            "stainedglass_min_identity": period["min_identity"],
            "period_method": period["method"],
            "period_median_bp": period["median"],
            "period_mad_bp": period["mad"],
            "unit_count": len(units),
            "repeat_rich_gap_count": sum(row["interval_type"] == "repeat_rich_gap" for row in gap_rows),
            "unassigned_interval_bp": sum(int(row["length_bp"]) for row in gap_rows),
            "whole_region_tiled": detection_info["whole_region_tiled"],
            "stainedglass_coherent_chains": detection_info["coherent_chains"],
            "stainedglass_repeat_hub_windows": detection_info["repeat_hub_windows"],
            "stainedglass_repeat_hub_runs": len(detection_info["repeat_hub_runs"]),
            "complete_unit_count": sum(unit.status == "complete" for unit in units),
            "refined_unit_length_basis": unit_length_stats["basis"],
            "refined_unit_median_bp": unit_length_stats["median_bp"],
            "boundary_refinement_high": sum(row["confidence"] == "high" for row in boundary_rows),
            "boundary_refinement_medium": sum(row["confidence"] == "medium" for row in boundary_rows),
            "boundary_refinement_low": sum(row["confidence"] == "low" for row in boundary_rows),
            "unit_orientation_status": "pass" if disconnected == 0 and orientation_agreement >= 0.95 else "downgrade",
            "unit_orientation_disconnected": disconnected,
            "unit_orientation_conflicts": conflicts,
            "unit_orientation_edges": orientation_info["orientation_edges"],
            "unit_orientation_weight_agreement": orientation_agreement,
            "sequence_reconstruction_mode": "stainedglass_reference_unit",
            "focal_repeat_family": focal_family_name or "",
            "focal_family_hits": len(focal_family_hits),
            **repeat_retention_qc,
            "repeat_projected_features": len(repeat_rows),
            "repeat_projection_zero_length_removed": zero_length_repeat_projections,
        }
    )
    if qc["repeat_retention_status"] != "pass":
        raise RuntimeError("repeat-retention QC failed")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    prefix = validate_output_prefix(str(args.prefix))
    plot_all_units_enabled = bool(getattr(args, "plot_all_units", False))
    outputs = {
        "input_qc": outdir / f"{prefix}.input_qc.tsv",
        "units": outdir / f"{prefix}.ampliconic_units.tsv",
        "unit_subtypes": outdir / f"{prefix}.unit_subtypes.tsv",
        "inter_unit_gaps": outdir / f"{prefix}.inter_unit_gaps.tsv",
        "gene_composition": outdir / f"{prefix}.unit_gene_composition.tsv",
        "repeat_components": outdir / f"{prefix}.repeat_components.tsv",
        "repeat_inventory": outdir / f"{prefix}.repeat_inventory.tsv",
        "boundary_refinement": outdir / f"{prefix}.boundary_refinement.tsv",
        "stainedglass_periods": outdir / f"{prefix}.stainedglass_periods.tsv",
        "oriented_units": outdir / f"{prefix}.units.oriented.fa",
        "unit_representative": outdir / f"{prefix}.unit_stainedglass_reference.fa",
        "main_pdf": outdir / f"{prefix}.structure.main.pdf",
        "main_png": outdir / f"{prefix}.structure.main.png",
        "subtypes_pdf": outdir / f"{prefix}.structure.unit_subtypes.pdf",
        "subtypes_png": outdir / f"{prefix}.structure.unit_subtypes.png",
        "period_pdf": outdir / f"{prefix}.stainedglass_period_evidence.pdf",
        "period_png": outdir / f"{prefix}.stainedglass_period_evidence.png",
        "summary": outdir / f"{prefix}.summary.md",
    }
    if plot_all_units_enabled:
        outputs.update(
            {
                "all_units_pdf": outdir / f"{prefix}.structure.all_units.pdf",
                "all_units_png": outdir / f"{prefix}.structure.all_units.png",
                "unit_plots_dir": outdir / f"{prefix}.structure.units",
            }
        )
    qc_rows = [
        {
            "metric": key,
            "value": value,
            "status": (
                "downgrade"
                if key == "unit_orientation_status" and value == "downgrade"
                else str(qc["coordinate_status"])
                if key in {"coordinate_status", "coordinate_reason"}
                else "pass"
            ),
        }
        for key, value in qc.items()
    ]
    _write_tsv(outputs["input_qc"], ["metric", "value", "status"], qc_rows)
    _write_tsv(outputs["units"], UNIT_COLUMNS, unit_rows)
    _write_tsv(outputs["unit_subtypes"], SUBTYPE_COLUMNS, subtype_rows)
    _write_tsv(outputs["inter_unit_gaps"], GAP_COLUMNS, gap_rows)
    _write_tsv(outputs["gene_composition"], GENE_COMPOSITION_COLUMNS, gene_rows)
    _write_tsv(outputs["repeat_components"], REPEAT_COMPONENT_COLUMNS, repeat_rows)
    _write_tsv(outputs["repeat_inventory"], REPEAT_INVENTORY_COLUMNS, repeat_inventory_rows)
    _write_tsv(outputs["boundary_refinement"], BOUNDARY_REFINEMENT_COLUMNS, boundary_rows)
    _write_tsv(outputs["stainedglass_periods"], STAINEDGLASS_PERIOD_COLUMNS, period_rows)
    write_fasta(unit_records, outputs["oriented_units"])
    write_fasta({f"stainedglass_reference_{reference_name}": reference_sequence}, outputs["unit_representative"])

    plot_data = {
        "region": target_region,
        "units": units,
        "features": [*repeats, *coding, *pseudogenes],
        "repeat_features": repeats,
        "projected": projected,
        "projected_repeats": projected_repeats,
        "subtypes": subtypes,
        "family_hits": focal_family_hits,
        "focal_family_name": focal_family_name,
        "analysis_mode": analysis_mode,
        "stainedglass_windows": stainedglass_window_support(stainedglass_hits, target_region),
        "coarse_boundaries": detection_info["candidate_boundaries"],
        "gaps": gap_rows,
        "phase_offsets": phase_offsets,
        "gene_colors": build_coding_gene_color_map(coding),
    }
    dpi = int(getattr(args, "dpi", 600))
    plot_region_overview(plot_data, outputs["main_pdf"], outputs["main_png"], dpi=dpi)
    plot_unit_subtypes(plot_data, outputs["subtypes_pdf"], outputs["subtypes_png"], dpi=dpi)
    individual_plot_paths: list[Path] = []
    if plot_all_units_enabled:
        individual_plot_paths = plot_all_units(
            plot_data,
            outputs["all_units_pdf"],
            outputs["all_units_png"],
            outputs["unit_plots_dir"],
            prefix,
            dpi=dpi,
        )
    plot_stainedglass_period_evidence(period, outputs["period_pdf"], outputs["period_png"], dpi=dpi)
    outputs["summary"].write_text(
        _family_free_summary_text(
            period,
            units,
            subtypes,
            gene_states,
            gap_rows,
            coding,
            pseudogenes,
            orientation_info,
            boundary_rows,
            analysis_mode,
            focal_family_name,
            len(focal_family_hits),
        )
    )
    missing_outputs = [
        path for path in [*outputs.values(), *individual_plot_paths] if not path.exists()
    ]
    if missing_outputs:
        raise RuntimeError(f"analysis did not materialize outputs: {missing_outputs}")
    return outputs


def run_analysis(args: argparse.Namespace) -> dict[str, Path]:
    """Run the complete portable analysis and return every materialized output path."""
    if getattr(args, "stainedglass_bed", None) is not None:
        return _run_stainedglass_analysis(args)
    target_region = parse_region_header(str(args.target_region).replace(",", ""))
    title, sequence, fasta_region = read_analysis_fasta(Path(args.fasta), target_region, getattr(args, "region_origin", None))
    repeat_source = read_gff(Path(args.repeat_gff))
    coding_source = read_gff(Path(args.coding_gff3))
    pseudogene_source = read_gff(Path(args.pseudogene_gff3))
    repeats = filter_features_to_region(repeat_source, target_region)
    coding = filter_features_to_region(coding_source, target_region)
    pseudogenes = filter_features_to_region(pseudogene_source, target_region)
    complete_models = _complete_coding_models(coding_source, target_region, fasta_region)
    qc = run_input_qc(sequence, fasta_region, repeats, complete_models)
    pseudogene_qc = summarize_pseudogene_qc(pseudogenes, fasta_region)
    qc.update(pseudogene_qc)
    qc["gene_models"] = int(qc["coding_models"]) + int(qc["pseudogene_models"])
    if qc["coordinate_status"] == "fail":
        raise RuntimeError(f"coding-coordinate QC failed: {qc['coordinate_reason']}")
    if qc["pseudogene_coordinate_status"] == "fail":
        raise RuntimeError("pseudogene-coordinate QC failed")

    hits = merge_family_fragments(select_repeat_hits(repeats, args.family))
    if not hits:
        raise ValueError(f"no valid repeat hits found for exact family name {args.family!r}")
    period = estimate_period(hits)
    if period["median"] is None:
        raise ValueError("at least two same-strand near-complete family hits are required to estimate the unit period")
    units = build_units(hits, period, target_region)
    if not units:
        raise ValueError("no ampliconic units could be constructed")

    phase_offsets = {
        unit.unit_id: unit_phase_offset(unit, hits) if unit.status == "complete" else 0
        for unit in units
    }
    unit_records = {
        unit.unit_id: circular_rephase(
            oriented_unit_sequence(unit, sequence, fasta_region), phase_offsets[unit.unit_id]
        )
        for unit in units
    }
    reference_candidates = _candidate_records(unit_records, units)
    hle_row, body_records = build_hle_evidence(
        args.family,
        hits,
        sequence,
        fasta_region,
        getattr(args, "heliano_bed", None),
        target_region,
    )

    reconstruction_mode = "kmer_medoid"
    mafft_status: str
    unit_representative_sequence: str
    body_representative_sequence: str
    requested_mafft = getattr(args, "mafft", None)
    mafft_executable = _resolve_executable(requested_mafft)
    if requested_mafft is None:
        mafft_status = "unused (not requested; kmer_medoid mode)"
    elif mafft_executable is None:
        mafft_status = "unavailable (requested executable not found; kmer_medoid mode)"
    else:
        try:
            unit_representative_sequence, _unit_alignment = run_mafft_consensus(reference_candidates, mafft_executable)
            body_representative_sequence, _body_alignment = run_mafft_consensus(body_records, mafft_executable)
        except (RuntimeError, ValueError) as exc:
            mafft_status = f"unavailable (alignment failed: {exc}; kmer_medoid mode)"
        else:
            reconstruction_mode = "mafft_consensus"
            mafft_status = "used (mafft_consensus mode)"
    if reconstruction_mode == "kmer_medoid":
        unit_medoid = choose_medoid(reference_candidates, k=min(21, min(map(len, reference_candidates.values()))))
        unit_representative_sequence = reference_candidates[unit_medoid]
        if body_records:
            body_medoid = choose_medoid(body_records, k=min(21, min(map(len, body_records.values()))))
            body_representative_sequence = body_records[body_medoid]
        else:
            body_medoid = "none"
            body_representative_sequence = ""

    projected: list[ProjectedFeature] = []
    projected_repeats: list[ProjectedFeature] = []
    for unit in units:
        for feature in [*repeats, *coding, *pseudogenes]:
            candidate = project_feature_to_unit(feature, unit)
            if candidate is not None:
                projected.append(candidate)
                if feature in repeats:
                    projected_repeats.append(candidate)
    phase_centers, subtypes = infer_body_relative_tspy_phases(units, projected, hits)
    projected = rephase_projected_features(projected, units, phase_offsets)
    projected_repeats = rephase_projected_features(projected_repeats, units, phase_offsets)

    reference_name = choose_medoid(reference_candidates, k=min(21, min(map(len, reference_candidates.values()))))
    reference_sequence = reference_candidates[reference_name]
    unit_rows: list[dict[str, object]] = []
    subtype_rows: list[dict[str, object]] = []
    orientation_supports: list[dict[str, int | float]] = []
    for unit in units:
        intersecting_hits = [hit for hit in hits if hit.start0 < unit.end0 and hit.end0 > unit.start0]
        family_intervals = [(max(hit.start0, unit.start0), min(hit.end0, unit.end0)) for hit in intersecting_hits]
        subtype = subtypes[unit.unit_id]
        oriented_support = oriented_kmer_support(
            reference_sequence,
            unit_records[unit.unit_id],
            k=min(21, len(reference_sequence), len(unit_records[unit.unit_id])),
        )
        support = oriented_support["fraction"]
        if unit.status == "complete":
            orientation_supports.append(oriented_support)
        unit_rows.append(
            {
                "unit_id": unit.unit_id,
                "chrom": unit.chrom,
                "start1": unit.start0 + 1,
                "end1": unit.end0,
                "strand": unit.strand,
                "arm": unit.arm,
                "status": unit.status,
                "length_bp": unit.end0 - unit.start0,
                "family_hit_count": len(intersecting_hits),
                "family_body_bp": _overlap_union_length(family_intervals),
                "period_deviation_bp": (unit.end0 - unit.start0) - float(period["median"]),
                "tspy_phase_A": _phase_state(subtype, "A"),
                "tspy_phase_B": _phase_state(subtype, "B"),
                "subtype": subtype,
                "anchor_support": support,
                "sequence_record": unit.unit_id,
            }
        )
        subtype_rows.append(
            {
                "unit_id": unit.unit_id,
                "tspy_phase_A": _phase_state(subtype, "A"),
                "tspy_phase_B": _phase_state(subtype, "B"),
                "subtype": subtype,
                "phase_A_center_bp": phase_centers.get("A", ""),
                "phase_B_center_bp": phase_centers.get("B", ""),
            }
        )

    if orientation_supports:
        support_values = [float(item["fraction"]) for item in orientation_supports]
        reverse_values = [float(item["reverse_fraction"]) for item in orientation_supports]
        failing_support = sum(value < 0.8 for value in support_values)
        reversed_support = sum(reverse > direct for direct, reverse in zip(support_values, reverse_values))
        orientation_status = "pass" if failing_support == 0 and reversed_support == 0 else "fail"
        qc.update(
            {
                "unit_orientation_status": orientation_status,
                "unit_orientation_min_support": min(support_values),
                "unit_orientation_median_support": statistics.median(support_values),
                "unit_orientation_below_0.8": failing_support,
                "unit_orientation_reverse_better": reversed_support,
            }
        )
        if orientation_status == "fail":
            raise RuntimeError(
                "oriented-unit QC failed: complete units are not homologous in one common direction"
            )
    else:
        qc.update(
            {
                "unit_orientation_status": "downgrade",
                "unit_orientation_min_support": "",
                "unit_orientation_median_support": "",
                "unit_orientation_below_0.8": "",
                "unit_orientation_reverse_better": "",
            }
        )

    period_median = float(period["median"])
    period_mad = float(period["mad"] or 0.0)
    outlier_threshold = max(3.0 * period_mad, 500.0)
    period_outliers = sum(
        unit.status == "complete" and abs((unit.end0 - unit.start0) - period_median) > outlier_threshold
        for unit in units
    )
    qc.update(
        {
            "period_median_bp": period["median"],
            "period_mad_bp": period["mad"],
            "period_min_bp": period["min"],
            "period_max_bp": period["max"],
            "period_plus_median_bp": period["plus_median"],
            "period_minus_median_bp": period["minus_median"],
            "period_intervals": period["n_intervals"],
            "period_outlier_unit_count": period_outliers,
        }
    )
    zero_length_repeat_projections = sum(
        feature.end0 <= feature.start0 for feature in projected_repeats
    )
    repeat_rows: list[dict[str, object]] = []
    for feature in projected_repeats:
        if feature.end0 <= feature.start0:
            continue
        repeat_rows.append(
            {
                "unit_id": feature.unit_id,
                "feature_type": feature.feature_type,
                "repeat_class": _repeat_class(feature) or feature.feature_type or "other_repeat",
                "relative_start0": feature.start0,
                "relative_end0": feature.end0,
                "strand": feature.strand,
                "genomic_start1": (feature.genomic_start0 + 1) if feature.genomic_start0 is not None else "",
                "genomic_end1": feature.genomic_end0 if feature.genomic_end0 is not None else "",
                "clipped": feature.clipped,
                "nested": _is_nested_repeat(feature),
                "attributes": _serialize_attributes(feature.attributes),
            }
        )
    repeat_inventory_rows = build_repeat_inventory(repeats, units, [], target_region)
    repeat_retention_qc = summarize_repeat_retention(
        repeats, repeat_inventory_rows, target_region
    )

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    prefix = validate_output_prefix(str(args.prefix))
    suffix = "consensus" if reconstruction_mode == "mafft_consensus" else "medoid"
    plot_all_units_enabled = bool(getattr(args, "plot_all_units", False))
    outputs = {
        "input_qc": outdir / f"{prefix}.input_qc.tsv",
        "units": outdir / f"{prefix}.ampliconic_units.tsv",
        "unit_subtypes": outdir / f"{prefix}.unit_subtypes.tsv",
        "repeat_components": outdir / f"{prefix}.repeat_components.tsv",
        "repeat_inventory": outdir / f"{prefix}.repeat_inventory.tsv",
        "hle_evidence": outdir / f"{prefix}.helitron_evidence.tsv",
        "oriented_units": outdir / f"{prefix}.units.oriented.fa",
        "unit_representative": outdir / f"{prefix}.unit_{suffix}.fa",
        "hle_representative": outdir / f"{prefix}.helitron_like_{suffix}.fa",
        "main_pdf": outdir / f"{prefix}.structure.main.pdf",
        "main_png": outdir / f"{prefix}.structure.main.png",
        "subtypes_pdf": outdir / f"{prefix}.structure.unit_subtypes.pdf",
        "subtypes_png": outdir / f"{prefix}.structure.unit_subtypes.png",
        "hle_pdf": outdir / f"{prefix}.helitron_boundary_evidence.pdf",
        "hle_png": outdir / f"{prefix}.helitron_boundary_evidence.png",
        "summary": outdir / f"{prefix}.summary.md",
    }
    if plot_all_units_enabled:
        outputs.update(
            {
                "all_units_pdf": outdir / f"{prefix}.structure.all_units.pdf",
                "all_units_png": outdir / f"{prefix}.structure.all_units.png",
                "unit_plots_dir": outdir / f"{prefix}.structure.units",
            }
        )

    dependency_statuses = _dependency_statuses(args, mafft_status)
    qc["target_region"] = f"{target_region.chrom}:{target_region.start1}-{target_region.end1}"
    qc["fasta_record"] = title
    qc["family_hits"] = len(hits)
    qc["unit_count"] = len(units)
    qc["sequence_reconstruction_mode"] = reconstruction_mode
    qc.update(repeat_retention_qc)
    qc["repeat_projected_features"] = len(repeat_rows)
    qc["repeat_projection_zero_length_removed"] = zero_length_repeat_projections
    if qc["repeat_retention_status"] != "pass":
        raise RuntimeError("repeat-retention QC failed")
    qc_rows = [
        {
            "metric": key,
            "value": value,
            "status": (
                str(qc["coordinate_status"])
                if key in {"coordinate_status", "coordinate_reason"}
                else "downgrade"
                if key == "genes_tested" and int(value) == 0
                else "pass"
            ),
        }
        for key, value in qc.items()
    ]
    qc_rows.extend(
        {"metric": f"dependency_{name}", "value": status, "status": "pass" if status.startswith("used") else "downgrade"}
        for name, status in dependency_statuses.items()
    )
    _write_tsv(outputs["input_qc"], ["metric", "value", "status"], qc_rows)
    _write_tsv(outputs["units"], UNIT_COLUMNS, unit_rows)
    _write_tsv(outputs["unit_subtypes"], SUBTYPE_COLUMNS, subtype_rows)
    _write_tsv(outputs["repeat_components"], REPEAT_COMPONENT_COLUMNS, repeat_rows)
    _write_tsv(outputs["repeat_inventory"], REPEAT_INVENTORY_COLUMNS, repeat_inventory_rows)
    _write_tsv(outputs["hle_evidence"], HLE_COLUMNS, [hle_row])
    write_fasta(unit_records, outputs["oriented_units"])
    write_fasta({f"unit_{suffix}": unit_representative_sequence}, outputs["unit_representative"])
    write_fasta({f"helitron_like_{suffix}": body_representative_sequence}, outputs["hle_representative"])

    plot_data = {
        "region": target_region,
        "units": units,
        "features": [*repeats, *coding, *pseudogenes],
        "repeat_features": repeats,
        "projected": projected,
        "projected_repeats": projected_repeats,
        "subtypes": subtypes,
        "family_hits": hits,
        "phase_offsets": phase_offsets,
        "hle_row": hle_row,
        "gene_colors": build_coding_gene_color_map(coding),
    }
    dpi = int(getattr(args, "dpi", 600))
    plot_region_overview(plot_data, outputs["main_pdf"], outputs["main_png"], dpi=dpi)
    plot_unit_subtypes(plot_data, outputs["subtypes_pdf"], outputs["subtypes_png"], dpi=dpi)
    individual_plot_paths: list[Path] = []
    if plot_all_units_enabled:
        individual_plot_paths = plot_all_units(
            plot_data,
            outputs["all_units_pdf"],
            outputs["all_units_png"],
            outputs["unit_plots_dir"],
            prefix,
            dpi=dpi,
        )
    plot_hle_boundary_evidence(plot_data, outputs["hle_pdf"], outputs["hle_png"], dpi=dpi)
    outputs["summary"].write_text(
        _summary_text(args.family, hle_row, period, units, subtypes, dependency_statuses, reconstruction_mode)
    )
    missing_outputs = [
        path for path in [*outputs.values(), *individual_plot_paths] if not path.exists()
    ]
    if missing_outputs:
        raise RuntimeError(f"analysis did not materialize outputs: {missing_outputs}")
    return outputs


def _positive_int(value: str) -> int:
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("value must be a positive integer")
    return parsed


def validate_output_prefix(value: str) -> str:
    parsed = str(value).strip()
    if (
        not parsed
        or parsed in {".", ".."}
        or Path(parsed).name != parsed
        or "/" in parsed
        or "\\" in parsed
        or any(ord(character) < 32 for character in parsed)
    ):
        raise ValueError("output prefix must be a non-empty filename stem without path separators")
    return parsed


def _output_prefix(value: str) -> str:
    try:
        return validate_output_prefix(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Resolve oriented ChrY ampliconic units from StainedGlass, an exact focal repeat "
            "family, or both; every repeat-GFF annotation is retained in all modes."
        )
    )
    parser.add_argument("--fasta", required=True, type=Path, help="unmasked single-record regional FASTA")
    parser.add_argument("--repeat-gff", required=True, type=Path, help="Earl Grey repeat GFF")
    parser.add_argument("--coding-gff3", required=True, type=Path, help="protein-coding gene GFF3")
    parser.add_argument("--pseudogene-gff3", required=True, type=Path, help="reliable pseudogene GFF3")
    parser.add_argument(
        "--family",
        help=(
            "exact repeat-family Name value; with --stainedglass-bed it is an overlay and "
            "does not filter other repeats"
        ),
    )
    parser.add_argument(
        "--stainedglass-bed",
        type=Path,
        help=(
            "StainedGlass eight-column BED or BED.gz; when supplied it controls unit discovery "
            "even if --family is also supplied"
        ),
    )
    parser.add_argument("--period-bp", type=_positive_int, help="optional manual period for StainedGlass mode")
    parser.add_argument("--target-region", required=True, help="target interval as chrom:start-end")
    parser.add_argument("--prefix", required=True, type=_output_prefix, help="output filename prefix")
    parser.add_argument("--outdir", required=True, type=Path, help="output directory")
    parser.add_argument("--region-origin", type=_positive_int, help="explicit 1-based origin for a coordinate-less FASTA")
    parser.add_argument("--mafft", help="optional MAFFT executable name or path")
    parser.add_argument("--minimap2", help="optional minimap2 executable name or path; status is recorded")
    parser.add_argument("--threads", type=_positive_int, default=1, help="worker threads for optional tools")
    parser.add_argument("--heliano-bed", type=Path, help="optional readable HELIANO RC.representative.bed")
    parser.add_argument("--dpi", type=_positive_int, default=600, help="PNG resolution (default: 600)")
    parser.add_argument(
        "--plot-all-units",
        action="store_true",
        help=(
            "write one combined all-unit PDF/PNG plus one PDF/PNG pair for every unit; "
            "disabled by default"
        ),
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if not args.family and args.stainedglass_bed is None:
        parser.error("either --family or --stainedglass-bed is required")
    if args.period_bp is not None and args.stainedglass_bed is None:
        parser.error("--period-bp requires --stainedglass-bed")
    run_analysis(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
