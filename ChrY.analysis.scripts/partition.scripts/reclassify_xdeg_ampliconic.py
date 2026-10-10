#!/usr/bin/env python3
"""Reclassify felid ChrY X-derived loci by ancestry and local structure.

The script deliberately keeps two independent labels:

1. ``evolutionary_class`` describes where a gene family came from.
2. ``structural_context`` describes the present-day sequence environment of
   each particular ChrY copy.

This avoids forcing amplified ancestral X/Y gametologs such as CUL4BY or SRY
into a single, biologically ambiguous "X-degenerate versus ampliconic" label.
"""

from __future__ import annotations

import argparse
import csv
import statistics
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable


TRACK_ORDER = ("FCA_ChrX", "FCA_ChrY", "PPA_ChrY", "PPA_ChrX")
REGION_ORDER = ("PAR", "X-degenerate")
KMER_SIZE = 31
KMER_MASK = (1 << (2 * KMER_SIZE)) - 1
BASE_TO_BITS = {"A": 0, "C": 1, "G": 2, "T": 3}

FOCAL_FAMILIES = {
    "CLDN34Y",
    "FELIDY_ORF1",
    "UTY",
    "DDX3Y",
    "USP9Y",
    "AP1S2Y",
    "EEF1A1Y",
    "UBE1Y",
    "KDM5D",
    "EIF1AY",
    "MBTPS2Y",
    "ZFY",
    "EIF2S3Y",
    "AMELY",
    "TXLNGY",
    "RBMY",
    "TMSB4Y",
    "OFD1Y",
    "RAB9AY",
    "SRY",
    "CUL4BY",
    "HSFY",
}

ANCESTRAL_XY_FAMILIES = {
    "UTY",
    "DDX3Y",
    "USP9Y",
    "AP1S2Y",
    "EEF1A1Y",
    "UBE1Y",
    "KDM5D",
    "EIF1AY",
    "MBTPS2Y",
    "ZFY",
    "EIF2S3Y",
    "AMELY",
    "TXLNGY",
    "RBMY",
    "TMSB4Y",
    "SRY",
    "CUL4BY",
    "HSFY",
}

EXPANDED_X_DERIVED_FAMILIES = {"TXLNGY", "SRY", "CUL4BY", "HSFY"}

# Exact conserved boundary copies.  They occur in the same order in both T2T
# assemblies: SHROOM2Y-[CLDN34Y-FELIDY_ORF1]-UTY.  Their coding sequence and
# distal flank can be duplicated, so local copy number alone is insufficient.
BOUNDARY_ANCHORS = {
    ("FCA", "FCA.ChrY.HC.CLDN34Y.001"),
    ("FCA", "FCA.ChrY.HC.FELIDY_ORF1.001"),
    ("PPA", "CLDN34.c2"),
    ("PPA", "PPA.YG.00046"),
}

LOW_COPY_ISLANDS = {
    ("FCA", "FCA.ChrY.HC.RBMY.001"),
    ("FCA", "FCA.ChrY.HC.TMSB4Y.001"),
    ("FCA", "FCA.ChrY.HC.OFD1Y.001"),
    ("PPA", "TMSB4X.c1"),
    ("PPA", "RAB9A.c1"),
    ("PPA", "RBMX.c1"),
    ("PPA", "KDM5D"),
    ("PPA", "UBA1.c1"),
}

MIXED_BOUNDARY_MODELS = {
    ("PPA", "PPA.ChrY.PR.OFD1Y.001"),
}

CORE_XDEG_FAMILIES = {
    "UTY",
    "DDX3Y",
    "USP9Y",
    "AP1S2Y",
    "EEF1A1Y",
    "UBE1Y",
    "KDM5D",
    "EIF1AY",
    "MBTPS2Y",
    "ZFY",
    "EIF2S3Y",
    "AMELY",
}

PLOT_OUTPUT_COLUMNS = (
    "track",
    "region",
    "order",
    "slot",
    "gene",
    "link_id",
    "homolog_group",
    "strand",
    "source_strand",
    "display_orientation",
    "status",
    "label",
    "label_show",
    "chrom",
    "start0",
    "end0",
    "source_id",
    "annotation_class",
    "evidence",
    "evolutionary_class",
    "structural_context",
    "family_state",
    "plot_role",
    "classification_confidence",
    "classification_note",
)

AUDIT_COLUMNS = (
    "species",
    "family",
    "gene_id",
    "start0",
    "end0",
    "strand",
    "annotation_class",
    "orf_status",
    "original_region_id",
    "original_region_class",
    "coding_or_reviewed_copy_number",
    "body_unique_31mer_fraction",
    "body_repeated_31mer_fraction",
    "body_median_31mer_copy",
    "flank_unique_31mer_fraction",
    "flank_repeated_31mer_fraction",
    "flank_median_31mer_copy",
    "palindrome_arm_ids",
    "palindrome_spacer_ids",
    "evolutionary_class",
    "structural_context",
    "family_state",
    "plot_role",
    "classification_confidence",
    "classification_reason",
)


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]


def canonical_kmer(sequence: str) -> int | None:
    """Return a 2-bit canonical code, or None if the k-mer contains ambiguity."""
    forward = 0
    reverse = 0
    k = len(sequence)
    for base in sequence.upper():
        value = BASE_TO_BITS.get(base)
        if value is None:
            return None
        forward = (forward << 2) | value
        reverse = (reverse >> 2) | ((3 - value) << (2 * (k - 1)))
    return min(forward, reverse)


def read_fasta(path: Path) -> str:
    with path.open() as handle:
        header = handle.readline()
        if not header.startswith(">"):
            raise ValueError(f"Not a FASTA file: {path}")
        return "".join(line.strip().upper() for line in handle if line.strip())


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[dict[str, object]], columns: Iterable[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(columns),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def sample_kmers(
    sequence: str,
    intervals: list[tuple[int, int]],
    max_samples: int = 500,
) -> list[int]:
    valid = [(max(0, start), min(len(sequence), end)) for start, end in intervals]
    total = sum(max(0, end - start - KMER_SIZE + 1) for start, end in valid)
    if total <= 0:
        return []
    step = max(1, total // max_samples)
    sampled: list[int] = []
    offset = 0
    for start, end in valid:
        last = end - KMER_SIZE
        if last < start:
            continue
        first = start + ((step - offset % step) % step)
        for position in range(first, last + 1, step):
            code = canonical_kmer(sequence[position : position + KMER_SIZE])
            if code is not None:
                sampled.append(code)
        offset += last - start + 1
    return sampled


def count_target_kmers(sequence: str, targets: set[int], cap: int = 99) -> dict[int, int]:
    counts = {target: 0 for target in targets}
    forward = 0
    reverse = 0
    valid_bases = 0
    reverse_shift = 2 * (KMER_SIZE - 1)
    for base in sequence:
        value = BASE_TO_BITS.get(base)
        if value is None:
            forward = 0
            reverse = 0
            valid_bases = 0
            continue
        forward = ((forward << 2) | value) & KMER_MASK
        reverse = (reverse >> 2) | ((3 - value) << reverse_shift)
        valid_bases += 1
        if valid_bases < KMER_SIZE:
            continue
        code = min(forward, reverse)
        if code in counts and counts[code] < cap:
            counts[code] += 1
    return counts


def kmer_metrics(kmers: list[int], counts: dict[int, int]) -> tuple[str, str, str]:
    if not kmers:
        return "NA", "NA", "NA"
    values = [counts[code] for code in kmers]
    unique = sum(value == 1 for value in values) / len(values)
    repeated = sum(value >= 2 for value in values) / len(values)
    median = statistics.median(values)
    return f"{unique:.4f}", f"{repeated:.4f}", f"{median:.2f}"


def overlap(start1: int, end1: int, start2: int, end2: int) -> bool:
    return start1 < end2 and end1 > start2


def palindrome_overlaps(
    row: dict[str, object], palindromes: list[dict[str, str]]
) -> tuple[str, str]:
    start = int(row["start0"])
    end = int(row["end0"])
    arms: list[str] = []
    spacers: list[str] = []
    for palindrome in palindromes:
        candidate = palindrome["candidate_id"]
        left_start = int(palindrome["left_start"])
        left_end = int(palindrome["left_end"])
        right_start = int(palindrome["right_start"])
        right_end = int(palindrome["right_end"])
        if overlap(start, end, left_start, left_end) or overlap(
            start, end, right_start, right_end
        ):
            arms.append(candidate)
        if overlap(start, end, left_end, right_start):
            spacers.append(candidate)
    return ",".join(sorted(set(arms))), ",".join(sorted(set(spacers)))


def evolutionary_class(family: str) -> str:
    if family in EXPANDED_X_DERIVED_FAMILIES:
        return "ancestral_XY_gametolog_expanded"
    if family in ANCESTRAL_XY_FAMILIES:
        return "ancestral_XY_gametolog"
    if family == "CLDN34Y":
        return "carnivore_X_derived_expanded"
    if family == "OFD1Y":
        return "X_derived_recurrent_transposition_candidate"
    if family == "RAB9AY":
        return "X_derived_locus_uncertain_route"
    if family == "FELIDY_ORF1":
        return "felid_Y_ORF_unresolved_origin"
    return "unresolved"


def context_call(species: str, gene_id: str, family: str) -> tuple[str, str, str, str]:
    key = (species, gene_id)
    if gene_id == "FCA_RAB9A_Y_like_fragment":
        return (
            "Xdeg_low_copy_island_remnant",
            "Xdeg_remnant",
            "moderate",
            "Reviewed 219-bp RAB9A/Y-like fragment at the distal edge of the FCA OFD1Y low-copy island; not promoted to a gene or reliable pseudogene.",
        )
    if key in BOUNDARY_ANCHORS:
        return (
            "Xdeg_boundary_anchor",
            "Xdeg_anchor",
            "moderate",
            "Conserved SHROOM2Y-[CLDN34Y-FELIDY_ORF1]-UTY order in both assemblies; duplicated sequence at the PAR/Xdeg boundary is recorded separately.",
        )
    if key in MIXED_BOUNDARY_MODELS:
        return (
            "Mixed_Xdeg_ampliconic_boundary",
            "Xdeg_candidate",
            "low",
            "Provisional OFD1Y projection spans the PPA ampliconic-to-low-copy transition; retain as a candidate locus rather than a clean regional anchor.",
        )
    if key in LOW_COPY_ISLANDS:
        return (
            "Xdeg_low_copy_island",
            "Xdeg_anchor",
            "moderate",
            "Single/low-copy locus is embedded between highly redundant ampliconic sequence; the same local-structure rule is used for FCA and PPA.",
        )
    if family in EXPANDED_X_DERIVED_FAMILIES:
        return (
            "Ampliconic",
            "excluded_ampliconic_X_derived",
            "high",
            "X-derived family, but every assembled copy lies in a high-redundancy ampliconic block; ancestry is not used as a structural label.",
        )
    if family in {"CLDN34Y", "FELIDY_ORF1"}:
        return (
            "Ampliconic",
            "excluded_ampliconic_copy",
            "high",
            "Additional family copy outside the conserved boundary anchor lies in an ampliconic block.",
        )
    if family in CORE_XDEG_FAMILIES:
        return (
            "Xdeg_core_backbone",
            "Xdeg_anchor",
            "high",
            "Low-copy X/Y gametolog in the conserved X-degenerate backbone.",
        )
    return (
        "Unresolved",
        "not_plotted",
        "low",
        "No rule matched this locus.",
    )


def load_gene_rows(project: Path, species: str) -> list[dict[str, object]]:
    path = project / f"results/03_landscape/{species}/gene_positions.tsv"
    rows: list[dict[str, object]] = []
    for row in read_tsv(path):
        if row["family"] not in FOCAL_FAMILIES:
            continue
        rows.append(
            {
                "species": species,
                "family": row["family"],
                "gene_id": row["gene_id"],
                "start0": int(row["start0"]),
                "end0": int(row["end0"]),
                "strand": row["strand"],
                "annotation_class": row["annotation_class"],
                "orf_status": row["orf_status"],
                "original_region_id": row["region_id"],
                "original_region_class": row["region_class"],
            }
        )
    if species == "FCA":
        review = read_tsv(project / "results/02_fca_pseudogenes/RAB9AY.review.tsv")
        if len(review) != 1:
            raise ValueError("Expected exactly one reviewed FCA RAB9A/Y-like fragment")
        item = review[0]
        rows.append(
            {
                "species": "FCA",
                "family": "RAB9AY",
                "gene_id": "FCA_RAB9A_Y_like_fragment",
                "start0": int(item["fca_start0"]),
                "end0": int(item["fca_end0"]),
                "strand": item["fca_strand"],
                "annotation_class": item["classification"],
                "orf_status": "short_fragment",
                "original_region_id": "AMP03",
                "original_region_class": "Ampliconic",
            }
        )
    return rows


def create_copy_audit(
    project: Path,
    fasta_paths: dict[str, Path],
) -> list[dict[str, object]]:
    audit: list[dict[str, object]] = []
    for species in ("FCA", "PPA"):
        sequence = read_fasta(fasta_paths[species])
        rows = load_gene_rows(project, species)
        family_counts = Counter(str(row["family"]) for row in rows)
        body_kmers: dict[str, list[int]] = {}
        flank_kmers: dict[str, list[int]] = {}
        targets: set[int] = set()
        for row in rows:
            gene_id = str(row["gene_id"])
            start = int(row["start0"])
            end = int(row["end0"])
            body = sample_kmers(sequence, [(start, end)])
            flank = sample_kmers(
                sequence,
                [(start - 20_000, start), (end, end + 20_000)],
            )
            body_kmers[gene_id] = body
            flank_kmers[gene_id] = flank
            targets.update(body)
            targets.update(flank)
        counts = count_target_kmers(sequence, targets)
        palindromes = read_tsv(
            project / f"results/05_palindromes/{species}/filtered_palindromes.tsv"
        )

        for row in rows:
            gene_id = str(row["gene_id"])
            family = str(row["family"])
            body_unique, body_repeated, body_median = kmer_metrics(
                body_kmers[gene_id], counts
            )
            flank_unique, flank_repeated, flank_median = kmer_metrics(
                flank_kmers[gene_id], counts
            )
            arms, spacers = palindrome_overlaps(row, palindromes)
            context, plot_role, confidence, reason = context_call(
                species, gene_id, family
            )
            audit.append(
                {
                    **row,
                    "coding_or_reviewed_copy_number": family_counts[family],
                    "body_unique_31mer_fraction": body_unique,
                    "body_repeated_31mer_fraction": body_repeated,
                    "body_median_31mer_copy": body_median,
                    "flank_unique_31mer_fraction": flank_unique,
                    "flank_repeated_31mer_fraction": flank_repeated,
                    "flank_median_31mer_copy": flank_median,
                    "palindrome_arm_ids": arms,
                    "palindrome_spacer_ids": spacers,
                    "evolutionary_class": evolutionary_class(family),
                    "structural_context": context,
                    "family_state": "pending",
                    "plot_role": plot_role,
                    "classification_confidence": confidence,
                    "classification_reason": reason,
                }
            )

    states = family_states(audit)
    for row in audit:
        row["family_state"] = states[(str(row["species"]), str(row["family"]))]
    return sorted(
        audit,
        key=lambda row: (
            str(row["species"]),
            int(row["start0"]),
            int(row["end0"]),
            str(row["gene_id"]),
        ),
    )


def family_states(rows: list[dict[str, object]]) -> dict[tuple[str, str], str]:
    grouped: dict[tuple[str, str], list[str]] = defaultdict(list)
    for row in rows:
        grouped[(str(row["species"]), str(row["family"]))].append(
            str(row["structural_context"])
        )
    result: dict[tuple[str, str], str] = {}
    for key, contexts in grouped.items():
        has_amp = "Ampliconic" in contexts
        has_anchor = any(context.startswith("Xdeg_") for context in contexts)
        has_mixed = "Mixed_Xdeg_ampliconic_boundary" in contexts
        if has_amp and has_anchor:
            state = "mixed_Xdeg_anchor_plus_ampliconic_copies"
        elif has_amp:
            state = "ampliconic_only_expansion"
        elif has_mixed and not has_anchor:
            state = "single_mixed_boundary_locus"
        elif has_anchor:
            state = "Xdeg_only_single_or_low_copy"
        else:
            state = "unresolved"
        result[key] = state
    return result


def summarize_families(audit: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in audit:
        grouped[(str(row["species"]), str(row["family"]))].append(row)
    summary: list[dict[str, object]] = []
    for (species, family), rows in sorted(grouped.items()):
        contexts = Counter(str(row["structural_context"]) for row in rows)
        roles = Counter(str(row["plot_role"]) for row in rows)
        summary.append(
            {
                "species": species,
                "family": family,
                "evolutionary_class": rows[0]["evolutionary_class"],
                "coding_or_reviewed_copy_count": len(rows),
                "Xdeg_anchor_copies": sum(
                    count for context, count in contexts.items() if context.startswith("Xdeg_")
                ),
                "mixed_boundary_copies": contexts["Mixed_Xdeg_ampliconic_boundary"],
                "ampliconic_copies": contexts["Ampliconic"],
                "plot_anchor_or_candidate_copies": roles["Xdeg_anchor"]
                + roles["Xdeg_candidate"]
                + roles["Xdeg_remnant"],
                "family_state": rows[0]["family_state"],
            }
        )
    return summary


def revised_partitions() -> dict[str, list[dict[str, object]]]:
    columns = ("chrom", "start0", "end0", "region_id", "class", "confidence", "subtype", "evidence")
    del columns  # documents the row layout; write order is specified by caller
    fca = [
        (0, 1_127_314, "RPT01", "Repeat-rich", "high", "terminal_repeat", "original_partition"),
        (1_127_314, 8_222_961, "PAR01", "PAR", "high", "PAR", "original_partition"),
        (8_222_961, 8_861_090, "XDG01", "X-degenerate", "high", "core_backbone", "original_partition_plus_conserved_boundary_anchor"),
        (8_861_090, 9_018_340, "AMP01", "Ampliconic", "high", "ampliconic", "original_partition"),
        (9_018_340, 9_182_865, "XDG02", "X-degenerate", "high", "core_backbone", "original_partition"),
        (9_182_865, 9_608_410, "AMP02", "Ampliconic", "high", "ampliconic", "original_partition"),
        (9_608_410, 10_006_000, "XDG03", "X-degenerate", "moderate", "core_backbone", "31mer_low_copy_transition;TXLNGY_excluded"),
        (10_006_000, 11_262_000, "AMP03A", "Ampliconic", "high", "ampliconic", "31mer_high_redundancy"),
        (11_262_000, 11_301_000, "XDG04", "X-degenerate", "moderate", "low_copy_island", "31mer_low_copy_island;contains_RBMY"),
        (11_301_000, 11_921_000, "AMP03B", "Ampliconic", "high", "ampliconic", "31mer_high_redundancy"),
        (11_921_000, 12_021_000, "XDG05", "X-degenerate", "moderate", "low_copy_island", "31mer_low_copy_island;contains_TMSB4Y_OFD1Y_RAB9A_fragment"),
        (12_021_000, 19_293_000, "AMP03C", "Ampliconic", "high", "ampliconic", "31mer_high_redundancy_plus_original_partition"),
        (19_293_000, 45_999_067, "CPLX01", "Repeat-rich_ampliconic_complex", "moderate", "complex", "original_partition"),
    ]
    ppa = [
        (0, 576_762, "RPT01", "Repeat-rich_satellite_cap", "high", "terminal_repeat", "original_partition"),
        (576_762, 7_572_927, "PAR01", "PAR", "high", "PAR", "original_partition"),
        (7_572_927, 7_593_358, "XDG00", "X-degenerate", "moderate", "boundary_anchor_mixed", "conserved_SHROOM2Y_CLDN34Y_FELIDY_ORF1_UTY_synteny;duplicated_sequence"),
        (7_593_358, 8_265_358, "XDG01", "X-degenerate", "high", "core_backbone", "original_partition"),
        (8_265_358, 8_352_507, "AMP02", "Ampliconic", "high", "ampliconic", "original_partition"),
        (8_352_507, 8_883_192, "XDG02", "X-degenerate", "high", "core_backbone", "original_partition"),
        (8_883_192, 9_558_848, "AMP03", "Ampliconic", "high", "ampliconic", "original_partition"),
        (9_558_848, 9_703_272, "XDG03", "X-degenerate", "high", "low_copy_island", "original_partition;contains_OFD1Y_tail_TMSB4Y_RAB9AY"),
        (9_703_272, 9_835_550, "AMP04", "Ampliconic", "high", "ampliconic", "original_partition"),
        (9_835_550, 9_897_584, "XDG04", "X-degenerate", "high", "low_copy_island", "original_partition;contains_RBMY"),
        (9_897_584, 10_008_683, "AMP05", "Ampliconic", "high", "ampliconic", "original_partition"),
        (10_008_683, 10_153_033, "XDG05", "X-degenerate", "high", "low_copy_island", "original_partition;contains_KDM5D_UBE1Y"),
        (10_153_033, 11_304_358, "AMP06", "Ampliconic", "high", "ampliconic", "original_partition"),
        (11_304_358, 11_399_358, "OMS01", "Other_MSY", "moderate", "other_MSY", "original_partition"),
        (11_399_358, 12_914_358, "AMP07", "Ampliconic", "high", "ampliconic", "original_partition"),
        (12_914_358, 44_117_358, "CPLX01", "Repeat-rich_ampliconic_complex", "moderate", "complex", "original_partition"),
    ]

    def materialize(items: list[tuple[object, ...]]) -> list[dict[str, object]]:
        return [
            {
                "chrom": "ChrY",
                "start0": start,
                "end0": end,
                "region_id": region_id,
                "class": class_name,
                "confidence": confidence,
                "subtype": subtype,
                "evidence": evidence,
            }
            for start, end, region_id, class_name, confidence, subtype, evidence in items
        ]

    return {"FCA": materialize(fca), "PPA": materialize(ppa)}


def plot_context_for_row(
    row: dict[str, str], audit_lookup: dict[tuple[str, str], dict[str, object]]
) -> tuple[str, str, str, str, str, str]:
    track = row["track"]
    region = row["region"]
    if region == "PAR":
        return (
            "pseudoautosomal",
            "PAR",
            "recombining_PAR",
            "PAR",
            "high",
            "Pseudoautosomal locus.",
        )
    if track.endswith("ChrX"):
        return (
            "X_reference",
            "X_reference",
            "X_reference",
            "X_reference",
            "high",
            "ChrX reference homolog used for order comparison.",
        )
    species = track.split("_", 1)[0]
    key = (species, row["source_id"])
    if key not in audit_lookup:
        raise ValueError(f"No copy-context audit row for plotting locus {key}")
    audit = audit_lookup[key]
    return (
        str(audit["evolutionary_class"]),
        str(audit["structural_context"]),
        str(audit["family_state"]),
        str(audit["plot_role"]),
        str(audit["classification_confidence"]),
        str(audit["classification_reason"]),
    )


def reclassify_plot_input(
    base_input: Path,
    project: Path,
    audit: list[dict[str, object]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = [dict(row) for row in read_tsv(base_input)]
    required = {"FCA_ChrX", "FCA_ChrY", "PPA_ChrY", "PPA_ChrX"}
    if {str(row["track"]) for row in rows} != required:
        raise ValueError("The revised plotting input must contain all four tracks")

    # TXLNGY is retained in both assemblies, but all current copies occupy
    # high-redundancy ampliconic sequence.  Excluding TXLNG from all four
    # strict-Xdeg tracks prevents an apparent PPA gene loss.
    rows = [
        row
        for row in rows
        if not (
            row["region"] == "X-degenerate" and row["homolog_group"] == "TXLNG"
        )
    ]

    # Correct the user's PPA CLDN34Y link typo so the same exact key is used
    # by both X chromosomes and both Y chromosomes.
    for row in rows:
        if row["region"] == "X-degenerate" and row["homolog_group"] == "CLDN34":
            row["link_id"] = "XDG_CLDN34_1"

    fca_gene_rows = {
        row["family"]: row
        for row in read_tsv(project / "results/03_landscape/FCA/gene_positions.tsv")
        if row["family"] in {"RBMY", "TMSB4Y", "OFD1Y"}
    }
    aliases = {
        "RBMY": ("RBMX", "XDG_RBMX_1"),
        "TMSB4Y": ("TMSB4X", "XDG_TMSB4X_1"),
        "OFD1Y": ("OFD1", "XDG_OFD1_1"),
    }
    for family in ("RBMY", "TMSB4Y", "OFD1Y"):
        source = fca_gene_rows[family]
        group, link_id = aliases[family]
        rows.append(
            {
                "track": "FCA_ChrY",
                "region": "X-degenerate",
                "order": 0,
                "slot": 0,
                "gene": family,
                "link_id": link_id,
                "homolog_group": group,
                "strand": source["strand"],
                "source_strand": source["strand"],
                "display_orientation": "forward",
                "status": "intact",
                "label": family,
                "label_show": 1,
                "chrom": "ChrY",
                "start0": source["start0"],
                "end0": source["end0"],
                "source_id": source["gene_id"],
                "annotation_class": source["annotation_class"],
                "evidence": "reclassified_31mer_low_copy_island",
            }
        )

    audit_lookup = {
        (str(row["species"]), str(row["gene_id"])): row for row in audit
    }
    for row in rows:
        (
            row["evolutionary_class"],
            row["structural_context"],
            row["family_state"],
            row["plot_role"],
            row["classification_confidence"],
            row["classification_note"],
        ) = plot_context_for_row(row, audit_lookup)

    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[(str(row["track"]), str(row["region"]))].append(row)
    finished: list[dict[str, object]] = []
    for track in TRACK_ORDER:
        for region in REGION_ORDER:
            group = grouped[(track, region)]
            orientations = {str(row["display_orientation"]) for row in group}
            if len(orientations) != 1:
                raise ValueError(f"Multiple display orientations for {(track, region)}")
            reverse = orientations == {"reverse"}
            group.sort(
                key=lambda row: (
                    int(row["start0"]),
                    int(row["end0"]),
                    str(row["source_id"]),
                ),
                reverse=reverse,
            )
            copy_counts = Counter(str(row["homolog_group"]) for row in group)
            copy_seen: Counter[str] = Counter()
            for order, row in enumerate(group, start=1):
                row["order"] = order
                row["slot"] = order
                homolog = str(row["homolog_group"])
                copy_seen[homolog] += 1
                if copy_counts[homolog] > 1 and not str(row["label"]).endswith(
                    f"-{copy_seen[homolog]}"
                ):
                    row["label"] = f"{row['gene']}-{copy_seen[homolog]}"
                finished.append(row)

    validate_plot_rows(finished)
    return finished


def validate_plot_rows(rows: list[dict[str, object]]) -> None:
    expected = {
        ("FCA_ChrX", "PAR"): 34,
        ("FCA_ChrY", "PAR"): 36,
        ("PPA_ChrY", "PAR"): 36,
        ("PPA_ChrX", "PAR"): 36,
        ("FCA_ChrX", "X-degenerate"): 14,
        ("FCA_ChrY", "X-degenerate"): 18,
        ("PPA_ChrY", "X-degenerate"): 18,
        ("PPA_ChrX", "X-degenerate"): 16,
    }
    observed = Counter((str(row["track"]), str(row["region"])) for row in rows)
    if observed != expected:
        raise ValueError(f"Unexpected reclassified plot counts: {observed}; expected {expected}")
    keys = Counter(
        (str(row["track"]), str(row["region"]), str(row["link_id"])) for row in rows
    )
    duplicates = [key for key, count in keys.items() if count > 1]
    if duplicates:
        raise ValueError(f"Duplicate exact link keys: {duplicates}")
    link_groups: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        link_groups[str(row["link_id"])].add(str(row["homolog_group"]))
    inconsistent = {key: value for key, value in link_groups.items() if len(value) > 1}
    if inconsistent:
        raise ValueError(f"Link keys map to multiple homolog groups: {inconsistent}")
    y_groups = {
        track: {
            str(row["homolog_group"])
            for row in rows
            if row["track"] == track and row["region"] == "X-degenerate"
        }
        for track in ("FCA_ChrY", "PPA_ChrY")
    }
    if y_groups["FCA_ChrY"] != y_groups["PPA_ChrY"]:
        raise ValueError(f"Asymmetric strict-Xdeg Y repertoires: {y_groups}")


def write_rules(path: Path) -> None:
    rows = [
        {
            "priority": 1,
            "criterion": "evolutionary_origin",
            "decision": "Record X/Y gametolog or X-derived ancestry independently of present structural context.",
            "purpose": "Allows an ancestral X-derived family to be structurally ampliconic after amplification.",
        },
        {
            "priority": 2,
            "criterion": "copy_level_local_structure",
            "decision": "Use local 31-mer redundancy, palindrome/duplicated-block overlap, and flanking context for each copy.",
            "purpose": "Prevents family-wide labels based only on total copy number or gene midpoint.",
        },
        {
            "priority": 3,
            "criterion": "conserved_boundary_synteny",
            "decision": "Retain conserved SHROOM2Y-[CLDN34Y-FELIDY_ORF1]-UTY copies as boundary Xdeg anchors even when the gene sequence is duplicated.",
            "purpose": "Distinguishes an ancestral/boundary anchor copy from additional ampliconic copies.",
        },
        {
            "priority": 4,
            "criterion": "low_copy_island",
            "decision": "Call a low-copy Xdeg island when a unique locus is bounded on both sides by high-redundancy ampliconic sequence.",
            "purpose": "Applies one structural rule to FCA and PPA RBMY/TMSB4Y/OFD1Y-like loci.",
        },
        {
            "priority": 5,
            "criterion": "strict_gene_order_plot",
            "decision": "Plot one low-copy/boundary anchor per family; exclude ampliconic-only X-derived families and list them explicitly as retained elsewhere.",
            "purpose": "Avoids false gene-loss inference and arbitrary choice among indistinguishable ampliconic copies.",
        },
    ]
    write_tsv(path, rows, ("priority", "criterion", "decision", "purpose"))


def main() -> None:
    project = Path(__file__).resolve().parents[1]
    default_base = project / "par_xdeg_4tracks/base_input.tsv"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--partitions-only", action="store_true",
                        help="Rebuild reviewed BED tables without genome input; boundaries are frozen decisions")
    parser.add_argument("--base-input", type=Path, default=default_base)
    parser.add_argument(
        "--fca-fasta",
        type=Path,
        default=project / "inputs/FCA.ChrY.fa",
    )
    parser.add_argument(
        "--ppa-fasta",
        type=Path,
        default=project / "inputs/PPA.ChrY.fa",
    )
    parser.add_argument(
        "--audit-output",
        type=Path,
        default=project / "results/04_comparative/Xderived_copy_context_audit.tsv",
    )
    parser.add_argument(
        "--summary-output",
        type=Path,
        default=project / "results/04_comparative/Xderived_family_context_summary.tsv",
    )
    parser.add_argument(
        "--plot-output",
        type=Path,
        default=project / "par_xdeg_4tracks/felid_PAR_Xdeg_input_reclassified.tsv",
    )
    args = parser.parse_args()

    partitions = revised_partitions()
    partition_columns = (
        "chrom",
        "start0",
        "end0",
        "region_id",
        "class",
        "confidence",
        "subtype",
        "evidence",
    )
    for species, rows in partitions.items():
        write_tsv(
            project
            / f"results/04_comparative/{species}.ChrY.structural_partition.revised.bed",
            rows,
            partition_columns,
        )

    write_rules(project / "results/04_comparative/Xdeg_ampliconic_classification_rules.tsv")
    if args.partitions_only:
        print("Wrote FCA/PPA reviewed partitions; no sequence-derived metrics were recomputed")
        return

    audit = create_copy_audit(
        project,
        {"FCA": args.fca_fasta, "PPA": args.ppa_fasta},
    )
    write_tsv(args.audit_output, audit, AUDIT_COLUMNS)
    summary = summarize_families(audit)
    write_tsv(
        args.summary_output,
        summary,
        (
            "species",
            "family",
            "evolutionary_class",
            "coding_or_reviewed_copy_count",
            "Xdeg_anchor_copies",
            "mixed_boundary_copies",
            "ampliconic_copies",
            "plot_anchor_or_candidate_copies",
            "family_state",
        ),
    )

    plot_rows = reclassify_plot_input(args.base_input, project, audit)
    write_tsv(args.plot_output, plot_rows, PLOT_OUTPUT_COLUMNS)
    print(f"Wrote {len(audit)} copy-audit rows to {args.audit_output}")
    print(f"Wrote {len(summary)} family summaries to {args.summary_output}")
    print(f"Wrote {len(plot_rows)} plotting rows to {args.plot_output}")


if __name__ == "__main__":
    main()
