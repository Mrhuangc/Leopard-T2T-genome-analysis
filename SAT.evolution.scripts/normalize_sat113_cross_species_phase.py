#!/usr/bin/env python3
"""Normalize circular SAT113 monomers to one cross-species reference phase.

The upstream satellite-array workflow phases copies to a species-specific global
consensus.  Different species can therefore use different linear starts for the
same circular monomer.  A fixed-coordinate MAFFT alignment hides that nuisance
for site-wise analyses, but raw exact-haplotype comparisons and VSEARCH still
see the artificial end break.

This script removes that discrepancy before exact haplotypes are defined for
downstream analyses.  It:

1. builds a copy-capped consensus for each species;
2. chooses the strand and circular rotation that best match one fixed reference;
3. refines each haplotype's rotation near that species-level anchor;
4. collapses identical reference-phased sequences and rebuilds every haplotype
   identifier and copy mapping.

Genomic coordinates and the original ``strand`` column are never modified.
Every transformation is written to explicit QC/mapping tables.
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, Optional, Tuple

import pandas as pd


CANONICAL = set("ACGT")
IUPAC = set("ACGTRYSWKMBDHVN")
REVCOMP = str.maketrans(
    "ACGTRYMKSWBDHVNacgtrymkswbdhvn",
    "TGCAYRKMSWVHDBNtgcayrkmswvhdbn",
)


def eprint(*args: object) -> None:
    print(*args, file=sys.stderr)


def normalize_fasta_id(name: str) -> str:
    return name.split()[0].split(";")[0]


def clean_sequence(sequence: str) -> str:
    return re.sub(r"[.\s-]", "", str(sequence)).upper()


def read_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    name: Optional[str] = None
    chunks: List[str] = []
    with open(path, encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield normalize_fasta_id(name), clean_sequence("".join(chunks))
                name = line[1:].strip()
                chunks = []
            else:
                if name is None:
                    raise ValueError(f"Sequence before FASTA header in {path}")
                chunks.append(line)
    if name is not None:
        yield normalize_fasta_id(name), clean_sequence("".join(chunks))


def read_fasta_dict(path: Path) -> Dict[str, str]:
    records: Dict[str, str] = {}
    for name, sequence in read_fasta(path):
        if name in records:
            raise ValueError(f"Duplicate FASTA identifier after normalization: {name}")
        if not sequence or any(base not in IUPAC for base in sequence):
            raise ValueError(f"Invalid nucleotide sequence for {name} in {path}")
        records[name] = sequence
    return records


def write_fasta(records: Iterable[Tuple[str, str]], path: Path, width: int = 80) -> None:
    with open(path, "w", encoding="utf-8") as out:
        for name, sequence in records:
            out.write(f">{name}\n")
            for start in range(0, len(sequence), width):
                out.write(sequence[start : start + width] + "\n")


def revcomp(sequence: str) -> str:
    return sequence.translate(REVCOMP)[::-1].upper()


def rotate_left(sequence: str, amount: int) -> str:
    if not sequence:
        return sequence
    amount %= len(sequence)
    return sequence[amount:] + sequence[:amount]


def signed_rotation(amount: int, period: int) -> int:
    """Return the shortest signed representation of a left rotation."""
    if period <= 0:
        return 0
    amount %= period
    return amount if amount <= period / 2 else amount - period


def circular_distance(first: int, second: int, period: int) -> int:
    if period <= 0:
        return 0
    delta = abs((first - second) % period)
    return int(min(delta, period - delta))


def substitution_cost(first: str, second: str) -> int:
    return 0 if first == second and first in CANONICAL else 1


def global_edit_distance(first: str, second: str) -> int:
    """Unit-cost global edit distance using O(len(second)) memory."""
    if len(first) < len(second):
        first, second = second, first
    previous = list(range(len(second) + 1))
    for i, base_first in enumerate(first, 1):
        current = [i] + [0] * len(second)
        for j, base_second in enumerate(second, 1):
            current[j] = min(
                previous[j] + 1,
                current[j - 1] + 1,
                previous[j - 1] + substitution_cost(base_first, base_second),
            )
        previous = current
    return previous[-1]


def banded_edit_distance(first: str, second: str, band: int) -> int:
    """Global edit distance in a diagonal band; falls back if the band is missed."""
    n_first, n_second = len(first), len(second)
    band = max(int(band), abs(n_first - n_second))
    infinity = n_first + n_second + 1
    previous = [infinity] * (n_second + 1)
    for j in range(0, min(n_second, band) + 1):
        previous[j] = j
    for i in range(1, n_first + 1):
        current = [infinity] * (n_second + 1)
        low = max(0, i - band)
        high = min(n_second, i + band)
        if low == 0:
            current[0] = i
        for j in range(max(1, low), high + 1):
            current[j] = min(
                previous[j] + 1,
                current[j - 1] + 1,
                previous[j - 1] + substitution_cost(first[i - 1], second[j - 1]),
            )
        previous = current
    result = previous[n_second]
    if result >= infinity:
        return global_edit_distance(first, second)
    return int(result)


def capped_species_haplotype_weights(
    metadata: pd.DataFrame,
    group_by: str,
    cap_per_unit: int,
) -> Dict[Tuple[str, str], float]:
    required = {"species", "haplotype_id", group_by}
    missing = sorted(required - set(metadata.columns))
    if missing:
        raise ValueError(f"Metadata lacks required columns: {', '.join(missing)}")
    frame = metadata.copy()
    frame["unit_uid"] = frame["species"].astype(str) + "|" + frame[group_by].astype(str)
    weights: Counter[Tuple[str, str]] = Counter()
    for _unit, subset in frame.groupby("unit_uid", sort=False):
        scale = 1.0
        if cap_per_unit > 0 and len(subset) > cap_per_unit:
            scale = cap_per_unit / len(subset)
        species = str(subset["species"].iloc[0])
        for haplotype, count in subset["haplotype_id"].astype(str).value_counts().items():
            weights[(species, str(haplotype))] += float(count) * scale
    return dict(weights)


def weighted_consensus(
    sequences: Mapping[str, str],
    weights: Mapping[str, float],
    length: int,
) -> Tuple[str, float]:
    counts = [{base: 0.0 for base in "ACGT"} for _ in range(length)]
    total_weight = 0.0
    for haplotype, weight in weights.items():
        sequence = sequences[haplotype]
        if len(sequence) != length or weight <= 0:
            continue
        total_weight += float(weight)
        for position, base in enumerate(sequence):
            if base in CANONICAL:
                counts[position][base] += float(weight)
    consensus: List[str] = []
    for position_counts in counts:
        maximum = max(position_counts.values())
        if maximum <= 0:
            consensus.append("N")
        else:
            consensus.append(
                sorted(base for base, count in position_counts.items() if count == maximum)[0]
            )
    return "".join(consensus), total_weight


def choose_species_anchor(
    species: str,
    sequences: Mapping[str, str],
    weights: Mapping[str, float],
    reference: str,
    orientation_mode: str,
) -> Tuple[Dict[str, object], str]:
    length_weights: Counter[int] = Counter()
    for haplotype, weight in weights.items():
        length_weights[len(sequences[haplotype])] += float(weight)
    if not length_weights:
        raise ValueError(f"{species}: no positively weighted haplotypes")
    reference_length = len(reference)
    if length_weights.get(reference_length, 0.0) > 0:
        anchor_length = reference_length
    else:
        anchor_length = min(
            length_weights,
            key=lambda length: (
                abs(length - reference_length),
                -length_weights[length],
                length,
            ),
        )
    consensus, anchor_weight = weighted_consensus(
        sequences, weights, anchor_length
    )
    total_weight = float(sum(length_weights.values()))
    orientations = [("forward", consensus)]
    if orientation_mode == "auto":
        orientations.append(("reverse_complement", revcomp(consensus)))

    scored: List[Tuple[int, int, int, str, str]] = []
    for orientation, oriented_consensus in orientations:
        orientation_rank = 0 if orientation == "forward" else 1
        for rotation in range(len(oriented_consensus)):
            phased = rotate_left(oriented_consensus, rotation)
            distance = global_edit_distance(phased, reference)
            scored.append(
                (distance, orientation_rank, rotation, orientation, phased)
            )
    scored.sort(key=lambda item: (item[0], item[1], item[2]))
    best = scored[0]
    second = scored[1] if len(scored) > 1 else best
    best_distance, _rank, rotation, orientation, phased_consensus = best
    denominator = max(len(consensus), len(reference), 1)
    similarity = 1.0 - best_distance / denominator
    second_similarity = 1.0 - second[0] / denominator
    row: Dict[str, object] = {
        "species": species,
        "reference_length": len(reference),
        "anchor_length": anchor_length,
        "anchor_capped_weight": anchor_weight,
        "total_capped_weight": total_weight,
        "anchor_weight_fraction": anchor_weight / max(total_weight, 1.0),
        "anchor_consensus_raw": consensus,
        "orientation_to_reference": orientation,
        "rotation_left_bp": rotation,
        "rotation_signed_bp": signed_rotation(rotation, anchor_length),
        "best_edit_distance": best_distance,
        "second_best_edit_distance": second[0],
        "edit_distance_margin": second[0] - best_distance,
        "reference_similarity": similarity,
        "second_best_similarity": second_similarity,
        "anchor_consensus_reference_phase": phased_consensus,
    }
    return row, phased_consensus


def circular_kmer_votes(sequence: str, reference: str, kmer_size: int) -> Counter[int]:
    votes: Counter[int] = Counter()
    if not sequence or not reference:
        return votes
    kmer_size = max(2, min(int(kmer_size), len(sequence), len(reference)))
    reference_circular = reference + reference[: kmer_size - 1]
    query_circular = sequence + sequence[: kmer_size - 1]
    reference_index: Dict[str, List[int]] = defaultdict(list)
    for position in range(len(reference)):
        kmer = reference_circular[position : position + kmer_size]
        if set(kmer) <= CANONICAL:
            reference_index[kmer].append(position)
    for query_position in range(len(sequence)):
        kmer = query_circular[query_position : query_position + kmer_size]
        for reference_position in reference_index.get(kmer, []):
            votes[(query_position - reference_position) % len(sequence)] += 1
    return votes


def ungapped_match_count(sequence: str, reference: str, rotation: int) -> int:
    phased = rotate_left(sequence, rotation)
    return sum(
        first == second and first in CANONICAL
        for first, second in zip(phased, reference)
    )


def choose_haplotype_rotation(
    sequence: str,
    reference: str,
    expected_rotation: int,
    search_radius: int,
    kmer_size: int,
    candidate_count: int,
    alignment_band: int,
) -> Dict[str, object]:
    length = len(sequence)
    expected_rotation %= length
    radius = min(max(0, int(search_radius)), max(0, length // 2))
    allowed = sorted({(expected_rotation + delta) % length for delta in range(-radius, radius + 1)})

    votes = circular_kmer_votes(sequence, reference, kmer_size)
    ranked = sorted(
        allowed,
        key=lambda rotation: (
            -votes.get(rotation, 0),
            -ungapped_match_count(sequence, reference, rotation),
            circular_distance(rotation, expected_rotation, length),
            rotation,
        ),
    )
    selected_candidates = ranked[: max(1, int(candidate_count))]
    if expected_rotation not in selected_candidates:
        selected_candidates.append(expected_rotation)

    scores: List[Tuple[int, int, int, int, int]] = []
    for rotation in sorted(set(selected_candidates)):
        phased = rotate_left(sequence, rotation)
        distance = banded_edit_distance(phased, reference, alignment_band)
        matches = ungapped_match_count(sequence, reference, rotation)
        scores.append(
            (
                distance,
                -votes.get(rotation, 0),
                -matches,
                circular_distance(rotation, expected_rotation, length),
                rotation,
            )
        )
    scores.sort()
    best = scores[0]
    second = scores[1] if len(scores) > 1 else best
    best_rotation = best[4]
    best_votes = -best[1]
    best_matches = -best[2]
    tied = sum(score[:3] == best[:3] for score in scores)
    return {
        "rotation_left_bp": best_rotation,
        "rotation_signed_bp": signed_rotation(best_rotation, length),
        "expected_rotation_left_bp": expected_rotation,
        "expected_rotation_signed_bp": signed_rotation(expected_rotation, length),
        "residual_rotation_bp": signed_rotation(best_rotation - expected_rotation, length),
        "best_edit_distance": best[0],
        "second_best_edit_distance": second[0],
        "edit_distance_margin": second[0] - best[0],
        "kmer_votes": best_votes,
        "ungapped_matches": best_matches,
        "candidate_rotations_tested": len(scores),
        "phase_tied": tied > 1,
    }


def validate_inputs(
    metadata: pd.DataFrame,
    copy_haplotypes: pd.DataFrame,
    haplotype_summary: pd.DataFrame,
    sequences: Mapping[str, str],
    group_by: str,
) -> None:
    required_metadata = {"seq_id", "species", "haplotype_id", group_by}
    missing = sorted(required_metadata - set(metadata.columns))
    if missing:
        raise ValueError(f"Metadata lacks required columns: {', '.join(missing)}")
    required_mapping = {"seq_id", "haplotype_id"}
    missing = sorted(required_mapping - set(copy_haplotypes.columns))
    if missing:
        raise ValueError(f"Copy-haplotype table lacks columns: {', '.join(missing)}")
    if metadata["seq_id"].duplicated().any():
        raise ValueError("Metadata seq_id values must be unique")
    if copy_haplotypes["seq_id"].duplicated().any():
        raise ValueError("Copy-haplotype seq_id values must be unique")
    check = metadata[["seq_id", "haplotype_id"]].merge(
        copy_haplotypes[["seq_id", "haplotype_id"]],
        on="seq_id",
        how="outer",
        suffixes=("_metadata", "_mapping"),
        indicator=True,
    )
    if not check["_merge"].eq("both").all() or not (
        check["haplotype_id_metadata"].astype(str)
        == check["haplotype_id_mapping"].astype(str)
    ).all():
        raise ValueError("Metadata and copy-haplotype assignments disagree")
    summary_ids = set(haplotype_summary["haplotype_id"].astype(str))
    fasta_ids = set(sequences)
    metadata_ids = set(metadata["haplotype_id"].astype(str))
    if summary_ids != fasta_ids or metadata_ids != fasta_ids:
        raise ValueError(
            "Haplotype identifiers disagree among metadata, summary, and FASTA: "
            f"summary={len(summary_ids)}, fasta={len(fasta_ids)}, metadata={len(metadata_ids)}"
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--copy-haplotypes", required=True, type=Path)
    parser.add_argument("--haplotype-summary", required=True, type=Path)
    parser.add_argument("--haplotype-fasta", required=True, type=Path)
    parser.add_argument("--reference-fasta", required=True, type=Path)
    parser.add_argument("--output-metadata", required=True, type=Path)
    parser.add_argument("--output-copy-haplotypes", required=True, type=Path)
    parser.add_argument("--output-haplotype-summary", required=True, type=Path)
    parser.add_argument("--output-haplotype-fasta", required=True, type=Path)
    parser.add_argument("--output-all-copies-fasta", required=True, type=Path)
    parser.add_argument("--output-phase-map", required=True, type=Path)
    parser.add_argument("--output-species-summary", required=True, type=Path)
    parser.add_argument("--output-summary", required=True, type=Path)
    parser.add_argument("--group-by", choices=["array_key", "region_key"], default="region_key")
    parser.add_argument("--cap-copies-per-unit", type=int, default=100)
    parser.add_argument("--orientation-mode", choices=["auto", "forward"], default="auto")
    parser.add_argument("--phase-search-radius", type=int, default=12)
    parser.add_argument("--kmer-size", type=int, default=7)
    parser.add_argument("--candidate-count", type=int, default=6)
    parser.add_argument("--alignment-band", type=int, default=16)
    parser.add_argument("--min-anchor-capped-weight", type=float, default=10.0)
    parser.add_argument("--min-anchor-similarity", type=float, default=0.70)
    parser.add_argument("--max-ambiguous-weight-fraction", type=float, default=0.25)
    parser.add_argument("--allow-poor-phase", action="store_true")
    args = parser.parse_args()

    for value, label in (
        (args.phase_search_radius, "--phase-search-radius"),
        (args.candidate_count, "--candidate-count"),
        (args.alignment_band, "--alignment-band"),
    ):
        if value < 1:
            raise ValueError(f"{label} must be >= 1")
    if not 0 <= args.max_ambiguous_weight_fraction <= 1:
        raise ValueError("--max-ambiguous-weight-fraction must be in [0, 1]")

    output_paths = (
        args.output_metadata,
        args.output_copy_haplotypes,
        args.output_haplotype_summary,
        args.output_haplotype_fasta,
        args.output_all_copies_fasta,
        args.output_phase_map,
        args.output_species_summary,
        args.output_summary,
    )
    for path in output_paths:
        path.parent.mkdir(parents=True, exist_ok=True)

    metadata = pd.read_csv(args.metadata, sep="\t", dtype=str, keep_default_na=False)
    copy_haplotypes = pd.read_csv(
        args.copy_haplotypes, sep="\t", dtype=str, keep_default_na=False
    )
    haplotype_summary = pd.read_csv(
        args.haplotype_summary, sep="\t", dtype=str, keep_default_na=False
    )
    sequences = read_fasta_dict(args.haplotype_fasta)
    references = list(read_fasta(args.reference_fasta))
    if len(references) != 1:
        raise ValueError("--reference-fasta must contain exactly one sequence")
    _reference_name, reference = references[0]
    if not reference or any(base not in IUPAC for base in reference):
        raise ValueError("Reference must be one non-empty nucleotide sequence")
    validate_inputs(
        metadata, copy_haplotypes, haplotype_summary, sequences, args.group_by
    )

    capped_weights = capped_species_haplotype_weights(
        metadata, args.group_by, args.cap_copies_per_unit
    )
    raw_counts = (
        metadata.groupby(["species", "haplotype_id"], sort=False)
        .size()
        .to_dict()
    )
    species_order = list(dict.fromkeys(metadata["species"].astype(str)))
    species_rows: List[Dict[str, object]] = []
    species_anchor: Dict[str, Dict[str, object]] = {}
    for species in species_order:
        species_haps = sorted(
            set(metadata.loc[metadata["species"] == species, "haplotype_id"].astype(str))
        )
        species_sequences = {hap: sequences[hap] for hap in species_haps}
        species_weights = {
            hap: capped_weights.get((species, hap), 0.0) for hap in species_haps
        }
        row, _phased_consensus = choose_species_anchor(
            species,
            species_sequences,
            species_weights,
            reference,
            args.orientation_mode,
        )
        anchor_status_reasons: List[str] = []
        if float(row["anchor_capped_weight"]) < args.min_anchor_capped_weight:
            anchor_status_reasons.append("LOW_ANCHOR_WEIGHT")
        if float(row["reference_similarity"]) < args.min_anchor_similarity:
            anchor_status_reasons.append("LOW_REFERENCE_SIMILARITY")
        if int(row["edit_distance_margin"]) <= 0:
            anchor_status_reasons.append("AMBIGUOUS_ROTATION_OR_ORIENTATION")
        row["anchor_status"] = "PASS" if not anchor_status_reasons else ";".join(anchor_status_reasons)
        species_rows.append(row)
        species_anchor[species] = row
        eprint(
            f"[{species}] orientation={row['orientation_to_reference']}; "
            f"rotation_left={row['rotation_left_bp']} "
            f"(signed {row['rotation_signed_bp']:+d}); "
            f"reference_similarity={float(row['reference_similarity']):.3f}"
        )

    phase_rows: List[Dict[str, object]] = []
    canonical_by_species_hap: Dict[Tuple[str, str], str] = {}
    for species in species_order:
        anchor = species_anchor[species]
        orientation = str(anchor["orientation_to_reference"])
        anchor_length = int(anchor["anchor_length"])
        anchor_signed_rotation = int(anchor["rotation_signed_bp"])
        species_haps = sorted(
            set(metadata.loc[metadata["species"] == species, "haplotype_id"].astype(str))
        )
        for haplotype in species_haps:
            raw_sequence = sequences[haplotype]
            oriented = revcomp(raw_sequence) if orientation == "reverse_complement" else raw_sequence
            expected_rotation = anchor_signed_rotation % len(oriented)
            rotation = choose_haplotype_rotation(
                oriented,
                reference,
                expected_rotation,
                args.phase_search_radius,
                args.kmer_size,
                args.candidate_count,
                args.alignment_band,
            )
            canonical = rotate_left(oriented, int(rotation["rotation_left_bp"]))
            canonical_by_species_hap[(species, haplotype)] = canonical
            raw_copy_count = int(raw_counts.get((species, haplotype), 0))
            capped_weight = float(capped_weights.get((species, haplotype), 0.0))
            phase_rows.append(
                {
                    "species": species,
                    "input_haplotype_id": haplotype,
                    "input_length": len(raw_sequence),
                    "orientation_to_reference": orientation,
                    "species_anchor_length": anchor_length,
                    "species_anchor_rotation_signed_bp": anchor_signed_rotation,
                    **rotation,
                    "raw_copy_count": raw_copy_count,
                    "capped_analysis_weight": capped_weight,
                    "canonical_sequence": canonical,
                }
            )

    phase_map = pd.DataFrame(phase_rows)
    canonical_counts: Counter[str] = Counter()
    for row in phase_map.itertuples(index=False):
        canonical_counts[str(row.canonical_sequence)] += int(row.raw_copy_count)
    sorted_sequences = sorted(canonical_counts, key=lambda seq: (-canonical_counts[seq], seq))
    sequence_to_haplotype = {
        sequence: f"HAP{i:07d}" for i, sequence in enumerate(sorted_sequences, 1)
    }
    phase_map["output_haplotype_id"] = phase_map["canonical_sequence"].map(
        sequence_to_haplotype
    )

    mapping_columns = [
        "species",
        "input_haplotype_id",
        "output_haplotype_id",
        "orientation_to_reference",
        "rotation_left_bp",
        "rotation_signed_bp",
        "expected_rotation_left_bp",
        "expected_rotation_signed_bp",
        "residual_rotation_bp",
        "phase_tied",
        "best_edit_distance",
        "second_best_edit_distance",
        "edit_distance_margin",
        "kmer_votes",
        "ungapped_matches",
        "candidate_rotations_tested",
        "raw_copy_count",
        "capped_analysis_weight",
        "canonical_sequence",
    ]
    phase_map[mapping_columns].to_csv(args.output_phase_map, sep="\t", index=False)

    map_for_metadata = phase_map[
        [
            "species",
            "input_haplotype_id",
            "output_haplotype_id",
            "orientation_to_reference",
            "rotation_left_bp",
            "rotation_signed_bp",
            "residual_rotation_bp",
            "phase_tied",
        ]
    ].copy()
    output_metadata = metadata.rename(
        columns={"haplotype_id": "species_phase_haplotype_id"}
    ).merge(
        map_for_metadata,
        left_on=["species", "species_phase_haplotype_id"],
        right_on=["species", "input_haplotype_id"],
        how="left",
        validate="many_to_one",
    )
    if output_metadata["output_haplotype_id"].isna().any():
        raise ValueError("Some copies could not be mapped to a reference-phased haplotype")
    output_metadata = output_metadata.drop(columns=["input_haplotype_id"])
    output_metadata = output_metadata.rename(
        columns={
            "output_haplotype_id": "haplotype_id",
            "orientation_to_reference": "cross_species_phase_orientation",
            "rotation_left_bp": "cross_species_rotation_left_bp",
            "rotation_signed_bp": "cross_species_rotation_signed_bp",
            "residual_rotation_bp": "cross_species_rotation_residual_bp",
            "phase_tied": "cross_species_phase_tied",
        }
    )
    original_columns = list(metadata.columns)
    haplotype_position = original_columns.index("haplotype_id")
    base_columns = [column for column in original_columns if column != "haplotype_id"]
    ordered_base = base_columns[:haplotype_position] + [
        "species_phase_haplotype_id",
        "haplotype_id",
        "cross_species_phase_orientation",
        "cross_species_rotation_left_bp",
        "cross_species_rotation_signed_bp",
        "cross_species_rotation_residual_bp",
        "cross_species_phase_tied",
    ] + base_columns[haplotype_position:]
    output_metadata = output_metadata[ordered_base]
    output_metadata.to_csv(args.output_metadata, sep="\t", index=False)
    output_metadata[["seq_id", "haplotype_id"]].to_csv(
        args.output_copy_haplotypes, sep="\t", index=False
    )

    write_fasta(
        (
            (
                str(row.seq_id),
                canonical_by_species_hap[
                    (str(row.species), str(row.species_phase_haplotype_id))
                ],
            )
            for row in output_metadata.itertuples(index=False)
        ),
        args.output_all_copies_fasta,
    )

    hap_rows: List[Dict[str, object]] = []
    for sequence in sorted_sequences:
        haplotype_id = sequence_to_haplotype[sequence]
        subset = output_metadata.loc[output_metadata["haplotype_id"] == haplotype_id]
        hap_rows.append(
            {
                "haplotype_id": haplotype_id,
                "n_copies": len(subset),
                "n_species": subset["species"].nunique(),
                "n_arrays": subset[["species", "array_key"]].drop_duplicates().shape[0],
                "n_regions": subset[["species", args.group_by]].drop_duplicates().shape[0],
                "species_list": ",".join(sorted(subset["species"].astype(str).unique())),
                "length": len(sequence),
                "n_fraction": sequence.count("N") / len(sequence),
                "sequence": sequence,
            }
        )
    output_haplotypes = pd.DataFrame(hap_rows)
    output_haplotypes.to_csv(args.output_haplotype_summary, sep="\t", index=False)
    write_fasta(
        (
            (f"{row.haplotype_id};size={row.n_copies};", row.sequence)
            for row in output_haplotypes.itertuples(index=False)
        ),
        args.output_haplotype_fasta,
    )

    species_summary = pd.DataFrame(species_rows)
    tied_weight_by_species = (
        phase_map.assign(
            tied_weight=lambda frame: frame["capped_analysis_weight"]
            * frame["phase_tied"].astype(bool).astype(float)
        )
        .groupby("species", sort=False)
        .agg(
            phase_tied_capped_weight=("tied_weight", "sum"),
            phase_total_capped_weight=("capped_analysis_weight", "sum"),
            max_abs_residual_rotation_bp=(
                "residual_rotation_bp",
                lambda values: max(abs(int(value)) for value in values),
            ),
        )
        .reset_index()
    )
    tied_weight_by_species["phase_tied_weight_fraction"] = (
        tied_weight_by_species["phase_tied_capped_weight"]
        / tied_weight_by_species["phase_total_capped_weight"].clip(lower=1.0)
    )
    species_summary = species_summary.merge(
        tied_weight_by_species, on="species", how="left", validate="one_to_one"
    )
    species_summary.to_csv(args.output_species_summary, sep="\t", index=False)

    input_hap_species = metadata[["species", "haplotype_id"]].drop_duplicates()
    input_species_counts = input_hap_species.groupby("haplotype_id")["species"].nunique()
    input_cross_species_haps = set(input_species_counts.loc[input_species_counts >= 2].index)
    input_cross_species_copies = int(
        metadata["haplotype_id"].isin(input_cross_species_haps).sum()
    )
    output_cross_species_haps = set(
        output_haplotypes.loc[
            output_haplotypes["n_species"] >= 2, "haplotype_id"
        ]
    )
    output_cross_species_copies = int(
        output_metadata["haplotype_id"].isin(output_cross_species_haps).sum()
    )
    total_capped_weight = float(phase_map["capped_analysis_weight"].sum())
    tied_capped_weight = float(
        phase_map.loc[
            phase_map["phase_tied"].astype(bool), "capped_analysis_weight"
        ].sum()
    )
    tied_fraction = tied_capped_weight / max(total_capped_weight, 1.0)
    anchor_failures = int((species_summary["anchor_status"] != "PASS").sum())
    phase_status = (
        "PASS"
        if anchor_failures == 0 and tied_fraction <= args.max_ambiguous_weight_fraction
        else "FAIL_PHASE_QC"
    )
    summary_rows = [
        {"metric": "phase_qc_status", "value": phase_status},
        {"metric": "reference_length", "value": len(reference)},
        {"metric": "species", "value": len(species_order)},
        {"metric": "copies", "value": len(output_metadata)},
        {"metric": "input_haplotypes", "value": len(sequences)},
        {"metric": "reference_phased_haplotypes", "value": len(output_haplotypes)},
        {
            "metric": "haplotypes_collapsed_after_reference_phasing",
            "value": len(sequences) - len(output_haplotypes),
        },
        {"metric": "species_anchor_failures", "value": anchor_failures},
        {"metric": "phase_tied_capped_weight_fraction", "value": tied_fraction},
        {
            "metric": "input_cross_species_exact_haplotype_copies",
            "value": input_cross_species_copies,
        },
        {
            "metric": "reference_phased_cross_species_exact_haplotype_copies",
            "value": output_cross_species_copies,
        },
        {
            "metric": "reference_phased_cross_species_exact_haplotype_fraction",
            "value": output_cross_species_copies / max(len(output_metadata), 1),
        },
    ]
    pd.DataFrame(summary_rows).to_csv(args.output_summary, sep="\t", index=False)
    eprint(
        f"Reference-phase normalization: {len(sequences):,} input haplotypes -> "
        f"{len(output_haplotypes):,}; cross-species exact copies "
        f"{input_cross_species_copies:,} -> {output_cross_species_copies:,}; "
        f"status={phase_status}"
    )
    if phase_status != "PASS" and not args.allow_poor_phase:
        raise RuntimeError(
            "Cross-species phase QC failed. Inspect cross_species_phase_species_summary.tsv "
            "and cross_species_phase_map.tsv. Use --allow-poor-phase only for diagnosis."
        )


if __name__ == "__main__":
    main()
