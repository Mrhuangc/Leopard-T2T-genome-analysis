#!/usr/bin/env python3
"""Build a fixed-coordinate SAT113 alignment against one monomer reference.

Ordinary joint multiple alignment can accumulate a different insertion column for
many highly divergent satellite monomers and inflate a ~113-bp repeat to hundreds
or thousands of columns.  This script instead adds every exact haplotype to one
ungapped reference with MAFFT ``--addfragments --keeplength``.  The resulting
alignment therefore has exactly one column per reference position.

If no reference is supplied, a high-support observed haplotype is used as a seed.
The seed is refined to a unit-capped, copy-weighted consensus and all haplotypes
are aligned once more to that consensus.  Insertions relative to the reference
are deliberately omitted from subtype fitting but their loss is reflected in the
per-haplotype QC table.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, Optional, Tuple

import numpy as np
import pandas as pd


CANONICAL = set("ACGT")


def eprint(*args: object) -> None:
    print(*args, file=sys.stderr)


def normalize_fasta_id(name: str) -> str:
    return name.split()[0].split(";")[0]


def clean_sequence(sequence: str) -> str:
    return "".join(sequence.split()).replace(".", "-").upper()


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
        records[name] = sequence
    return records


def write_fasta(records: Iterable[Tuple[str, str]], path: Path, width: int = 80) -> None:
    with open(path, "w", encoding="utf-8") as out:
        for name, sequence in records:
            out.write(f">{name}\n")
            for start in range(0, len(sequence), width):
                out.write(sequence[start : start + width] + "\n")


def weighted_mode(values: pd.Series) -> int:
    parsed = pd.to_numeric(values, errors="coerce").dropna().round().astype(int)
    parsed = parsed.loc[parsed > 0]
    if parsed.empty:
        raise ValueError("No positive monomer lengths were available")
    counts = parsed.value_counts()
    maximum = counts.max()
    return int(sorted(counts.loc[counts == maximum].index)[0])


def choose_seed_reference(
    haplotypes: Mapping[str, str],
    summary: pd.DataFrame,
    metadata: pd.DataFrame,
) -> Tuple[str, str, int]:
    if "local_consensus_length" in metadata.columns:
        target_length = weighted_mode(metadata["local_consensus_length"])
    else:
        expanded_lengths: List[int] = []
        for row in summary.itertuples(index=False):
            count = max(1, int(getattr(row, "n_copies", 1)))
            expanded_lengths.extend([int(getattr(row, "length"))] * count)
        target_length = int(np.median(expanded_lengths))

    rows: List[Tuple[Tuple[object, ...], str]] = []
    for row in summary.itertuples(index=False):
        hap = str(row.haplotype_id)
        if hap not in haplotypes:
            continue
        sequence = haplotypes[hap].replace("-", "")
        ambiguous = sum(base not in CANONICAL for base in sequence) / max(len(sequence), 1)
        score = (
            abs(len(sequence) - target_length),
            ambiguous,
            -int(getattr(row, "n_species", 1)),
            -int(getattr(row, "n_arrays", 1)),
            -int(getattr(row, "n_copies", 1)),
            hap,
        )
        rows.append((score, hap))
    if not rows:
        raise ValueError("No haplotype in the summary was found in the FASTA file")
    rows.sort(key=lambda item: item[0])
    seed_id = rows[0][1]
    return seed_id, haplotypes[seed_id].replace("-", ""), target_length


def capped_haplotype_weights(
    metadata: pd.DataFrame,
    group_by: str,
    cap_per_unit: int,
) -> Dict[str, float]:
    if group_by not in metadata.columns:
        raise ValueError(f"Metadata lacks grouping column: {group_by}")
    required = {"species", "haplotype_id", group_by}
    missing = sorted(required - set(metadata.columns))
    if missing:
        raise ValueError(f"Metadata lacks required columns: {', '.join(missing)}")
    frame = metadata.copy()
    frame["unit_uid"] = frame["species"].astype(str) + "|" + frame[group_by].astype(str)
    weights: Counter[str] = Counter()
    for _unit, subset in frame.groupby("unit_uid", sort=False):
        scale = 1.0
        if cap_per_unit > 0 and len(subset) > cap_per_unit:
            scale = cap_per_unit / len(subset)
        for hap, count in subset["haplotype_id"].astype(str).value_counts().items():
            weights[hap] += float(count) * scale
    return dict(weights)


def run_mafft_add(
    mafft: str,
    query_fasta: Path,
    reference_sequence: str,
    threads: int,
    strategy: str,
    workdir: Path,
    label: str,
    log_handle,
) -> Tuple[Dict[str, str], str]:
    reference_path = workdir / f"{label}.reference.fa"
    output_path = workdir / f"{label}.alignment.fa"
    write_fasta([("SAT113_REFERENCE", reference_sequence)], reference_path)
    command = [
        mafft,
        "--thread",
        str(threads),
        {"auto": "--auto", "multipair": "--multipair", "6merpair": "--6merpair"}[
            strategy
        ],
        "--addfragments",
        str(query_fasta),
        "--keeplength",
        str(reference_path),
    ]
    log_handle.write("COMMAND\t" + " ".join(command) + "\n")
    log_handle.flush()
    with open(output_path, "w", encoding="utf-8") as out:
        result = subprocess.run(
            command,
            stdout=out,
            stderr=log_handle,
            text=True,
            check=False,
        )
    if result.returncode != 0:
        raise RuntimeError(f"MAFFT failed in {label} with exit code {result.returncode}")
    aligned = read_fasta_dict(output_path)
    if "SAT113_REFERENCE" not in aligned:
        raise ValueError("MAFFT output does not contain SAT113_REFERENCE")
    aligned_reference = aligned.pop("SAT113_REFERENCE")
    expected_length = len(reference_sequence)
    if len(aligned_reference) != expected_length or "-" in aligned_reference:
        raise ValueError(
            "MAFFT did not preserve the ungapped reference length; "
            f"expected {expected_length}, observed {len(aligned_reference)}"
        )
    bad_lengths = {name: len(seq) for name, seq in aligned.items() if len(seq) != expected_length}
    if bad_lengths:
        example = next(iter(bad_lengths.items()))
        raise ValueError(
            f"MAFFT --keeplength output has inconsistent query lengths; example={example}"
        )
    return aligned, aligned_reference


def consensus_from_alignment(
    alignment: Mapping[str, str],
    reference: str,
    weights: Mapping[str, float],
) -> Tuple[str, np.ndarray]:
    counts = np.zeros((len(reference), 4), dtype=float)
    bases = "ACGT"
    for hap, sequence in alignment.items():
        weight = float(weights.get(hap, 0.0))
        if weight <= 0:
            continue
        chars = np.frombuffer(sequence.encode("ascii"), dtype="S1")
        for code, base in enumerate(bases):
            counts[:, code] += weight * (chars == base.encode("ascii"))
    totals = counts.sum(axis=1)
    consensus: List[str] = []
    for position, total in enumerate(totals):
        if total <= 0:
            fallback = reference[position]
            consensus.append(fallback if fallback in CANONICAL else "N")
        else:
            consensus.append(bases[int(np.argmax(counts[position]))])
    return "".join(consensus), totals


def alignment_qc_table(
    alignment: Mapping[str, str],
    haplotypes: Mapping[str, str],
    reference: str,
    raw_counts: Mapping[str, float],
    capped_weights: Mapping[str, float],
    min_coverage: float,
) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    for hap, raw_sequence in haplotypes.items():
        sequence = alignment[hap]
        non_gap = sum(base != "-" for base in sequence)
        canonical_positions = [i for i, base in enumerate(sequence) if base in CANONICAL]
        comparable = [i for i in canonical_positions if reference[i] in CANONICAL]
        mismatches = sum(sequence[i] != reference[i] for i in comparable)
        raw_length = len(raw_sequence.replace("-", ""))
        canonical_coverage = len(canonical_positions) / len(reference)
        retained_fraction = non_gap / max(raw_length, 1)
        rows.append(
            {
                "haplotype_id": hap,
                "raw_length": raw_length,
                "reference_length": len(reference),
                "mapped_non_gap_bases": non_gap,
                "mapped_canonical_bases": len(canonical_positions),
                "reference_coverage": non_gap / len(reference),
                "canonical_reference_coverage": canonical_coverage,
                "query_fraction_retained": retained_fraction,
                "mismatch_fraction_comparable": mismatches / max(len(comparable), 1),
                "raw_copy_count": float(raw_counts.get(hap, 0.0)),
                "capped_analysis_weight": float(capped_weights.get(hap, 0.0)),
                "passes_alignment_qc": bool(
                    canonical_coverage >= min_coverage and retained_fraction >= min_coverage
                ),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--haplotype-fasta", required=True, type=Path)
    parser.add_argument("--haplotype-summary", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--output-alignment", required=True, type=Path)
    parser.add_argument("--output-reference", required=True, type=Path)
    parser.add_argument("--output-qc", required=True, type=Path)
    parser.add_argument("--output-summary", required=True, type=Path)
    parser.add_argument("--output-coordinate-map", required=True, type=Path)
    parser.add_argument("--log", required=True, type=Path)
    parser.add_argument("--reference-fasta", type=Path)
    parser.add_argument("--mafft", default="mafft")
    parser.add_argument(
        "--mafft-strategy",
        choices=["multipair", "auto", "6merpair"],
        default="multipair",
        help="multipair is recommended for short but highly divergent SAT113 copies.",
    )
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--group-by", choices=["array_key", "region_key"], default="region_key")
    parser.add_argument("--cap-copies-per-unit", type=int, default=100)
    parser.add_argument("--refinement-rounds", type=int, default=2)
    parser.add_argument("--min-alignment-coverage", type=float, default=0.80)
    parser.add_argument("--max-failed-weight-fraction", type=float, default=0.05)
    parser.add_argument("--allow-poor-alignment", action="store_true")
    args = parser.parse_args()

    if shutil.which(args.mafft) is None:
        raise RuntimeError(f"MAFFT executable not found: {args.mafft}")
    if args.refinement_rounds < 0:
        raise ValueError("--refinement-rounds must be >= 0")
    if not 0 < args.min_alignment_coverage <= 1:
        raise ValueError("--min-alignment-coverage must be in (0, 1]")

    for path in (
        args.output_alignment,
        args.output_reference,
        args.output_qc,
        args.output_summary,
        args.output_coordinate_map,
        args.log,
    ):
        path.parent.mkdir(parents=True, exist_ok=True)

    haplotypes = read_fasta_dict(args.haplotype_fasta)
    summary = pd.read_csv(args.haplotype_summary, sep="\t", dtype={"haplotype_id": str})
    metadata = pd.read_csv(args.metadata, sep="\t", dtype={"haplotype_id": str})
    hap_order = list(summary["haplotype_id"].astype(str))
    missing = sorted(set(hap_order) - set(haplotypes))
    extra = sorted(set(haplotypes) - set(hap_order))
    if missing or extra:
        raise ValueError(
            f"Haplotype FASTA/summary mismatch: missing={len(missing)}, extra={len(extra)}"
        )
    haplotypes = {hap: haplotypes[hap].replace("-", "") for hap in hap_order}
    raw_counts = metadata["haplotype_id"].astype(str).value_counts().to_dict()
    capped_weights = capped_haplotype_weights(
        metadata, args.group_by, args.cap_copies_per_unit
    )

    seed_id = "user_reference"
    target_length: Optional[int] = None
    if args.reference_fasta:
        references = list(read_fasta(args.reference_fasta))
        if len(references) != 1:
            raise ValueError("--reference-fasta must contain exactly one sequence")
        _name, reference = references[0]
        reference = reference.replace("-", "")
    else:
        seed_id, reference, target_length = choose_seed_reference(
            haplotypes, summary, metadata
        )
    if not reference or any(base not in set("ACGTRYSWKMBDHVN") for base in reference):
        raise ValueError("Reference must be a non-empty nucleotide sequence")

    with tempfile.TemporaryDirectory(prefix="sat113_reference_alignment_") as tmp_name:
        tmp = Path(tmp_name)
        clean_queries = tmp / "haplotypes.clean.fa"
        write_fasta(((hap, haplotypes[hap]) for hap in hap_order), clean_queries)
        with open(args.log, "w", encoding="utf-8") as log_handle:
            log_handle.write(f"seed_id\t{seed_id}\n")
            log_handle.write(f"seed_length\t{len(reference)}\n")
            if target_length is not None:
                log_handle.write(f"target_length_mode\t{target_length}\n")
            last_alignment: Optional[Dict[str, str]] = None
            rounds_completed = 0
            reference_changed_after_last_alignment = False
            for round_index in range(args.refinement_rounds):
                aligned, _aligned_reference = run_mafft_add(
                    args.mafft,
                    clean_queries,
                    reference,
                    args.threads,
                    args.mafft_strategy,
                    tmp,
                    f"refine_{round_index + 1}",
                    log_handle,
                )
                if set(aligned) != set(hap_order):
                    raise ValueError(
                        "MAFFT output query identifiers do not match input haplotypes"
                    )
                new_reference, _totals = consensus_from_alignment(
                    aligned, reference, capped_weights
                )
                changes = sum(a != b for a, b in zip(reference, new_reference))
                log_handle.write(
                    f"refinement_round\t{round_index + 1}\treference_changes\t{changes}\n"
                )
                log_handle.flush()
                rounds_completed += 1
                last_alignment = aligned
                if new_reference == reference:
                    reference_changed_after_last_alignment = False
                    break
                reference = new_reference
                reference_changed_after_last_alignment = True

            # If the last refinement changed the reference, align once to the
            # final consensus.  Otherwise the last alignment can be reused.
            if (
                last_alignment is None
                or args.refinement_rounds == 0
                or reference_changed_after_last_alignment
            ):
                final_alignment, _ = run_mafft_add(
                    args.mafft,
                    clean_queries,
                    reference,
                    args.threads,
                    args.mafft_strategy,
                    tmp,
                    "final",
                    log_handle,
                )
            else:
                final_alignment = last_alignment

    if set(final_alignment) != set(hap_order):
        raise ValueError("Final alignment is missing or adding haplotypes")
    write_fasta(((hap, final_alignment[hap]) for hap in hap_order), args.output_alignment)
    write_fasta([("SAT113_FIXED_REFERENCE", reference)], args.output_reference)

    qc = alignment_qc_table(
        final_alignment,
        haplotypes,
        reference,
        raw_counts,
        capped_weights,
        args.min_alignment_coverage,
    )
    qc.to_csv(args.output_qc, sep="\t", index=False)
    failed = ~qc["passes_alignment_qc"]
    failed_raw_fraction = float(
        qc.loc[failed, "raw_copy_count"].sum() / max(qc["raw_copy_count"].sum(), 1.0)
    )
    failed_weight_fraction = float(
        qc.loc[failed, "capped_analysis_weight"].sum()
        / max(qc["capped_analysis_weight"].sum(), 1.0)
    )
    occupancy = np.zeros(len(reference), dtype=float)
    total_weight = sum(capped_weights.values())
    for hap, sequence in final_alignment.items():
        weight = capped_weights.get(hap, 0.0)
        occupancy += weight * np.array([base in CANONICAL for base in sequence], dtype=float)
    occupancy /= max(total_weight, 1.0)
    pd.DataFrame(
        {
            "reference_position_1based": np.arange(1, len(reference) + 1),
            "reference_base": list(reference),
            "weighted_canonical_occupancy": occupancy,
        }
    ).to_csv(args.output_coordinate_map, sep="\t", index=False)

    status = (
        "PASS"
        if failed_weight_fraction <= args.max_failed_weight_fraction
        else "FAIL_TOO_MANY_POORLY_ALIGNED_HAPLOTYPES"
    )
    summary_rows = [
        {"metric": "alignment_status", "value": status},
        {"metric": "reference_seed_id", "value": seed_id},
        {"metric": "reference_length", "value": len(reference)},
        {"metric": "refinement_rounds_completed", "value": rounds_completed},
        {"metric": "haplotypes", "value": len(hap_order)},
        {"metric": "raw_copies", "value": float(qc["raw_copy_count"].sum())},
        {"metric": "failed_haplotypes", "value": int(failed.sum())},
        {"metric": "failed_raw_copy_fraction", "value": failed_raw_fraction},
        {"metric": "failed_capped_weight_fraction", "value": failed_weight_fraction},
        {"metric": "median_canonical_reference_coverage", "value": qc["canonical_reference_coverage"].median()},
        {"metric": "minimum_weighted_site_occupancy", "value": float(occupancy.min())},
        {"metric": "median_weighted_site_occupancy", "value": float(np.median(occupancy))},
    ]
    pd.DataFrame(summary_rows).to_csv(args.output_summary, sep="\t", index=False)
    eprint(
        f"Fixed reference length={len(reference)}; failed capped weight="
        f"{failed_weight_fraction:.3%}; status={status}"
    )
    if status != "PASS" and not args.allow_poor_alignment:
        raise RuntimeError(
            "Reference alignment QC failed. Inspect the alignment QC table; this can indicate "
            "incorrect circular phasing/orientation or an unsuitable reference. Re-run with a "
            "validated --reference-fasta. Use --allow-poor-alignment only for diagnosis."
        )


if __name__ == "__main__":
    main()
