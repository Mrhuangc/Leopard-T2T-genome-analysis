#!/usr/bin/env python3
"""Prepare jointly phased SAT113 monomers from satellite_array_pipeline_v2.1 outputs.

The script deliberately reuses the two phase transformations already recorded by
satellite_array_pipeline_v2.1.py:

1. final.phase_report.tsv: raw copy -> local-array phase
2. array_consensus.phase_report.tsv (global_consensus_phase): local-array phase
   -> global SAT phase

This avoids treating different circular starts of the same monomer as biological
variants.  Input sequences can come from the original TRASH2 CSV, a copy FASTA,
or the assembly FASTA.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, Optional, Sequence, Tuple

import pandas as pd


DNA_RE = re.compile(r"^[ACGTRYSWKMBDHVNacgtryswkmbdhvn.-]+$")
COPY_ID_RE = re.compile(r"^(?P<chrom>.+):(?P<start>\d+)-(?P<end>\d+)\((?P<strand>[+-])\)")
REVCOMP = str.maketrans(
    "ACGTRYMKSWBDHVNacgtrymkswbdhvn",
    "TGCAYRKMSWVHDBNtgcayrkmswvhdbn",
)


def eprint(*args: object) -> None:
    print(*args, file=sys.stderr)


def open_text(path: Path, mode: str = "rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode, encoding="utf-8", newline="")


def read_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    name: Optional[str] = None
    chunks: List[str] = []
    with open_text(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                name = line[1:].strip().split()[0]
                chunks = []
            else:
                if name is None:
                    raise ValueError(f"FASTA sequence before header in {path}")
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks).upper()


def write_fasta(records: Iterable[Tuple[str, str]], path: Path, width: int = 80) -> None:
    with open(path, "w", encoding="utf-8") as out:
        for name, seq in records:
            out.write(f">{name}\n")
            for i in range(0, len(seq), width):
                out.write(seq[i : i + width] + "\n")


def revcomp(seq: str) -> str:
    return seq.translate(REVCOMP)[::-1]


def clean_sequence(value: str) -> str:
    return re.sub(r"[.\s-]", "", value).upper()


def parse_integer(value: object) -> Optional[int]:
    """Parse ordinary or scientific-notation integer text without truncation."""
    text = str(value).strip().replace(",", "")
    if not text:
        return None
    try:
        number = float(text)
    except ValueError:
        return None
    if not math.isfinite(number) or abs(number - round(number)) > 1e-6:
        return None
    return int(round(number))


def resolve_path(value: str, manifest_dir: Path) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = manifest_dir / path
    return path.resolve()


def normalize_orientation(value: object) -> str:
    text = str(value).strip().lower().replace("-", "_")
    if text in {"reverse", "reverse_complement", "revcomp", "rc", "minus"}:
        return "reverse_complement"
    if text in {"forward", "plus", "+", "none", "nan", ""}:
        return "forward"
    raise ValueError(f"Unrecognized phase orientation: {value!r}")


def transform_sequence(seq: str, orientation: object, rotation_left_bp: object) -> str:
    seq = clean_sequence(seq)
    if normalize_orientation(orientation) == "reverse_complement":
        seq = revcomp(seq)
    if not seq:
        return seq
    rotation = parse_integer(rotation_left_bp)
    if rotation is None:
        rotation = 0
    rotation %= len(seq)
    return seq[rotation:] + seq[:rotation]


def load_copy_table(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    required = {
        "copy_id",
        "chrom",
        "array_key",
        "start",
        "end",
        "strand",
        "phase_status",
        "quality_tier",
    }
    missing = sorted(required - set(table.columns))
    if missing:
        raise ValueError(f"{path}: missing columns: {', '.join(missing)}")
    for column in ("start", "end", "raw_length", "local_consensus_length"):
        if column in table.columns:
            table[column] = table[column].map(parse_integer)
    return table


def select_phase_rows(path: Path, key_column: str, preferred_stage: str) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    required = {key_column, "orientation", "rotation_left_bp"}
    missing = sorted(required - set(table.columns))
    if missing:
        raise ValueError(f"{path}: missing columns: {', '.join(missing)}")
    if "stage" in table.columns and (table["stage"] == preferred_stage).any():
        table = table.loc[table["stage"] == preferred_stage].copy()
    # If a legacy report lacks the expected stage, the last row for each key is
    # the least surprising behavior and is explicitly reported in QC.
    table = table.drop_duplicates(subset=[key_column], keep="last")
    return table[[key_column, "orientation", "rotation_left_bp"]].rename(
        columns={
            "orientation": f"{preferred_stage}_orientation",
            "rotation_left_bp": f"{preferred_stage}_rotation_left_bp",
        }
    )


def extract_copy_tuple(text: str) -> Optional[Tuple[str, int, int, str]]:
    match = COPY_ID_RE.match(text)
    if not match:
        return None
    return (
        match.group("chrom"),
        int(match.group("start")),
        int(match.group("end")),
        match.group("strand"),
    )


def best_sequence_cell(row: Sequence[str]) -> Optional[str]:
    candidates: List[str] = []
    for cell in row:
        compact = re.sub(r"\s+", "", str(cell))
        if len(compact) >= 20 and DNA_RE.fullmatch(compact):
            sequence = clean_sequence(compact)
            if sequence:
                candidates.append(sequence)
    return max(candidates, key=len) if candidates else None


def sniff_rows(path: Path) -> Iterator[List[str]]:
    with open_text(path) as handle:
        sample = handle.read(65536)
        handle.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
        except csv.Error:
            dialect = csv.excel
        reader = csv.reader(handle, dialect)
        for row in reader:
            if row and any(str(cell).strip() for cell in row):
                yield [str(cell).strip() for cell in row]


def load_sequences_from_trash_csv(path: Path, copies: pd.DataFrame) -> Tuple[Dict[str, str], Dict[str, int]]:
    """Map a headered or headerless TRASH2 CSV to copy_id by coordinates.

    The parser does not hard-code the column positions.  It scans numeric fields
    and accepts scientific notation (for example 5e+05), then chooses a target
    coordinate whose span best matches the detected DNA sequence length.
    """
    target_by_chrom: Dict[str, Dict[Tuple[int, int], List[str]]] = defaultdict(lambda: defaultdict(list))
    target_global: Dict[Tuple[int, int], List[Tuple[str, str]]] = defaultdict(list)
    chroms = set(copies["chrom"].astype(str))
    for row in copies.itertuples(index=False):
        chrom = str(row.chrom)
        start, end = int(row.start), int(row.end)
        target_by_chrom[chrom][(start, end)].append(str(row.copy_id))
        target_global[(start, end)].append((chrom, str(row.copy_id)))

    found: Dict[str, str] = {}
    stats = Counter()
    for row in sniff_rows(path):
        stats["csv_rows"] += 1
        sequence = best_sequence_cell(row)
        if sequence is None:
            stats["rows_without_sequence"] += 1
            continue

        embedded_ids: List[Tuple[str, int, int, str]] = []
        for cell in row:
            parsed = extract_copy_tuple(cell)
            if parsed:
                embedded_ids.append(parsed)

        chrom_candidates = [cell for cell in row if cell in chroms]
        numeric = [(i, parse_integer(cell)) for i, cell in enumerate(row)]
        numeric = [(i, value) for i, value in numeric if value is not None]
        coordinate_pairs: set[Tuple[int, int]] = set()
        for idx, (i, first) in enumerate(numeric):
            for j, second in numeric[idx + 1 :]:
                if j - i > 5:
                    break
                if first is not None and second is not None and second >= first:
                    coordinate_pairs.add((first, second))

        candidates: List[Tuple[int, str]] = []
        for chrom, start, end, _strand in embedded_ids:
            for copy_id in target_by_chrom.get(chrom, {}).get((start, end), []):
                candidates.append((abs((end - start + 1) - len(sequence)), copy_id))
        for chrom in chrom_candidates:
            for start, end in coordinate_pairs:
                for copy_id in target_by_chrom[chrom].get((start, end), []):
                    candidates.append((abs((end - start + 1) - len(sequence)), copy_id))
        if not chrom_candidates:
            for start, end in coordinate_pairs:
                global_hits = target_global.get((start, end), [])
                if len(global_hits) == 1:
                    chrom, copy_id = global_hits[0]
                    candidates.append((abs((end - start + 1) - len(sequence)), copy_id))

        if not candidates:
            stats["rows_without_coordinate_match"] += 1
            continue
        candidates.sort(key=lambda item: (item[0], item[1]))
        best_delta = candidates[0][0]
        best_ids = sorted({copy_id for delta, copy_id in candidates if delta == best_delta})
        if len(best_ids) > 1:
            # Duplicate coordinate records intentionally receive the same input
            # sequence; exact coordinate duplicates are removed before this step.
            stats["ambiguous_rows"] += 1
        for copy_id in best_ids:
            if copy_id in found and found[copy_id] != sequence:
                stats["conflicting_duplicate_rows"] += 1
                continue
            found[copy_id] = sequence
            stats["mapped_sequences"] += 1
    return found, dict(stats)


def load_sequences_from_copy_fasta(path: Path, copies: pd.DataFrame) -> Tuple[Dict[str, str], Dict[str, int]]:
    by_copy_id = set(copies["copy_id"].astype(str))
    by_tuple: Dict[Tuple[str, int, int, str], str] = {}
    for row in copies.itertuples(index=False):
        by_tuple[(str(row.chrom), int(row.start), int(row.end), str(row.strand))] = str(row.copy_id)
    found: Dict[str, str] = {}
    stats = Counter()
    for name, sequence in read_fasta(path):
        stats["fasta_records"] += 1
        copy_id: Optional[str] = name if name in by_copy_id else None
        if copy_id is None:
            parsed = extract_copy_tuple(name)
            if parsed:
                copy_id = by_tuple.get(parsed)
        if copy_id is None:
            stats["unmapped_records"] += 1
            continue
        found[copy_id] = clean_sequence(sequence)
        stats["mapped_sequences"] += 1
    return found, dict(stats)


def load_sequences_from_genome(path: Path, copies: pd.DataFrame, samtools: str) -> Tuple[Dict[str, str], Dict[str, int]]:
    if shutil.which(samtools) is None:
        raise RuntimeError(f"samtools executable not found: {samtools}")
    coord_to_copy: Dict[Tuple[str, int, int], str] = {}
    regions: List[str] = []
    for row in copies.itertuples(index=False):
        key = (str(row.chrom), int(row.start), int(row.end))
        coord_to_copy[key] = str(row.copy_id)
        regions.append(f"{key[0]}:{key[1]}-{key[2]}")
    found: Dict[str, str] = {}
    stats = Counter()
    with tempfile.TemporaryDirectory(prefix="sat113_faidx_") as tmp:
        region_path = Path(tmp) / "regions.txt"
        output_path = Path(tmp) / "copies.fa"
        region_path.write_text("\n".join(regions) + "\n", encoding="utf-8")
        with open(output_path, "w", encoding="utf-8") as out:
            command = [samtools, "faidx", "-r", str(region_path), str(path)]
            result = subprocess.run(command, stdout=out, stderr=subprocess.PIPE, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"samtools faidx failed for {path}:\n{result.stderr}")
        for name, sequence in read_fasta(output_path):
            stats["fasta_records"] += 1
            parsed = re.match(r"^(.+):(\d+)-(\d+)$", name)
            if not parsed:
                stats["unmapped_records"] += 1
                continue
            key = (parsed.group(1), int(parsed.group(2)), int(parsed.group(3)))
            copy_id = coord_to_copy.get(key)
            if copy_id is None:
                stats["unmapped_records"] += 1
                continue
            found[copy_id] = clean_sequence(sequence)
            stats["mapped_sequences"] += 1
    return found, dict(stats)


def load_source_sequences(
    source_type: str,
    source_path: Path,
    copies: pd.DataFrame,
    samtools: str,
) -> Tuple[Dict[str, str], Dict[str, int]]:
    normalized = source_type.strip().lower().replace("-", "_")
    if normalized in {"trash_csv", "csv", "trash2_csv"}:
        return load_sequences_from_trash_csv(source_path, copies)
    if normalized in {"copy_fasta", "copies_fasta", "fasta", "fa"}:
        return load_sequences_from_copy_fasta(source_path, copies)
    if normalized in {"genome_fasta", "genome", "assembly_fasta"}:
        return load_sequences_from_genome(source_path, copies, samtools)
    raise ValueError(f"Unsupported source_type={source_type!r}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Prepare globally phased SAT113 copies from v2.1 output reports."
    )
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--min-length-ratio", type=float, default=0.80)
    parser.add_argument("--max-length-ratio", type=float, default=1.20)
    parser.add_argument("--max-n-fraction", type=float, default=0.05)
    parser.add_argument(
        "--include-non-full-length",
        action="store_true",
        help="Do not restrict quality_tier to full_length (not recommended for the main test).",
    )
    args = parser.parse_args()

    manifest_path = args.manifest.resolve()
    manifest_dir = manifest_path.parent
    manifest = pd.read_csv(manifest_path, sep="\t", dtype=str, keep_default_na=False)
    required_manifest = {
        "species",
        "copy_table",
        "sequence_source",
        "source_type",
        "final_phase_report",
        "array_phase_report",
    }
    missing = sorted(required_manifest - set(manifest.columns))
    if missing:
        raise ValueError(f"Manifest is missing columns: {', '.join(missing)}")
    if manifest["species"].duplicated().any():
        duplicates = ", ".join(manifest.loc[manifest["species"].duplicated(), "species"])
        raise ValueError(f"Each species must occur once in the manifest; duplicates: {duplicates}")

    args.outdir.mkdir(parents=True, exist_ok=True)
    all_records: List[Dict[str, object]] = []
    qc_rows: List[Dict[str, object]] = []

    for dataset in manifest.itertuples(index=False):
        species = str(dataset.species).strip()
        copy_path = resolve_path(str(dataset.copy_table), manifest_dir)
        source_path = resolve_path(str(dataset.sequence_source), manifest_dir)
        final_phase_path = resolve_path(str(dataset.final_phase_report), manifest_dir)
        array_phase_path = resolve_path(str(dataset.array_phase_report), manifest_dir)
        source_type = str(dataset.source_type).strip()
        sequence_orientation = str(getattr(dataset, "sequence_orientation", "as_is") or "as_is").strip().lower()
        if sequence_orientation not in {"as_is", "feature", "genomic_plus"}:
            raise ValueError(
                f"{species}: sequence_orientation must be as_is, feature, or genomic_plus; "
                f"got {sequence_orientation!r}"
            )

        eprint(f"[{species}] reading {copy_path.name}")
        copies = load_copy_table(copy_path)
        initial_n = len(copies)
        copies = copies.loc[copies["phase_status"].str.upper().eq("PASS")].copy()
        after_pass = len(copies)
        if not args.include_non_full_length:
            copies = copies.loc[copies["quality_tier"].eq("full_length")].copy()
        after_quality = len(copies)

        # PTI-SAT113 contains known duplicate coordinates.  Coordinate-level
        # deduplication is done for every species to prevent pseudoreplication.
        coord_cols = ["chrom", "start", "end", "strand"]
        duplicate_mask = copies.duplicated(subset=coord_cols, keep="first")
        duplicate_coords = int(duplicate_mask.sum())
        copies = copies.loc[~duplicate_mask].copy()

        final_phase = select_phase_rows(final_phase_path, "copy_id", "final_local_consensus")
        array_phase = select_phase_rows(array_phase_path, "array_key", "global_consensus_phase")
        copies = copies.merge(final_phase, on="copy_id", how="left", validate="one_to_one")
        copies = copies.merge(array_phase, on="array_key", how="left", validate="many_to_one")
        missing_phase = copies[
            "final_local_consensus_orientation"
        ].eq("") | copies["global_consensus_phase_orientation"].eq("")
        missing_phase |= copies["final_local_consensus_orientation"].isna()
        missing_phase |= copies["global_consensus_phase_orientation"].isna()
        n_missing_phase = int(missing_phase.sum())
        copies = copies.loc[~missing_phase].copy()

        source_sequences, source_stats = load_source_sequences(
            source_type, source_path, copies, args.samtools
        )
        kept = 0
        missing_sequence = 0
        length_filtered = 0
        n_filtered = 0
        invalid_filtered = 0
        for row in copies.itertuples(index=False):
            raw_seq = source_sequences.get(str(row.copy_id))
            if not raw_seq:
                missing_sequence += 1
                continue
            sequence = clean_sequence(raw_seq)
            if sequence_orientation == "feature" and str(row.strand) == "-":
                sequence = revcomp(sequence)
            # genomic_plus and as_is both leave the source unchanged.  The
            # distinction is retained in metadata for reproducibility.
            sequence = transform_sequence(
                sequence,
                getattr(row, "final_local_consensus_orientation"),
                getattr(row, "final_local_consensus_rotation_left_bp"),
            )
            sequence = transform_sequence(
                sequence,
                getattr(row, "global_consensus_phase_orientation"),
                getattr(row, "global_consensus_phase_rotation_left_bp"),
            )
            if not sequence or not re.fullmatch(r"[ACGTRYSWKMBDHVN]+", sequence):
                invalid_filtered += 1
                continue
            local_length = getattr(row, "local_consensus_length", None)
            if local_length is None or pd.isna(local_length) or int(local_length) <= 0:
                local_length = len(sequence)
            ratio = len(sequence) / int(local_length)
            if ratio < args.min_length_ratio or ratio > args.max_length_ratio:
                length_filtered += 1
                continue
            n_fraction = sequence.count("N") / len(sequence)
            if n_fraction > args.max_n_fraction:
                n_filtered += 1
                continue
            record: Dict[str, object] = {
                "species": species,
                "copy_id": str(row.copy_id),
                "chrom": str(row.chrom),
                "array_key": str(row.array_key),
                "region_key": str(getattr(row, "region_key", "")),
                "start": int(row.start),
                "end": int(row.end),
                "strand": str(row.strand),
                "raw_length": int(getattr(row, "raw_length", len(raw_seq))),
                "local_consensus_length": int(local_length),
                "phased_length": len(sequence),
                "n_fraction": n_fraction,
                "source_type": source_type,
                "sequence_orientation": sequence_orientation,
                "sequence": sequence,
            }
            all_records.append(record)
            kept += 1

        qc_base = {
            "species": species,
            "input_copy_rows": initial_n,
            "pass_rows": after_pass,
            "quality_filtered_rows": after_quality,
            "duplicate_coordinate_rows_removed": duplicate_coords,
            "missing_phase_rows": n_missing_phase,
            "source_sequences_mapped": len(source_sequences),
            "missing_source_sequence_rows": missing_sequence,
            "length_filtered_rows": length_filtered,
            "n_filtered_rows": n_filtered,
            "invalid_sequence_rows": invalid_filtered,
            "final_kept_copies": kept,
        }
        qc_base.update({f"source_{key}": value for key, value in source_stats.items()})
        qc_rows.append(qc_base)
        eprint(
            f"[{species}] kept {kept:,}; removed duplicate coordinates={duplicate_coords:,}, "
            f"missing sequences={missing_sequence:,}"
        )

    if not all_records:
        raise RuntimeError("No monomer copies passed preparation; inspect preparation_qc.tsv and manifest paths.")

    records = pd.DataFrame(all_records)
    records = records.sort_values(["species", "chrom", "array_key", "start", "end", "copy_id"]).reset_index(drop=True)
    records.insert(0, "seq_id", [f"COPY{i:08d}" for i in range(1, len(records) + 1)])

    sequence_counts = Counter(records["sequence"])
    sorted_sequences = sorted(sequence_counts, key=lambda seq: (-sequence_counts[seq], seq))
    sequence_to_haplotype = {seq: f"HAP{i:07d}" for i, seq in enumerate(sorted_sequences, 1)}
    records["haplotype_id"] = records["sequence"].map(sequence_to_haplotype)

    metadata_columns = [column for column in records.columns if column != "sequence"]
    records[metadata_columns].to_csv(args.outdir / "copy_metadata.tsv", sep="\t", index=False)
    records[["seq_id", "haplotype_id"]].to_csv(
        args.outdir / "copy_haplotype.tsv", sep="\t", index=False
    )
    write_fasta(
        ((row.seq_id, row.sequence) for row in records.itertuples(index=False)),
        args.outdir / "all_copies.phased.fa",
    )

    hap_rows: List[Dict[str, object]] = []
    for sequence in sorted_sequences:
        haplotype_id = sequence_to_haplotype[sequence]
        subset = records.loc[records["haplotype_id"] == haplotype_id]
        hap_rows.append(
            {
                "haplotype_id": haplotype_id,
                "n_copies": len(subset),
                "n_species": subset["species"].nunique(),
                "n_arrays": subset[["species", "array_key"]].drop_duplicates().shape[0],
                "species_list": ",".join(sorted(subset["species"].unique())),
                "length": len(sequence),
                "n_fraction": sequence.count("N") / len(sequence),
                "sequence": sequence,
            }
        )
    haplotypes = pd.DataFrame(hap_rows)
    haplotypes.to_csv(args.outdir / "exact_haplotype_summary.tsv", sep="\t", index=False)
    write_fasta(
        (
            (f"{row.haplotype_id};size={row.n_copies};", row.sequence)
            for row in haplotypes.itertuples(index=False)
        ),
        args.outdir / "exact_haplotypes.fa",
    )
    pd.DataFrame(qc_rows).fillna(0).to_csv(
        args.outdir / "preparation_qc.tsv", sep="\t", index=False
    )

    summary = pd.DataFrame(
        [
            {"metric": "kept_copies", "value": len(records)},
            {"metric": "exact_haplotypes", "value": len(haplotypes)},
            {"metric": "species", "value": records["species"].nunique()},
            {
                "metric": "arrays",
                "value": records[["species", "array_key"]].drop_duplicates().shape[0],
            },
            {
                "metric": "copies_in_cross_array_exact_haplotypes",
                "value": int(
                    haplotypes.loc[haplotypes["n_arrays"] >= 2, "n_copies"].sum()
                ),
            },
            {
                "metric": "copies_in_cross_species_exact_haplotypes",
                "value": int(
                    haplotypes.loc[haplotypes["n_species"] >= 2, "n_copies"].sum()
                ),
            },
        ]
    )
    summary.to_csv(args.outdir / "preparation_summary.tsv", sep="\t", index=False)
    eprint(f"Done: {len(records):,} copies and {len(haplotypes):,} exact haplotypes")


if __name__ == "__main__":
    main()
