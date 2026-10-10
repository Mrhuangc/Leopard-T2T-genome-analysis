#!/usr/bin/env python3
"""Calculate all-versus-all nucleotide identity for repeat-unit FASTA sequences.

The script uses minimap2 for pairwise alignment and reports both alignment identity
and a coverage-penalized whole-unit identity.  It does not cluster the units or
infer a phylogeny.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import math
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, Sequence


@dataclass(frozen=True)
class FastaRecord:
    identifier: str
    description: str
    sequence: str


@dataclass(frozen=True)
class PafHit:
    query: str
    query_length: int
    query_start: int
    query_end: int
    strand: str
    target: str
    target_length: int
    target_start: int
    target_end: int
    matching_bases: int
    alignment_block_length: int
    mapq: int


@dataclass(frozen=True)
class PairRow:
    unit_a: str
    unit_b: str
    length_a: int
    length_b: int
    aligned_identity_pct: Optional[float]
    coverage_a_pct: Optional[float]
    coverage_b_pct: Optional[float]
    min_coverage_pct: Optional[float]
    whole_unit_identity_pct: Optional[float]
    matching_bases: Optional[int]
    alignment_block_length: Optional[int]
    strand: str
    mapq: Optional[int]
    status: str


VALID_IUPAC = set("ACGTRYSWKMBDHVN")


def read_fasta(path: Path) -> list[FastaRecord]:
    records: list[FastaRecord] = []
    header: Optional[str] = None
    chunks: list[str] = []

    def finish_record() -> None:
        if header is None:
            return
        identifier = header.split()[0]
        sequence = "".join(chunks).replace(" ", "").upper()
        if not identifier:
            raise ValueError("empty FASTA identifier")
        if not sequence:
            raise ValueError(f"empty FASTA sequence: {identifier}")
        invalid = sorted(set(sequence) - VALID_IUPAC)
        if invalid:
            raise ValueError(
                f"invalid FASTA character(s) in {identifier}: {''.join(invalid)}"
            )
        records.append(FastaRecord(identifier, header, sequence))

    with path.open() as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                finish_record()
                header = line[1:].strip()
                chunks = []
            elif header is None:
                raise ValueError(
                    f"sequence before first FASTA header at line {line_number}"
                )
            else:
                chunks.append(line)
    finish_record()

    if len(records) < 2:
        raise ValueError("FASTA must contain at least two sequences")
    identifiers = [record.identifier for record in records]
    duplicates = sorted({name for name in identifiers if identifiers.count(name) > 1})
    if duplicates:
        raise ValueError("duplicate FASTA ID(s): " + ", ".join(duplicates))
    return records


def parse_paf_line(line: str) -> PafHit:
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 12:
        raise ValueError(f"invalid PAF line: expected >=12 fields, got {len(fields)}")
    try:
        hit = PafHit(
            query=fields[0],
            query_length=int(fields[1]),
            query_start=int(fields[2]),
            query_end=int(fields[3]),
            strand=fields[4],
            target=fields[5],
            target_length=int(fields[6]),
            target_start=int(fields[7]),
            target_end=int(fields[8]),
            matching_bases=int(fields[9]),
            alignment_block_length=int(fields[10]),
            mapq=int(fields[11]),
        )
    except ValueError as exc:
        raise ValueError(f"invalid numeric value in PAF line: {line.rstrip()}") from exc
    if hit.strand not in {"+", "-"}:
        raise ValueError(f"invalid PAF strand: {hit.strand}")
    if hit.alignment_block_length <= 0:
        raise ValueError("PAF alignment block length must be positive")
    return hit


def hit_to_pair(hit: PafHit, unit_a: str, unit_b: str) -> PairRow:
    if (hit.query, hit.target) == (unit_a, unit_b):
        length_a, length_b = hit.query_length, hit.target_length
        span_a = hit.query_end - hit.query_start
        span_b = hit.target_end - hit.target_start
    elif (hit.query, hit.target) == (unit_b, unit_a):
        length_a, length_b = hit.target_length, hit.query_length
        span_a = hit.target_end - hit.target_start
        span_b = hit.query_end - hit.query_start
    else:
        raise ValueError("PAF hit does not match requested unit pair")

    coverage_a = 100.0 * span_a / length_a
    coverage_b = 100.0 * span_b / length_b
    return PairRow(
        unit_a=unit_a,
        unit_b=unit_b,
        length_a=length_a,
        length_b=length_b,
        aligned_identity_pct=100.0
        * hit.matching_bases
        / hit.alignment_block_length,
        coverage_a_pct=coverage_a,
        coverage_b_pct=coverage_b,
        min_coverage_pct=min(coverage_a, coverage_b),
        whole_unit_identity_pct=100.0
        * hit.matching_bases
        / max(length_a, length_b),
        matching_bases=hit.matching_bases,
        alignment_block_length=hit.alignment_block_length,
        strand=hit.strand,
        mapq=hit.mapq,
        status="",
    )


def _row_rank(row: PairRow) -> tuple[float, float, float, int, int]:
    return (
        row.whole_unit_identity_pct or 0.0,
        row.min_coverage_pct or 0.0,
        row.aligned_identity_pct or 0.0,
        row.mapq or 0,
        row.alignment_block_length or 0,
    )


def build_pairwise_rows(
    identifiers: Sequence[str],
    lengths: dict[str, int],
    paf_lines: Iterable[str],
    min_coverage_pct: float,
) -> list[PairRow]:
    order = {name: index for index, name in enumerate(identifiers)}
    best: dict[tuple[str, str], PairRow] = {}
    for line in paf_lines:
        if not line.strip() or line.startswith("#"):
            continue
        hit = parse_paf_line(line)
        if hit.query == hit.target:
            continue
        if hit.query not in order or hit.target not in order:
            raise ValueError(
                f"PAF contains an unknown sequence: {hit.query} or {hit.target}"
            )
        if order[hit.query] < order[hit.target]:
            key = (hit.query, hit.target)
        else:
            key = (hit.target, hit.query)
        row = hit_to_pair(hit, *key)
        if key not in best or _row_rank(row) > _row_rank(best[key]):
            best[key] = row

    rows: list[PairRow] = []
    for unit_a, unit_b in itertools.combinations(identifiers, 2):
        key = (unit_a, unit_b)
        if key not in best:
            rows.append(
                PairRow(
                    unit_a,
                    unit_b,
                    lengths[unit_a],
                    lengths[unit_b],
                    None,
                    None,
                    None,
                    None,
                    None,
                    None,
                    None,
                    ".",
                    None,
                    "NO_ALIGNMENT",
                )
            )
            continue
        row = best[key]
        status = (
            "PASS"
            if (row.min_coverage_pct or 0.0) >= min_coverage_pct
            else "LOW_COVERAGE"
        )
        rows.append(PairRow(**{**row.__dict__, "status": status}))
    return rows


def run_minimap2(
    executable: str,
    preset: str,
    threads: int,
    fasta: Path,
    sequence_count: int,
    paf_path: Path,
) -> None:
    resolved = shutil.which(executable)
    if resolved is None:
        candidate = Path(executable)
        if not candidate.is_file():
            raise RuntimeError(
                f"minimap2 executable not found: {executable}. Install minimap2 or "
                "provide its path with --minimap2."
            )
        resolved = str(candidate)
    command = [
        resolved,
        "-x",
        preset,
        "-c",
        "--eqx",
        "--secondary=yes",
        "-N",
        str(max(sequence_count * 4, 50)),
        "-p",
        "0",
        "-f",
        "0",
        "-t",
        str(threads),
        str(fasta),
        str(fasta),
    ]
    with paf_path.open("w") as output:
        result = subprocess.run(command, stdout=output, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"minimap2 failed with exit code {result.returncode}:\n{result.stderr.strip()}"
        )


def _format_value(value: Optional[float | int]) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float):
        return f"{value:.6f}"
    return str(value)


def write_pairwise_table(path: Path, rows: Sequence[PairRow]) -> None:
    columns = list(PairRow.__dataclass_fields__)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        for row in rows:
            writer.writerow([_format_value(getattr(row, column)) for column in columns])


def make_matrix(
    identifiers: Sequence[str], rows: Sequence[PairRow], attribute: str
) -> list[list[float]]:
    index = {name: i for i, name in enumerate(identifiers)}
    matrix = [[math.nan for _ in identifiers] for _ in identifiers]
    for i in range(len(identifiers)):
        matrix[i][i] = 100.0
    for row in rows:
        i, j = index[row.unit_a], index[row.unit_b]
        value = getattr(row, attribute)
        if value is not None:
            matrix[i][j] = matrix[j][i] = float(value)
    return matrix


def write_matrix(path: Path, identifiers: Sequence[str], matrix: Sequence[Sequence[float]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["unit_id", *identifiers])
        for identifier, values in zip(identifiers, matrix):
            writer.writerow(
                [identifier]
                + ["NA" if math.isnan(value) else f"{value:.6f}" for value in values]
            )


def plot_heatmap(
    pdf_path: Path,
    png_path: Path,
    identifiers: Sequence[str],
    matrix: Sequence[Sequence[float]],
    dpi: int,
) -> None:
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError as exc:
        raise RuntimeError(
            "plotting requires matplotlib and numpy; install them or use --no-plot"
        ) from exc

    data = np.asarray(matrix, dtype=float)
    masked = np.ma.masked_invalid(data)
    finite = data[np.isfinite(data)]
    vmin = max(0.0, math.floor(float(finite.min()))) if finite.size else 0.0
    size = max(8.0, min(18.0, 0.28 * len(identifiers) + 3.0))
    fig, axis = plt.subplots(figsize=(size, size))
    colour_map = plt.get_cmap("viridis").copy()
    colour_map.set_bad("#d9d9d9")
    image = axis.imshow(masked, cmap=colour_map, vmin=vmin, vmax=100.0)
    ticks = range(len(identifiers))
    font_size = max(4.0, min(8.0, 180.0 / max(len(identifiers), 1)))
    axis.set_xticks(ticks, labels=identifiers, rotation=90, fontsize=font_size)
    axis.set_yticks(ticks, labels=identifiers, fontsize=font_size)
    axis.set_xlabel("Unit")
    axis.set_ylabel("Unit")
    axis.set_title("Pairwise whole-unit nucleotide identity (%)")
    colour_bar = fig.colorbar(image, ax=axis, fraction=0.046, pad=0.04)
    colour_bar.set_label("Whole-unit identity (%)")
    fig.tight_layout()
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate all-versus-all unit-sequence identity with minimap2. "
            "No clustering or phylogeny is performed."
        )
    )
    parser.add_argument("--fasta", required=True, type=Path, help="multi-FASTA of units")
    parser.add_argument("--outdir", type=Path, default=Path("."))
    parser.add_argument("--prefix", help="output prefix (default: FASTA stem)")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument(
        "--preset",
        choices=("asm5", "asm10", "asm20"),
        default="asm5",
        help="minimap2 preset; asm5 is appropriate for highly similar copies",
    )
    parser.add_argument("--minimap2", default="minimap2", help="minimap2 executable")
    parser.add_argument(
        "--min-coverage",
        type=float,
        default=80.0,
        help="minimum coverage of both units for PASS status (percent; default 80)",
    )
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--no-plot", action="store_true")
    parser.add_argument(
        "--keep-paf", action="store_true", help="retain raw minimap2 PAF output"
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if args.threads < 1:
        raise ValueError("--threads must be >= 1")
    if not 0.0 <= args.min_coverage <= 100.0:
        raise ValueError("--min-coverage must be between 0 and 100")
    if args.dpi < 72:
        raise ValueError("--dpi must be >= 72")
    if not args.fasta.is_file():
        raise FileNotFoundError(f"FASTA not found: {args.fasta}")

    records = read_fasta(args.fasta)
    identifiers = [record.identifier for record in records]
    lengths = {record.identifier: len(record.sequence) for record in records}
    prefix = args.prefix or args.fasta.stem
    args.outdir.mkdir(parents=True, exist_ok=True)

    final_paf = args.outdir / f"{prefix}.minimap2.paf"
    with tempfile.TemporaryDirectory(prefix="unit_identity_", dir=args.outdir) as temp_dir:
        temporary_paf = Path(temp_dir) / "alignments.paf"
        run_minimap2(
            args.minimap2,
            args.preset,
            args.threads,
            args.fasta,
            len(records),
            temporary_paf,
        )
        with temporary_paf.open() as handle:
            rows = build_pairwise_rows(
                identifiers, lengths, handle, min_coverage_pct=args.min_coverage
            )
        if args.keep_paf:
            shutil.copyfile(temporary_paf, final_paf)

    pairwise_path = args.outdir / f"{prefix}.pairwise_identity.tsv"
    aligned_matrix_path = args.outdir / f"{prefix}.aligned_identity.matrix.tsv"
    whole_matrix_path = args.outdir / f"{prefix}.whole_unit_identity.matrix.tsv"
    write_pairwise_table(pairwise_path, rows)
    aligned_matrix = make_matrix(identifiers, rows, "aligned_identity_pct")
    whole_matrix = make_matrix(identifiers, rows, "whole_unit_identity_pct")
    write_matrix(aligned_matrix_path, identifiers, aligned_matrix)
    write_matrix(whole_matrix_path, identifiers, whole_matrix)

    if not args.no_plot:
        plot_heatmap(
            args.outdir / f"{prefix}.identity_heatmap.pdf",
            args.outdir / f"{prefix}.identity_heatmap.png",
            identifiers,
            whole_matrix,
            args.dpi,
        )

    expected = len(records) * (len(records) - 1) // 2
    pass_count = sum(row.status == "PASS" for row in rows)
    low_count = sum(row.status == "LOW_COVERAGE" for row in rows)
    missing_count = sum(row.status == "NO_ALIGNMENT" for row in rows)
    print(
        f"Sequences: {len(records)}; pairs: {expected}; PASS: {pass_count}; "
        f"LOW_COVERAGE: {low_count}; NO_ALIGNMENT: {missing_count}",
        file=sys.stderr,
    )
    print(pairwise_path)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
