#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Split a multi-FASTA file into multiple FASTA files by sequence count.

Example:
    python split_fasta_by_records.py \
        -i PPA.centromere.fa \
        -o PPA.centromere.split \
        -p PPA_cen \
        -n 5
"""

import argparse
import os
from pathlib import Path


def read_fasta(path):
    """Yield (header, sequence) records from FASTA."""
    header = None
    seq_lines = []

    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_lines)
                header = line[1:].strip()
                seq_lines = []
            else:
                if header is None:
                    raise ValueError("FASTA sequence encountered before the first header.")
                seq_lines.append(line.strip())

    if header is not None:
        yield header, "".join(seq_lines)


def write_fasta(records, output_file, line_width=60):
    """Write FASTA records."""
    with open(output_file, "w") as out:
        for header, seq in records:
            out.write(f">{header}\n")
            for i in range(0, len(seq), line_width):
                out.write(seq[i:i + line_width] + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Split a multi-FASTA file into groups containing N sequences."
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input multi-FASTA file."
    )
    parser.add_argument(
        "-o", "--outdir",
        required=True,
        help="Output directory."
    )
    parser.add_argument(
        "-p", "--prefix",
        default="split",
        help="Output FASTA prefix. Default: split"
    )
    parser.add_argument(
        "-n", "--number",
        type=int,
        default=5,
        help="Number of sequences per output FASTA. Default: 5"
    )
    parser.add_argument(
        "--manifest",
        default=None,
        help="Optional output manifest TSV file. Default: <outdir>/<prefix>.manifest.tsv"
    )

    args = parser.parse_args()

    if args.number < 1:
        raise ValueError("--number must be >= 1")

    input_fasta = Path(args.input)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if not input_fasta.exists():
        raise FileNotFoundError(f"Input FASTA not found: {input_fasta}")

    records = list(read_fasta(input_fasta))

    if not records:
        raise ValueError(f"No FASTA records found in: {input_fasta}")

    manifest_file = (
        Path(args.manifest)
        if args.manifest
        else outdir / f"{args.prefix}.manifest.tsv"
    )

    total_seq = len(records)
    total_bp = sum(len(seq) for _, seq in records)
    n_files = (total_seq + args.number - 1) // args.number

    with open(manifest_file, "w") as manifest:
        manifest.write(
            "group\tfasta_file\tsequence_number\ttotal_length_bp\tsequence_ids\n"
        )

        for group_index, start in enumerate(
            range(0, total_seq, args.number), start=1
        ):
            batch = records[start:start + args.number]

            output_fasta = outdir / f"{args.prefix}.group{group_index:02d}.fa"
            write_fasta(batch, output_fasta)

            seq_ids = [header.split()[0] for header, _ in batch]
            batch_bp = sum(len(seq) for _, seq in batch)

            manifest.write(
                f"group{group_index:02d}\t"
                f"{output_fasta}\t"
                f"{len(batch)}\t"
                f"{batch_bp}\t"
                f"{','.join(seq_ids)}\n"
            )

            length_mb = batch_bp / 1e6
            status = "OK"
            if length_mb > 40:
                status = "WARNING: total length > 40 Mb"

            print(
                f"[INFO] group{group_index:02d}: "
                f"{len(batch)} sequences, {length_mb:.2f} Mb, "
                f"{output_fasta} [{status}]"
            )

    print()
    print(f"[DONE] Input sequences: {total_seq}")
    print(f"[DONE] Input total length: {total_bp / 1e6:.2f} Mb")
    print(f"[DONE] Output FASTA files: {n_files}")
    print(f"[DONE] Manifest: {manifest_file}")


if __name__ == "__main__":
    main()
