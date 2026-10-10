#!/usr/bin/env python3
import argparse
import gzip
from collections import defaultdict


def smart_open(path):
    """Open plain text or gzip-compressed files automatically."""
    with open(path, "rb") as f:
        magic = f.read(2)

    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt")
    return open(path, "rt")


def read_sizes(path):
    """Read a two-column genome.sizes file."""
    sizes = {}
    with smart_open(path) as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue
            parts = line.rstrip().split()
            if len(parts) < 2:
                continue
            sizes[parts[0]] = int(parts[1])
    return sizes


def merge_intervals(intervals):
    intervals = sorted(intervals)
    merged = []

    for start, end in intervals:
        if end <= start:
            continue

        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)

    return merged


def parse_tags(parts):
    """Return PAF optional tags as a set, e.g. {'tp:A:P', ...}."""
    if len(parts) <= 12:
        return set()
    return set(parts[12:])


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Generate PUR BED from old genome vs T2T genome PAF. "
            "Covered regions can be defined by strict MAPQ/length filtering, "
            "plus optional rescue of long primary alignments with lower MAPQ."
        )
    )

    parser.add_argument("-i", "--paf", required=True, help="Input PAF file; plain text or gzip-compressed")
    parser.add_argument("-o", "--out", required=True, help="Output PUR BED")
    parser.add_argument(
        "-s", "--sizes",
        default=None,
        help=(
            "Optional T2T genome sizes file. Recommended. "
            "If provided, chromosomes with no retained alignments will still be output as PUR."
        )
    )

    parser.add_argument(
        "--min_mapq",
        type=int,
        default=20,
        help="Minimum MAPQ for the strict filter, default: 20"
    )
    parser.add_argument(
        "--min_len",
        type=int,
        default=1000,
        help="Minimum alignment length on T2T for the strict filter, default: 1000"
    )
    parser.add_argument(
        "--primary_only",
        action="store_true",
        help="Only use primary alignments with tp:A:P"
    )

    parser.add_argument(
        "--long_primary_min_len",
        type=int,
        default=1000000,
        help=(
            "Rescue low-MAPQ but long primary alignments if their T2T span is at least this length. "
            "Default: 1000000 bp. Set to 0 to disable this rescue rule."
        )
    )

    parser.add_argument(
        "--kept_paf",
        default=None,
        help=(
            "Optional output file listing retained PAF alignments with filter_reason "
            "(strict or long_primary_rescue). Useful for checking."
        )
    )

    args = parser.parse_args()

    # t2t_len stores target chromosome lengths.
    # If --sizes is provided, use it as the authoritative source.
    # Otherwise, collect target lengths from all PAF records before filtering.
    t2t_len = read_sizes(args.sizes) if args.sizes else {}
    covered = defaultdict(list)

    n_total = 0
    n_kept = 0
    n_strict = 0
    n_rescue = 0
    n_non_primary_skipped = 0

    kept_handle = open(args.kept_paf, "w") if args.kept_paf else None
    if kept_handle:
        kept_handle.write(
            "filter_reason\tqname\tqlen\tqstart\tqend\tstrand\t"
            "tname\ttlen\ttstart\ttend\tnmatch\talen\tmapq\n"
        )

    with smart_open(args.paf) as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue

            parts = line.rstrip().split("\t")
            if len(parts) < 12:
                continue

            n_total += 1

            # Standard PAF fields
            q_name = parts[0]
            q_len = int(parts[1])
            q_start = int(parts[2])
            q_end = int(parts[3])
            strand = parts[4]

            t_chr = parts[5]
            t_len = int(parts[6])
            t_start = int(parts[7])
            t_end = int(parts[8])

            n_match = int(parts[9])
            aln_block_len = int(parts[10])
            mapq = int(parts[11])

            tags = parse_tags(parts)
            is_primary = "tp:A:P" in tags

            # Collect target lengths from all PAF lines, not only retained lines.
            # This avoids losing a chromosome when all of its alignments are filtered.
            if t_chr not in t2t_len:
                t2t_len[t_chr] = t_len

            t_span = t_end - t_start
            if t_span <= 0:
                continue

            if args.primary_only and not is_primary:
                n_non_primary_skipped += 1
                continue

            strict_pass = (mapq >= args.min_mapq and t_span >= args.min_len)

            rescue_pass = (
                args.long_primary_min_len > 0
                and is_primary
                and t_span >= args.long_primary_min_len
            )

            if not (strict_pass or rescue_pass):
                continue

            if strict_pass:
                reason = "strict"
                n_strict += 1
            else:
                reason = "long_primary_rescue"
                n_rescue += 1

            n_kept += 1
            covered[t_chr].append((t_start, t_end))

            if kept_handle:
                kept_handle.write(
                    f"{reason}\t{q_name}\t{q_len}\t{q_start}\t{q_end}\t{strand}\t"
                    f"{t_chr}\t{t_len}\t{t_start}\t{t_end}\t"
                    f"{n_match}\t{aln_block_len}\t{mapq}\n"
                )

    if kept_handle:
        kept_handle.close()

    total_pur = 0

    with open(args.out, "w") as out:
        for chrom in sorted(t2t_len.keys()):
            chrom_len = t2t_len[chrom]
            merged = merge_intervals(covered[chrom])

            prev_end = 0

            for start, end in merged:
                if start > prev_end:
                    pur_start = prev_end
                    pur_end = start
                    pur_len = pur_end - pur_start
                    total_pur += pur_len
                    out.write(f"{chrom}\t{pur_start}\t{pur_end}\tPUR\t{pur_len}\n")

                prev_end = max(prev_end, end)

            if prev_end < chrom_len:
                pur_start = prev_end
                pur_end = chrom_len
                pur_len = pur_end - pur_start
                total_pur += pur_len
                out.write(f"{chrom}\t{pur_start}\t{pur_end}\tPUR\t{pur_len}\n")

    print(f"Total PAF records: {n_total}")
    print(f"Kept alignments: {n_kept}")
    print(f"  strict filter kept: {n_strict}")
    print(f"  long primary rescue kept: {n_rescue}")
    if args.primary_only:
        print(f"Non-primary skipped by --primary_only: {n_non_primary_skipped}")

    print(f"Total PUR length: {total_pur} bp")
    print(f"Total PUR length: {total_pur / 1_000_000:.6f} Mb")


if __name__ == "__main__":
    main()
