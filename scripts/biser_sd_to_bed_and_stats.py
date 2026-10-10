#!/usr/bin/env python3
"""
Convert BISER BEDPE calls into BED/TSV files and summarize SD length/distribution.

BISER coordinates are already 0-based, half-open, therefore they can be written
to BED directly without coordinate conversion.

Typical use:
  python biser_sd_to_bed_and_stats.py \
      -i FCA.biser \
      -o FCA.BISER.SD1kb \
      --min-len 1000 \
      --length-field max_len \
      --fai FCA.genome.softmask.fai \
      --elem FCA.biser.elem.txt

Outputs:
  <prefix>.pairs.bedpe                 BEDPE: one SD relationship per line
  <prefix>.pairs.annotated.tsv         BEDPE plus BISER metrics/annotations
  <prefix>.mates.bed                   BED6: unique genomic SD-copy intervals
  <prefix>.mates.annotated.tsv         metadata for unique SD-copy intervals
  <prefix>.merged_loci.bed             union of mate intervals on each chromosome
  <prefix>.summary.tsv                 overall summary
  <prefix>.length_distribution.tsv     number and reported length by length bin
  <prefix>.chrom_distribution.tsv      chromosome-level distribution
  <prefix>.orientation_summary.tsv     intra/inter-chromosomal and orientation classes
  <prefix>.elem_*.{bed,tsv}            optional elementary/core SD decomposition outputs

Notes:
  * Default filtering uses BISER's `max_len` field (the longer mate length).
    This avoids using start-end spans that may contain large internal gaps in
    BISER CIGAR strings.
  * `merged_loci.bed` reports genomic span coverage, not ungapped aligned bases.
    Interpret it as an SD-span annotation layer.
  * `--max-total-error` is optional. Do not apply it merely because an interval
    lies in an SD; choose a cutoff according to your biological definition.
"""

from __future__ import annotations
import argparse
import csv
import re
import sys
from pathlib import Path
from collections import defaultdict, Counter
from statistics import median

PAIR_COLS = [
    "chr1", "start1", "end1", "chr2", "start2", "end2",
    "reference", "total_error_pct", "strand1", "strand2",
    "max_len", "aln_len", "cigar", "optional"
]

LENGTH_BINS = [
    (1_000, 2_000, "1-2kb"),
    (2_000, 5_000, "2-5kb"),
    (5_000, 10_000, "5-10kb"),
    (10_000, 50_000, "10-50kb"),
    (50_000, 100_000, "50-100kb"),
    (100_000, 500_000, "100-500kb"),
    (500_000, None, ">=500kb"),
]

def natural_key(s: str):
    """Natural-ish chromosome ordering, e.g. ChrA2 before ChrA10."""
    return [int(x) if x.isdigit() else x.lower() for x in re.split(r"(\d+)", s)]

def parse_optional(text: str, key: str):
    m = re.search(rf"(?:^|;){re.escape(key)}=([0-9.]+)", text or "")
    return float(m.group(1)) if m else None

def row_to_pair(fields, line_no: int):
    if len(fields) < 14:
        raise ValueError(f"line {line_no}: expected >=14 columns, got {len(fields)}")
    d = dict(zip(PAIR_COLS, fields[:14]))
    for key in ("start1", "end1", "start2", "end2", "max_len", "aln_len"):
        d[key] = int(d[key])
    d["total_error_pct"] = float(d["total_error_pct"])
    d["mismatch_pct"] = parse_optional(d["optional"], "X")
    d["gap_pct"] = parse_optional(d["optional"], "ID")
    d["len1_span"] = d["end1"] - d["start1"]
    d["len2_span"] = d["end2"] - d["start2"]
    d["line_no"] = line_no
    return d

def metric_value(d, name: str):
    if name == "max_len":
        return d["max_len"]
    if name == "aln_len":
        return d["aln_len"]
    if name == "min_mate_span":
        return min(d["len1_span"], d["len2_span"])
    if name == "max_mate_span":
        return max(d["len1_span"], d["len2_span"])
    raise ValueError(f"Unsupported length field: {name}")

def mate_tuple(d, side: int):
    return (d[f"chr{side}"], d[f"start{side}"], d[f"end{side}"], d[f"strand{side}"])

def canonical_pair_key(d):
    a, b = mate_tuple(d, 1), mate_tuple(d, 2)
    return tuple(sorted((a, b)))

def choose_better(old, new):
    """Retain the lower-error representative of exact reciprocal calls."""
    old_key = (old["total_error_pct"], -old["max_len"], old["line_no"])
    new_key = (new["total_error_pct"], -new["max_len"], new["line_no"])
    return new if new_key < old_key else old

def write_tsv(path: Path, header, rows):
    with path.open("w", newline="") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)

def out_path(prefix: Path, tail: str) -> Path:
    """Append a filename tail without treating dots in the prefix as suffixes."""
    return Path(str(prefix) + tail)

def read_fai(fai_path: Path):
    sizes = {}
    with fai_path.open() as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            sizes[f[0]] = int(f[1])
    return sizes

def merge_intervals(intervals):
    """Input iterable of (start, end); merge overlapping or touching intervals."""
    intervals = sorted(intervals)
    if not intervals:
        return []
    merged = [[intervals[0][0], intervals[0][1]]]
    for s, e in intervals[1:]:
        if s <= merged[-1][1]:
            if e > merged[-1][1]:
                merged[-1][1] = e
        else:
            merged.append([s, e])
    return [(s, e) for s, e in merged]

def length_bin_name(x):
    for lo, hi, label in LENGTH_BINS:
        if x >= lo and (hi is None or x < hi):
            return label
    return "<1kb"

def read_biser(path: Path, min_len: int, length_field: str,
               max_total_error, keep_reciprocal: bool):
    raw_n = 0
    pass_length_n = 0
    pass_error_n = 0
    pair_map = {} if not keep_reciprocal else None
    kept = []
    with path.open() as fh:
        for line_no, line in enumerate(fh, 1):
            if not line.strip() or line.startswith("#"):
                continue
            raw_n += 1
            d = row_to_pair(line.rstrip("\n").split("\t"), line_no)
            d["length_metric"] = metric_value(d, length_field)
            if d["length_metric"] < min_len:
                continue
            pass_length_n += 1
            if max_total_error is not None and d["total_error_pct"] > max_total_error:
                continue
            pass_error_n += 1
            if keep_reciprocal:
                kept.append(d)
            else:
                key = canonical_pair_key(d)
                if key not in pair_map:
                    pair_map[key] = d
                else:
                    pair_map[key] = choose_better(pair_map[key], d)

    if pair_map is not None:
        kept = list(pair_map.values())
    kept.sort(key=lambda d: (natural_key(d["chr1"]), d["start1"], natural_key(d["chr2"]), d["start2"]))
    for i, d in enumerate(kept, 1):
        d["sd_id"] = f"SD{i:07d}"
    return kept, {
        "raw_biser_records": raw_n,
        "pass_length_filter_before_reciprocal_dedup": pass_length_n,
        "pass_error_filter_before_reciprocal_dedup": pass_error_n,
        "retained_sd_pair_calls": len(kept),
    }

def write_pair_outputs(prefix: Path, pairs):
    # BEDPE without header
    bedpe = out_path(prefix, ".pairs.bedpe")
    annot = out_path(prefix, ".pairs.annotated.tsv")
    header = [
        "sd_id", "chr1", "start1", "end1", "chr2", "start2", "end2",
        "reference", "total_error_pct", "mismatch_pct", "gap_pct",
        "strand1", "strand2", "max_len", "aln_len", "length_metric_bp",
        "intra_or_interchromosomal", "direct_or_inverted", "cigar", "optional"
    ]
    annot_rows = []
    with bedpe.open("w") as out:
        for d in pairs:
            out.write("\t".join(map(str, [
                d["chr1"], d["start1"], d["end1"], d["chr2"], d["start2"], d["end2"],
                d["sd_id"], d["max_len"], d["strand1"], d["strand2"]
            ])) + "\n")
            chrom_class = "intra" if d["chr1"] == d["chr2"] else "inter"
            ori_class = "direct" if d["strand1"] == d["strand2"] else "inverted"
            annot_rows.append([
                d["sd_id"], d["chr1"], d["start1"], d["end1"], d["chr2"], d["start2"], d["end2"],
                d["reference"], d["total_error_pct"],
                "" if d["mismatch_pct"] is None else d["mismatch_pct"],
                "" if d["gap_pct"] is None else d["gap_pct"],
                d["strand1"], d["strand2"], d["max_len"], d["aln_len"], d["length_metric"],
                chrom_class, ori_class, d["cigar"], d["optional"]
            ])
    write_tsv(annot, header, annot_rows)
    return bedpe, annot

def make_mates(pairs):
    """Collapse exact same loci while retaining all pair memberships."""
    loci = {}
    for d in pairs:
        for side, label in ((1, "mate1"), (2, "mate2")):
            k = mate_tuple(d, side)
            entry = loci.setdefault(k, {
                "chrom": k[0], "start": k[1], "end": k[2], "strand": k[3],
                "sd_ids": set(), "mate_labels": set(), "max_len_values": [],
                "errors": []
            })
            entry["sd_ids"].add(d["sd_id"])
            entry["mate_labels"].add(label)
            entry["max_len_values"].append(d["max_len"])
            entry["errors"].append(d["total_error_pct"])
    entries = list(loci.values())
    entries.sort(key=lambda x: (natural_key(x["chrom"]), x["start"], x["end"], x["strand"]))
    for i, e in enumerate(entries, 1):
        e["copy_id"] = f"SDcopy{i:07d}"
    return entries

def write_mate_outputs(prefix: Path, mates):
    bed = out_path(prefix, ".mates.bed")
    annot = out_path(prefix, ".mates.annotated.tsv")
    with bed.open("w") as out:
        for e in mates:
            out.write(f'{e["chrom"]}\t{e["start"]}\t{e["end"]}\t{e["copy_id"]}\t0\t{e["strand"]}\n')
    header = [
        "copy_id", "chrom", "start", "end", "span_bp", "strand",
        "n_partner_sd_pairs", "sd_ids", "max_pair_length_bp",
        "median_pair_total_error_pct", "mate_roles"
    ]
    rows = []
    for e in mates:
        rows.append([
            e["copy_id"], e["chrom"], e["start"], e["end"], e["end"]-e["start"], e["strand"],
            len(e["sd_ids"]), ",".join(sorted(e["sd_ids"])), max(e["max_len_values"]),
            round(median(e["errors"]), 4), ",".join(sorted(e["mate_labels"]))
        ])
    write_tsv(annot, header, rows)
    return bed, annot

def write_merged_loci(prefix: Path, mates):
    by_chr = defaultdict(list)
    for e in mates:
        by_chr[e["chrom"]].append((e["start"], e["end"]))
    merged_by_chr = {chrom: merge_intervals(iv) for chrom, iv in by_chr.items()}
    bed = out_path(prefix, ".merged_loci.bed")
    with bed.open("w") as out:
        i = 1
        for chrom in sorted(merged_by_chr, key=natural_key):
            for s, e in merged_by_chr[chrom]:
                out.write(f"{chrom}\t{s}\t{e}\tSDlocus{i:07d}\n")
                i += 1
    return bed, merged_by_chr

def write_stats(prefix: Path, pairs, mates, merged_by_chr, summary_base, args, fai_sizes=None):
    # Length bins: report pair calls, using selected length metric.
    by_bin = defaultdict(list)
    for d in pairs:
        by_bin[length_bin_name(d["length_metric"])].append(d)
    length_rows = []
    for lo, hi, label in LENGTH_BINS:
        xs = by_bin.get(label, [])
        lengths = [d["length_metric"] for d in xs]
        length_rows.append([
            label, lo, "" if hi is None else hi, len(xs), sum(lengths),
            "" if not lengths else min(lengths), "" if not lengths else round(median(lengths), 2),
            "" if not lengths else max(lengths)
        ])
    write_tsv(out_path(prefix, ".length_distribution.tsv"),
              ["length_bin", "lower_inclusive_bp", "upper_exclusive_bp", "n_sd_pairs",
               "sum_reported_sd_length_bp", "min_bp", "median_bp", "max_bp"],
              length_rows)

    # Orientation classes at pair level.
    orient = Counter()
    for d in pairs:
        chrom_class = "intra" if d["chr1"] == d["chr2"] else "inter"
        ori_class = "direct" if d["strand1"] == d["strand2"] else "inverted"
        orient[(chrom_class, ori_class)] += 1
    orient_rows = []
    for chrom_class in ("intra", "inter"):
        for ori_class in ("direct", "inverted"):
            n = orient[(chrom_class, ori_class)]
            orient_rows.append([chrom_class, ori_class, n,
                                round(100*n/len(pairs), 4) if pairs else 0])
    write_tsv(out_path(prefix, ".orientation_summary.tsv"),
              ["chromosomal_class", "orientation_class", "n_sd_pairs", "fraction_pct"],
              orient_rows)

    # Chromosome distribution from unique mate loci and their merged span.
    mate_by_chr = defaultdict(list)
    endpoints_by_chr = Counter()
    for e in mates:
        mate_by_chr[e["chrom"]].append(e)
    for d in pairs:
        endpoints_by_chr[d["chr1"]] += 1
        endpoints_by_chr[d["chr2"]] += 1

    all_chroms = set(mate_by_chr) | set(merged_by_chr)
    if fai_sizes:
        all_chroms |= set(fai_sizes)
    chrom_rows = []
    for chrom in sorted(all_chroms, key=natural_key):
        merged = merged_by_chr.get(chrom, [])
        merged_bp = sum(e-s for s, e in merged)
        n_mates = len(mate_by_chr.get(chrom, []))
        chrom_len = fai_sizes.get(chrom) if fai_sizes else None
        chrom_rows.append([
            chrom,
            endpoints_by_chr.get(chrom, 0),
            n_mates,
            len(merged),
            merged_bp,
            "" if chrom_len is None else chrom_len,
            "" if chrom_len is None else round(merged_bp / chrom_len * 100, 6),
            "" if chrom_len is None else round(merged_bp / chrom_len * 1_000_000, 6),
        ])
    write_tsv(out_path(prefix, ".chrom_distribution.tsv"),
              ["chrom", "n_sd_pair_endpoints", "n_unique_sd_copy_loci", "n_merged_sd_loci",
               "merged_sd_span_bp", "chrom_length_bp", "merged_sd_coverage_pct",
               "merged_sd_span_bp_per_Mb"],
              chrom_rows)

    # Overall summary.
    total_merged_bp = sum(sum(e-s for s, e in ivs) for ivs in merged_by_chr.values())
    summary_rows = [
        ["input", "biser_file", args.input],
        ["filter", "min_sd_len_bp", args.min_len],
        ["filter", "length_field", args.length_field],
        ["filter", "max_total_error_pct", "" if args.max_total_error is None else args.max_total_error],
        ["filter", "reciprocal_handling", "kept" if args.keep_reciprocal else "exact reciprocal pairs deduplicated"],
    ]
    for k, v in summary_base.items():
        summary_rows.append(["count", k, v])
    summary_rows.extend([
        ["count", "unique_sd_copy_loci", len(mates)],
        ["length", "sum_reported_pair_lengths_bp", sum(d["length_metric"] for d in pairs)],
        ["length", "median_reported_pair_length_bp", round(median([d["length_metric"] for d in pairs]), 2) if pairs else 0],
        ["span", "merged_genomic_sd_span_bp", total_merged_bp],
        ["span", "merged_genomic_sd_span_Mb", round(total_merged_bp / 1_000_000, 6)],
    ])
    if fai_sizes:
        genome_bp = sum(fai_sizes.values())
        summary_rows.extend([
            ["reference", "genome_length_bp_from_fai", genome_bp],
            ["reference", "merged_sd_coverage_pct_of_genome", round(total_merged_bp / genome_bp * 100, 6)],
        ])
    write_tsv(out_path(prefix, ".summary.tsv"), ["section", "metric", "value"], summary_rows)

def process_elem(elem_path: Path, prefix: Path, min_len: int):
    """Optional: export BISER elementary/core decomposition. Does not alter SD pairs."""
    rows = []
    with elem_path.open() as fh:
        for line_no, line in enumerate(fh, 1):
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                raise ValueError(f"{elem_path}: line {line_no} has <8 columns")
            ref, chrom, start, end, elem_id, ln, score, strand = f[:8]
            extra = f[8] if len(f) > 8 else ""
            start, end, ln, score = int(start), int(end), int(ln), int(score)
            if ln >= min_len:
                rows.append((ref, chrom, start, end, elem_id, ln, score, strand, extra))
    rows.sort(key=lambda r: (natural_key(r[1]), r[2], r[3], r[4]))
    bed = out_path(prefix, ".elem_ge_minlen.bed")
    with bed.open("w") as out:
        for r in rows:
            # BED6: extra trailing columns stay in separate TSV
            out.write(f"{r[1]}\t{r[2]}\t{r[3]}\t{r[4]}\t0\t{r[7]}\n")
    write_tsv(out_path(prefix, ".elem_ge_minlen.annotated.tsv"),
              ["reference", "chrom", "start", "end", "element_id", "length_bp",
               "internal_score", "strand", "extra_flag"],
              rows)
    grouped = defaultdict(list)
    for r in rows:
        grouped[r[4]].append(r)
    summary_rows = []
    for elem_id, rs in sorted(grouped.items(), key=lambda kv: (-len(kv[1]), kv[0])):
        lengths = [r[5] for r in rs]
        summary_rows.append([
            elem_id, len(rs), len(set(r[1] for r in rs)), sum(lengths),
            min(lengths), round(median(lengths), 2), max(lengths),
            sum(1 for r in rs if r[8] == "CORE")
        ])
    write_tsv(out_path(prefix, ".elem_ge_minlen.by_id.tsv"),
              ["element_id", "n_instances", "n_chromosomes", "sum_length_bp",
               "min_length_bp", "median_length_bp", "max_length_bp", "n_extra_CORE_flag"],
              summary_rows)
    return len(rows), len(grouped)

def main():
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Convert BISER BEDPE to SD BED tracks and summary tables."
    )
    p.add_argument("-i", "--input", required=True, help="BISER main BEDPE-like output file.")
    p.add_argument("-o", "--output-prefix", required=True,
                   help="Prefix for all output files; parent directory is created automatically.")
    p.add_argument("--min-len", type=int, default=1000,
                   help="Minimum SD length for retention.")
    p.add_argument("--length-field", choices=["max_len", "aln_len", "min_mate_span", "max_mate_span"],
                   default="max_len",
                   help="Field used for the minimum-length filter and length-bin statistics.")
    p.add_argument("--max-total-error", type=float, default=None,
                   help="Optional maximum BISER total error percentage. Disabled by default.")
    p.add_argument("--keep-reciprocal", action="store_true",
                   help="Do not collapse exact reciprocal/redundant BEDPE calls.")
    p.add_argument("--fai", default=None,
                   help="Optional FASTA index (.fai) for chromosome-normalized SD distribution.")
    p.add_argument("--elem", default=None,
                   help="Optional BISER .elem output for elementary/core SD decomposition export.")
    args = p.parse_args()

    if args.min_len < 1:
        p.error("--min-len must be >= 1")
    if args.max_total_error is not None and not (0 <= args.max_total_error <= 100):
        p.error("--max-total-error must be between 0 and 100")

    prefix = Path(args.output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    input_path = Path(args.input)

    try:
        pairs, summary = read_biser(input_path, args.min_len, args.length_field,
                                    args.max_total_error, args.keep_reciprocal)
    except Exception as e:
        sys.exit(f"ERROR while reading {input_path}: {e}")

    fai_sizes = read_fai(Path(args.fai)) if args.fai else None
    write_pair_outputs(prefix, pairs)
    mates = make_mates(pairs)
    write_mate_outputs(prefix, mates)
    _, merged_by_chr = write_merged_loci(prefix, mates)
    write_stats(prefix, pairs, mates, merged_by_chr, summary, args, fai_sizes)

    if args.elem:
        n_elem, n_ids = process_elem(Path(args.elem), prefix, args.min_len)
        with out_path(prefix, ".summary.tsv").open("a") as out:
            out.write(f"elem\tretained_elementary_regions_ge_minlen\t{n_elem}\n")
            out.write(f"elem\tunique_element_ids_ge_minlen\t{n_ids}\n")

    print(f"Done. Retained {len(pairs):,} SD-pair calls and {len(mates):,} unique SD-copy loci.")
    print(f"Main BED for overlap annotation: {out_path(prefix, '.merged_loci.bed')}")
    print(f"Pairwise SD relationships:       {out_path(prefix, '.pairs.bedpe')}")

if __name__ == "__main__":
    main()
