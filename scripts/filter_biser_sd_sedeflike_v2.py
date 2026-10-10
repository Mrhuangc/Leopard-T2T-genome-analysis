#!/usr/bin/env python3
"""Filter BISER SD calls with the same *decision logic* used by the CHM13 SEDEF workflow.

Important implementation notes
------------------------------
This script does not treat BISER's N/S CIGAR operations as alignment gaps.
BISER defines ``aln_len`` as the alignment span and reports mismatch (X) and
indel (ID) rates for that alignment; N/S are not included in those quantities.

The retained-call criteria mirror sedef_to_bed.py in Vollger et al.:
  1) aln_len >= 1,000 bp
  2) fracMatch >= 0.90 (estimated from BISER X and CIGAR M)
  3) fracMatchIndel >= 0.50 (estimated from BISER X and aln_len)
  4) satellite fraction <=0.70 for at least one mate by default.

The original SEDEF workflow filters Satellite coverage only on the input-side
mate before making a symmetric track. For orientation-invariant BISER calls,
``--satellite-mode either`` implements the intended behavior: retain the pair
when at least one mate is <= the satellite threshold. ``both`` is deliberately
more stringent and should not be used to reproduce the CHM13 workflow.

The result is a BISER-defined SD callset filtered with SEDEF-like thresholds;
it is not expected to exactly reproduce a SEDEF callset because BISER and
SEDEF use different seed/chaining/refinement algorithms.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
ERR_RE = re.compile(r"X=([0-9.]+);ID=([0-9.]+)")


@dataclass
class Call:
    idx: int
    chrom1: str
    start1: int
    end1: int
    chrom2: str
    start2: int
    end2: int
    name: str
    score_pct: float
    strand1: str
    strand2: str
    max_len: int
    aln_len: int
    cigar: str
    x_pct: float
    id_pct: float
    m_bases: int
    i_bases: int
    d_bases: int
    id_runs: int
    mismatch_est: float
    frac_match: float
    frac_match_indel: float
    sat1: float = 0.0
    sat2: float = 0.0
    nfrac1: float = 0.0
    nfrac2: float = 0.0

    @property
    def max_mate_span(self) -> int:
        return max(self.end1 - self.start1, self.end2 - self.start2)


def open_text(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def read_merged_bed(path: str) -> Dict[str, List[Tuple[int, int]]]:
    by_chr: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    with open_text(path) as fh:
        for raw in fh:
            if not raw.strip() or raw.startswith("#"):
                continue
            f = raw.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            try:
                s, e = int(f[1]), int(f[2])
            except ValueError:
                continue
            if e > s:
                by_chr[f[0]].append((s, e))
    merged: Dict[str, List[Tuple[int, int]]] = {}
    for chrom, xs in by_chr.items():
        xs.sort()
        out: List[Tuple[int, int]] = []
        for s, e in xs:
            if not out or s > out[-1][1]:
                out.append((s, e))
            else:
                out[-1] = (out[-1][0], max(out[-1][1], e))
        merged[chrom] = out
    return merged


def overlap_bp(index: Dict[str, List[Tuple[int, int]]], chrom: str, start: int, end: int) -> int:
    total = 0
    for s, e in index.get(chrom, []):
        if e <= start:
            continue
        if s >= end:
            break
        total += max(0, min(end, e) - max(start, s))
    return total


def parse_biser(path: str) -> Iterable[Call]:
    with open_text(path) as fh:
        for idx, raw in enumerate(fh, start=1):
            if not raw.strip() or raw.startswith("#"):
                continue
            f = raw.rstrip("\n").split("\t")
            if len(f) < 14:
                raise ValueError(f"{path}: line {idx} has {len(f)} columns; expected 14")
            ops = CIGAR_RE.findall(f[12])
            if not ops:
                raise ValueError(f"{path}: line {idx} has no parseable CIGAR")
            m_bases = i_bases = d_bases = id_runs = 0
            for n, op in ops:
                n = int(n)
                if op in ("M", "=", "X"):
                    m_bases += n
                elif op == "I":
                    i_bases += n
                    id_runs += 1
                elif op == "D":
                    d_bases += n
                    id_runs += 1
                # Critically: BISER N/S are omitted from aln_len/score/ID. Do not
                # reinterpret them as indel bases for the SEDEF-like criteria.
            em = ERR_RE.search(f[13])
            if em is None:
                raise ValueError(f"{path}: line {idx} lacks X=<pct>;ID=<pct>")
            x_pct, id_pct = map(float, em.groups())
            aln_len = int(f[11])
            mismatch_est = x_pct * aln_len / 100.0
            match_est = max(0.0, m_bases - mismatch_est)
            frac_match = match_est / m_bases if m_bases else 0.0
            frac_match_indel = match_est / aln_len if aln_len else 0.0
            yield Call(
                idx=idx,
                chrom1=f[0], start1=int(f[1]), end1=int(f[2]),
                chrom2=f[3], start2=int(f[4]), end2=int(f[5]),
                name=f[6], score_pct=float(f[7]), strand1=f[8], strand2=f[9],
                max_len=int(f[10]), aln_len=aln_len, cigar=f[12],
                x_pct=x_pct, id_pct=id_pct, m_bases=m_bases,
                i_bases=i_bases, d_bases=d_bases, id_runs=id_runs,
                mismatch_est=mismatch_est, frac_match=frac_match,
                frac_match_indel=frac_match_indel,
            )


def canonical_key(c: Call):
    a = (c.chrom1, c.start1, c.end1)
    b = (c.chrom2, c.start2, c.end2)
    return (a, b) if a <= b else (b, a)


def canonicalize(c: Call) -> Call:
    if (c.chrom1, c.start1, c.end1) <= (c.chrom2, c.start2, c.end2):
        return c
    return Call(
        idx=c.idx, chrom1=c.chrom2, start1=c.start2, end1=c.end2,
        chrom2=c.chrom1, start2=c.start1, end2=c.end1, name=c.name,
        score_pct=c.score_pct, strand1=c.strand2, strand2=c.strand1,
        max_len=c.max_len, aln_len=c.aln_len, cigar=c.cigar, x_pct=c.x_pct,
        id_pct=c.id_pct, m_bases=c.m_bases, i_bases=c.i_bases, d_bases=c.d_bases,
        id_runs=c.id_runs, mismatch_est=c.mismatch_est, frac_match=c.frac_match,
        frac_match_indel=c.frac_match_indel, sat1=c.sat2, sat2=c.sat1,
        nfrac1=c.nfrac2, nfrac2=c.nfrac1,
    )


def exact_pair_dedup(calls: Iterable[Call]) -> List[Call]:
    """Match the original helper's intent: only collapse exact mate-coordinate duplicates.

    Deliberately does NOT use a 50%/80% reciprocal-overlap clustering step;
    the Vollger workflow does not apply that kind of overlap deduplication.
    """
    best: Dict[tuple, Call] = {}
    for raw in calls:
        c = canonicalize(raw)
        k = canonical_key(c)
        old = best.get(k)
        if old is None or (c.frac_match, c.frac_match_indel, c.aln_len, -c.score_pct) > (
            old.frac_match, old.frac_match_indel, old.aln_len, -old.score_pct
        ):
            best[k] = c
    return sorted(best.values(), key=lambda c: (c.chrom1, c.start1, c.end1, c.chrom2, c.start2, c.end2))


def write_tsv(path: str, calls: List[Call], prefix: str) -> None:
    header = [
        "chrom1", "start1", "end1", "chrom2", "start2", "end2", "sd_id",
        "score_pct", "strand1", "strand2", "aligned_bp", "aln_len", "chain_span",
        "gap_bases", "gap_runs", "mismatch_est", "gap_compressed_identity", "gap_fraction",
        "satellite_fraction1", "satellite_fraction2", "N_fraction1", "N_fraction2",
        "biser_name", "cigar", "ungapped_aligned_bp", "fracMatchIndel_est",
        "max_mate_span_bp", "X_pct", "ID_pct",
    ]
    with open(path, "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        for i, c in enumerate(calls, start=1):
            # ``aligned_bp`` is deliberately aln_len here: it is the field used by
            # the SEDEF workflow's >=1 kb threshold, not the previous M-only proxy.
            w.writerow([
                c.chrom1, c.start1, c.end1, c.chrom2, c.start2, c.end2, f"{prefix}_SD{i:07d}",
                f"{c.score_pct:.4f}", c.strand1, c.strand2, c.aln_len, c.aln_len,
                c.max_mate_span, c.i_bases + c.d_bases, c.id_runs,
                f"{c.mismatch_est:.4f}", f"{c.frac_match:.7f}", f"{c.id_pct/100.0:.7f}",
                f"{c.sat1:.7f}", f"{c.sat2:.7f}", f"{c.nfrac1:.7f}", f"{c.nfrac2:.7f}",
                c.name, c.cigar, c.m_bases, f"{c.frac_match_indel:.7f}",
                c.max_mate_span, f"{c.x_pct:.4f}", f"{c.id_pct:.4f}",
            ])


def write_mates(path: str, calls: List[Call], prefix: str) -> None:
    with open(path, "wt") as out:
        out.write("#chrom\tstart\tend\tsd_id\tmate\tstrand\taln_len\tfracMatch\tfracMatchIndel\tsatellite_fraction\n")
        for i, c in enumerate(calls, start=1):
            sid = f"{prefix}_SD{i:07d}"
            out.write(f"{c.chrom1}\t{c.start1}\t{c.end1}\t{sid}\t1\t{c.strand1}\t{c.aln_len}\t{c.frac_match:.7f}\t{c.frac_match_indel:.7f}\t{c.sat1:.7f}\n")
            out.write(f"{c.chrom2}\t{c.start2}\t{c.end2}\t{sid}\t2\t{c.strand2}\t{c.aln_len}\t{c.frac_match:.7f}\t{c.frac_match_indel:.7f}\t{c.sat2:.7f}\n")


def main() -> None:
    ap = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
    ap.add_argument("-i", "--input", required=True, help="BISER .biser output (14 columns)")
    ap.add_argument("-o", "--out-prefix", required=True, help="Output prefix")
    ap.add_argument("--satellite-bed", help="RepeatMasker Satellite BED (BED3+); may be unmerged")
    ap.add_argument("--satellite-mode", choices=("either", "both", "first"), default="either",
                    help="either: at least one mate <= threshold [recommended]; both: both mates <= threshold; first: literal input-side only")
    ap.add_argument("--max-satellite-fraction", type=float, default=0.70)
    ap.add_argument("--gap-bed", help="Optional assembly N/gap BED; calls are annotated, not filtered unless --max-N-fraction is supplied")
    ap.add_argument("--max-N-fraction", type=float, default=None)
    ap.add_argument("--min-aln-len", type=int, default=1000, help="Minimum BISER aln_len, including I/D gaps [1000]")
    ap.add_argument("--min-frac-match", type=float, default=0.90, help="Minimum estimated BISER/SEDEF-like fracMatch [0.90]")
    ap.add_argument("--min-frac-match-indel", type=float, default=0.50, help="Minimum estimated BISER/SEDEF-like fracMatchIndel [0.50]")
    ap.add_argument("--prefix", help="SD ID prefix; default: output basename")
    args = ap.parse_args()

    if not (0 <= args.max_satellite_fraction <= 1 and 0 < args.min_frac_match <= 1 and 0 < args.min_frac_match_indel <= 1):
        ap.error("Fractions must be in (0,1], except --max-satellite-fraction may be zero.")
    if args.max_N_fraction is not None and not (0 <= args.max_N_fraction <= 1):
        ap.error("--max-N-fraction must be within [0,1].")
    if args.max_N_fraction is not None and not args.gap_bed:
        ap.error("--max-N-fraction requires --gap-bed.")

    sat = read_merged_bed(args.satellite_bed) if args.satellite_bed else {}
    gaps = read_merged_bed(args.gap_bed) if args.gap_bed else {}
    stage = defaultdict(int)
    passed: List[Call] = []

    for c in parse_biser(args.input):
        stage["raw"] += 1
        if c.aln_len < args.min_aln_len:
            continue
        stage["aln_len_ge_min"] += 1
        if c.frac_match < args.min_frac_match:
            continue
        stage["fracMatch_ge_min"] += 1
        if c.frac_match_indel < args.min_frac_match_indel:
            continue
        stage["fracMatchIndel_ge_min"] += 1

        if sat:
            c.sat1 = overlap_bp(sat, c.chrom1, c.start1, c.end1) / max(1, c.end1 - c.start1)
            c.sat2 = overlap_bp(sat, c.chrom2, c.start2, c.end2) / max(1, c.end2 - c.start2)
            if args.satellite_mode == "either":
                keep_sat = min(c.sat1, c.sat2) <= args.max_satellite_fraction
            elif args.satellite_mode == "both":
                keep_sat = max(c.sat1, c.sat2) <= args.max_satellite_fraction
            else:
                keep_sat = c.sat1 <= args.max_satellite_fraction
            if not keep_sat:
                continue
        stage["satellite_pass"] += 1

        if gaps:
            c.nfrac1 = overlap_bp(gaps, c.chrom1, c.start1, c.end1) / max(1, c.end1 - c.start1)
            c.nfrac2 = overlap_bp(gaps, c.chrom2, c.start2, c.end2) / max(1, c.end2 - c.start2)
            if args.max_N_fraction is not None and max(c.nfrac1, c.nfrac2) > args.max_N_fraction:
                continue
        stage["gap_mask_pass"] += 1
        passed.append(c)

    final = exact_pair_dedup(passed)
    prefix = args.prefix or os.path.basename(args.out_prefix)
    raw_file = args.out_prefix + ".pass_before_exact_dedup.bedpe.tsv"
    final_file = args.out_prefix + ".final.bedpe.tsv"
    mate_file = args.out_prefix + ".final.mates.bed"
    write_tsv(raw_file, passed, prefix)
    write_tsv(final_file, final, prefix)
    write_mates(mate_file, final, prefix)

    with open(args.out_prefix + ".summary.tsv", "wt") as out:
        out.write("metric\tvalue\n")
        for key in ("raw", "aln_len_ge_min", "fracMatch_ge_min", "fracMatchIndel_ge_min", "satellite_pass", "gap_mask_pass"):
            out.write(f"{key}\t{stage[key]}\n")
        out.write(f"pass_before_exact_coordinate_dedup\t{len(passed)}\n")
        out.write(f"final_exact_coordinate_dedup\t{len(final)}\n")
        out.write(f"min_aln_len\t{args.min_aln_len}\n")
        out.write(f"min_fracMatch\t{args.min_frac_match}\n")
        out.write(f"min_fracMatchIndel\t{args.min_frac_match_indel}\n")
        out.write(f"satellite_mode\t{args.satellite_mode}\n")
        out.write(f"max_satellite_fraction\t{args.max_satellite_fraction}\n")

    print(f"Final exact-coordinate-deduplicated SD pairs: {len(final):,}")
    print(f"Final table: {final_file}")
    print(f"Copy-level BED: {mate_file}")
    if not args.satellite_bed:
        print("WARNING: no Satellite BED supplied; the <=70% Satellite step was not applied.", file=sys.stderr)


if __name__ == "__main__":
    main()
