#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cdr_merge_rank.v1.py

Post-process fragmented CDR-Finder calls for heterogeneous non-human centromeres.

Design
------
* CDR-Finder calls are discovery seeds, not final biological CDRs.
* Seeds are clipped to the curated centromere BED.
* Adjacent seeds are merged in two stages:
    1) gap <= 25 kb: merge unconditionally;
    2) 25-50 kb: merge only when methylation across the gap remains depleted.
* 5-kb windows are the discovery/scoring scale.
* 10-kb windows provide broad-domain validation.
* 1-kb windows refine outer boundaries.
* Candidate ranking is driven by methylation depletion AUC, not length alone.
* SAT113 is OPTIONAL ANNOTATION ONLY. It is never used to rank or select a CDR.

Coordinate conventions
----------------------
BED / bedMethyl / outputs: 0-based, half-open.
GFF: 1-based, inclusive and converted internally to BED coordinates.

Recommended defaults were chosen for felid centromeres after inspecting PPA/PTI:
  hard gap 25 kb; conditional gap 50 kb; minimum seed support 50 kb;
  minimum seed density 0.60; minimum depletion 10 percentage points.
"""

from __future__ import annotations

import argparse
import gzip
import math
import os
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


@dataclass
class Interval:
    chrom: str
    start: int
    end: int

    @property
    def length(self) -> int:
        return max(0, self.end - self.start)


@dataclass
class Window:
    chrom: str
    start: int
    end: int
    size: int
    cpg_sites: int = 0
    valid_calls: float = 0.0
    modified_calls: float = 0.0
    methyl_pct: Optional[float] = None
    mean_site_cov: float = 0.0
    coverage_pass: bool = False

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass
class MergedSeed:
    chrom: str
    start: int
    end: int
    members: List[Interval]
    hard_links: int = 0
    conditional_links: int = 0


@dataclass
class Candidate:
    chrom: str
    seed_start: int
    seed_end: int
    start: int
    end: int
    seed_support_bp: int
    seed_density: float
    n_seed_intervals: int
    hard_links: int
    conditional_links: int
    methyl_pct: float
    baseline_pct: float
    depletion_pp: float
    depletion_auc_pp_kb: float
    coverage_fraction_5kb: float
    broad_support_fraction_10kb: float
    fine_support_fraction_1kb: float
    sat113_overlap_bp: int = 0
    sat113_overlap_fraction: float = 0.0
    sat113_distance_bp: Optional[int] = None
    confidence: str = "FAIL"
    role: str = "REJECTED"
    rank_chr: int = 0
    candidate_id: str = ""

    @property
    def length(self) -> int:
        return max(0, self.end - self.start)


# ------------------------------- utilities -----------------------------------

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def open_text(path: str):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def as_int(x: str, default: int = 0) -> int:
    try:
        return int(float(x))
    except Exception:
        return default


def as_float(x: str, default: float = float("nan")) -> float:
    try:
        return float(x)
    except Exception:
        return default


def overlap_bp(a0: int, a1: int, b0: int, b1: int) -> int:
    return max(0, min(a1, b1) - max(a0, b0))


def union_intervals(xs: Iterable[Tuple[int, int]]) -> List[Tuple[int, int]]:
    vals = sorted((int(s), int(e)) for s, e in xs if e > s)
    if not vals:
        return []
    out: List[List[int]] = [[vals[0][0], vals[0][1]]]
    for s, e in vals[1:]:
        if s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def union_bp(xs: Iterable[Tuple[int, int]]) -> int:
    return sum(e - s for s, e in union_intervals(xs))


def union_overlap_bp(xs: Iterable[Tuple[int, int]], start: int, end: int) -> int:
    clipped = []
    for s, e in xs:
        cs, ce = max(start, s), min(end, e)
        if ce > cs:
            clipped.append((cs, ce))
    return union_bp(clipped)


def quantile(values: List[float], q: float) -> Optional[float]:
    vals = sorted(v for v in values if v is not None and math.isfinite(v))
    if not vals:
        return None
    if len(vals) == 1:
        return vals[0]
    q = max(0.0, min(1.0, q))
    p = (len(vals) - 1) * q
    lo = int(math.floor(p))
    hi = int(math.ceil(p))
    if lo == hi:
        return vals[lo]
    f = p - lo
    return vals[lo] * (1 - f) + vals[hi] * f


def fmt(x, nd=4):
    if x is None:
        return "NA"
    if isinstance(x, float) and not math.isfinite(x):
        return "NA"
    if isinstance(x, float):
        return f"{x:.{nd}f}"
    return str(x)


# ------------------------------- input ---------------------------------------

def read_centromere_bed(path: str) -> Dict[str, Tuple[int, int]]:
    regions: Dict[str, Tuple[int, int]] = {}
    with open_text(path) as fh:
        for ln, line in enumerate(fh, 1):
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                f = line.split()
            if len(f) < 3:
                raise ValueError(f"{path}: line {ln} has <3 columns")
            chrom, s, e = f[0], as_int(f[1], -1), as_int(f[2], -1)
            if s < 0 or e <= s:
                if ln == 1:
                    continue
                raise ValueError(f"{path}: invalid interval at line {ln}: {line.rstrip()}")
            if chrom in regions:
                raise ValueError(
                    f"{path}: multiple centromere intervals for {chrom}. "
                    "This postprocessor expects one curated centromere interval per chromosome."
                )
            regions[chrom] = (s, e)
    if not regions:
        raise ValueError(f"No valid centromere intervals found in {path}")
    return regions


def read_raw_cdr(path: str, regions: Dict[str, Tuple[int, int]]) -> Dict[str, List[Interval]]:
    out: Dict[str, List[Interval]] = defaultdict(list)
    with open_text(path) as fh:
        for ln, line in enumerate(fh, 1):
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                f = line.split()
            if len(f) < 3:
                continue
            chrom, s, e = f[0], as_int(f[1], -1), as_int(f[2], -1)
            if chrom not in regions or s < 0 or e <= s:
                continue
            rs, re = regions[chrom]
            s, e = max(s, rs), min(e, re)
            if e > s:
                out[chrom].append(Interval(chrom, s, e))

    # Collapse overlapping/touching raw calls first; preserve only biological gaps.
    collapsed: Dict[str, List[Interval]] = defaultdict(list)
    for chrom, xs in out.items():
        for s, e in union_intervals((x.start, x.end) for x in xs):
            collapsed[chrom].append(Interval(chrom, s, e))
    return collapsed


def read_sat113(path: Optional[str], fmt_mode: str = "auto") -> Dict[str, List[Tuple[int, int]]]:
    """Read SAT113 BED or GFF. Entire supplied file is treated as SAT113."""
    out: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    if not path:
        return out

    mode = fmt_mode
    if mode == "auto":
        low = str(path).lower()
        if low.endswith(('.gff', '.gff3', '.gff.gz', '.gff3.gz')):
            mode = "gff"
        else:
            mode = "bed"

    with open_text(path) as fh:
        for ln, line in enumerate(fh, 1):
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if mode == "gff":
                if len(f) < 5:
                    continue
                chrom = f[0]
                gs, ge = as_int(f[3], -1), as_int(f[4], -1)
                if gs < 1 or ge < gs:
                    continue
                s, e = gs - 1, ge
            else:
                if len(f) < 3:
                    continue
                chrom = f[0]
                s, e = as_int(f[1], -1), as_int(f[2], -1)
                if s < 0 or e <= s:
                    continue
            out[chrom].append((s, e))

    for chrom in list(out):
        out[chrom] = union_intervals(out[chrom])
    return out


# ----------------------------- methyl windows --------------------------------

def build_windows(
    methyl_bed: str,
    regions: Dict[str, Tuple[int, int]],
    scales: List[int],
    modified_mode: str,
    min_site_valid_cov: int,
    min_cpg_per_kb: float,
    min_cpg_absolute: int,
    min_mean_site_cov: float,
    min_valid_calls: int,
) -> Dict[int, Dict[str, List[Window]]]:
    aggs = {
        scale: {
            chrom: [[0, 0.0, 0.0] for _ in range(int(math.ceil((re-rs)/scale)))]
            for chrom, (rs, re) in regions.items()
        }
        for scale in scales
    }

    seen = set()
    n_rows = 0
    n_kept = 0
    n_dup = 0
    with open_text(methyl_bed) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            n_rows += 1
            f = line.rstrip("\n").split("\t")
            if len(f) < 14:
                continue
            chrom = f[0]
            if chrom not in regions:
                continue
            start = as_int(f[1], -1)
            if start < 0:
                continue
            rs, re = regions[chrom]
            if not (rs <= start < re):
                continue
            key = (chrom, start)
            if key in seen:
                n_dup += 1
                continue
            seen.add(key)

            valid = as_int(f[9], 0)
            if valid < min_site_valid_cov:
                continue
            listed = as_int(f[11], 0)
            other = as_int(f[13], 0)
            modified = listed + other if modified_mode == "total" else listed
            modified = min(max(0, modified), max(0, valid))
            n_kept += 1

            for scale in scales:
                idx = (start - rs) // scale
                if 0 <= idx < len(aggs[scale][chrom]):
                    aggs[scale][chrom][idx][0] += 1
                    aggs[scale][chrom][idx][1] += valid
                    aggs[scale][chrom][idx][2] += modified

    eprint(f"[methyl] rows={n_rows:,}; unique retained CpGs={n_kept:,}; duplicate positions skipped={n_dup:,}")

    out: Dict[int, Dict[str, List[Window]]] = {s: {} for s in scales}
    for scale in scales:
        for chrom, (rs, re) in regions.items():
            ws: List[Window] = []
            for idx, (cpg, valid, mod) in enumerate(aggs[scale][chrom]):
                s = rs + idx * scale
                e = min(re, s + scale)
                wlen = e - s
                methyl = 100.0 * mod / valid if valid > 0 else None
                mean_cov = valid / cpg if cpg > 0 else 0.0
                min_cpg = max(min_cpg_absolute, int(math.ceil(min_cpg_per_kb * wlen / 1000.0)))
                cov_pass = (
                    cpg >= min_cpg and mean_cov >= min_mean_site_cov and valid >= min_valid_calls
                )
                ws.append(Window(
                    chrom=chrom, start=s, end=e, size=scale,
                    cpg_sites=cpg, valid_calls=valid, modified_calls=mod,
                    methyl_pct=methyl, mean_site_cov=mean_cov,
                    coverage_pass=cov_pass,
                ))
            out[scale][chrom] = ws
    return out


def raw_overlap_fraction(w: Window, raw: List[Interval]) -> float:
    if not raw or w.length <= 0:
        return 0.0
    ov = union_overlap_bp(((x.start, x.end) for x in raw), w.start, w.end)
    return ov / w.length


def chromosome_baselines(
    windows: Dict[int, Dict[str, List[Window]]],
    raw: Dict[str, List[Interval]],
    q: float,
) -> Dict[int, Dict[str, float]]:
    out: Dict[int, Dict[str, float]] = defaultdict(dict)
    for scale, by_chr in windows.items():
        for chrom, ws in by_chr.items():
            pool = [
                w.methyl_pct for w in ws
                if w.coverage_pass and w.methyl_pct is not None
                and raw_overlap_fraction(w, raw.get(chrom, [])) < 0.20
            ]
            if len(pool) < 5:
                pool = [w.methyl_pct for w in ws if w.coverage_pass and w.methyl_pct is not None]
            b = quantile(pool, q)
            out[scale][chrom] = float("nan") if b is None else b
    return out


def low_window(w: Window, baseline: float, min_drop: float, max_ratio: float, max_methyl: float) -> bool:
    if not w.coverage_pass or w.methyl_pct is None or not math.isfinite(baseline) or baseline <= 0:
        return False
    return (
        baseline - w.methyl_pct >= min_drop
        and w.methyl_pct <= baseline * max_ratio
        and w.methyl_pct <= max_methyl
    )


def weighted_fraction_in_interval(
    ws: List[Window], start: int, end: int, predicate
) -> Tuple[float, int, int]:
    elig = 0
    positive = 0
    for w in ws:
        ov = overlap_bp(start, end, w.start, w.end)
        if ov <= 0 or not w.coverage_pass or w.methyl_pct is None:
            continue
        elig += ov
        if predicate(w):
            positive += ov
    return (positive / elig if elig > 0 else 0.0, positive, elig)


# ------------------------------- merging -------------------------------------

def merge_raw_seeds(
    raw: Dict[str, List[Interval]],
    windows_5k: Dict[str, List[Window]],
    baselines_5k: Dict[str, float],
    hard_gap: int,
    conditional_gap: int,
    bridge_min_fraction: float,
    bridge_min_drop: float,
    bridge_max_ratio: float,
    bridge_max_methylation: float,
) -> List[MergedSeed]:
    merged: List[MergedSeed] = []
    for chrom, xs in sorted(raw.items()):
        if not xs:
            continue
        xs = sorted(xs, key=lambda x: (x.start, x.end))
        cur = MergedSeed(chrom, xs[0].start, xs[0].end, [xs[0]])
        for nxt in xs[1:]:
            gap = nxt.start - cur.end
            do_merge = False
            conditional = False
            if gap <= hard_gap:
                do_merge = True
            elif gap <= conditional_gap:
                baseline = baselines_5k.get(chrom, float("nan"))
                frac, _pos, elig = weighted_fraction_in_interval(
                    windows_5k.get(chrom, []), cur.end, nxt.start,
                    lambda w: low_window(
                        w, baseline, bridge_min_drop, bridge_max_ratio, bridge_max_methylation
                    )
                )
                if elig > 0 and frac >= bridge_min_fraction:
                    do_merge = True
                    conditional = True

            if do_merge:
                cur.end = max(cur.end, nxt.end)
                cur.members.append(nxt)
                if conditional:
                    cur.conditional_links += 1
                else:
                    cur.hard_links += 1
            else:
                merged.append(cur)
                cur = MergedSeed(chrom, nxt.start, nxt.end, [nxt])
        merged.append(cur)
    return merged


# ------------------------------ candidate QC ---------------------------------

def local_baseline(
    chrom: str,
    start: int,
    end: int,
    ws5: Dict[str, List[Window]],
    chrom_baseline: Dict[str, float],
    flank: int,
    min_windows: int,
    q: float,
) -> float:
    vals = []
    for w in ws5.get(chrom, []):
        if not w.coverage_pass or w.methyl_pct is None:
            continue
        in_left = overlap_bp(max(0, start-flank), start, w.start, w.end) > 0
        in_right = overlap_bp(end, end+flank, w.start, w.end) > 0
        if in_left or in_right:
            vals.append(w.methyl_pct)
    b = quantile(vals, q) if len(vals) >= min_windows else None
    if b is None:
        return chrom_baseline.get(chrom, float("nan"))
    return b


def interval_methyl(ws: List[Window], start: int, end: int) -> Tuple[float, float, int]:
    valid = 0.0
    mod = 0.0
    covered = 0
    for w in ws:
        ov = overlap_bp(start, end, w.start, w.end)
        if ov <= 0 or w.valid_calls <= 0:
            continue
        frac = ov / w.length
        valid += w.valid_calls * frac
        mod += w.modified_calls * frac
        covered += ov
    methyl = 100.0 * mod / valid if valid > 0 else float("nan")
    return methyl, valid, covered


def depletion_auc(ws: List[Window], start: int, end: int, baseline: float) -> Tuple[float, float]:
    auc = 0.0
    eligible_bp = 0
    span = max(1, end - start)
    for w in ws:
        ov = overlap_bp(start, end, w.start, w.end)
        if ov <= 0 or not w.coverage_pass or w.methyl_pct is None:
            continue
        eligible_bp += ov
        auc += max(0.0, baseline - w.methyl_pct) * ov / 1000.0
    return auc, eligible_bp / span


def refine_boundaries(
    chrom: str,
    start: int,
    end: int,
    region: Tuple[int, int],
    ws1: List[Window],
    baseline: float,
    pad: int,
    gap_windows: int,
    edge_min_drop: float,
    edge_max_ratio: float,
    edge_max_methyl: float,
    min_support_fraction: float,
) -> Tuple[int, int, float]:
    def is_low(w: Window) -> bool:
        return low_window(w, baseline, edge_min_drop, edge_max_ratio, edge_max_methyl)

    elig = [w for w in ws1 if overlap_bp(start, end, w.start, w.end) > 0 and w.coverage_pass and w.methyl_pct is not None]
    if not elig:
        return start, end, 0.0
    low_bp = sum(overlap_bp(start, end, w.start, w.end) for w in elig if is_low(w))
    elig_bp = sum(overlap_bp(start, end, w.start, w.end) for w in elig)
    support_fraction = low_bp / elig_bp if elig_bp else 0.0
    if support_fraction < min_support_fraction:
        return start, end, support_fraction

    rs, re = region
    search_s = max(rs, start - pad)
    search_e = min(re, end + pad)
    search = [w for w in ws1 if overlap_bp(search_s, search_e, w.start, w.end) > 0]
    low_idx = [i for i, w in enumerate(search) if is_low(w) and overlap_bp(start, end, w.start, w.end) > 0]
    if not low_idx:
        return start, end, support_fraction

    li, ri = min(low_idx), max(low_idx)
    left_low = li
    rebound = 0
    i = li - 1
    while i >= 0:
        w = search[i]
        if not w.coverage_pass or w.methyl_pct is None:
            rebound += 1
        elif is_low(w):
            left_low = i
            rebound = 0
        else:
            rebound += 1
        if rebound > gap_windows:
            break
        i -= 1

    right_low = ri
    rebound = 0
    i = ri + 1
    while i < len(search):
        w = search[i]
        if not w.coverage_pass or w.methyl_pct is None:
            rebound += 1
        elif is_low(w):
            right_low = i
            rebound = 0
        else:
            rebound += 1
        if rebound > gap_windows:
            break
        i += 1

    new_s = search[left_low].start
    new_e = search[right_low].end
    if new_e <= new_s:
        return start, end, support_fraction
    return new_s, new_e, support_fraction


def sat_annotation(chrom: str, start: int, end: int, sat: Dict[str, List[Tuple[int, int]]]) -> Tuple[int, float, Optional[int]]:
    xs = sat.get(chrom, [])
    if not xs:
        return 0, 0.0, None
    ov = union_overlap_bp(xs, start, end)
    if ov > 0:
        return ov, ov / max(1, end-start), 0
    d = min(
        (start - e if e <= start else s - end)
        for s, e in xs
        if e <= start or s >= end
    )
    return 0, 0.0, max(0, d)


def classify(c: Candidate, args) -> str:
    if (
        c.seed_support_bp >= args.min_supported_bp
        and c.seed_density >= args.min_density
        and c.depletion_pp >= args.min_depletion
        and c.broad_support_fraction_10kb >= args.min_broad_support
        and c.coverage_fraction_5kb >= args.min_coverage_fraction
    ):
        return "HIGH"
    if (
        c.seed_support_bp >= args.min_supported_bp
        and c.seed_density >= args.moderate_min_density
        and c.depletion_pp >= args.moderate_min_depletion
        and c.broad_support_fraction_10kb >= args.moderate_min_broad_support
        and c.coverage_fraction_5kb >= args.moderate_min_coverage_fraction
    ):
        return "MODERATE"
    if (
        c.seed_support_bp >= args.low_min_supported_bp
        and c.seed_density >= args.low_min_density
        and c.depletion_pp >= args.low_min_depletion
        and c.coverage_fraction_5kb >= args.low_min_coverage_fraction
    ):
        return "LOW"
    return "FAIL"


def grade_roles(cands: List[Candidate], secondary_ratio: float):
    by_chr: Dict[str, List[Candidate]] = defaultdict(list)
    for c in cands:
        by_chr[c.chrom].append(c)

    for chrom, xs in by_chr.items():
        usable = [c for c in xs if c.confidence in {"HIGH", "MODERATE"}]
        uncertain = False
        if not usable:
            usable = [c for c in xs if c.confidence == "LOW"]
            uncertain = bool(usable)

        for c in xs:
            if c.confidence == "FAIL":
                c.role = "REJECTED"
            elif c.confidence == "LOW":
                c.role = "LOW_CONFIDENCE"
            else:
                c.role = "CANDIDATE"

        # Always assign stable chromosome-level IDs, even when every candidate
        # fails the biological filters. This makes QC tables unambiguous.
        ranked = sorted(xs, key=lambda c: (-c.depletion_auc_pp_kb, -c.depletion_pp, c.start))
        for rank, c in enumerate(ranked, 1):
            c.rank_chr = rank
            c.candidate_id = f"{chrom}.CDR{rank:02d}"

        if not usable:
            continue

        usable.sort(key=lambda c: (-c.depletion_auc_pp_kb, -c.depletion_pp, -c.seed_support_bp, c.start))
        primary = usable[0]
        primary.role = "PRIMARY_UNCERTAIN" if uncertain else "PRIMARY"
        p_auc = max(primary.depletion_auc_pp_kb, 1e-12)
        for c in usable[1:]:
            if uncertain:
                c.role = "LOW_CONFIDENCE"
            elif c.depletion_auc_pp_kb >= secondary_ratio * p_auc:
                c.role = "SECONDARY"
            else:
                c.role = "CANDIDATE"


# ------------------------------- outputs -------------------------------------

def write_windows(outdir: Path, prefix: str, windows: Dict[int, Dict[str, List[Window]]]):
    for scale, by_chr in windows.items():
        path = outdir / f"{prefix}.windows.{scale//1000 if scale%1000==0 else scale}kb.tsv.gz"
        with gzip.open(path, "wt") as out:
            out.write("chrom\tstart\tend\twindow_bp\tcpg_sites\tvalid_calls\tmodified_calls\tmean_site_cov\tmethylation_pct\tcoverage_pass\n")
            for chrom in sorted(by_chr):
                for w in by_chr[chrom]:
                    out.write("\t".join(map(str, [
                        w.chrom, w.start, w.end, w.size, w.cpg_sites,
                        int(w.valid_calls), int(w.modified_calls), fmt(w.mean_site_cov),
                        fmt(w.methyl_pct), int(w.coverage_pass)
                    ])) + "\n")


def write_outputs(outdir: Path, prefix: str, raw, merged, cands, regions, args):
    outdir.mkdir(parents=True, exist_ok=True)

    with (outdir / f"{prefix}.raw_CDR.clipped.bed").open("w") as out:
        for chrom in sorted(raw):
            for x in raw[chrom]:
                out.write(f"{chrom}\t{x.start}\t{x.end}\n")

    with (outdir / f"{prefix}.merged_seed_domains.tsv").open("w") as out:
        out.write("chrom\tstart\tend\tspan_bp\tn_seed_intervals\tseed_support_bp\tseed_density\thard_links\tconditional_links\n")
        for m in sorted(merged, key=lambda x: (x.chrom, x.start, x.end)):
            support = union_bp((x.start, x.end) for x in m.members)
            density = support / max(1, m.end-m.start)
            out.write(f"{m.chrom}\t{m.start}\t{m.end}\t{m.end-m.start}\t{len(m.members)}\t{support}\t{density:.6f}\t{m.hard_links}\t{m.conditional_links}\n")

    header = [
        "chrom","start","end","candidate_id","rank_chr","role","confidence","length_bp",
        "seed_start","seed_end","seed_support_bp","seed_density","n_seed_intervals",
        "hard_links","conditional_links","methylation_pct","baseline_pct","depletion_pp",
        "depletion_auc_pp_kb","coverage_fraction_5kb","broad_support_fraction_10kb",
        "fine_support_fraction_1kb","SAT113_overlap_bp","SAT113_overlap_fraction","SAT113_distance_bp"
    ]
    with (outdir / f"{prefix}.CDR.graded.tsv").open("w") as out:
        out.write("\t".join(header)+"\n")
        for c in sorted(cands, key=lambda x: (x.chrom, x.rank_chr or 999, x.start)):
            row = [
                c.chrom,c.start,c.end,c.candidate_id,c.rank_chr,c.role,c.confidence,c.length,
                c.seed_start,c.seed_end,c.seed_support_bp,fmt(c.seed_density),c.n_seed_intervals,
                c.hard_links,c.conditional_links,fmt(c.methyl_pct),fmt(c.baseline_pct),fmt(c.depletion_pp),
                fmt(c.depletion_auc_pp_kb),fmt(c.coverage_fraction_5kb),fmt(c.broad_support_fraction_10kb),
                fmt(c.fine_support_fraction_1kb),c.sat113_overlap_bp,fmt(c.sat113_overlap_fraction),
                "NA" if c.sat113_distance_bp is None else c.sat113_distance_bp,
            ]
            out.write("\t".join(map(str,row))+"\n")

    primary_roles = {"PRIMARY", "PRIMARY_UNCERTAIN"}
    with (outdir / f"{prefix}.primary_CDR.bed").open("w") as out:
        for c in sorted((x for x in cands if x.role in primary_roles), key=lambda x: (x.chrom,x.start)):
            out.write(f"{c.chrom}\t{c.start}\t{c.end}\t{c.candidate_id}\t{c.role}\t{c.confidence}\n")

    with (outdir / f"{prefix}.primary_secondary_CDR.bed").open("w") as out:
        for c in sorted((x for x in cands if x.role in primary_roles | {"SECONDARY"}), key=lambda x: (x.chrom,x.start)):
            out.write(f"{c.chrom}\t{c.start}\t{c.end}\t{c.candidate_id}\t{c.role}\t{c.confidence}\n")

    with (outdir / f"{prefix}.chromosome_summary.tsv").open("w") as out:
        out.write("chrom\tcentromere_start\tcentromere_end\tn_candidates\tprimary_id\tprimary_role\tprimary_confidence\tprimary_start\tprimary_end\tprimary_auc_pp_kb\tprimary_SAT113_distance_bp\n")
        by_chr = defaultdict(list)
        for c in cands:
            by_chr[c.chrom].append(c)
        for chrom in sorted(regions):
            rs,re = regions[chrom]
            xs = by_chr.get(chrom, [])
            ps = [x for x in xs if x.role in primary_roles]
            if ps:
                p = ps[0]
                out.write(f"{chrom}\t{rs}\t{re}\t{len(xs)}\t{p.candidate_id}\t{p.role}\t{p.confidence}\t{p.start}\t{p.end}\t{fmt(p.depletion_auc_pp_kb)}\t{p.sat113_distance_bp if p.sat113_distance_bp is not None else 'NA'}\n")
            else:
                out.write(f"{chrom}\t{rs}\t{re}\t{len(xs)}\tNA\tNO_PRIMARY\tNA\tNA\tNA\tNA\tNA\n")

    params = vars(args).copy()
    with (outdir / f"{prefix}.parameters.tsv").open("w") as out:
        out.write("parameter\tvalue\n")
        for k in sorted(params):
            out.write(f"{k}\t{params[k]}\n")

    methods = f"""Primary CDR post-processing method\n\nCDR-Finder calls were treated only as discovery seeds. Calls were clipped to the curated centromere interval and overlapping calls were first collapsed. Adjacent seed intervals separated by <= {args.hard_gap:,} bp were merged directly. Gaps between {args.hard_gap:,} and {args.conditional_gap:,} bp were bridged only when at least {args.bridge_min_fraction:.2f} of eligible 5-kb windows across the gap remained hypomethylated under relaxed criteria (decrease >= {args.bridge_min_drop:g} percentage points and methylation <= {args.bridge_max_ratio:g} of the chromosome-specific centromere baseline).\n\nCpG methylation was summarized at 1-, 5- and 10-kb resolution. The 5-kb scale was used for candidate scoring, 10-kb windows for broad-domain validation, and 1-kb windows for outer-boundary refinement. The local baseline for each candidate was the {args.baseline_quantile:.2f} quantile of well-covered 5-kb windows within {args.flank:,} bp on either side of the candidate, with a chromosome-wide centromeric fallback. Candidate strength was quantified by methylation depletion and by depletion AUC (the summed positive difference between baseline and window methylation, in percentage-point x kb units).\n\nHigh-confidence candidates required >= {args.min_supported_bp:,} bp of CDR-Finder seed support, seed density >= {args.min_density:.2f}, methylation depletion >= {args.min_depletion:g} percentage points, >= {args.min_broad_support:.2f} support at 10-kb resolution, and >= {args.min_coverage_fraction:.2f} well-covered 5-kb sequence. Candidates were ranked within each chromosome primarily by depletion AUC. The strongest HIGH/MODERATE candidate was designated Primary CDR; additional candidates with depletion AUC >= {args.secondary_ratio:.2f} of the Primary were designated Secondary CDRs. If no HIGH/MODERATE candidate was present, the strongest LOW-confidence candidate was reported as PRIMARY_UNCERTAIN rather than forcing a high-confidence call. SAT113 overlap and distance, when supplied, were recorded only as independent annotations and were not used for ranking or CDR selection.\n"""
    (outdir / f"{prefix}.METHODS.txt").write_text(methods)


# -------------------------------- main ---------------------------------------

def get_args():
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Merge fragmented CDR-Finder calls and select graded Primary/Secondary CDRs."
    )
    p.add_argument("--raw-cdr", required=True, help="CDR-Finder *_CDR.bed (discovery seeds).")
    p.add_argument("--methyl-bed", required=True, help="CDR-Finder/Modkit *_methyl.bed[.gz].")
    p.add_argument("--centromere-bed", required=True, help="Curated centromere BED; one interval/chromosome.")
    p.add_argument("--sat113", default=None, help="Optional SAT113-only BED/GFF. Used only for annotation.")
    p.add_argument("--sat-format", choices=["auto","bed","gff"], default="auto")
    p.add_argument("--outdir", required=True)
    p.add_argument("--prefix", required=True)

    p.add_argument("--fine-window", type=int, default=1000)
    p.add_argument("--primary-window", type=int, default=5000)
    p.add_argument("--broad-window", type=int, default=10000)
    p.add_argument("--modified-mode", choices=["total","listed"], default="total")
    p.add_argument("--min-site-valid-cov", type=int, default=1)
    p.add_argument("--min-cpg-per-kb", type=float, default=1.0)
    p.add_argument("--min-cpg-absolute", type=int, default=3)
    p.add_argument("--min-mean-site-cov", type=float, default=5.0)
    p.add_argument("--min-valid-calls", type=int, default=0)

    p.add_argument("--hard-gap", type=int, default=25000, help="Always merge seed gaps <= this bp.")
    p.add_argument("--conditional-gap", type=int, default=50000, help="Maximum gap eligible for methylation-aware bridging.")
    p.add_argument("--bridge-min-fraction", type=float, default=0.50)
    p.add_argument("--bridge-min-drop", type=float, default=8.0)
    p.add_argument("--bridge-max-ratio", type=float, default=0.90)
    p.add_argument("--bridge-max-methylation", type=float, default=70.0)

    p.add_argument("--baseline-quantile", type=float, default=0.75)
    p.add_argument("--flank", type=int, default=100000)
    p.add_argument("--min-baseline-windows", type=int, default=4)

    p.add_argument("--broad-min-drop", type=float, default=8.0)
    p.add_argument("--broad-max-ratio", type=float, default=0.90)
    p.add_argument("--broad-max-methylation", type=float, default=70.0)

    p.add_argument("--fine-search-pad", type=int, default=25000)
    p.add_argument("--refine-edge-gap-windows", type=int, default=1)
    p.add_argument("--edge-min-drop", type=float, default=5.0)
    p.add_argument("--edge-max-ratio", type=float, default=0.95)
    p.add_argument("--edge-max-methylation", type=float, default=75.0)
    p.add_argument("--fine-min-support", type=float, default=0.30)

    # HIGH confidence
    p.add_argument("--min-supported-bp", type=int, default=50000)
    p.add_argument("--min-density", type=float, default=0.60)
    p.add_argument("--min-depletion", type=float, default=10.0)
    p.add_argument("--min-broad-support", type=float, default=0.50)
    p.add_argument("--min-coverage-fraction", type=float, default=0.50)
    # MODERATE
    p.add_argument("--moderate-min-density", type=float, default=0.50)
    p.add_argument("--moderate-min-depletion", type=float, default=8.0)
    p.add_argument("--moderate-min-broad-support", type=float, default=0.30)
    p.add_argument("--moderate-min-coverage-fraction", type=float, default=0.40)
    # LOW
    p.add_argument("--low-min-supported-bp", type=int, default=30000)
    p.add_argument("--low-min-density", type=float, default=0.35)
    p.add_argument("--low-min-depletion", type=float, default=5.0)
    p.add_argument("--low-min-coverage-fraction", type=float, default=0.30)
    p.add_argument("--secondary-ratio", type=float, default=0.50,
                   help="Secondary must have AUC >= this fraction of Primary AUC.")
    p.add_argument("--write-windows", action="store_true", help="Write 1/5/10-kb window tables (gzipped).")

    a = p.parse_args()
    if not (a.fine_window < a.primary_window < a.broad_window):
        p.error("Expected fine-window < primary-window < broad-window")
    if a.hard_gap < 0 or a.conditional_gap < a.hard_gap:
        p.error("Expected 0 <= hard-gap <= conditional-gap")
    for name in ["bridge_min_fraction","min_density","min_broad_support","min_coverage_fraction",
                 "moderate_min_density","moderate_min_broad_support","moderate_min_coverage_fraction",
                 "low_min_density","low_min_coverage_fraction","secondary_ratio","fine_min_support"]:
        v = getattr(a, name)
        if not 0 <= v <= 1:
            p.error(f"--{name.replace('_','-')} must be between 0 and 1")
    return a


def main():
    args = get_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    eprint("[1/8] Reading curated centromere BED...")
    regions = read_centromere_bed(args.centromere_bed)
    eprint(f"  chromosomes={len(regions)}")

    eprint("[2/8] Reading/clipping CDR-Finder discovery calls...")
    raw = read_raw_cdr(args.raw_cdr, regions)
    eprint(f"  clipped seed intervals={sum(len(v) for v in raw.values())}")

    eprint("[3/8] Summarizing methylation at 1/5/10-kb scales...")
    scales = [args.fine_window, args.primary_window, args.broad_window]
    windows = build_windows(
        args.methyl_bed, regions, scales, args.modified_mode,
        args.min_site_valid_cov, args.min_cpg_per_kb, args.min_cpg_absolute,
        args.min_mean_site_cov, args.min_valid_calls,
    )
    baselines = chromosome_baselines(windows, raw, args.baseline_quantile)

    eprint("[4/8] Merging raw CDR seeds with 25-kb hard / 50-kb methylation-aware rules...")
    merged = merge_raw_seeds(
        raw, windows[args.primary_window], baselines[args.primary_window],
        args.hard_gap, args.conditional_gap,
        args.bridge_min_fraction, args.bridge_min_drop,
        args.bridge_max_ratio, args.bridge_max_methylation,
    )
    eprint(f"  merged seed domains={len(merged)}")

    eprint("[5/8] Refining/scoring merged candidates...")
    sat = read_sat113(args.sat113, args.sat_format)
    cands: List[Candidate] = []
    for m in merged:
        seed_support = union_bp((x.start, x.end) for x in m.members)
        seed_density0 = seed_support / max(1, m.end-m.start)
        baseline = local_baseline(
            m.chrom, m.start, m.end,
            windows[args.primary_window], baselines[args.primary_window],
            args.flank, args.min_baseline_windows, args.baseline_quantile,
        )
        if not math.isfinite(baseline):
            baseline = baselines[args.primary_window].get(m.chrom, float("nan"))

        rs, re = regions[m.chrom]
        rstart, rend, fine_frac = refine_boundaries(
            m.chrom, m.start, m.end, (rs,re), windows[args.fine_window].get(m.chrom, []),
            baseline, args.fine_search_pad, args.refine_edge_gap_windows,
            args.edge_min_drop, args.edge_max_ratio, args.edge_max_methylation,
            args.fine_min_support,
        )

        # Safety: do not let fine-scale refinement discard >50% of discovery support.
        refined_support = union_overlap_bp(((x.start,x.end) for x in m.members), rstart, rend)
        if refined_support < 0.50 * seed_support:
            rstart, rend = m.start, m.end
            refined_support = seed_support

        density = refined_support / max(1, rend-rstart)
        methyl, _valid, _covered = interval_methyl(windows[args.fine_window].get(m.chrom, []), rstart, rend)
        depletion = baseline - methyl if math.isfinite(baseline) and math.isfinite(methyl) else float("nan")
        auc, cov_frac = depletion_auc(windows[args.primary_window].get(m.chrom, []), rstart, rend, baseline)
        broad_frac, _bpos, _belig = weighted_fraction_in_interval(
            windows[args.broad_window].get(m.chrom, []), rstart, rend,
            lambda w, b=baseline: low_window(w, b, args.broad_min_drop, args.broad_max_ratio, args.broad_max_methylation)
        )
        sat_bp, sat_frac, sat_dist = sat_annotation(m.chrom, rstart, rend, sat)

        c = Candidate(
            chrom=m.chrom, seed_start=m.start, seed_end=m.end,
            start=rstart, end=rend,
            seed_support_bp=refined_support, seed_density=density,
            n_seed_intervals=len(m.members), hard_links=m.hard_links,
            conditional_links=m.conditional_links,
            methyl_pct=methyl, baseline_pct=baseline, depletion_pp=depletion,
            depletion_auc_pp_kb=auc, coverage_fraction_5kb=cov_frac,
            broad_support_fraction_10kb=broad_frac,
            fine_support_fraction_1kb=fine_frac,
            sat113_overlap_bp=sat_bp, sat113_overlap_fraction=sat_frac,
            sat113_distance_bp=sat_dist,
        )
        c.confidence = classify(c, args)
        cands.append(c)

    eprint("[6/8] Ranking Primary/Secondary/uncertain candidates by depletion AUC...")
    grade_roles(cands, args.secondary_ratio)

    eprint("[7/8] Writing outputs...")
    write_outputs(outdir, args.prefix, raw, merged, cands, regions, args)
    if args.write_windows:
        write_windows(outdir, args.prefix, windows)

    prim = [c for c in cands if c.role in {"PRIMARY","PRIMARY_UNCERTAIN"}]
    eprint(f"[8/8] Done. Primary/uncertain primary calls={len(prim)} -> {outdir / (args.prefix + '.primary_CDR.bed')}")


if __name__ == "__main__":
    main()
