#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Predict candidate centromere regions using:
1. CL1_113 repeat monomer positions
2. TRASH satellite annotation
3. EarlGrey repeat annotation
4. intra-chromosomal Hi-C interaction density
5. synteny breakpoint regions

Main logic:
- If breakpoint region exists, use it as a prior search interval.
  Then refine inward according to CL1_113 and satellite/repeat density.
- If no breakpoint region exists, infer centromere from CL1_113 + satellite/repeat density.
- Hi-C is used only as weak auxiliary evidence, not as a hard filter.

Author: ChatGPT
"""

import argparse
import pandas as pd
import numpy as np
from collections import defaultdict


def read_fai(fai):
    chrom_len = {}
    with open(fai) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip().split("\t")
            chrom_len[parts[0]] = int(parts[1])
    return chrom_len


def read_bed3(path, allow_na=False):
    """
    Read BED-like file. Only first 3 columns are used.
    If allow_na=True, rows containing NA in start/end will be kept as NaN.
    """
    df = pd.read_csv(path, sep="\t", header=None, comment="#", dtype=str)

    if df.shape[1] < 3:
        raise ValueError(f"{path} should contain at least 3 columns.")

    df = df.iloc[:, :3]
    df.columns = ["chrom", "start", "end"]

    if allow_na:
        df["start"] = pd.to_numeric(df["start"], errors="coerce")
        df["end"] = pd.to_numeric(df["end"], errors="coerce")
    else:
        df["start"] = pd.to_numeric(df["start"], errors="coerce")
        df["end"] = pd.to_numeric(df["end"], errors="coerce")
        df = df.dropna(subset=["start", "end"])

    df["start"] = df["start"].astype("Int64")
    df["end"] = df["end"].astype("Int64")
    return df


def read_hic(path):
    """
    Expected columns:
    chrom start end value
    With or without header.
    """
    df = pd.read_csv(path, sep="\t", comment="#")

    if list(df.columns[:4]) == ["chrom", "start", "end", "value"]:
        pass
    else:
        df = pd.read_csv(path, sep="\t", header=None, comment="#")
        if df.shape[1] < 4:
            raise ValueError("Hi-C file should contain at least 4 columns: chrom start end value")
        df = df.iloc[:, :4]
        df.columns = ["chrom", "start", "end", "value"]

    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df["value"] = pd.to_numeric(df["value"], errors="coerce")
    df = df.dropna(subset=["chrom", "start", "end", "value"])
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    return df


def overlap_len(a_start, a_end, b_start, b_end):
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def compute_window_signal(chrom, chrom_len, bed_df, window_size, step_size, region_start=None, region_end=None):
    """
    Calculate bp coverage density of BED intervals in sliding windows.
    Return dataframe:
    chrom start end bp density count
    """
    if region_start is None:
        region_start = 0
    if region_end is None:
        region_end = chrom_len

    region_start = max(0, int(region_start))
    region_end = min(chrom_len, int(region_end))

    sub = bed_df[bed_df["chrom"] == chrom].copy()
    sub = sub.dropna(subset=["start", "end"])

    windows = []
    pos = region_start
    while pos < region_end:
        w_start = pos
        w_end = min(pos + window_size, region_end)
        if w_end <= w_start:
            break

        bp = 0
        count = 0
        if not sub.empty:
            candidate = sub[(sub["start"] < w_end) & (sub["end"] > w_start)]
            count = candidate.shape[0]
            for _, row in candidate.iterrows():
                bp += overlap_len(w_start, w_end, int(row["start"]), int(row["end"]))

        density = bp / (w_end - w_start)
        windows.append([chrom, w_start, w_end, bp, density, count])
        pos += step_size

    return pd.DataFrame(windows, columns=["chrom", "start", "end", "bp", "density", "count"])


def get_cluster_from_windows(win_df, score_col, min_score, merge_gap):
    """
    Merge consecutive/high-scoring windows into clusters.
    """
    high = win_df[win_df[score_col] >= min_score].copy()
    if high.empty:
        return []

    high = high.sort_values(["chrom", "start", "end"])
    clusters = []

    cur_start = None
    cur_end = None
    scores = []

    for _, row in high.iterrows():
        s = int(row["start"])
        e = int(row["end"])
        sc = float(row[score_col])

        if cur_start is None:
            cur_start, cur_end = s, e
            scores = [sc]
        else:
            if s <= cur_end + merge_gap:
                cur_end = max(cur_end, e)
                scores.append(sc)
            else:
                clusters.append((cur_start, cur_end, np.mean(scores), max(scores)))
                cur_start, cur_end = s, e
                scores = [sc]

    if cur_start is not None:
        clusters.append((cur_start, cur_end, np.mean(scores), max(scores)))

    return clusters


def robust_zscore(series):
    """
    Robust z-score using median and MAD.
    """
    x = np.array(series, dtype=float)
    med = np.nanmedian(x)
    mad = np.nanmedian(np.abs(x - med))
    if mad == 0 or np.isnan(mad):
        sd = np.nanstd(x)
        if sd == 0 or np.isnan(sd):
            return np.zeros_like(x)
        return (x - np.nanmean(x)) / sd
    return 0.6745 * (x - med) / mad


def normalize_0_1(series):
    x = np.array(series, dtype=float)
    mn = np.nanmin(x)
    mx = np.nanmax(x)
    if mx == mn:
        return np.zeros_like(x)
    return (x - mn) / (mx - mn)


def build_evidence_windows(
    chrom,
    chrom_len,
    cl113,
    trash,
    repeats,
    hic,
    window_size,
    step_size,
    search_start,
    search_end,
    w_cl113,
    w_trash,
    w_repeat,
    w_hic
):
    """
    Build integrated evidence score in sliding windows.
    Hi-C is treated as weak evidence:
      lower intra-chromosomal interaction value gets slightly higher score.
    """
    cl_win = compute_window_signal(chrom, chrom_len, cl113, window_size, step_size, search_start, search_end)
    tr_win = compute_window_signal(chrom, chrom_len, trash, window_size, step_size, search_start, search_end)
    rp_win = compute_window_signal(chrom, chrom_len, repeats, window_size, step_size, search_start, search_end)

    win = cl_win[["chrom", "start", "end"]].copy()
    win["cl113_bp"] = cl_win["bp"]
    win["cl113_density"] = cl_win["density"]
    win["cl113_count"] = cl_win["count"]

    win["trash_bp"] = tr_win["bp"]
    win["trash_density"] = tr_win["density"]
    win["trash_count"] = tr_win["count"]

    win["repeat_bp"] = rp_win["bp"]
    win["repeat_density"] = rp_win["density"]
    win["repeat_count"] = rp_win["count"]

    # Hi-C value per window: weighted mean of overlapping Hi-C bins.
    hic_sub = hic[hic["chrom"] == chrom].copy()
    hic_values = []
    for _, row in win.iterrows():
        w_start, w_end = int(row["start"]), int(row["end"])
        cand = hic_sub[(hic_sub["start"] < w_end) & (hic_sub["end"] > w_start)]
        if cand.empty:
            hic_values.append(np.nan)
        else:
            weights = []
            values = []
            for _, h in cand.iterrows():
                ov = overlap_len(w_start, w_end, int(h["start"]), int(h["end"]))
                if ov > 0:
                    weights.append(ov)
                    values.append(float(h["value"]))
            if len(weights) == 0:
                hic_values.append(np.nan)
            else:
                hic_values.append(np.average(values, weights=weights))

    win["hic_value"] = hic_values

    # normalize signals
    win["cl113_score"] = normalize_0_1(win["cl113_density"])
    win["trash_score"] = normalize_0_1(win["trash_density"])
    win["repeat_score"] = normalize_0_1(win["repeat_density"])

    # lower Hi-C value as weak positive evidence
    if win["hic_value"].notna().sum() > 2:
        hic_z = robust_zscore(win["hic_value"])
        # lower value -> higher score
        win["hic_low_score"] = normalize_0_1(-hic_z)
    else:
        win["hic_low_score"] = 0.0

    win["total_score"] = (
        w_cl113 * win["cl113_score"] +
        w_trash * win["trash_score"] +
        w_repeat * win["repeat_score"] +
        w_hic * win["hic_low_score"]
    )

    return win


def refine_region_from_evidence(
    chrom,
    chrom_len,
    windows,
    breakpoint_start,
    breakpoint_end,
    min_total_score,
    min_cl113_density,
    min_trash_density,
    boundary_pad,
    merge_gap,
    allow_expand_outside_breakpoint,
    max_expand_outside_breakpoint
):
    """
    Decide candidate centromere region.

    Key rule:
    - If breakpoint exists, refine inward using evidence-rich windows.
    - Do not blindly expand to breakpoint boundary.
    - Only expand outside breakpoint if high-confidence evidence lies outside.
    """

    has_bp = breakpoint_start is not None and breakpoint_end is not None

    # High confidence windows:
    # 1. CL1_113-positive windows
    # 2. strong satellite windows
    # 3. high integrated score windows
    evidence = windows[
        (windows["cl113_density"] >= min_cl113_density) |
        (windows["trash_density"] >= min_trash_density) |
        (windows["total_score"] >= min_total_score)
    ].copy()

    if evidence.empty:
        # fallback: use best total score window
        best = windows.sort_values("total_score", ascending=False).head(1)
        if best.empty:
            return None
        start = int(best["start"].iloc[0])
        end = int(best["end"].iloc[0])
        note = "low_confidence_best_window_only"
    else:
        # Prefer cluster containing CL1_113.
        cl_windows = evidence[evidence["cl113_density"] > 0]
        if not cl_windows.empty:
            start = int(cl_windows["start"].min())
            end = int(cl_windows["end"].max())

            # add nearby high-evidence windows that touch/close to CL1_113 cluster
            nearby = evidence[
                (evidence["start"] <= end + merge_gap) &
                (evidence["end"] >= start - merge_gap)
            ]
            start = int(nearby["start"].min())
            end = int(nearby["end"].max())
            note = "CL1_113_supported"
        else:
            # no CL1_113 inside search interval; use strongest total-score cluster
            clusters = get_cluster_from_windows(evidence, "total_score", min_total_score, merge_gap)
            if clusters:
                clusters = sorted(clusters, key=lambda x: (x[3], x[2], x[1] - x[0]), reverse=True)
                start, end, _, _ = clusters[0]
                note = "satellite_repeat_supported_no_CL1_113"
            else:
                best = evidence.sort_values("total_score", ascending=False).head(1)
                start = int(best["start"].iloc[0])
                end = int(best["end"].iloc[0])
                note = "best_evidence_window"

    # Apply padding around evidence-supported region.
    start = max(0, start - boundary_pad)
    end = min(chrom_len, end + boundary_pad)

    if has_bp:
        bp_start, bp_end = int(breakpoint_start), int(breakpoint_end)

        if allow_expand_outside_breakpoint:
            # Allow limited expansion only if evidence region goes outside breakpoint.
            allowed_start = max(0, bp_start - max_expand_outside_breakpoint)
            allowed_end = min(chrom_len, bp_end + max_expand_outside_breakpoint)
            start = max(start, allowed_start)
            end = min(end, allowed_end)
        else:
            # Strict inward refinement: candidate must be inside breakpoint.
            # If evidence extends outside, cut it at breakpoint boundary.
            start = max(start, bp_start)
            end = min(end, bp_end)

        # If strict cutting makes invalid region, fallback to breakpoint-overlapping best evidence.
        if start >= end:
            overlap_win = windows[
                (windows["start"] < bp_end) &
                (windows["end"] > bp_start)
            ].copy()
            if not overlap_win.empty:
                best = overlap_win.sort_values("total_score", ascending=False).head(1)
                start = max(bp_start, int(best["start"].iloc[0]) - boundary_pad)
                end = min(bp_end, int(best["end"].iloc[0]) + boundary_pad)
                note += ";fallback_best_inside_breakpoint"
            else:
                start, end = bp_start, bp_end
                note += ";fallback_whole_breakpoint"

    if start < 0:
        start = 0
    if end > chrom_len:
        end = chrom_len

    if start >= end:
        return None

    return start, end, note


def summarize_region(chrom, start, end, cl113, trash, repeats, hic):
    def cov_bp(df):
        sub = df[(df["chrom"] == chrom) & (df["start"] < end) & (df["end"] > start)]
        total = 0
        for _, r in sub.iterrows():
            total += overlap_len(start, end, int(r["start"]), int(r["end"]))
        return total, sub.shape[0]

    cl_bp, cl_n = cov_bp(cl113)
    tr_bp, tr_n = cov_bp(trash)
    rp_bp, rp_n = cov_bp(repeats)

    hic_sub = hic[(hic["chrom"] == chrom) & (hic["start"] < end) & (hic["end"] > start)]
    hic_mean = np.nan
    if not hic_sub.empty:
        hic_mean = hic_sub["value"].mean()

    length = end - start
    return {
        "region_len": length,
        "cl113_bp": cl_bp,
        "cl113_count": cl_n,
        "cl113_density": cl_bp / length if length > 0 else 0,
        "trash_bp": tr_bp,
        "trash_count": tr_n,
        "trash_density": tr_bp / length if length > 0 else 0,
        "repeat_bp": rp_bp,
        "repeat_count": rp_n,
        "repeat_density": rp_bp / length if length > 0 else 0,
        "hic_mean": hic_mean
    }


def main():
    parser = argparse.ArgumentParser(
        description="Predict candidate centromere regions using CL1_113, TRASH satellite, EarlGrey repeats, Hi-C and synteny breakpoints."
    )

    parser.add_argument("--fai", required=True, help="Genome FASTA .fai file")
    parser.add_argument("--repeat-bed", required=True, help="EarlGrey repeat annotation BED")
    parser.add_argument("--trash-bed", required=True, help="TRASH satellite annotation BED")
    parser.add_argument("--cl113-bed", required=True, help="113 bp repeat monomer BED")
    parser.add_argument("--hic", required=True, help="Hi-C intra-chromosomal interaction density file: chrom start end value")
    parser.add_argument("--breakpoint-bed", required=True, help="Synteny breakpoint BED. Rows with NA are allowed")

    parser.add_argument("-o", "--output-prefix", default="predicted_centromere", help="Output prefix")

    parser.add_argument("--window-size", type=int, default=200000, help="Sliding window size. Default: 200000")
    parser.add_argument("--step-size", type=int, default=None, help="Sliding step size. Default: same as window-size")

    parser.add_argument("--search-flank", type=int, default=5000000,
                        help="If no breakpoint exists, search around CL1_113 merged region with this flank. Default: 5 Mb")

    parser.add_argument("--boundary-pad", type=int, default=200000,
                        help="Padding added to evidence-supported region. Default: 200 kb")

    parser.add_argument("--merge-gap", type=int, default=400000,
                        help="Merge evidence windows separated by <= this gap. Default: 400 kb")

    parser.add_argument("--min-total-score", type=float, default=0.35,
                        help="Minimum integrated score for high-evidence windows. Default: 0.35")

    parser.add_argument("--min-cl113-density", type=float, default=0.0001,
                        help="Minimum CL1_113 bp density to define CL1_113-supported window. Default: 0.0001")

    parser.add_argument("--min-trash-density", type=float, default=0.02,
                        help="Minimum TRASH satellite bp density for satellite-supported window. Default: 0.02")

    parser.add_argument("--w-cl113", type=float, default=0.50, help="Weight for CL1_113 evidence. Default: 0.50")
    parser.add_argument("--w-trash", type=float, default=0.30, help="Weight for TRASH satellite evidence. Default: 0.30")
    parser.add_argument("--w-repeat", type=float, default=0.15, help="Weight for EarlGrey repeat density. Default: 0.15")
    parser.add_argument("--w-hic", type=float, default=0.05, help="Weight for Hi-C low interaction evidence. Default: 0.05")

    parser.add_argument("--allow-expand-outside-breakpoint", action="store_true",
                        help="Allow limited expansion outside synteny breakpoint if evidence supports it.")

    parser.add_argument("--max-expand-outside-breakpoint", type=int, default=1000000,
                        help="Max allowed expansion outside breakpoint when --allow-expand-outside-breakpoint is used. Default: 1 Mb")

    parser.add_argument("--whole-chrom-if-no-cl113", action="store_true",
                        help="If a chromosome lacks CL1_113 and breakpoint, scan whole chromosome using TRASH/repeat density.")

    args = parser.parse_args()

    if args.step_size is None:
        args.step_size = args.window_size

    chrom_len = read_fai(args.fai)
    repeats = read_bed3(args.repeat_bed)
    trash = read_bed3(args.trash_bed)
    cl113 = read_bed3(args.cl113_bed)
    hic = read_hic(args.hic)
    breakpoints = read_bed3(args.breakpoint_bed, allow_na=True)

    results = []
    detail_windows = []

    for chrom, clen in chrom_len.items():
        bp_sub = breakpoints[breakpoints["chrom"] == chrom].copy()
        bp_sub_valid = bp_sub.dropna(subset=["start", "end"])

        has_bp = not bp_sub_valid.empty

        if has_bp:
            # If multiple breakpoint intervals exist on one chromosome,
            # use their merged span as prior search interval.
            bp_start = int(bp_sub_valid["start"].min())
            bp_end = int(bp_sub_valid["end"].max())

            search_start = max(0, bp_start)
            search_end = min(clen, bp_end)
            search_mode = "breakpoint_prior_refine_inward"

        else:
            bp_start = None
            bp_end = None

            cl_sub = cl113[cl113["chrom"] == chrom].dropna(subset=["start", "end"])

            if not cl_sub.empty:
                cl_start = int(cl_sub["start"].min())
                cl_end = int(cl_sub["end"].max())
                search_start = max(0, cl_start - args.search_flank)
                search_end = min(clen, cl_end + args.search_flank)
                search_mode = "no_breakpoint_CL1_113_centered"
            else:
                if args.whole_chrom_if_no_cl113:
                    search_start = 0
                    search_end = clen
                    search_mode = "no_breakpoint_no_CL1_113_whole_chrom_scan"
                else:
                    # No usable evidence. Output NA-like record.
                    results.append({
                        "chrom": chrom,
                        "start": "NA",
                        "end": "NA",
                        "source": "no_breakpoint_no_CL1_113_skipped",
                        "breakpoint_start": "NA",
                        "breakpoint_end": "NA",
                        "search_start": "NA",
                        "search_end": "NA",
                        "region_len": "NA",
                        "cl113_bp": 0,
                        "cl113_count": 0,
                        "cl113_density": 0,
                        "trash_bp": 0,
                        "trash_count": 0,
                        "trash_density": 0,
                        "repeat_bp": 0,
                        "repeat_count": 0,
                        "repeat_density": 0,
                        "hic_mean": "NA"
                    })
                    continue

        windows = build_evidence_windows(
            chrom=chrom,
            chrom_len=clen,
            cl113=cl113,
            trash=trash,
            repeats=repeats,
            hic=hic,
            window_size=args.window_size,
            step_size=args.step_size,
            search_start=search_start,
            search_end=search_end,
            w_cl113=args.w_cl113,
            w_trash=args.w_trash,
            w_repeat=args.w_repeat,
            w_hic=args.w_hic
        )

        windows["search_mode"] = search_mode
        windows["breakpoint_start"] = bp_start if bp_start is not None else "NA"
        windows["breakpoint_end"] = bp_end if bp_end is not None else "NA"
        detail_windows.append(windows)

        refined = refine_region_from_evidence(
            chrom=chrom,
            chrom_len=clen,
            windows=windows,
            breakpoint_start=bp_start,
            breakpoint_end=bp_end,
            min_total_score=args.min_total_score,
            min_cl113_density=args.min_cl113_density,
            min_trash_density=args.min_trash_density,
            boundary_pad=args.boundary_pad,
            merge_gap=args.merge_gap,
            allow_expand_outside_breakpoint=args.allow_expand_outside_breakpoint,
            max_expand_outside_breakpoint=args.max_expand_outside_breakpoint
        )

        if refined is None:
            results.append({
                "chrom": chrom,
                "start": "NA",
                "end": "NA",
                "source": "no_valid_region",
                "breakpoint_start": bp_start if bp_start is not None else "NA",
                "breakpoint_end": bp_end if bp_end is not None else "NA",
                "search_start": search_start,
                "search_end": search_end,
                "region_len": "NA",
                "cl113_bp": 0,
                "cl113_count": 0,
                "cl113_density": 0,
                "trash_bp": 0,
                "trash_count": 0,
                "trash_density": 0,
                "repeat_bp": 0,
                "repeat_count": 0,
                "repeat_density": 0,
                "hic_mean": "NA"
            })
            continue

        start, end, note = refined
        summary = summarize_region(chrom, start, end, cl113, trash, repeats, hic)

        rec = {
            "chrom": chrom,
            "start": start,
            "end": end,
            "source": f"{search_mode};{note}",
            "breakpoint_start": bp_start if bp_start is not None else "NA",
            "breakpoint_end": bp_end if bp_end is not None else "NA",
            "search_start": search_start,
            "search_end": search_end
        }
        rec.update(summary)
        results.append(rec)

    result_df = pd.DataFrame(results)

    # Main BED output
    bed_out = f"{args.output_prefix}.candidate_centromere.bed"
    bed_df = result_df[result_df["start"] != "NA"][["chrom", "start", "end", "source"]]
    bed_df.to_csv(bed_out, sep="\t", index=False, header=False)

    # Detailed summary
    summary_out = f"{args.output_prefix}.candidate_centromere.summary.tsv"
    result_df.to_csv(summary_out, sep="\t", index=False)

    # Window evidence table
    if detail_windows:
        win_out = f"{args.output_prefix}.window_evidence.tsv"
        pd.concat(detail_windows, ignore_index=True).to_csv(win_out, sep="\t", index=False)

    print(f"[Done] Candidate centromere BED: {bed_out}")
    print(f"[Done] Summary table: {summary_out}")
    if detail_windows:
        print(f"[Done] Window evidence table: {win_out}")


if __name__ == "__main__":
    main()
