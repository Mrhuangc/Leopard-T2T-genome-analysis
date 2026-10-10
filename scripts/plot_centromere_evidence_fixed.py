#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
plot_centromere_evidence_fixed.py

Visualize chromosome-wise centromere evidence tracks:
  - CL1_113
  - FA_SAT / CL2
  - TRASH satellite
  - EarlGrey repeats
  - optional intra-chromosomal Hi-C interaction values and low-value regions
  - synteny breakpoint prior regions
  - final candidate centromere regions

Key fixes in this version:
  1. merge_intervals() is fully defined before get_low_hic_spans() uses it.
  2. Window coverage is calculated from interval UNION, so density cannot exceed 1.
     Annotation count is still calculated from raw intervals.
  3. --hic is optional. Without Hi-C, plots are still generated normally.
  4. Dense filled signal tracks can be rasterized in PDF/SVG output to reduce file size.
  5. The overview x-axis is displayed on top and always includes chromosome end.
"""

import os
import argparse
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Patch
from matplotlib.lines import Line2D


# =========================
# Matplotlib export settings
# =========================
mpl.rcParams["pdf.compression"] = 9
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["svg.fonttype"] = "none"


# =========================
# I/O
# =========================
def read_fai(fai_file):
    chrom_len = {}
    with open(fai_file, "r") as f:
        for line in f:
            if not line.strip():
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 2:
                raise ValueError(f"Invalid FAI line: {line.rstrip()}")
            chrom_len[fields[0]] = int(fields[1])
    return chrom_len


def read_bed_generic(path, allow_na=False):
    """Read first three columns of BED-like file: chrom, start, end."""
    df = pd.read_csv(path, sep="\t", header=None, comment="#", dtype=str)
    if df.shape[1] < 3:
        raise ValueError(f"{path} 至少需要 3 列: chrom start end")

    df = df.iloc[:, :3].copy()
    df.columns = ["chrom", "start", "end"]
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")

    if allow_na:
        return df

    df = df.dropna(subset=["start", "end"]).copy()
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df = df[df["end"] > df["start"]].copy()
    return df


def read_breakpoint_bed(path):
    """Breakpoint BED permits NA in start/end fields."""
    return read_bed_generic(path, allow_na=True)


def read_hic(path):
    """
    Read Hi-C track with four columns: chrom, start, end, value.
    Accepts header lines such as '#chrom\tstart\tend\tvalue', a normal header,
    or no header.
    """
    raw = pd.read_csv(path, sep="\t", header=None, dtype=str, comment=None)
    if raw.shape[1] < 4:
        raise ValueError("Hi-C 文件至少需要4列: chrom start end value")

    df = raw.iloc[:, :4].copy()
    df.columns = ["chrom", "start", "end", "value"]
    df["chrom"] = df["chrom"].astype(str).str.lstrip("#")
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df["value"] = pd.to_numeric(df["value"], errors="coerce")
    df = df.dropna(subset=["start", "end", "value"]).copy()
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df = df[df["end"] > df["start"]].copy()
    return df


def load_optional_hic(path):
    if path is None:
        print("[INFO] No Hi-C file provided. Hi-C track will be skipped.")
        return None
    if not os.path.exists(path):
        print(f"[WARNING] Hi-C file not found: {path}. Hi-C track will be skipped.")
        return None
    if os.path.getsize(path) == 0:
        print(f"[WARNING] Hi-C file is empty: {path}. Hi-C track will be skipped.")
        return None
    try:
        hic = read_hic(path)
    except Exception as exc:
        print(f"[WARNING] Failed to read Hi-C file: {path}. Hi-C track will be skipped.")
        print(f"[WARNING] {exc}")
        return None
    if hic.empty:
        print(f"[WARNING] No valid Hi-C records found in: {path}. Hi-C track will be skipped.")
        return None
    print(f"[INFO] Hi-C track enabled: {path}")
    return hic


# =========================
# Utility
# =========================
def ensure_outdir(path):
    os.makedirs(path, exist_ok=True)


def overlap_len(a_start, a_end, b_start, b_end):
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def merge_intervals(intervals, max_gap=0):
    """
    Merge intervals if overlapping or separated by <= max_gap.
    This function is used both for plotting regions and Hi-C low spans.
    """
    if not intervals:
        return []

    cleaned = []
    for s, e in intervals:
        if pd.isna(s) or pd.isna(e):
            continue
        s, e = int(s), int(e)
        if e > s:
            cleaned.append((s, e))
    if not cleaned:
        return []

    cleaned.sort(key=lambda x: (x[0], x[1]))
    merged = [[cleaned[0][0], cleaned[0][1]]]
    for s, e in cleaned[1:]:
        if s <= merged[-1][1] + max_gap:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [(s, e) for s, e in merged]


def merge_intervals_for_coverage(intervals):
    """
    Merge only overlapping or directly touching intervals for union coverage.
    The returned union makes bp/window density bounded by 0..1.
    """
    return merge_intervals(intervals, max_gap=0)


def minmax_norm(values):
    x = np.asarray(values, dtype=float)
    if len(x) == 0 or np.isfinite(x).sum() == 0:
        return np.zeros_like(x)
    xmin = np.nanmin(x)
    xmax = np.nanmax(x)
    if xmax == xmin:
        return np.zeros_like(x)
    return (x - xmin) / (xmax - xmin)


def mb(value):
    return value / 1e6


def get_midpoints(df):
    if df.empty:
        return np.array([])
    return (df["start"].to_numpy() + df["end"].to_numpy()) / 2.0


def make_mb_ticks(xmax):
    """Generate readable Mb ticks and always show the chromosome end."""
    if xmax <= 80:
        step = 10
    elif xmax <= 160:
        step = 20
    elif xmax <= 300:
        step = 50
    else:
        step = 100

    ticks = list(np.arange(0, xmax, step, dtype=float))
    if len(ticks) == 0 or ticks[0] != 0:
        ticks.insert(0, 0.0)
    if abs(ticks[-1] - xmax) > 1e-8:
        ticks.append(float(xmax))
    return ticks


def format_mb_tick(tick, xmax):
    if abs(tick - xmax) < 1e-8:
        return f"{tick:.2f}"
    return f"{tick:.0f}"


# =========================
# Window signal calculation
# =========================
def compute_window_coverage(intervals_df, chrom, chrom_len, window_size, step_size):
    """
    Calculate sliding-window evidence efficiently.

    - bp/density: based on UNION coverage, so density is always in [0, 1].
    - count: number of raw annotation intervals overlapping each window,
      computed with numpy.searchsorted instead of a Python nested loop.
    """
    starts = np.arange(0, chrom_len, step_size, dtype=np.int64)
    ends = np.minimum(starts + window_size, chrom_len)
    mids = (starts + ends) / 2.0

    sub = intervals_df.loc[
        intervals_df["chrom"].eq(chrom),
        ["start", "end"]
    ].dropna()

    if sub.empty:
        zeros_f = np.zeros(len(starts), dtype=float)
        zeros_i = np.zeros(len(starts), dtype=np.int64)
        return pd.DataFrame({
            "chrom": chrom,
            "start": starts,
            "end": ends,
            "mid": mids,
            "bp": zeros_f,
            "density": zeros_f.copy(),
            "count": zeros_i
        })

    arr = sub.to_numpy(dtype=np.int64, copy=False)
    arr = arr[arr[:, 1] > arr[:, 0]]
    if len(arr) == 0:
        zeros_f = np.zeros(len(starts), dtype=float)
        zeros_i = np.zeros(len(starts), dtype=np.int64)
        return pd.DataFrame({
            "chrom": chrom,
            "start": starts,
            "end": ends,
            "mid": mids,
            "bp": zeros_f,
            "density": zeros_f.copy(),
            "count": zeros_i
        })

    # Sort raw intervals once.
    order = np.lexsort((arr[:, 1], arr[:, 0]))
    arr = arr[order]

    # Efficient overlap counts:
    # count(window) = number(start < window_end) - number(end <= window_start)
    raw_starts = np.sort(arr[:, 0])
    raw_ends = np.sort(arr[:, 1])
    counts = (
        np.searchsorted(raw_starts, ends, side="left")
        - np.searchsorted(raw_ends, starts, side="right")
    ).astype(np.int64)

    # Build interval UNION in one pass.
    union_intervals = []
    cur_s, cur_e = int(arr[0, 0]), int(arr[0, 1])
    for s, e in arr[1:]:
        s, e = int(s), int(e)
        if s <= cur_e:
            if e > cur_e:
                cur_e = e
        else:
            union_intervals.append((cur_s, cur_e))
            cur_s, cur_e = s, e
    union_intervals.append((cur_s, cur_e))

    # Sliding-window UNION coverage with moving pointers.
    bp_cov = np.zeros(len(starts), dtype=float)
    j = 0
    n_union = len(union_intervals)

    for wi, (ws, we) in enumerate(zip(starts, ends)):
        while j < n_union and union_intervals[j][1] <= ws:
            j += 1

        k = j
        cov = 0
        while k < n_union and union_intervals[k][0] < we:
            s, e = union_intervals[k]
            cov += min(we, e) - max(ws, s)
            k += 1
        bp_cov[wi] = cov

    window_lengths = ends - starts
    density = np.divide(
        bp_cov,
        window_lengths,
        out=np.zeros_like(bp_cov),
        where=window_lengths > 0
    )
    np.clip(density, 0.0, 1.0, out=density)

    return pd.DataFrame({
        "chrom": chrom,
        "start": starts,
        "end": ends,
        "mid": mids,
        "bp": bp_cov,
        "density": density,
        "count": counts
    })


def compute_hic_track(hic_df, chrom, chrom_len, window_size, step_size):
    """
    Compute weighted Hi-C signal per sliding window using moving pointers.
    This avoids restarting from the first Hi-C interval for every window.
    """
    starts = np.arange(0, chrom_len, step_size, dtype=np.int64)
    ends = np.minimum(starts + window_size, chrom_len)
    mids = (starts + ends) / 2.0
    values = np.full(len(starts), np.nan, dtype=float)

    if hic_df is None or hic_df.empty:
        return pd.DataFrame({
            "chrom": chrom, "start": starts, "end": ends, "mid": mids, "value": values
        })

    sub = hic_df.loc[
        hic_df["chrom"].eq(chrom),
        ["start", "end", "value"]
    ].dropna()

    if sub.empty:
        return pd.DataFrame({
            "chrom": chrom, "start": starts, "end": ends, "mid": mids, "value": values
        })

    arr = sub.to_numpy(dtype=float, copy=False)
    order = np.lexsort((arr[:, 1], arr[:, 0]))
    arr = arr[order]

    j = 0
    n = len(arr)

    for wi, (ws, we) in enumerate(zip(starts, ends)):
        while j < n and arr[j, 1] <= ws:
            j += 1

        k = j
        weighted_sum = 0.0
        total_weight = 0

        while k < n and arr[k, 0] < we:
            hs, he, hv = arr[k]
            ov = min(we, int(he)) - max(ws, int(hs))
            if ov > 0:
                weighted_sum += float(hv) * ov
                total_weight += ov
            k += 1

        if total_weight > 0:
            values[wi] = weighted_sum / total_weight

    return pd.DataFrame({
        "chrom": chrom, "start": starts, "end": ends, "mid": mids, "value": values
    })


def get_low_hic_spans(hic_raw_df, chrom, q=0.10, merge_gap=0):
    """Return merged regions occupied by the bottom-q Hi-C windows."""
    if hic_raw_df is None or hic_raw_df.empty:
        return [], np.nan

    sub = hic_raw_df[hic_raw_df["chrom"] == chrom].copy()
    if sub.empty:
        return [], np.nan

    threshold = sub["value"].quantile(q)
    low = sub[sub["value"] <= threshold].copy()
    if low.empty:
        return [], threshold

    spans = low[["start", "end"]].sort_values(["start", "end"]).values.tolist()
    return merge_intervals(spans, max_gap=merge_gap), threshold


# =========================
# Region collection
# =========================
def get_regions_for_chrom(df, chrom):
    sub = df[df["chrom"] == chrom].dropna(subset=["start", "end"]).copy()
    if sub.empty:
        return []
    spans = sub[["start", "end"]].values.tolist()
    return merge_intervals(spans, max_gap=0)


# =========================
# Plot one chromosome
# =========================
def plot_one_chrom(
    chrom,
    chrom_len,
    cl1_df,
    cl2_df,
    trash_df,
    repeat_df,
    hic_df,
    breakpoint_df,
    cen_df,
    window_size,
    step_size,
    hic_low_q,
    out_png,
    out_svg,
    pdf_handle=None,
    rasterize_tracks=True,
    png_dpi=300,
    pdf_dpi=300
):
    has_hic_for_chrom = hic_df is not None and not hic_df[hic_df["chrom"] == chrom].empty

    print(f"    [TRACK] {chrom}: CL1")
    cl1_win = compute_window_coverage(cl1_df, chrom, chrom_len, window_size, step_size)
    print(f"    [TRACK] {chrom}: CL2")
    cl2_win = compute_window_coverage(cl2_df, chrom, chrom_len, window_size, step_size)
    print(f"    [TRACK] {chrom}: TRASH")
    trash_win = compute_window_coverage(trash_df, chrom, chrom_len, window_size, step_size)
    print(f"    [TRACK] {chrom}: EarlGrey repeat")
    repeat_win = compute_window_coverage(repeat_df, chrom, chrom_len, window_size, step_size)
    print(f"    [TRACK] {chrom}: Hi-C")
    hic_win = compute_hic_track(hic_df, chrom, chrom_len, window_size, step_size)

    cl1_norm = minmax_norm(cl1_win["density"].values)
    cl2_norm = minmax_norm(cl2_win["density"].values)
    trash_norm = minmax_norm(trash_win["density"].values)
    repeat_norm = minmax_norm(repeat_win["density"].values)

    low_hic_spans, hic_thr = get_low_hic_spans(hic_df, chrom, q=hic_low_q, merge_gap=0)
    bp_regions = get_regions_for_chrom(breakpoint_df, chrom)
    cen_regions = get_regions_for_chrom(cen_df, chrom)

    cl1_mid = get_midpoints(cl1_df[cl1_df["chrom"] == chrom]) / 1e6
    cl2_mid = get_midpoints(cl2_df[cl2_df["chrom"] == chrom]) / 1e6
    trash_mid = get_midpoints(trash_df[trash_df["chrom"] == chrom]) / 1e6

    x_mb = cl1_win["mid"].values / 1e6
    xmax = chrom_len / 1e6
    xticks = make_mb_ticks(xmax)
    xtick_labels = [format_mb_tick(t, xmax) for t in xticks]

    C = {
        "cl1": "#d73027",
        "cl2": "#e6ab02",
        "trash": "#7570b3",
        "repeat": "#4c78a8",
        "hic": "#666666",
        "hic_low": "#b2182b",
        "bp_fill": "#92c5de",
        "bp_edge": "#2166ac",
        "cen_fill": "#fdae61",
        "cen_edge": "#e6550d",
        "chrom_bar": "#9e9e9e"
    }

    fig = plt.figure(figsize=(15, 11), constrained_layout=False)
    gs = fig.add_gridspec(
        nrows=6,
        ncols=1,
        height_ratios=[1.8, 1, 1, 1, 1, 1.25],
        hspace=0.20
    )
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[1, 0], sharex=ax0)
    ax2 = fig.add_subplot(gs[2, 0], sharex=ax0)
    ax3 = fig.add_subplot(gs[3, 0], sharex=ax0)
    ax4 = fig.add_subplot(gs[4, 0], sharex=ax0)
    ax5 = fig.add_subplot(gs[5, 0], sharex=ax0)
    all_axes = [ax0, ax1, ax2, ax3, ax4, ax5]

    for ax in all_axes:
        ax.set_xlim(0, xmax)
        ax.set_xticks(xticks)

    # ---------- Overview ----------
    ax0.set_ylim(0, 1.15)
    for s, e in bp_regions:
        ax0.axvspan(mb(s), mb(e), color=C["bp_fill"], alpha=0.30, lw=0, zorder=0)
        ax0.axvline(mb(s), color=C["bp_edge"], lw=1.4, ls="--", zorder=1)
        ax0.axvline(mb(e), color=C["bp_edge"], lw=1.4, ls="--", zorder=1)
    if has_hic_for_chrom:
        for s, e in low_hic_spans:
            ax0.axvspan(mb(s), mb(e), color=C["hic_low"], alpha=0.10, lw=0, zorder=0)
    for s, e in cen_regions:
        ax0.axvspan(mb(s), mb(e), color=C["cen_fill"], alpha=0.35, lw=0, zorder=0)
        ax0.axvline(mb(s), color=C["cen_edge"], lw=1.8, zorder=2)
        ax0.axvline(mb(e), color=C["cen_edge"], lw=1.8, zorder=2)

    ax0.hlines(1.02, 0, xmax, color=C["chrom_bar"], lw=6, zorder=3)
    overview_rows = [("SAT113", cl1_mid, 0.75, C["cl1"]),
                     ("SAT483(FA_SAT)", cl2_mid, 0.50, C["cl2"]),
                     ("TRASH", trash_mid, 0.25, C["trash"])]
    label_x = -0.01 * xmax
    for label, mids, ypos, color in overview_rows:
        if len(mids) > 0:
            ax0.scatter(mids, np.full_like(mids, ypos), s=20, color=color, marker="|", linewidths=1.0)
        ax0.text(label_x, ypos, label, ha="right", va="center", fontsize=10, color=color)
    ax0.text(label_x, 1.02, "Chr", ha="right", va="center", fontsize=10)
    ax0.set_yticks([])
    ax0.set_ylabel("Overview", fontsize=11)
    for spine in ["top", "right", "left"]:
        ax0.spines[spine].set_visible(False)
    ax0.tick_params(axis="x", labelbottom=False, bottom=False, labeltop=True, top=True, labelsize=9, pad=2)
    ax0.set_xticklabels(xtick_labels)

    # ---------- Density tracks ----------
    def draw_density_track(ax, x, y_norm, color, label, raw_max):
        for s, e in bp_regions:
            ax.axvspan(mb(s), mb(e), color=C["bp_fill"], alpha=0.18, lw=0, zorder=0)
        for s, e in cen_regions:
            ax.axvspan(mb(s), mb(e), color=C["cen_fill"], alpha=0.18, lw=0, zorder=0)
        ax.fill_between(x, 0, y_norm, color=color, alpha=0.60, linewidth=0, zorder=2,
                        rasterized=rasterize_tracks)
        ax.plot(x, y_norm, color=color, lw=1.1, zorder=3, rasterized=rasterize_tracks)
        ax.set_ylim(0, 1.05)
        ax.set_yticks([0, 1])
        ax.set_yticklabels(["0", "1"], fontsize=8)
        ax.set_ylabel(label, rotation=0, labelpad=42, fontsize=10, va="center")
        ax.text(1.002, 0.88, f"max={raw_max:.4g}", transform=ax.transAxes,
                ha="left", va="center", fontsize=8, color=color)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(axis="x", labelbottom=False, bottom=False)

    draw_density_track(ax1, x_mb, cl1_norm, C["cl1"], "CL1", cl1_win["density"].max())
    draw_density_track(ax2, x_mb, cl2_norm, C["cl2"], "CL2", cl2_win["density"].max())
    draw_density_track(ax3, x_mb, trash_norm, C["trash"], "TRASH", trash_win["density"].max())
    draw_density_track(ax4, x_mb, repeat_norm, C["repeat"], "Repeat", repeat_win["density"].max())

    # ---------- Hi-C track ----------
    for s, e in bp_regions:
        ax5.axvspan(mb(s), mb(e), color=C["bp_fill"], alpha=0.18, lw=0, zorder=0)
    for s, e in cen_regions:
        ax5.axvspan(mb(s), mb(e), color=C["cen_fill"], alpha=0.18, lw=0, zorder=0)

    if has_hic_for_chrom:
        hic_x = hic_win["mid"].values / 1e6
        hic_y = hic_win["value"].values
        ax5.plot(hic_x, hic_y, color=C["hic"], lw=1.0, zorder=2, rasterized=rasterize_tracks)
        low_mask = np.isfinite(hic_y) & (hic_y <= hic_thr)
        if low_mask.any():
            ax5.scatter(hic_x[low_mask], hic_y[low_mask], s=18, color=C["hic_low"], zorder=3,
                        rasterized=rasterize_tracks)
        ax5.axhline(hic_thr, color=C["hic_low"], lw=1.2, ls="--", alpha=0.9)
        ax5.text(1.002, 0.88, f"Q{int(hic_low_q * 100)}={hic_thr:.4g}",
                 transform=ax5.transAxes, ha="left", va="center", fontsize=8, color=C["hic_low"])
    else:
        ax5.text(0.5, 0.5, "Hi-C data not provided for this chromosome",
                 transform=ax5.transAxes, ha="center", va="center", fontsize=10, color="#666666")
        ax5.set_ylim(0, 1)

    ax5.set_ylabel("Hi-C", rotation=0, labelpad=42, fontsize=10, va="center")
    ax5.spines["top"].set_visible(False)
    ax5.spines["right"].set_visible(False)
    ax5.tick_params(axis="x", labelbottom=True, bottom=True, labelsize=9)
    ax5.set_xticklabels(xtick_labels)
    ax5.set_xlabel("Genomic position (Mb)", fontsize=11)

    bp_text = "None" if not bp_regions else "; ".join(f"{mb(s):.2f}-{mb(e):.2f} Mb" for s, e in bp_regions)
    cen_text = "None" if not cen_regions else "; ".join(f"{mb(s):.2f}-{mb(e):.2f} Mb" for s, e in cen_regions)
    hic_title = f"Hi-C low threshold = bottom {int(hic_low_q * 100)}%" if has_hic_for_chrom else "Hi-C = not provided"
    fig.suptitle(
        f"{chrom}  |  length = {xmax:.2f} Mb\n"
        f"Predicted centromere: {cen_text}   |   Breakpoint region: {bp_text}\n"
        f"Window = {window_size:,} bp, Step = {step_size:,} bp, {hic_title}",
        fontsize=13, y=0.995
    )

    legend_handles = [
        Patch(facecolor=C["cen_fill"], edgecolor=C["cen_edge"], alpha=0.35, label="Predicted centromere region"),
        Patch(facecolor=C["bp_fill"], edgecolor=C["bp_edge"], alpha=0.30, label="Synteny breakpoint region"),
        Line2D([0], [0], color=C["cl1"], lw=2, label="SAT113 density"),
        Line2D([0], [0], color=C["cl2"], lw=2, label="SAT483(FASAT) density"),
        Line2D([0], [0], color=C["trash"], lw=2, label="TRASH density"),
        Line2D([0], [0], color=C["repeat"], lw=2, label="EarlGrey repeat density")
    ]
    if has_hic_for_chrom:
        legend_handles.insert(2, Patch(facecolor=C["hic_low"], edgecolor=C["hic_low"], alpha=0.10,
                                       label="Low Hi-C region"))
        legend_handles.append(Line2D([0], [0], color=C["hic"], lw=1.5, label="Hi-C interaction value"))

    fig.legend(handles=legend_handles, loc="upper center", bbox_to_anchor=(0.5, 0.925),
               ncol=3, frameon=False, fontsize=9)
    fig.subplots_adjust(top=0.84, left=0.10, right=0.92, bottom=0.07, hspace=0.20)

    fig.savefig(out_png, dpi=png_dpi)
    fig.savefig(out_svg, dpi=pdf_dpi)
    if pdf_handle is not None:
        pdf_handle.savefig(fig, dpi=pdf_dpi)
    plt.close(fig)


# =========================
# Main
# =========================
def main():
    parser = argparse.ArgumentParser(
        description="Visualize centromere evidence tracks for each chromosome. Hi-C is optional."
    )
    parser.add_argument("--hic", default=None, help="Optional: intra-chromosomal Hi-C track: chrom start end value")
    parser.add_argument("--repeat-bed", required=True, help="EarlGrey repeat annotation BED")
    parser.add_argument("--cl1-bed", required=True, help="CL1_113 BED")
    parser.add_argument("--cl2-bed", required=True, help="FA_SAT / CL2 BED")
    parser.add_argument("--breakpoint-bed", required=True, help="Synteny breakpoint BED; NA allowed")
    parser.add_argument("--trash-bed", required=True, help="TRASH satellite BED")
    parser.add_argument("--fai", required=True, help="Genome FASTA .fai file")
    parser.add_argument("--cen-bed", required=True, help="Final candidate centromere BED")
    parser.add_argument("--window-size", type=int, default=200000, help="Sliding window size. Default: 200000")
    parser.add_argument("--step-size", type=int, default=50000, help="Sliding step size. Default: 50000")
    parser.add_argument("--hic-low-q", type=float, default=0.10, help="Low Hi-C quantile. Default: 0.10")
    parser.add_argument("--chroms", nargs="*", default=None, help="Optional chromosome IDs to plot")
    parser.add_argument("--outdir", default="centromere_evidence_plots", help="Output directory")
    parser.add_argument("--png-dpi", type=int, default=300, help="PNG resolution. Default: 300")
    parser.add_argument("--pdf-dpi", type=int, default=300, help="Raster layer resolution in PDF. Default: 300")
    parser.add_argument("--no-rasterize", action="store_true", help="Keep dense tracks as vectors (larger PDF)")
    args = parser.parse_args()

    if args.window_size <= 0 or args.step_size <= 0:
        raise ValueError("--window-size and --step-size must be positive integers.")
    if not (0 < args.hic_low_q < 1):
        raise ValueError("--hic-low-q must be between 0 and 1.")

    ensure_outdir(args.outdir)
    png_dir = os.path.join(args.outdir, "png")
    svg_dir = os.path.join(args.outdir, "svg")
    ensure_outdir(png_dir)
    ensure_outdir(svg_dir)

    print("[INFO] Reading input files ...")
    chrom_len = read_fai(args.fai)
    repeat_df = read_bed_generic(args.repeat_bed)
    cl1_df = read_bed_generic(args.cl1_bed)
    cl2_df = read_bed_generic(args.cl2_bed)
    trash_df = read_bed_generic(args.trash_bed)
    breakpoint_df = read_breakpoint_bed(args.breakpoint_bed)
    cen_df = read_bed_generic(args.cen_bed)
    hic_df = load_optional_hic(args.hic)

    if args.chroms:
        unknown = [c for c in args.chroms if c not in chrom_len]
        if unknown:
            print(f"[WARNING] Chromosomes absent from FAI and skipped: {', '.join(unknown)}")
        plot_chroms = [c for c in args.chroms if c in chrom_len]
    else:
        plot_chroms = list(chrom_len.keys())
    if not plot_chroms:
        raise ValueError("No valid chromosomes to plot.")

    pdf_path = os.path.join(args.outdir, "all_chromosomes.centromere_evidence.pdf")
    print(f"[INFO] Plotting {len(plot_chroms)} chromosomes ...")
    with PdfPages(pdf_path) as pdf:
        for chrom in plot_chroms:
            print(f"[INFO] Plotting {chrom} ...")
            plot_one_chrom(
                chrom=chrom,
                chrom_len=chrom_len[chrom],
                cl1_df=cl1_df,
                cl2_df=cl2_df,
                trash_df=trash_df,
                repeat_df=repeat_df,
                hic_df=hic_df,
                breakpoint_df=breakpoint_df,
                cen_df=cen_df,
                window_size=args.window_size,
                step_size=args.step_size,
                hic_low_q=args.hic_low_q,
                out_png=os.path.join(png_dir, f"{chrom}.centromere_evidence.png"),
                out_svg=os.path.join(svg_dir, f"{chrom}.centromere_evidence.svg"),
                pdf_handle=pdf,
                rasterize_tracks=not args.no_rasterize,
                png_dpi=args.png_dpi,
                pdf_dpi=args.pdf_dpi
            )

    print("[DONE]")
    print(f"PDF: {pdf_path}")
    print(f"PNG dir: {png_dir}")
    print(f"SVG dir: {svg_dir}")
    print(f"[INFO] Dense tracks rasterized: {'No' if args.no_rasterize else 'Yes'}")


if __name__ == "__main__":
    main()
