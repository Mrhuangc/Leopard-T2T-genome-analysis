#!/usr/bin/env python3
"""
Plot a local SyRI-supported assembly-discrepancy panel.

This version is designed for local INV/TRANS visualisations of the Panthera
pardus T2T assembly and addresses four practical issues:
  1) neighbouring SYN blocks can be fetched from a genome-wide SyRI table;
  2) the plotted T2T interval is clipped to the chromosome end inferred from
     depth data (or supplied explicitly), preventing an artificial blank tail;
  3) earlier-assembly coordinates are locally rescaled by default so that the
     upper and lower chromosome bars have the same displayed width while all
     structural ribbons remain visible;
  4) optional PUR masking can suppress SYN ribbons that pass through
     alignment-ambiguous regions, without reclassifying them as SyRI SVs.

Input SyRI columns (tab-delimited, no header):
ref_chr ref_start ref_end ... qry_chr qry_start qry_end event_id ... event_type ...

Outputs: PDF, SVG and PNG.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle

# Retain editable text for Illustrator.
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["font.family"] = "DejaVu Sans"

COLORS = {
    "syn": "#B7BEC8",
    "trans": "#4CAF50",
    "inv": "#F0A339",
    "dup": "#6ABF96",
    "breakpoint": "#B2182B",
    "repeat": "#DAB6D9",
    "sd": "#6C757D",
    "pur": "#C85D8B",
    "hifi": "#6FB4E7",
    "ont": "#72C5AE",
    "gene": "#333333",
}

# Panthera pardus T2T chromosome names -> chromosome names in the GCI depth files.
DEFAULT_COVERAGE_CHR = {
    "ChrA1": "Chr1", "ChrC1": "Chr2", "ChrB1": "Chr3", "ChrA2": "Chr4",
    "ChrC2": "Chr5", "ChrB2": "Chr6", "ChrB3": "Chr7", "ChrB4": "Chr8",
    "ChrA3": "Chr9", "ChrD1": "Chr10", "ChrD4": "Chr11", "ChrD3": "Chr12",
    "ChrD2": "Chr13", "ChrF2": "Chr14", "ChrF3": "Chr15", "ChrE1": "Chr16",
    "ChrE2": "Chr17", "ChrE3": "Chr18", "ChrX": "ChrX", "ChrY": "ChrY",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Plot a SyRI-supported local assembly discrepancy panel.",
    )
    p.add_argument("--syri", required=True,
                   help="Small local SyRI context file used to identify the focal event.")
    p.add_argument("--syri-all", default=None,
                   help="Genome-wide SyRI table. Overlapping SYN/SV records are added as local context.")
    p.add_argument("--t2t-chr", required=True, help="Reference/T2T chromosome name in SyRI.")
    p.add_argument("--query-chr", default=None,
                   help="Earlier-assembly/query chromosome name. Defaults to --t2t-chr.")
    p.add_argument("--event-id", default=None,
                   help="Focal SyRI event, e.g. INV89. Defaults to first non-SYN event in --syri.")
    p.add_argument("--hifi-depth", required=True, help="10-kb HiFi depth: chrom start end mean_depth.")
    p.add_argument("--ont-depth", required=True, help="10-kb ONT depth: chrom start end mean_depth.")
    p.add_argument("--coverage-chr", default=None,
                   help="Chromosome name in depth files. Defaults to built-in PPA mapping.")
    p.add_argument("--t2t-chrom-length", type=int, default=None,
                   help="Optional exact T2T chromosome length (bp). Overrides length inferred from depth.")
    p.add_argument("--t2t-sizes", default=None,
                   help="Optional 2-column T2T chrom.sizes file (chrom length). Used if --t2t-chrom-length is absent.")
    p.add_argument("--repeats", default=None,
                   help="Complete RepeatMasker BED (chrom start end class ...). Omit to skip repeat-density track.")
    p.add_argument("--repeat-window", type=int, default=50_000,
                   help="Window size for repeat occupancy.")
    p.add_argument("--sd", default=None, help="Segmental duplication BED (chrom start end).")
    p.add_argument("--pur", default=None, help="PUR BED (chrom start end ...).")
    p.add_argument("--mask-pur-synteny", action="store_true",
                   help=("Mask the portions of SYN ribbons and chromosome bars that overlap retained PUR intervals. "
                         "This is a visual ambiguity mask; it does not change SyRI calls."))
    p.add_argument("--pur-mask-min-len", type=int, default=10_000,
                   help="Only PUR intervals at least this long are used for the synteny ambiguity mask (bp).")
    p.add_argument("--pur-mask-pad", type=int, default=0,
                   help="Expand each retained PUR interval by this many bp on both sides before masking.")
    p.add_argument("--pur-mask-style", choices=("blank", "hatch"), default="blank",
                   help="Rendering style for PUR-masked parts of the two synteny bars.")
    p.add_argument("--genes", default=None, help="Gene BED (chrom start end gene_name).")
    p.add_argument("--flank", type=int, default=1_000_000,
                   help="T2T reference flank added around the focal local SyRI records (bp).")
    p.add_argument("--start", type=int, default=None, help="Manual T2T/reference plot start (bp).")
    p.add_argument("--end", type=int, default=None, help="Manual T2T/reference plot end (bp).")
    p.add_argument("--query-display", choices=("scaled", "native"), default="scaled",
                   help=("How to position the earlier-assembly bar. 'scaled' rescales its local span "
                         "to the T2T panel width; 'native' uses raw query coordinates."))
    p.add_argument("--depth-cap", type=float, default=3.0,
                   help="Coverage normalized by chromosome median and capped at this fold value.")
    p.add_argument("--out", required=True, help="Output prefix (without file extension).")
    return p.parse_args()


def read_syri(path: str) -> pd.DataFrame:
    rows: List[dict] = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line_no, line in enumerate(fh, 1):
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 11:
                raise ValueError(f"{path}: line {line_no} has <11 columns; not a supported SyRI table.")
            try:
                rows.append({
                    "ref_chr": f[0], "ref_start": int(float(f[1])), "ref_end": int(float(f[2])),
                    "ref_orient": f[3], "qry_orient": f[4], "qry_chr": f[5],
                    "qry_start": int(float(f[6])), "qry_end": int(float(f[7])),
                    "event_id": f[8], "event_type": f[10].upper(),
                })
            except ValueError as exc:
                raise ValueError(f"{path}: cannot parse coordinates at line {line_no}: {line.rstrip()}") from exc
    df = pd.DataFrame(rows)
    if df.empty:
        raise ValueError(f"{path}: no valid SyRI records found.")
    # Local panels require ascending intervals; event direction is represented graphically for INV.
    df[["ref_start", "ref_end"]] = np.sort(df[["ref_start", "ref_end"]].to_numpy(), axis=1)
    df[["qry_start", "qry_end"]] = np.sort(df[["qry_start", "qry_end"]].to_numpy(), axis=1)
    return df


def read_intervals(path: Optional[str], chrom: str, start: int, end: int) -> pd.DataFrame:
    if path is None:
        return pd.DataFrame(columns=["chrom", "start", "end", "name"])
    rows: List[Tuple[str, int, int, str]] = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3 or f[0] != chrom:
                continue
            try:
                s, e = int(float(f[1])), int(float(f[2]))
            except ValueError:
                continue
            if e < start or s > end:
                continue
            rows.append((f[0], s, e, f[3] if len(f) > 3 else ""))
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "name"])


def chromosome_length_from_sizes(path: str, chrom: str) -> Optional[int]:
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split()
            if len(f) >= 2 and f[0] == chrom:
                try:
                    return int(float(f[1]))
                except ValueError:
                    return None
    return None


def chromosome_length_from_depth(path: str, chrom: str) -> Optional[int]:
    max_end = None
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3 or f[0] != chrom:
                continue
            try:
                end = int(float(f[2]))
            except ValueError:
                continue
            max_end = end if max_end is None else max(max_end, end)
    return max_end


def read_depth(path: str, chrom: str, start: int, end: int) -> Tuple[pd.DataFrame, float]:
    all_depth: List[float] = []
    region: List[Tuple[int, int, float]] = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4 or f[0] != chrom:
                continue
            try:
                s, e, d = int(float(f[1])), int(float(f[2])), float(f[3])
            except ValueError:
                continue
            all_depth.append(d)
            if e >= start and s <= end:
                region.append((s, e, d))
    if not all_depth:
        raise ValueError(f"No depth rows for {chrom} in {path}. Check --coverage-chr.")
    return pd.DataFrame(region, columns=["start", "end", "depth"]), float(np.nanmedian(np.asarray(all_depth)))


def clip_syri_to_ref(df: pd.DataFrame, start: int, end: int) -> pd.DataFrame:
    """Clip reference intervals to panel bounds and interpolate the matching query edges."""
    out: List[dict] = []
    for _, row in df.iterrows():
        rs, re = int(row.ref_start), int(row.ref_end)
        qs, qe = int(row.qry_start), int(row.qry_end)
        left, right = max(rs, start), min(re, end)
        if right <= left or re <= rs:
            continue
        frac_left = (left - rs) / (re - rs)
        frac_right = (right - rs) / (re - rs)
        q_left = qs + frac_left * (qe - qs)
        q_right = qs + frac_right * (qe - qs)
        d = row.to_dict()
        d.update({
            "ref_start": int(left), "ref_end": int(right),
            "qry_start": int(round(min(q_left, q_right))),
            "qry_end": int(round(max(q_left, q_right))),
        })
        out.append(d)
    return pd.DataFrame(out, columns=df.columns)



def merge_bp_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Merge a list of half-open-ish bp intervals, keeping only positive spans."""
    ordered = sorted((int(s), int(e)) for s, e in intervals if int(e) > int(s))
    merged: List[List[int]] = []
    for s, e in ordered:
        if not merged or s > merged[-1][1]:
            merged.append([s, e])
        else:
            merged[-1][1] = max(merged[-1][1], e)
    return [(s, e) for s, e in merged]


def pur_mask_intervals(pur: pd.DataFrame, start: int, end: int,
                       min_len: int, pad: int) -> List[Tuple[int, int]]:
    """Return clipped/merged PUR intervals retained as an ambiguity mask."""
    if pur.empty:
        return []
    raw: List[Tuple[int, int]] = []
    for s, e in pur[["start", "end"]].to_numpy():
        s, e = int(s), int(e)
        if e - s < min_len:
            continue
        raw.append((max(start, s - pad), min(end, e + pad)))
    return merge_bp_intervals(raw)


def interval_intersections(a: List[Tuple[int, int]], b: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Intersections of two interval lists; output is merged."""
    hits: List[Tuple[int, int]] = []
    for s1, e1 in a:
        for s2, e2 in b:
            s, e = max(s1, s2), min(e1, e2)
            if e > s:
                hits.append((s, e))
    return merge_bp_intervals(hits)


def subtract_intervals(start: int, end: int, masks: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Return portions of [start,end] not covered by merged masks."""
    pieces: List[Tuple[int, int]] = []
    cursor = int(start)
    for ms, me in masks:
        if me <= cursor:
            continue
        if ms >= end:
            break
        if ms > cursor:
            pieces.append((cursor, min(ms, end)))
        cursor = max(cursor, me)
        if cursor >= end:
            break
    if cursor < end:
        pieces.append((cursor, end))
    return [(s, e) for s, e in pieces if e > s]


def ref_piece_to_query(row: pd.Series, ref_left: int, ref_right: int) -> Tuple[int, int]:
    """Linearly interpolate query edges for a clipped reference subinterval."""
    rs, re = int(row.ref_start), int(row.ref_end)
    qs, qe = int(row.qry_start), int(row.qry_end)
    if re <= rs:
        return min(qs, qe), max(qs, qe)
    f1 = (ref_left - rs) / (re - rs)
    f2 = (ref_right - rs) / (re - rs)
    q1 = qs + f1 * (qe - qs)
    q2 = qs + f2 * (qe - qs)
    return int(round(min(q1, q2))), int(round(max(q1, q2)))


def syn_ref_intervals(df: pd.DataFrame) -> List[Tuple[int, int]]:
    """Merged reference spans labelled SYN in the local context."""
    if df.empty:
        return []
    syn = df[df["event_type"].str.upper() == "SYN"]
    return merge_bp_intervals([(int(r.ref_start), int(r.ref_end)) for _, r in syn.iterrows()])


def split_synteny_by_masks(df: pd.DataFrame, masks: List[Tuple[int, int]]) -> pd.DataFrame:
    """Split SYN rows around PUR masks; non-SYN rows are retained unchanged."""
    if df.empty or not masks:
        return df.copy()
    rows: List[dict] = []
    for _, row in df.iterrows():
        if row.event_type.upper() != "SYN":
            rows.append(row.to_dict())
            continue
        for left, right in subtract_intervals(int(row.ref_start), int(row.ref_end), masks):
            qleft, qright = ref_piece_to_query(row, left, right)
            d = row.to_dict()
            d.update({"ref_start": left, "ref_end": right, "qry_start": qleft, "qry_end": qright})
            rows.append(d)
    return pd.DataFrame(rows, columns=df.columns)


def query_mask_intervals(df: pd.DataFrame, ref_masks: List[Tuple[int, int],], mapper) -> List[Tuple[float, float]]:
    """Project masked SYN portions into the displayed query coordinate system."""
    if df.empty or not ref_masks:
        return []
    projected: List[Tuple[float, float]] = []
    syn = df[df["event_type"].str.upper() == "SYN"]
    for _, row in syn.iterrows():
        for ms, me in ref_masks:
            left, right = max(int(row.ref_start), ms), min(int(row.ref_end), me)
            if right <= left:
                continue
            qleft, qright = ref_piece_to_query(row, left, right)
            x1, x2 = map(float, mapper([qleft, qright]))
            projected.append((min(x1, x2), max(x1, x2)))
    # Work in Mb/display units; merge with a tiny tolerance to prevent adjoining seam lines.
    projected = sorted((s, e) for s, e in projected if e > s)
    merged: List[List[float]] = []
    for s, e in projected:
        if not merged or s > merged[-1][1] + 1e-10:
            merged.append([s, e])
        else:
            merged[-1][1] = max(merged[-1][1], e)
    return [(s, e) for s, e in merged]


def draw_mask_rectangles(ax: plt.Axes, ref_masks: List[Tuple[int, int]], query_masks: List[Tuple[float, float]],
                         style: str) -> None:
    """Blank or hatch PUR-overlapping portions of both chromosome bars."""
    hatch = "////" if style == "hatch" else None
    edge = "#999999" if style == "hatch" else "#D8D8D8"
    for s, e in ref_masks:
        ax.add_patch(Rectangle((s / 1e6, 0.17), (e - s) / 1e6, 0.13,
                               facecolor="white", edgecolor=edge, lw=0.32, hatch=hatch, zorder=8))
    for s, e in query_masks:
        ax.add_patch(Rectangle((s, 1.04), e - s, 0.13,
                               facecolor="white", edgecolor=edge, lw=0.32, hatch=hatch, zorder=8))


def union_intervals(intervals: pd.DataFrame, start: int, end: int) -> List[Tuple[int, int]]:
    if intervals.empty:
        return []
    ordered = sorted((max(start, int(s)), min(end, int(e)))
                     for s, e in intervals[["start", "end"]].to_numpy())
    merged: List[List[int]] = []
    for s, e in ordered:
        if e <= s:
            continue
        if not merged or s > merged[-1][1]:
            merged.append([s, e])
        else:
            merged[-1][1] = max(merged[-1][1], e)
    return [(s, e) for s, e in merged]


def interval_occupancy(intervals: pd.DataFrame, start: int, end: int, window: int) -> Tuple[np.ndarray, np.ndarray]:
    left = (start // window) * window
    edges = np.arange(left, end + window, window, dtype=np.int64)
    if edges[-1] < end:
        edges = np.append(edges, end)
    cov = np.zeros(len(edges) - 1, dtype=float)
    for s, e in union_intervals(intervals, start, end):
        first = max(0, np.searchsorted(edges, s, side="right") - 1)
        last = min(len(cov) - 1, np.searchsorted(edges, e - 1, side="right") - 1)
        for i in range(first, last + 1):
            overlap = max(0, min(e, int(edges[i + 1])) - max(s, int(edges[i])))
            cov[i] += overlap
    widths = np.diff(edges).astype(float)
    mids = (edges[:-1] + edges[1:]) / 2.0
    mask = (edges[1:] > start) & (edges[:-1] < end)
    return mids[mask], (cov / widths)[mask]


def event_style(event_type: str) -> Tuple[str, str, str]:
    event_type = event_type.upper()
    if event_type == "INV":
        return COLORS["inv"], "inverted block", "Intra-chromosomal inversion"
    if event_type == "TRANS":
        return COLORS["trans"], "misordered block", "Intra-chromosomal block mis-ordering"
    if event_type.startswith("DUP"):
        return COLORS["dup"], "duplicated block", "Intra-chromosomal duplication"
    return COLORS["trans"], "rearranged block", f"SyRI {event_type} event"


def query_mapper(local: pd.DataFrame, ref_start: int, ref_end: int, mode: str):
    qmin = int(local[["qry_start", "qry_end"]].min().min())
    qmax = int(local[["qry_start", "qry_end"]].max().max())
    rleft, rright = ref_start / 1e6, ref_end / 1e6
    if mode == "scaled" or qmax <= qmin:
        qspan = max(1.0, float(qmax - qmin))

        def mapper(q):
            return rleft + (np.asarray(q, dtype=float) - qmin) / qspan * (rright - rleft)

        return mapper, (rleft, rright), qmin, qmax

    def mapper(q):
        return np.asarray(q, dtype=float) / 1e6

    return mapper, (qmin / 1e6, qmax / 1e6), qmin, qmax


def add_breakpoints(ax: plt.Axes, event_start: int, event_end: int) -> None:
    for x in (event_start / 1e6, event_end / 1e6):
        ax.axvline(x, color=COLORS["breakpoint"], lw=0.8, ls=(0, (3, 3)), alpha=0.9, zorder=30)


def style_track(ax: plt.Axes, label: str, xlim: Tuple[float, float], show_x: bool = False) -> None:
    ax.set_xlim(*xlim)
    ax.set_facecolor("white")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", left=False, labelleft=False)
    if not show_x:
        ax.tick_params(axis="x", bottom=False, labelbottom=False)
        ax.spines["bottom"].set_visible(False)
    else:
        ax.spines["bottom"].set_linewidth(0.6)
    ax.text(-0.025, 0.5, label, transform=ax.transAxes, ha="right", va="center", fontsize=8)


def draw_synteny(ax: plt.Axes, local_full: pd.DataFrame, local_draw: pd.DataFrame,
                 ref_start: int, ref_end: int, event: pd.Series, t2t_chr: str,
                 query_chr: str, query_display: str, ref_masks: List[Tuple[int, int]],
                 mask_style: str) -> None:
    xlim = (ref_start / 1e6, ref_end / 1e6)
    map_query, qbar, _, _ = query_mapper(local_full, ref_start, ref_end, query_display)
    query_masks = query_mask_intervals(local_full, ref_masks, map_query)
    event_color, block_word, _ = event_style(event.event_type)

    # Ribbons. Inversions are drawn as crossed ribbons to make the orientation explicit.
    for _, row in local_draw.iterrows():
        typ = row.event_type.upper()
        color = COLORS["syn"]
        if typ == "TRANS":
            color = COLORS["trans"]
        elif typ == "INV":
            color = COLORS["inv"]
        elif typ.startswith("DUP"):
            color = COLORS["dup"]
        q1, q2 = float(map_query(row.qry_start)), float(map_query(row.qry_end))
        r1, r2 = row.ref_start / 1e6, row.ref_end / 1e6
        if typ == "INV":
            points = [(q1, 1.04), (q2, 1.04), (r1, 0.30), (r2, 0.30)]
        else:
            points = [(q1, 1.04), (q2, 1.04), (r2, 0.30), (r1, 0.30)]
        ax.add_patch(Polygon(points, closed=True, facecolor=color, edgecolor="none", alpha=0.80, zorder=1))

    # Equal-width bars in scaled mode; native local query span in native mode.
    ax.add_patch(Rectangle((qbar[0], 1.04), qbar[1] - qbar[0], 0.13,
                           facecolor="#F1F1F1", edgecolor="#7F7F7F", lw=0.6, zorder=3))
    ax.add_patch(Rectangle((ref_start / 1e6, 0.17), (ref_end - ref_start) / 1e6, 0.13,
                           facecolor="#F1F1F1", edgecolor="#7F7F7F", lw=0.6, zorder=3))

    # Highlight focal intervals in both assemblies.
    q_event_left, q_event_right = sorted(map(float, map_query([event.qry_start, event.qry_end])))
    ax.add_patch(Rectangle((q_event_left, 1.04), max(q_event_right - q_event_left, 0.002), 0.13,
                           facecolor=event_color, edgecolor="none", alpha=0.90, zorder=4))
    ax.add_patch(Rectangle((event.ref_start / 1e6, 0.17),
                           (event.ref_end - event.ref_start) / 1e6, 0.13,
                           facecolor=event_color, edgecolor="none", alpha=0.90, zorder=4))

    # PUR overlaps in SYN blocks are visually removed after all ribbons/highlights are drawn.
    # This deliberately indicates an ambiguous local alignment, rather than converting it
    # into a SyRI rearrangement call.
    if ref_masks:
        draw_mask_rectangles(ax, ref_masks, query_masks, mask_style)

    ax.text(-0.025, 0.75, f"Earlier assembly\n{query_chr}", transform=ax.transAxes,
            ha="right", va="center", fontsize=8)
    ax.text(-0.025, 0.14, f"T2T reference\n{t2t_chr}", transform=ax.transAxes,
            ha="right", va="center", fontsize=8)
    ax.text((event.ref_start + event.ref_end) / 2e6, 0.42,
            f"{event.event_id}: {((event.ref_end - event.ref_start) / 1e6):.2f} Mb {block_word}",
            ha="center", va="bottom", color="#7A376D" if event.event_type == "INV" else "#9A4F00",
            fontsize=8, fontweight="bold")
    ax.set_ylim(0, 1.38)
    ax.set_xlim(*xlim)
    ax.axis("off")


def draw_event_track(ax: plt.Axes, event: pd.Series, xlim: Tuple[float, float]) -> None:
    color, _, descriptor = event_style(event.event_type)
    s, e = event.ref_start / 1e6, event.ref_end / 1e6
    ax.add_patch(Rectangle((s, 0.27), e - s, 0.46, facecolor=color, edgecolor="#7F4A70", lw=0.45))

    # Avoid collision between the event text and the left-side track label for events
    # located close to the beginning of the panel (e.g. ChrE1/INV106).
    panel_width = xlim[1] - xlim[0]
    midpoint = (s + e) / 2.0
    rel_midpoint = (midpoint - xlim[0]) / panel_width if panel_width > 0 else 0.5
    if rel_midpoint < 0.22:
        tx, ha = e + panel_width * 0.012, "left"
    elif rel_midpoint > 0.78:
        tx, ha = s - panel_width * 0.012, "right"
    else:
        tx, ha = midpoint, "center"
    ax.text(tx, 0.5, f"{descriptor} ({event.event_id})",
            ha=ha, va="center", fontsize=8, color="#3D2A00")
    add_breakpoints(ax, event.ref_start, event.ref_end)
    ax.set_ylim(0, 1)
    style_track(ax, "SyRI event", xlim)


def draw_interval_track(ax: plt.Axes, df: pd.DataFrame, label: str, color: str,
                        event: pd.Series, xlim: Tuple[float, float], min_display_mb: float = 0.0) -> None:
    if not df.empty:
        for _, row in df.iterrows():
            s, e = row.start / 1e6, row.end / 1e6
            width = max(e - s, min_display_mb)
            x0 = (s + e) / 2.0 - width / 2.0
            ax.add_patch(Rectangle((x0, 0.18), width, 0.64,
                                   facecolor=color, edgecolor="none", alpha=0.90))
    add_breakpoints(ax, event.ref_start, event.ref_end)
    ax.set_ylim(0, 1)
    style_track(ax, label, xlim)


def draw_repeat_track(ax: plt.Axes, repeats: pd.DataFrame, start: int, end: int, window: int,
                      event: pd.Series, xlim: Tuple[float, float]) -> None:
    if repeats.empty:
        ax.text(0.01, 0.5, "no repeat intervals loaded", transform=ax.transAxes,
                ha="left", va="center", fontsize=7, color="#777777")
    else:
        mids, occ = interval_occupancy(repeats, start, end, window)
        ax.fill_between(mids / 1e6, 0, occ, color=COLORS["repeat"], linewidth=0)
        ax.plot(mids / 1e6, occ, color="#A982A9", lw=0.45)
    add_breakpoints(ax, event.ref_start, event.ref_end)
    ax.set_ylim(0, 1.0)
    style_track(ax, f"Repeat\noccupancy\n({window // 1000} kb)", xlim)


def draw_depth_track(ax: plt.Axes, depth: pd.DataFrame, median_depth: float, label: str, color: str,
                     event: pd.Series, xlim: Tuple[float, float], cap: float) -> None:
    if depth.empty:
        ax.text(0.01, 0.5, "no depth values in window", transform=ax.transAxes,
                ha="left", va="center", fontsize=7, color="#777777")
    else:
        x = (depth["start"].to_numpy() + depth["end"].to_numpy()) / 2e6
        yraw = depth["depth"].to_numpy(dtype=float) / median_depth
        y = np.minimum(yraw, cap)
        ax.fill_between(x, 0, y, color=color, alpha=0.70, linewidth=0)
        ax.plot(x, y, color=color, lw=0.45)
        capped = yraw > cap
        if capped.any():
            ax.scatter(x[capped], np.repeat(cap * 0.965, capped.sum()), marker="^", s=8,
                       color="#A30000", linewidths=0, zorder=5)
    ax.axhline(1.0, color="#777777", lw=0.55, ls=(0, (2, 2)))
    add_breakpoints(ax, event.ref_start, event.ref_end)
    ax.set_ylim(0, cap)
    style_track(ax, f"{label}\n(depth / chr median)", xlim)
    ax.text(1.005, 1.0 / cap, "1x", transform=ax.transAxes, va="center", ha="left", fontsize=6, color="#666666")


def pack_genes(genes: pd.DataFrame, max_lanes: int = 5, gap: int = 30_000) -> pd.DataFrame:
    if genes.empty:
        genes["lane"] = []
        return genes
    g = genes.copy().sort_values(["start", "end"], kind="mergesort")
    lane_ends = [-np.inf] * max_lanes
    lanes: List[int] = []
    keep: List[bool] = []
    for _, row in g.iterrows():
        placed = False
        for lane in range(max_lanes):
            if row.start >= lane_ends[lane] + gap:
                lane_ends[lane] = row.end
                lanes.append(lane)
                keep.append(True)
                placed = True
                break
        if not placed:
            lanes.append(-1)
            keep.append(False)
    g["lane"] = lanes
    return g[np.asarray(keep, dtype=bool)].copy()


def draw_gene_track(ax: plt.Axes, genes: pd.DataFrame, event: pd.Series,
                    xlim: Tuple[float, float], region_width: int) -> None:
    g = pack_genes(genes)
    if g.empty:
        ax.text(0.01, 0.5, "no genes in window", transform=ax.transAxes,
                ha="left", va="center", fontsize=7, color="#777777")
        nlanes = 1
    else:
        nlanes = max(1, int(g["lane"].max()) + 1)
        labels = g.assign(length=g.end - g.start).sort_values("length", ascending=False).head(12)
        label_ids = set(labels.index)
        for idx, row in g.iterrows():
            y = nlanes - row.lane - 0.5
            x1, x2 = row.start / 1e6, row.end / 1e6
            ax.plot([x1, x2], [y, y], color=COLORS["gene"], lw=0.75, solid_capstyle="butt")
            ax.add_patch(Rectangle((x1, y - 0.12), max(x2 - x1, 0.002), 0.24,
                                   facecolor=COLORS["gene"], edgecolor="none"))
            if idx in label_ids and (row.end - row.start) >= max(10_000, region_width * 0.005):
                ax.text((x1 + x2) / 2, y + 0.20, str(row["name"]), ha="center", va="bottom", fontsize=6.5)
    add_breakpoints(ax, event.ref_start, event.ref_end)
    ax.set_ylim(0, nlanes + 0.65)
    style_track(ax, "Gene annotation", xlim, show_x=True)
    ax.set_xlabel("T2T reference coordinate (Mb)", fontsize=9)
    ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(nbins=7))
    ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%.1f"))
    ax.tick_params(axis="x", labelsize=8)


def main() -> None:
    args = parse_args()
    focal_source = read_syri(args.syri)
    query_chr = args.query_chr or args.t2t_chr
    focal = focal_source[(focal_source.ref_chr == args.t2t_chr) &
                         (focal_source.qry_chr == query_chr)].copy()
    if focal.empty:
        raise ValueError(f"No focal SyRI rows found for reference={args.t2t_chr}, query={query_chr}.")

    sv = focal[focal.event_type != "SYN"].copy()
    if sv.empty:
        raise ValueError("The focal SyRI file contains no non-SYN event to highlight.")
    if args.event_id is not None:
        sv = sv[sv.event_id == args.event_id]
        if sv.empty:
            raise ValueError(f"Event ID {args.event_id!r} was not found in the focal SyRI file.")
    event = sv.iloc[0].copy()

    coverage_chr = args.coverage_chr or DEFAULT_COVERAGE_CHR.get(args.t2t_chr, args.t2t_chr)
    if args.t2t_chrom_length is not None:
        t2t_len = args.t2t_chrom_length
    elif args.t2t_sizes is not None:
        t2t_len = chromosome_length_from_sizes(args.t2t_sizes, args.t2t_chr)
    else:
        t2t_len = chromosome_length_from_depth(args.hifi_depth, coverage_chr)
    if t2t_len is None:
        raise ValueError("Could not determine the T2T chromosome length. Provide --t2t-chrom-length or --t2t-sizes.")

    # Use the whole focal input extent (event + immediately supplied local SYN rows) as the seed.
    if args.start is not None and args.end is not None:
        region_start, region_end = args.start, args.end
    else:
        region_start = int(focal.ref_start.min()) - args.flank
        region_end = int(focal.ref_end.max()) + args.flank
    region_start = max(0, region_start)
    region_end = min(int(t2t_len), region_end)
    if region_end <= region_start:
        raise ValueError("Invalid plotting region after clipping to the T2T chromosome length.")

    context_source = read_syri(args.syri_all) if args.syri_all else focal_source
    context = context_source[(context_source.ref_chr == args.t2t_chr) &
                             (context_source.qry_chr == query_chr)].copy()
    if context.empty:
        context = focal.copy()
    context = clip_syri_to_ref(context, region_start, region_end)
    if context.empty:
        raise ValueError("No SyRI blocks overlap the selected panel interval.")

    hifi, hifi_median = read_depth(args.hifi_depth, coverage_chr, region_start, region_end)
    ont, ont_median = read_depth(args.ont_depth, coverage_chr, region_start, region_end)
    repeats = read_intervals(args.repeats, args.t2t_chr, region_start, region_end)
    sd = read_intervals(args.sd, args.t2t_chr, region_start, region_end)
    pur = read_intervals(args.pur, args.t2t_chr, region_start, region_end)
    genes = read_intervals(args.genes, args.t2t_chr, region_start, region_end)

    # A PUR interval is an independent mappability/assembly-complexity annotation.
    # Do not edit SyRI calls. Instead, when requested, mask only its overlap with SYN
    # blocks in the figure and retain SV ribbons/calls unchanged.
    syn_pur_masks: List[Tuple[int, int]] = []
    context_draw = context
    if args.mask_pur_synteny:
        retained_pur = pur_mask_intervals(pur, region_start, region_end,
                                          args.pur_mask_min_len, args.pur_mask_pad)
        syn_pur_masks = interval_intersections(retained_pur, syn_ref_intervals(context))
        context_draw = split_synteny_by_masks(context, syn_pur_masks)

    xlim = (region_start / 1e6, region_end / 1e6)

    print(f"Selected event: {event.event_id} ({event.event_type})")
    print(f"T2T block: {event.ref_chr}:{event.ref_start}-{event.ref_end} ({(event.ref_end - event.ref_start) / 1e6:.3f} Mb)")
    print(f"Query block: {event.qry_chr}:{event.qry_start}-{event.qry_end}")
    print(f"Panel: {args.t2t_chr}:{region_start}-{region_end} (chromosome length={t2t_len})")
    print(f"Context blocks: {len(context)}; displayed ribbons: {len(context_draw)}; query display={args.query_display}")
    if args.mask_pur_synteny:
        masked_bp = sum(e - s for s, e in syn_pur_masks)
        print(f"PUR ambiguity mask: {len(syn_pur_masks)} SYN-overlapping interval(s), {masked_bp:,} bp "
              f"(min PUR length={args.pur_mask_min_len:,} bp; pad={args.pur_mask_pad:,} bp)")
    print(f"Coverage chromosome: {coverage_chr}; HiFi median={hifi_median:.2f}x; ONT median={ont_median:.2f}x")
    if args.repeats and repeats.empty:
        print("WARNING: no repeat intervals overlap this region. Check chromosome naming and full RepeatMasker BED.", file=sys.stderr)

    fig = plt.figure(figsize=(11.5, 8.1), constrained_layout=False)
    gs = fig.add_gridspec(nrows=8, ncols=1,
                          height_ratios=[2.1, 0.65, 0.45, 0.45, 0.9, 1.0, 1.0, 1.35],
                          hspace=0.13)
    ax_syn = fig.add_subplot(gs[0])
    ax_event = fig.add_subplot(gs[1])
    ax_sd = fig.add_subplot(gs[2], sharex=ax_event)
    ax_pur = fig.add_subplot(gs[3], sharex=ax_event)
    ax_rep = fig.add_subplot(gs[4], sharex=ax_event)
    ax_hifi = fig.add_subplot(gs[5], sharex=ax_event)
    ax_ont = fig.add_subplot(gs[6], sharex=ax_event)
    ax_genes = fig.add_subplot(gs[7], sharex=ax_event)

    draw_synteny(ax_syn, context, context_draw, region_start, region_end, event, args.t2t_chr,
                 query_chr, args.query_display, syn_pur_masks, args.pur_mask_style)
    draw_event_track(ax_event, event, xlim)
    draw_interval_track(ax_sd, sd, "Segmental\nduplications", COLORS["sd"], event, xlim, min_display_mb=0.002)
    draw_interval_track(ax_pur, pur, "PUR", COLORS["pur"], event, xlim, min_display_mb=0.025)
    draw_repeat_track(ax_rep, repeats, region_start, region_end, args.repeat_window, event, xlim)
    draw_depth_track(ax_hifi, hifi, hifi_median, "HiFi coverage", COLORS["hifi"], event, xlim, args.depth_cap)
    draw_depth_track(ax_ont, ont, ont_median, "ONT coverage", COLORS["ont"], event, xlim, args.depth_cap)
    draw_gene_track(ax_genes, genes, event, xlim, region_end - region_start)

    fig.suptitle(
        f"Putative assembly error supported by SyRI and local sequence features: {args.t2t_chr} / {event.event_id}",
        x=0.52, y=0.987, fontsize=11, fontweight="bold"
    )
    note = "Vertical dashed lines: T2T breakpoints. Depth is normalized to the chromosome median; red triangles indicate values above the display cap."
    if args.query_display == "scaled":
        note += " Earlier-assembly local coordinates are rescaled only for visual alignment of chromosome bars."
    if args.mask_pur_synteny:
        note += (f" White bar/ribbon gaps mark PUR-overlapping SYN segments (PUR ≥ {args.pur_mask_min_len // 1000} kb); "
                 "they are treated as alignment-ambiguous, not as SyRI SVs.")
    fig.text(0.115, 0.012, note, fontsize=6.25, ha="left", color="#555555")
    fig.subplots_adjust(left=0.16, right=0.985, top=0.95, bottom=0.065)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    for ext, dpi in [("pdf", None), ("svg", None), ("png", 400)]:
        target = str(out) + "." + ext
        fig.savefig(target, dpi=dpi, bbox_inches="tight")
        print(f"Wrote {target}")
    plt.close(fig)


if __name__ == "__main__":
    main()
