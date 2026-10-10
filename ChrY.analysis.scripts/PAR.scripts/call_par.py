#!/usr/bin/env python3
"""Call terminal X–Y alignment chains as PAR candidates from a PAF file."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from dataclasses import dataclass

from common import PafRecord, parse_paf


@dataclass
class Chain:
    qname: str
    qlen: int
    qstart: int
    qend: int
    tname: str
    tlen: int
    tstart: int
    tend: int
    strand: str
    aligned_bp: int
    weighted_matches: int
    n_alignments: int
    oq_end: int
    ot_end: int

    @property
    def identity(self) -> float:
        return self.weighted_matches / self.aligned_bp if self.aligned_bp else 0.0


def terminal_side(start: int, end: int, length: int, end_window: int) -> str:
    left = start <= end_window
    right = length - end <= end_window
    if left and right:
        return "both"
    if left:
        return "left"
    if right:
        return "right"
    return "internal"


def oriented_target(rec: PafRecord) -> tuple[int, int]:
    if rec.strand == "+":
        return rec.tstart, rec.tend
    return rec.tlen - rec.tend, rec.tlen - rec.tstart


def chain_records(records: list[PafRecord], merge_gap: int) -> list[Chain]:
    grouped: dict[tuple[str, str, str], list[tuple[PafRecord, int, int]]] = defaultdict(list)
    for rec in records:
        ot_start, ot_end = oriented_target(rec)
        grouped[(rec.qname, rec.tname, rec.strand)].append((rec, ot_start, ot_end))

    output: list[Chain] = []
    for group in grouped.values():
        group.sort(key=lambda item: (item[0].qstart, item[1]))
        current: Chain | None = None
        for rec, ot_start, ot_end in group:
            if current is None:
                current = Chain(
                    rec.qname, rec.qlen, rec.qstart, rec.qend,
                    rec.tname, rec.tlen, rec.tstart, rec.tend, rec.strand,
                    rec.aln_len, rec.matches, 1, rec.qend, ot_end,
                )
                continue
            q_gap = rec.qstart - current.oq_end
            t_gap = ot_start - current.ot_end
            if q_gap <= merge_gap and t_gap <= merge_gap and q_gap >= -merge_gap and t_gap >= -merge_gap:
                current.qstart = min(current.qstart, rec.qstart)
                current.qend = max(current.qend, rec.qend)
                current.tstart = min(current.tstart, rec.tstart)
                current.tend = max(current.tend, rec.tend)
                current.aligned_bp += rec.aln_len
                current.weighted_matches += rec.matches
                current.n_alignments += 1
                current.oq_end = max(current.oq_end, rec.qend)
                current.ot_end = max(current.ot_end, ot_end)
            else:
                output.append(current)
                current = Chain(
                    rec.qname, rec.qlen, rec.qstart, rec.qend,
                    rec.tname, rec.tlen, rec.tstart, rec.tend, rec.strand,
                    rec.aln_len, rec.matches, 1, rec.qend, ot_end,
                )
        if current is not None:
            output.append(current)
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--paf", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--min-aln", type=int, default=10_000)
    parser.add_argument("--min-identity", type=float, default=0.95)
    parser.add_argument("--end-window", type=int, default=10_000_000)
    parser.add_argument("--merge-gap", type=int, default=100_000)
    args = parser.parse_args()

    passing = [
        rec for rec in parse_paf(args.paf)
        if rec.aln_len >= args.min_aln and rec.identity >= args.min_identity
    ]
    chains = chain_records(passing, args.merge_gap)
    terminal = []
    for chain in chains:
        q_side = terminal_side(chain.qstart, chain.qend, chain.qlen, args.end_window)
        t_side = terminal_side(chain.tstart, chain.tend, chain.tlen, args.end_window)
        if q_side != "internal" and t_side != "internal":
            terminal.append((chain, q_side, t_side))
    terminal.sort(key=lambda item: (item[0].aligned_bp, item[0].identity), reverse=True)

    fields = [
        "candidate_id", "y_name", "y_start", "y_end", "y_length", "y_side",
        "x_name", "x_start", "x_end", "x_length", "x_side", "strand",
        "aligned_bp", "mean_identity", "n_alignments", "status",
    ]
    with open(args.out, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for index, (chain, q_side, t_side) in enumerate(terminal, 1):
            writer.writerow({
                "candidate_id": f"PAR{index}",
                "y_name": chain.qname,
                "y_start": chain.qstart,
                "y_end": chain.qend,
                "y_length": chain.qlen,
                "y_side": q_side,
                "x_name": chain.tname,
                "x_start": chain.tstart,
                "x_end": chain.tend,
                "x_length": chain.tlen,
                "x_side": t_side,
                "strand": chain.strand,
                "aligned_bp": chain.aligned_bp,
                "mean_identity": f"{chain.identity:.6f}",
                "n_alignments": chain.n_alignments,
                "status": "best_terminal_candidate" if index == 1 else "secondary_terminal_candidate",
            })


if __name__ == "__main__":
    main()
