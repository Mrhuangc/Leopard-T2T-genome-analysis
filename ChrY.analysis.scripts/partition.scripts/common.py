#!/usr/bin/env python3
"""Shared parsers and coordinate helpers for the leopard Y workflow."""

from __future__ import annotations

import gzip
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, TextIO

import yaml


def load_config(path: str | Path) -> dict:
    with open(path, "rt", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    if not isinstance(data, dict):
        raise ValueError(f"Configuration is not a mapping: {path}")
    return data


def open_text(path: str | Path) -> TextIO:
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "rt", encoding="utf-8")


def is_real_path(value: object) -> bool:
    if value is None:
        return False
    text = str(value).strip()
    if not text or text.startswith("/path/to/"):
        return False
    return Path(text).exists()


def read_fai(path: str | Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    with open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            lengths[fields[0]] = int(fields[1])
    return lengths


@dataclass(frozen=True)
class PafRecord:
    qname: str
    qlen: int
    qstart: int
    qend: int
    strand: str
    tname: str
    tlen: int
    tstart: int
    tend: int
    matches: int
    aln_len: int
    mapq: int

    @property
    def identity(self) -> float:
        return self.matches / self.aln_len if self.aln_len else math.nan


def parse_paf(path: str | Path) -> Iterator[PafRecord]:
    with open_text(path) as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise ValueError(f"PAF line {line_no} has fewer than 12 columns: {path}")
            yield PafRecord(
                qname=fields[0],
                qlen=int(fields[1]),
                qstart=int(fields[2]),
                qend=int(fields[3]),
                strand=fields[4],
                tname=fields[5],
                tlen=int(fields[6]),
                tstart=int(fields[7]),
                tend=int(fields[8]),
                matches=int(fields[9]),
                aln_len=int(fields[10]),
                mapq=int(fields[11]),
            )


def merge_intervals(intervals: Iterable[tuple[int, int]], gap: int = 0) -> list[tuple[int, int]]:
    cleaned = sorted((max(0, int(s)), max(0, int(e))) for s, e in intervals if e > s)
    if not cleaned:
        return []
    merged: list[list[int]] = [[cleaned[0][0], cleaned[0][1]]]
    for start, end in cleaned[1:]:
        if start <= merged[-1][1] + gap:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(start, end) for start, end in merged]


def interval_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def add_interval_coverage(
    buckets: list[list[tuple[int, int]]],
    start: int,
    end: int,
    window: int,
) -> None:
    """Append clipped intervals to per-window buckets using 0-based half-open coordinates."""
    if end <= start or not buckets:
        return
    first = max(0, start // window)
    last = min(len(buckets) - 1, (end - 1) // window)
    for index in range(first, last + 1):
        w_start = index * window
        w_end = w_start + window
        clipped_start = max(start, w_start)
        clipped_end = min(end, w_end)
        if clipped_end > clipped_start:
            buckets[index].append((clipped_start, clipped_end))


def covered_bases(intervals: Iterable[tuple[int, int]]) -> int:
    return sum(end - start for start, end in merge_intervals(intervals))
