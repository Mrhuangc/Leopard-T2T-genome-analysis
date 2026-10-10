#!/usr/bin/env python3
"""Parse a multi-segment Ppa/Ppa_or minimap2 PAF into structural events."""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path


TOKEN_RE = re.compile(
    r":\d+|=[A-Za-z]+|\*[A-Za-z]{2}|\+[A-Za-z]+|-[A-Za-z]+|~[A-Za-z]{2}\d+[A-Za-z]{2}"
)
REGION_RE = re.compile(r"^[^|]+\|[^:]+:(\d+)-(\d+)\|")


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("paf", type=Path)
    parser.add_argument("events_tsv", type=Path)
    parser.add_argument("summary_json", type=Path)
    parser.add_argument("--minimum-event-bp", type=int, default=50)
    parser.add_argument("--minimum-mapq", type=int, default=30)
    return parser.parse_args()


def parse_region(name: str) -> tuple[int, int]:
    match = REGION_RE.match(name)
    if not match:
        raise ValueError(f"source interval is absent from FASTA header: {name}")
    return int(match.group(1)), int(match.group(2))


def parse_tags(fields: list[str]) -> dict[str, str]:
    tags: dict[str, str] = {}
    for tag in fields[12:]:
        if len(tag) >= 6 and tag[2:5] in {":A:", ":i:", ":f:", ":Z:"}:
            tags[tag[:2]] = tag[5:]
    return tags


def interval(start0: int, length: int, source_start1: int) -> tuple[str | int, str | int]:
    if length == 0:
        return "", ""
    return source_start1 + start0, source_start1 + start0 + length - 1


def main() -> int:
    args = arguments()
    if args.minimum_event_bp < 1:
        raise ValueError("minimum event size must be positive")

    records: list[dict[str, object]] = []
    with args.paf.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise ValueError(f"PAF line {line_number} has fewer than 12 columns")
            tags = parse_tags(fields)
            records.append(
                {
                    "line": line_number,
                    "qname": fields[0],
                    "qlen": int(fields[1]),
                    "qstart": int(fields[2]),
                    "qend": int(fields[3]),
                    "strand": fields[4],
                    "tname": fields[5],
                    "tlen": int(fields[6]),
                    "tstart": int(fields[7]),
                    "tend": int(fields[8]),
                    "mapq": int(fields[11]),
                    "tags": tags,
                }
            )

    if not records:
        raise ValueError("PAF is empty")
    qnames = {str(record["qname"]) for record in records}
    tnames = {str(record["tname"]) for record in records}
    if len(qnames) != 1 or len(tnames) != 1:
        raise ValueError("expected exactly one query region and one target region")

    q_source_start, q_source_end = parse_region(next(iter(qnames)))
    t_source_start, t_source_end = parse_region(next(iter(tnames)))
    qlen = int(records[0]["qlen"])
    tlen = int(records[0]["tlen"])
    if q_source_end - q_source_start + 1 != qlen:
        raise ValueError("query FASTA header interval disagrees with PAF query length")
    if t_source_end - t_source_start + 1 != tlen:
        raise ValueError("target FASTA header interval disagrees with PAF target length")

    forward = [
        record
        for record in records
        if record["strand"] == "+"
        and int(record["mapq"]) >= args.minimum_mapq
        and "cg" in record["tags"]
        and "cs" in record["tags"]
    ]
    forward.sort(key=lambda record: (int(record["qstart"]), int(record["tstart"])))
    if not forward:
        raise ValueError("no qualified plus-strand alignment segments with cg and cs tags")

    for previous, current in zip(forward, forward[1:]):
        if int(current["qstart"]) < int(previous["qend"]):
            raise ValueError("qualified plus-strand segments overlap on Ppa")
        if int(current["tstart"]) < int(previous["tend"]):
            raise ValueError("qualified plus-strand segments are not collinear on Ppa_or")

    events: list[dict[str, object]] = []
    ppa_extra_all = 0
    ppa_or_extra_all = 0
    substitutions = 0

    def add_event(
        segment: str,
        event_type: str,
        qstart0: int,
        qbp: int,
        tstart0: int,
        tbp: int,
    ) -> None:
        if max(qbp, tbp) < args.minimum_event_bp:
            return
        qabs_start, qabs_end = interval(qstart0, qbp, q_source_start)
        tabs_start, tabs_end = interval(tstart0, tbp, t_source_start)
        events.append(
            {
                "segment": segment,
                "event_type": event_type,
                "Ppa_length_bp": qbp,
                "Ppa_or_length_bp": tbp,
                "net_Ppa_minus_Ppa_or_bp": qbp - tbp,
                "Ppa_start_1": qabs_start,
                "Ppa_end_1": qabs_end,
                "Ppa_or_start_1": tabs_start,
                "Ppa_or_end_1": tabs_end,
                "Ppa_breakpoint_1": q_source_start + qstart0,
                "Ppa_or_breakpoint_1": t_source_start + tstart0,
            }
        )

    for segment_index, record in enumerate(forward, start=1):
        qpos = int(record["qstart"])
        tpos = int(record["tstart"])
        cs = str(record["tags"]["cs"])
        tokens = TOKEN_RE.findall(cs)
        if "".join(tokens) != cs:
            raise ValueError(f"unsupported cs operation in PAF line {record['line']}")
        for token in tokens:
            operation = token[0]
            if operation == ":":
                length = int(token[1:])
                qpos += length
                tpos += length
            elif operation == "=":
                length = len(token) - 1
                qpos += length
                tpos += length
            elif operation == "*":
                substitutions += 1
                qpos += 1
                tpos += 1
            elif operation == "+":
                length = len(token) - 1
                ppa_extra_all += length
                add_event(str(segment_index), "Ppa_insertion", qpos, length, tpos, 0)
                qpos += length
            elif operation == "-":
                length = len(token) - 1
                ppa_or_extra_all += length
                add_event(str(segment_index), "Ppa_or_insertion", qpos, 0, tpos, length)
                tpos += length
            elif operation == "~":
                match = re.search(r"\d+", token)
                if not match:
                    raise ValueError(f"malformed cs long-gap operation: {token}")
                length = int(match.group())
                ppa_or_extra_all += length
                add_event(str(segment_index), "Ppa_or_long_gap", qpos, 0, tpos, length)
                tpos += length
        if qpos != int(record["qend"]) or tpos != int(record["tend"]):
            raise ValueError(
                f"cs consumption disagrees with PAF line {record['line']}: "
                f"query {qpos}/{record['qend']}, target {tpos}/{record['tend']}"
            )

    for gap_index, (previous, current) in enumerate(zip(forward, forward[1:]), start=1):
        qstart = int(previous["qend"])
        tstart = int(previous["tend"])
        qgap = int(current["qstart"]) - qstart
        tgap = int(current["tstart"]) - tstart
        ppa_extra_all += qgap
        ppa_or_extra_all += tgap
        add_event(f"gap_{gap_index}", "unaligned_chain_gap", qstart, qgap, tstart, tgap)

    fieldnames = [
        "segment",
        "event_type",
        "Ppa_length_bp",
        "Ppa_or_length_bp",
        "net_Ppa_minus_Ppa_or_bp",
        "Ppa_start_1",
        "Ppa_end_1",
        "Ppa_or_start_1",
        "Ppa_or_end_1",
        "Ppa_breakpoint_1",
        "Ppa_or_breakpoint_1",
    ]
    args.events_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.events_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(events)

    forward_q_aligned = sum(int(record["qend"]) - int(record["qstart"]) for record in forward)
    forward_t_aligned = sum(int(record["tend"]) - int(record["tstart"]) for record in forward)
    reverse = [record for record in records if record["strand"] == "-"]
    summary = {
        "schema_version": "2.0",
        "status": "pass",
        "method": "collinear plus-strand cg/cs chain; reverse overlaps audited separately",
        "coordinate_system": "reported genomic intervals are 1-based inclusive",
        "parameters": {
            "minimum_event_bp": args.minimum_event_bp,
            "minimum_mapq": args.minimum_mapq,
        },
        "counts": {
            "input_paf_records": len(records),
            "qualified_forward_segments": len(forward),
            "reverse_segments": len(reverse),
            "reported_structural_events": len(events),
            "Ppa_extra_bp_including_chain_gaps": ppa_extra_all,
            "Ppa_or_extra_bp_including_chain_gaps": ppa_or_extra_all,
            "net_Ppa_minus_Ppa_or_bp_within_chain": ppa_extra_all - ppa_or_extra_all,
            "single_base_substitutions_within_forward_segments": substitutions,
            "Ppa_forward_aligned_bp": forward_q_aligned,
            "Ppa_or_forward_aligned_bp": forward_t_aligned,
        },
        "fractions": {
            "Ppa_forward_aligned_fraction": forward_q_aligned / qlen,
            "Ppa_or_forward_aligned_fraction": forward_t_aligned / tlen,
        },
        "reverse_segment_spans": [
            {
                "Ppa_bp": int(record["qend"]) - int(record["qstart"]),
                "Ppa_or_bp": int(record["tend"]) - int(record["tstart"]),
                "mapq": int(record["mapq"]),
            }
            for record in reverse
        ],
    }
    args.summary_json.parent.mkdir(parents=True, exist_ok=True)
    args.summary_json.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
