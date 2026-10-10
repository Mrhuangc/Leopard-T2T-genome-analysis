"""Keep only complete Exonerate alignment blocks after a controlled timeout."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


def _complete_block(lines: list[str]) -> bool:
    return (
        any(line.startswith("vulgar:") for line in lines)
        and any(line.startswith("# --- START OF GFF DUMP ---") for line in lines)
        and any(line.startswith("# --- END OF GFF DUMP ---") for line in lines)
    )


def finalize_chunk(
    raw_path: Path,
    output_path: Path,
    summary_path: Path,
    exit_code: int,
    timeout_seconds: int,
    exhaustive: bool = True,
    suboptimal: bool = False,
    timeout_is_fatal: bool = True,
    elapsed_seconds: int | None = None,
    bestn: int | None = None,
) -> dict[str, Any]:
    """Remove a possibly truncated final alignment while retaining completed blocks."""

    if exit_code not in {0, 124}:
        raise ValueError(f"unsupported Exonerate exit code: {exit_code}")
    if timeout_seconds < 1:
        raise ValueError("timeout_seconds must be positive")
    if elapsed_seconds is not None and elapsed_seconds < 0:
        raise ValueError("elapsed_seconds must be non-negative")
    if bestn is not None and bestn < 1:
        raise ValueError("bestn must be positive when specified")

    preamble: list[str] = []
    current: list[str] = []
    retained_blocks: list[list[str]] = []
    raw_blocks = 0
    global_completion = False

    with raw_path.open("rt", encoding="utf-8", errors="replace") as source:
        for line in source:
            global_completion |= "completed exonerate analysis" in line
            if line.startswith("C4 Alignment:"):
                if current and _complete_block(current):
                    retained_blocks.append(current)
                current = [line]
                raw_blocks += 1
            elif current:
                current.append(line)
            else:
                preamble.append(line)

    if exit_code == 0:
        if not global_completion:
            raise ValueError("Exonerate exited zero without its completion marker")
        if current and _complete_block(current):
            retained_blocks.append(current)

    vulgar_lines = 0
    gff_dumps = 0
    with output_path.open("wt", encoding="utf-8", newline="\n") as output:
        output.writelines(preamble)
        for block in retained_blocks:
            vulgar_lines += sum(line.startswith("vulgar:") for line in block)
            gff_dumps += sum(
                line.startswith("# --- START OF GFF DUMP ---") for line in block
            )
            output.writelines(block)
        output.write(
            "# FELIDAE_OR_CHUNK_STATUS "
            f"{'completed' if exit_code == 0 else 'soft_timeout'}\n"
        )

    summary: dict[str, Any] = {
        "schema_version": "1.2",
        "status": "pass",
        "termination": "completed" if exit_code == 0 else "soft_timeout",
        "command_exit_code": exit_code,
        "timeout_seconds": timeout_seconds,
        "elapsed_seconds": elapsed_seconds,
        "output_file": output_path.name,
        "search_parameters": {
            "model": "protein2genome",
            "exhaustive": exhaustive,
            "suboptimal": suboptimal,
            "bestn": bestn,
            "timeout_is_fatal": timeout_is_fatal,
        },
        "counts": {
            "raw_alignment_blocks": raw_blocks,
            "retained_complete_alignment_blocks": len(retained_blocks),
            "dropped_incomplete_alignment_blocks": raw_blocks - len(retained_blocks),
            "retained_vulgar_alignments": vulgar_lines,
            "retained_gff_dumps": gff_dumps,
            "raw_bytes": raw_path.stat().st_size,
            "retained_bytes": output_path.stat().st_size,
        },
    }
    with summary_path.open("wt", encoding="utf-8", newline="\n") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("raw", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--exit-code", required=True, type=int)
    parser.add_argument("--timeout-seconds", required=True, type=int)
    parser.add_argument("--exhaustive", required=True, choices=("TRUE", "FALSE"))
    parser.add_argument("--suboptimal", required=True, choices=("yes", "no"))
    parser.add_argument("--timeout-is-fatal", required=True, choices=("true", "false"))
    parser.add_argument("--elapsed-seconds", required=True, type=int)
    parser.add_argument("--bestn", type=int)
    args = parser.parse_args()
    finalize_chunk(
        args.raw,
        args.output,
        args.summary,
        args.exit_code,
        args.timeout_seconds,
        args.exhaustive == "TRUE",
        args.suboptimal == "yes",
        args.timeout_is_fatal == "true",
        args.elapsed_seconds,
        args.bestn,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
