"""CLI for generating, validating, and gating Exonerate smoke-test artifacts."""

from __future__ import annotations

import argparse
from pathlib import Path

from or_mining.exonerate_smoke import check_smoke_gate, validate_smoke_result, write_smoke_inputs


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    generate = subparsers.add_parser("generate")
    generate.add_argument("--query", required=True, type=Path)
    generate.add_argument("--target", required=True, type=Path)
    generate.add_argument("--summary", required=True, type=Path)

    validate = subparsers.add_parser("validate")
    validate.add_argument("result", type=Path)
    validate.add_argument("chunk_summary", type=Path)
    validate.add_argument("--output", required=True, type=Path)

    gate = subparsers.add_parser("gate")
    gate.add_argument("summary", type=Path)
    gate.add_argument("--output", required=True, type=Path)

    args = parser.parse_args()
    if args.command == "generate":
        write_smoke_inputs(args.query, args.target, args.summary)
    elif args.command == "validate":
        validate_smoke_result(args.result, args.chunk_summary, args.output)
    else:
        check_smoke_gate(args.summary, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
