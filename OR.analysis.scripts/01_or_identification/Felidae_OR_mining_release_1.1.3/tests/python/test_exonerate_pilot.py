import json
from pathlib import Path

import pytest

from or_mining.exonerate_pilot import validate_pilot


def _pilot_record(
    path: Path, result_name: str, termination: str, alignments: int
) -> None:
    path.write_text(
        json.dumps(
            {
                "status": "pass",
                "termination": termination,
                "output_file": result_name,
                "search_parameters": {
                    "model": "protein2genome",
                    "exhaustive": False,
                    "suboptimal": False,
                    "bestn": None,
                    "timeout_is_fatal": True,
                },
                "counts": {
                    "retained_complete_alignment_blocks": alignments,
                    "retained_vulgar_alignments": alignments,
                    "retained_gff_dumps": alignments,
                },
            }
        ),
        encoding="utf-8",
    )


def test_real_data_pilot_opens_gate_only_with_complete_output(tmp_path: Path) -> None:
    results = [tmp_path / f"r{i}.txt" for i in range(2)]
    statuses = [tmp_path / f"s{i}.json" for i in range(2)]
    for result in results:
        result.write_text("alignment\n", encoding="utf-8")
    _pilot_record(statuses[0], results[0].name, "completed", 1)
    _pilot_record(statuses[1], results[1].name, "completed", 0)
    summary = validate_pilot(
        results,
        statuses,
        tmp_path / "pilot.json",
        tmp_path / "gate.ok",
        expected_pairs=2,
    )
    assert summary["status"] == "pass"
    assert (tmp_path / "gate.ok").read_text(encoding="utf-8") == "status=pass\n"


def test_real_data_pilot_rejects_any_timeout(tmp_path: Path) -> None:
    result = tmp_path / "result.txt"
    status = tmp_path / "status.json"
    result.write_text("", encoding="utf-8")
    _pilot_record(status, result.name, "soft_timeout", 0)
    with pytest.raises(ValueError, match="timed out"):
        validate_pilot(
            [result],
            [status],
            tmp_path / "pilot.json",
            tmp_path / "gate.ok",
            expected_pairs=1,
        )


def test_real_data_pilot_requires_multiple_positive_pairs(tmp_path: Path) -> None:
    results = [tmp_path / f"r{i}.txt" for i in range(2)]
    statuses = [tmp_path / f"s{i}.json" for i in range(2)]
    for result in results:
        result.write_text("alignment\n", encoding="utf-8")
    _pilot_record(statuses[0], results[0].name, "completed", 1)
    _pilot_record(statuses[1], results[1].name, "completed", 0)
    with pytest.raises(ValueError, match="positive pairs"):
        validate_pilot(
            results,
            statuses,
            tmp_path / "pilot.json",
            tmp_path / "gate.ok",
            expected_pairs=2,
            minimum_positive_pairs=2,
        )
