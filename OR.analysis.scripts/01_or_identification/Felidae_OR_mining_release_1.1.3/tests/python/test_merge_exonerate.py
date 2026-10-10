import json
from pathlib import Path

import pytest

from or_mining.merge_exonerate import merge_outputs


def _write_cell(
    tmp_path: Path, stem: str, *, suboptimal: bool, alignment: bool = True
) -> tuple[Path, Path]:
    result = tmp_path / f"{stem}.exonerate.txt"
    result.write_text(
        (
            "C4 Alignment:\n"
            "vulgar: q 0 10 + t 0 30 + 50 M 10 30\n"
            "# --- START OF GFF DUMP ---\n"
            "# --- END OF GFF DUMP ---\n"
            if alignment
            else ""
        )
        + "# FELIDAE_OR_CHUNK_STATUS completed\n",
        encoding="utf-8",
    )
    status = tmp_path / f"{stem}.exonerate.summary.json"
    status.write_text(
        json.dumps(
            {
                "status": "pass",
                "termination": "completed",
                "output_file": result.name,
                "search_parameters": {
                    "model": "protein2genome",
                    "exhaustive": True,
                    "suboptimal": suboptimal,
                    "bestn": None,
                    "timeout_is_fatal": True,
                },
                "counts": {
                    "retained_complete_alignment_blocks": int(alignment),
                    "dropped_incomplete_alignment_blocks": 0,
                },
            }
        ),
        encoding="utf-8",
    )
    return result, status


def test_merge_records_consistent_search_parameters(tmp_path: Path) -> None:
    result, status = _write_cell(tmp_path, "cell", suboptimal=False)
    summary = merge_outputs(
        [result],
        [status],
        tmp_path / "merged.txt",
        tmp_path / "merged.json",
    )
    assert summary["search_parameters"]["suboptimal"] is False
    assert summary["search_parameters"]["bestn"] is None


def test_merge_rejects_mixed_search_parameters(tmp_path: Path) -> None:
    result_a, status_a = _write_cell(tmp_path, "a", suboptimal=False)
    result_b, status_b = _write_cell(tmp_path, "b", suboptimal=True)
    with pytest.raises(ValueError, match="inconsistent search parameters"):
        merge_outputs(
            [result_a, result_b],
            [status_a, status_b],
            tmp_path / "merged.txt",
            tmp_path / "merged.json",
        )


def test_merge_rejects_legacy_status_without_search_parameters(tmp_path: Path) -> None:
    result, status = _write_cell(tmp_path, "legacy", suboptimal=False)
    record = json.loads(status.read_text(encoding="utf-8"))
    del record["search_parameters"]
    status.write_text(json.dumps(record), encoding="utf-8")
    with pytest.raises(ValueError, match="missing auditable"):
        merge_outputs(
            [result],
            [status],
            tmp_path / "merged.txt",
            tmp_path / "merged.json",
        )


def test_merge_rejects_completed_but_empty_search(tmp_path: Path) -> None:
    result, status = _write_cell(tmp_path, "empty", suboptimal=False, alignment=False)
    with pytest.raises(ValueError, match="no biological alignment"):
        merge_outputs(
            [result],
            [status],
            tmp_path / "merged.txt",
            tmp_path / "merged.json",
        )
