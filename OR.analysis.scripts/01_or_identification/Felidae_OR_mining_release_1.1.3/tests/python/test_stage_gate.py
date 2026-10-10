import json
from pathlib import Path

import pytest
import yaml

from or_mining.stage_gate import SUMMARY_PATTERNS, check_stage


def _params(tmp_path: Path, species: str = "Ppa") -> Path:
    path = tmp_path / "params.yml"
    path.write_text(
        yaml.safe_dump({"species": species, "outdir": str(tmp_path / "results")}),
        encoding="utf-8",
    )
    return path


def test_stage_gate_accepts_warning_summary(tmp_path: Path) -> None:
    params = _params(tmp_path)
    summary = tmp_path / "results" / SUMMARY_PATTERNS[18].format(s="Ppa")
    summary.parent.mkdir(parents=True)
    summary.write_text(json.dumps({"status": "pass_with_warnings"}), encoding="utf-8")
    actual, payload = check_stage(params, 18)
    assert actual == summary
    assert payload["status"] == "pass_with_warnings"


def test_exonerate_gate_rejects_partial_search(tmp_path: Path) -> None:
    params = _params(tmp_path)
    summary = tmp_path / "results" / SUMMARY_PATTERNS[10].format(s="Ppa")
    summary.parent.mkdir(parents=True)
    summary.write_text(
        json.dumps(
            {
                "status": "pass",
                "search_completeness": "partial",
                "counts": {"soft_timed_out_chunks": 1},
            }
        ),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="not complete"):
        check_stage(params, 10)


def test_initial_search_gate_requires_all_shards_and_valid_rows(tmp_path: Path) -> None:
    params = _params(tmp_path)
    summary = tmp_path / "results" / SUMMARY_PATTERNS[0].format(s="Ppa")
    summary.parent.mkdir(parents=True)
    summary.write_text(
        json.dumps(
            {
                "status": "pass",
                "search_completeness": "complete",
                "counts": {"input_query_shards": 2, "completed_query_shards": 2},
                "checks": {
                    "invalid_rows": 0,
                    "all_query_shards_completed": True,
                    "query_manifest_matches_parts": True,
                    "merged_output_nonempty": True,
                    "canonical_16_column_contract": True,
                },
            }
        ),
        encoding="utf-8",
    )
    actual, _ = check_stage(params, 0)
    assert actual == summary

    payload = json.loads(summary.read_text())
    payload["counts"]["completed_query_shards"] = 1
    summary.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="every initial TBLASTN query shard"):
        check_stage(params, 0)
