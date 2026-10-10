import csv
import json
from pathlib import Path

from or_mining.strict_refinement import (
    STRICT_SEARCH_PARAMETERS,
    merge_strict_results,
    prepare_strict_pairs,
    validate_strict_pilot,
)


def test_prepare_strict_pairs_uses_one_query_and_local_window_per_candidate(
    tmp_path: Path,
) -> None:
    reference = tmp_path / "reference.fa"
    species = tmp_path / "species.fa"
    targets = tmp_path / "targets.fa"
    regions = tmp_path / "regions.tsv"
    candidates = tmp_path / "candidates.tsv"
    reference.write_text(">q1\n" + "M" * 300 + "\n", encoding="utf-8")
    species.write_text(">q2\n" + "A" * 320 + "\n", encoding="utf-8")
    targets.write_text(">r1\n" + "A" * 1000 + "\n", encoding="utf-8")
    regions.write_text(
        "region_id\tsseqid\texpanded_start_1\texpanded_end_1\n"
        "r1\tchr1\t1001\t2000\n",
        encoding="utf-8",
    )
    candidates.write_text(
        "candidate_id\tregion_id\tquery_id\tsseqid\tstrand\tgenomic_start_1\t"
        "genomic_end_1\trefinement_start_1\trefinement_end_1\tquery_span_aa\t"
        "exonerate_score\n"
        "c1\tr1\tq1\tchr1\t+\t1301\t1500\t1201\t1600\t250\t500\n"
        "c2\tr1\tq2\tchr1\t-\t1701\t1800\t1601\t1900\t300\t600\n",
        encoding="utf-8",
    )
    manifest = tmp_path / "manifest.tsv"
    summary = prepare_strict_pairs(
        reference,
        species,
        targets,
        regions,
        candidates,
        tmp_path / "pairs",
        manifest,
        tmp_path / "summary.json",
    )
    with manifest.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["candidate_id"] for row in rows] == ["c2", "c1"]
    assert rows[0]["target_bases"] == "300"
    assert len(list((tmp_path / "pairs").glob("*.query.prot.fa"))) == 2
    assert len(list((tmp_path / "pairs").glob("*.target.nucl.fa"))) == 2
    assert summary["counts"]["strict_pairs"] == 2


def _strict_cell(tmp_path: Path, stem: str) -> tuple[Path, Path]:
    result = tmp_path / f"{stem}.exonerate.txt"
    result.write_text(
        "C4 Alignment:\n"
        "vulgar: q 0 10 . t 0 30 + 50 M 10 30\n"
        "# --- START OF GFF DUMP ---\n"
        "# --- END OF GFF DUMP ---\n",
        encoding="utf-8",
    )
    status = tmp_path / f"{stem}.exonerate.summary.json"
    status.write_text(
        json.dumps(
            {
                "status": "pass",
                "termination": "completed",
                "elapsed_seconds": 1,
                "output_file": result.name,
                "search_parameters": STRICT_SEARCH_PARAMETERS,
                "counts": {
                    "retained_complete_alignment_blocks": 1,
                    "retained_vulgar_alignments": 1,
                    "retained_gff_dumps": 1,
                    "dropped_incomplete_alignment_blocks": 0,
                },
            }
        ),
        encoding="utf-8",
    )
    return result, status


def test_strict_pilot_and_merge_require_one_bestn1_alignment_per_candidate(
    tmp_path: Path,
) -> None:
    assert STRICT_SEARCH_PARAMETERS["exhaustive"] is False
    first_result, first_status = _strict_cell(tmp_path, "a")
    second_result, second_status = _strict_cell(tmp_path, "b")
    pilot = validate_strict_pilot(
        [first_result, second_result],
        [first_status, second_status],
        tmp_path / "pilot.json",
        tmp_path / "gate.ok",
        expected_pairs=2,
    )
    assert pilot["counts"]["positive_pairs"] == 2
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "strict_pair_id\tcandidate_id\n"
        "strict_1\tc1\n"
        "strict_2\tc2\n",
        encoding="utf-8",
    )
    merged = merge_strict_results(
        [first_result, second_result],
        [first_status, second_status],
        manifest,
        tmp_path / "merged.txt",
        tmp_path / "merged.json",
    )
    assert merged["status"] == "pass"
    assert merged["counts"]["expected_candidates"] == 2
    assert merged["counts"]["positive_candidates"] == 2


def test_strict_pilot_rejects_a_zero_alignment_task(tmp_path: Path) -> None:
    result, status = _strict_cell(tmp_path, "empty")
    record = json.loads(status.read_text(encoding="utf-8"))
    record["counts"]["retained_complete_alignment_blocks"] = 0
    status.write_text(json.dumps(record), encoding="utf-8")

    try:
        validate_strict_pilot(
            [result],
            [status],
            tmp_path / "pilot.json",
            tmp_path / "gate.ok",
            expected_pairs=1,
        )
    except ValueError as error:
        assert "retained 0 alignment blocks" in str(error)
    else:
        raise AssertionError("a zero-alignment pilot task must fail the gate")

    assert not (tmp_path / "gate.ok").exists()
