"""Tests for the merged TBLASTN data contract."""

from pathlib import Path

import pytest

from or_mining.blast import BlastFormatError, parse_blast_row, validate_merged_file


PLUS = (
    "query1\tChrA1\t25.185\t270\t175\t9\t68\t321\t100\t900\t4.60e-06\t51.2"
    "\t0\t1\t333\t1000\n"
)
MINUS = (
    "query2\tChrA1\t26.636\t214\t125\t6\t105\t301\t900\t300\t3.38e-07\t54.3"
    "\t0\t-3\t309\t1000\n"
)


def test_parse_coordinates_on_both_strands() -> None:
    plus = parse_blast_row(PLUS)
    minus = parse_blast_row(MINUS)

    assert (plus.strand, plus.genome_start_0, plus.genome_end_0) == ("+", 99, 900)
    assert (minus.strand, minus.genome_start_0, minus.genome_end_0) == ("-", 299, 900)


def test_rejects_strand_coordinate_disagreement() -> None:
    bad = PLUS.replace("\t100\t900\t", "\t900\t100\t")

    with pytest.raises(BlastFormatError, match="positive sframe"):
        parse_blast_row(bad)


def test_validation_report_counts_rows(tmp_path: Path) -> None:
    merged = tmp_path / "merged.tsv"
    merged.write_text(PLUS + MINUS, encoding="utf-8")

    report = validate_merged_file(merged)

    assert report["status"] == "pass"
    assert report["counts"] == {
        "total_rows": 2,
        "valid_rows": 2,
        "invalid_rows": 0,
        "plus_strand_rows": 1,
        "minus_strand_rows": 1,
        "unique_query_ids": 2,
        "unique_subject_ids": 1,
    }


def test_validation_report_records_bad_column_count(tmp_path: Path) -> None:
    merged = tmp_path / "bad.tsv"
    merged.write_text("too\tfew\tcolumns\n", encoding="utf-8")

    report = validate_merged_file(merged)

    assert report["status"] == "fail"
    assert report["counts"]["invalid_rows"] == 1
    assert report["errors"][0]["line"] == 1

