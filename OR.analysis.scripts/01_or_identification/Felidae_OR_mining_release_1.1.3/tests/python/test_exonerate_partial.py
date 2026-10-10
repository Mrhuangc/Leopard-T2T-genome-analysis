from pathlib import Path

from or_mining.finalize_exonerate_chunk import finalize_chunk


def test_soft_timeout_drops_only_trailing_incomplete_alignment(tmp_path: Path) -> None:
    raw = tmp_path / "raw.txt"
    raw.write_text(
        "Command line: exonerate\n"
        "C4 Alignment:\n"
        "vulgar: q1 0 10 + t1 0 30 + 50 M 10 30\n"
        "# --- START OF GFF DUMP ---\n"
        "t1\texonerate\tgene\t1\t30\t.\t+\t.\tID=g1\n"
        "# --- END OF GFF DUMP ---\n"
        "#\nATG\n--\n"
        "C4 Alignment:\n"
        "vulgar: q2 0 10 + t2 0 30 + 50 M 10 30\n",
        encoding="utf-8",
    )
    output = tmp_path / "clean.txt"
    summary = finalize_chunk(raw, output, tmp_path / "summary.json", 124, 21600)
    text = output.read_text(encoding="utf-8")
    assert "vulgar: q1" in text
    assert "vulgar: q2" not in text
    assert summary["termination"] == "soft_timeout"
    assert summary["search_parameters"]["exhaustive"] is True
    assert summary["search_parameters"]["suboptimal"] is False
    assert summary["search_parameters"]["bestn"] is None
    assert summary["search_parameters"]["timeout_is_fatal"] is True
    assert summary["counts"]["retained_complete_alignment_blocks"] == 1
    assert summary["counts"]["dropped_incomplete_alignment_blocks"] == 1


def test_completed_run_retains_final_complete_alignment(tmp_path: Path) -> None:
    raw = tmp_path / "raw.txt"
    raw.write_text(
        "C4 Alignment:\n"
        "vulgar: q1 0 10 + t1 0 30 + 50 M 10 30\n"
        "# --- START OF GFF DUMP ---\n"
        "t1\texonerate\tgene\t1\t30\t.\t+\t.\tID=g1\n"
        "# --- END OF GFF DUMP ---\n"
        "-- completed exonerate analysis\n",
        encoding="utf-8",
    )
    output = tmp_path / "clean.txt"
    summary = finalize_chunk(raw, output, tmp_path / "summary.json", 0, 21600)
    assert "vulgar: q1" in output.read_text(encoding="utf-8")
    assert summary["termination"] == "completed"
    assert summary["counts"]["retained_complete_alignment_blocks"] == 1
    assert summary["counts"]["dropped_incomplete_alignment_blocks"] == 0
