import csv
from pathlib import Path

from or_mining.exonerate_model_selection import select_exonerate_candidates


def _gff_block(
    query: str,
    target: str,
    qstart: int,
    qend: int,
    exon_intervals: list[tuple[int, int]],
    score: int,
) -> str:
    start = min(interval[0] for interval in exon_intervals)
    end = max(interval[1] for interval in exon_intervals)
    target_bases = sum(right - left + 1 for left, right in exon_intervals)
    operations: list[str] = []
    query_left = qend - qstart
    if len(exon_intervals) == 1:
        operations = ["M", str(query_left), str(target_bases)]
    else:
        intron_start = exon_intervals[0][1] + 1
        intron_end = exon_intervals[1][0] - 1
        left_aa = (exon_intervals[0][1] - exon_intervals[0][0] + 1) // 3
        right_aa = query_left - left_aa
        operations = [
            "M",
            str(left_aa),
            str(exon_intervals[0][1] - exon_intervals[0][0] + 1),
            "I",
            "0",
            str(intron_end - intron_start + 1),
            "M",
            str(right_aa),
            str(exon_intervals[1][1] - exon_intervals[1][0] + 1),
        ]
    lines = [
        "C4 Alignment:\n",
        f"vulgar: {query} {qstart} {qend} . {target} {start - 1} {end} + {score} "
        + " ".join(operations)
        + "\n",
        "# --- START OF GFF DUMP ---\n",
        f"{target}\ttest\tgene\t{start}\t{end}\t{score}\t+\t.\t"
        "gene_id 1 ; identity 80.0 ; similarity 90.0\n",
    ]
    for index, (left, right) in enumerate(exon_intervals):
        lines.append(f"{target}\ttest\tcds\t{left}\t{right}\t.\t+\t.\t\n")
        lines.append(
            f"{target}\ttest\texon\t{left}\t{right}\t.\t+\t.\t"
            "insertions 0 ; deletions 0\n"
        )
        if index + 1 < len(exon_intervals):
            intron_start = right + 1
            intron_end = exon_intervals[index + 1][0] - 1
            lines.append(
                f"{target}\ttest\tintron\t{intron_start}\t{intron_end}\t.\t+\t.\t"
                f"intron_id {index + 1}\n"
            )
    lines.extend(
        [
            "# --- END OF GFF DUMP ---\n",
            "#\n",
            "A" * target_bases + "\n",
            "\n",
        ]
    )
    return "".join(lines)


def _blast_row(query: str, start: int, end: int) -> str:
    length = end - start + 1
    return "\t".join(
        map(
            str,
            (
                query,
                "chr1",
                80,
                length,
                0,
                0,
                1,
                length,
                start,
                end,
                1e-20,
                100,
                0,
                1,
                400,
                5000,
            ),
        )
    )


def _write_contract_files(tmp_path: Path) -> tuple[Path, Path, Path]:
    regions = tmp_path / "regions.tsv"
    regions.write_text(
        "region_id\tsseqid\texpanded_start_1\texpanded_end_1\tslen\n"
        "r1\tchr1\t1001\t3000\t5000\n",
        encoding="utf-8",
    )
    targets = tmp_path / "targets.fa"
    targets.write_text(">r1\n" + "A" * 2000 + "\n", encoding="utf-8")
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "region_id\tselected_query_ids\n"
        "r1\tq_single,q_long,q_short\n",
        encoding="utf-8",
    )
    return regions, targets, manifest


def test_nc_tiered_selection_keeps_shortest_overlapping_spliced_model(
    tmp_path: Path,
) -> None:
    regions, targets, manifest = _write_contract_files(tmp_path)
    raw = tmp_path / "raw.txt"
    raw.write_text(
        _gff_block("q_single", "r1", 0, 200, [(101, 700)], 900)
        + _gff_block("q_long", "r1", 0, 200, [(201, 500), (601, 900)], 500)
        + _gff_block("q_short", "r1", 0, 190, [(221, 500), (601, 890)], 600),
        encoding="utf-8",
    )
    hits = tmp_path / "hits.tsv"
    hits.write_text(
        _blast_row("q_long", 1201, 1500)
        + "\n"
        + _blast_row("q_long", 1601, 1900)
        + "\n",
        encoding="utf-8",
    )
    prefix = tmp_path / "test"
    summary = select_exonerate_candidates(
        raw,
        regions,
        targets,
        manifest,
        hits,
        prefix,
    )
    with Path(f"{prefix}.multiexon_refinement_candidates.tsv").open(
        encoding="utf-8", newline=""
    ) as handle:
        candidates = list(csv.DictReader(handle, delimiter="\t"))
    assert len(candidates) == 1
    assert candidates[0]["query_id"] == "q_short"
    assert candidates[0]["genomic_start_1"] == "1221"
    assert candidates[0]["genomic_end_1"] == "1890"
    assert candidates[0]["selection_tier_aa"] == "180"
    assert summary["counts"]["input_alignment_blocks"] == 3
    assert summary["counts"]["models_with_introns"] == 2
    assert summary["counts"]["selected_refinement_candidates"] == 1
    assert summary["counts"]["rejected_single_exon"] == 1


def test_chimera_rule_rejects_more_tblastn_segments_than_exons(tmp_path: Path) -> None:
    regions, targets, manifest = _write_contract_files(tmp_path)
    raw = tmp_path / "raw.txt"
    raw.write_text(
        _gff_block("q_long", "r1", 0, 200, [(201, 500), (801, 1100)], 500),
        encoding="utf-8",
    )
    hits = tmp_path / "hits.tsv"
    hits.write_text(
        "\n".join(
            (
                _blast_row("q_long", 1201, 1300),
                _blast_row("q_long", 1401, 1500),
                _blast_row("q_long", 1801, 1900),
            )
        )
        + "\n",
        encoding="utf-8",
    )
    try:
        select_exonerate_candidates(
            raw,
            regions,
            targets,
            manifest,
            hits,
            tmp_path / "test",
        )
    except ValueError as exc:
        assert "no multi-exon candidates" in str(exc)
    else:
        raise AssertionError("probable chimera was not rejected")
