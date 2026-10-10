"""Tests for candidate FASTA renaming and getorf coordinate normalization."""

from pathlib import Path

from or_mining.extract_candidates import rename_extracted_regions
from or_mining.fasta import read_fasta
from or_mining.orfs import normalize_getorf, parse_getorf_header


REGION_HEADER = (
    "region_id\tsseqid\tstrand\thit_start_1\thit_end_1\t"
    "expanded_start_1\texpanded_end_1\thit_width_bp\texpanded_width_bp\t"
    "hit_count\tquery_count\tbest_evalue\tmax_bitscore\tslen\tfaidx_region\n"
)


def _region(
    region_id: str, expanded_start: int, expanded_end: int, hit_strand: str = "+"
) -> str:
    width = expanded_end - expanded_start + 1
    return (
        f"{region_id}\tChr1\t{hit_strand}\t200\t300\t{expanded_start}\t"
        f"{expanded_end}\t101\t{width}\t2\t2\t1e-20\t100\t2000\t"
        f"Chr1:{expanded_start}-{expanded_end}\n"
    )


def test_parse_getorf_forward_and_reverse_headers() -> None:
    forward = parse_getorf_header("Ppa_ORC000001_1 [10 - 309]")
    reverse = parse_getorf_header("Ppa_ORC000001_2 [900 - 601] (REVERSE SENSE)")

    assert (forward.region_id, forward.strand) == ("Ppa_ORC000001", "+")
    assert (reverse.local_start_1, reverse.local_end_1, reverse.strand) == (900, 601, "-")


def test_rename_extracted_regions_checks_order_and_length(tmp_path: Path) -> None:
    regions = tmp_path / "regions.tsv"
    raw = tmp_path / "raw.fa"
    output = tmp_path / "candidate.fa"
    summary = tmp_path / "summary.json"
    regions.write_text(REGION_HEADER + _region("Ppa_ORC000001", 101, 1100), encoding="utf-8")
    raw.write_text(">Chr1:101-1100\n" + "A" * 1000 + "\n", encoding="utf-8")

    result = rename_extracted_regions(raw, regions, output, summary)

    assert result["extracted_sequences"] == 1
    assert list(read_fasta(output)) == [("Ppa_ORC000001", "A" * 1000)]


def test_normalize_getorf_maps_strands_and_removes_duplicate_loci(tmp_path: Path) -> None:
    regions = tmp_path / "regions.tsv"
    raw = tmp_path / "getorf.fa"
    output_fasta = tmp_path / "orfs.fa"
    output_tsv = tmp_path / "orfs.tsv"
    summary = tmp_path / "summary.json"
    regions.write_text(
        REGION_HEADER
        + _region("Ppa_ORC000001", 101, 1100)
        + _region("Ppa_ORC000002", 1, 1000),
        encoding="utf-8",
    )
    sequence_a = "ATG" * 100
    sequence_b = "CAT" * 100
    raw.write_text(
        ">Ppa_ORC000001_1 [10 - 309]\n"
        + sequence_a
        + "\n>Ppa_ORC000001_2 [900 - 601] (REVERSE SENSE)\n"
        + sequence_b
        + "\n>Ppa_ORC000002_1 [110 - 409]\n"
        + sequence_a
        + "\n",
        encoding="utf-8",
    )

    result = normalize_getorf(
        raw, regions, output_fasta, output_tsv, summary, species="Ppa"
    )

    assert result["counts"] == {
        "raw_getorf_records": 3,
        "unique_orfs": 2,
        "duplicate_records_removed": 1,
        "plus_strand_orfs": 1,
        "minus_strand_orfs": 1,
        "orfs_with_ambiguous_bases": 0,
    }
    rows = output_tsv.read_text(encoding="utf-8").splitlines()
    assert "\tChr1\t110\t409\t+\t300\t" in rows[1]
    assert "\tChr1\t701\t1000\t-\t300\t" in rows[2]

