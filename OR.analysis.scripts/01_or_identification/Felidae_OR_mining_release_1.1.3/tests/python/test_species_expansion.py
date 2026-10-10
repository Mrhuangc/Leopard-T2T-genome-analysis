"""Tests for species-specific OR query preparation and TBLASTN expansion."""

from pathlib import Path

from or_mining.expansion_hits import merge_species_expansion_hits
from or_mining.species_queries import prepare_species_queries


def _hit(query: str, start: int = 100) -> str:
    return (
        f"{query}\tChr1\t90.0\t100\t10\t0\t1\t100\t{start}\t{start + 299}"
        "\t1e-30\t200\t0\t1\t100\t2000\n"
    )


def test_prepare_queries_removes_n_without_clustering(tmp_path: Path) -> None:
    proteins = tmp_path / "or.prot.fa"
    nucleotides = tmp_path / "or.nt.fa"
    prefix = tmp_path / "sp"
    proteins.write_text(">a\nMAAA\n>b\nMBBB\n>c\nMCCC\n", encoding="utf-8")
    nucleotides.write_text(
        ">a\nATGGCTGCTGCT\n>b\nATGNNNGCTGCT\n>c\nATGTGTTGTTGT\n",
        encoding="utf-8",
    )

    summary = prepare_species_queries(
        proteins,
        nucleotides,
        prefix,
        cluster_threshold=500,
    )

    assert summary["cd_hit_applied"] is False
    assert summary["counts"] == {
        "input_candidates": 3,
        "clear_candidates": 2,
        "ambiguous_candidates": 1,
        "query_representatives": 2,
    }
    assert Path(f"{prefix}.ambiguous_OR.ids.txt").read_text(encoding="utf-8") == "b\n"


def test_merge_species_hits_filters_queries_at_threshold(tmp_path: Path) -> None:
    initial = tmp_path / "initial.tsv"
    species = tmp_path / "species.tsv"
    prefix = tmp_path / "sp"
    initial.write_text(_hit("known", 10), encoding="utf-8")
    species.write_text(
        _hit("too_many", 400) + _hit("too_many", 800) + _hit("kept", 1200),
        encoding="utf-8",
    )

    summary = merge_species_expansion_hits(
        initial,
        species,
        prefix,
        max_hits_per_query=2,
    )

    assert summary["counts"]["outlier_species_queries"] == 1
    assert summary["counts"]["filtered_species_rows"] == 1
    assert summary["counts"]["merged_rows"] == 2
    assert Path(f"{prefix}.species.outlier_query_ids.txt").read_text(
        encoding="utf-8"
    ) == "too_many\n"
