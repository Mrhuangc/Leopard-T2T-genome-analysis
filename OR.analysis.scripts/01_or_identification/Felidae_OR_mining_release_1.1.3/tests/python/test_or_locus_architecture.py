from pathlib import Path

from or_mining.or_locus_architecture import analyze_or_locus_architecture


def test_locus_architecture_builds_coordinate_clusters(tmp_path: Path) -> None:
    fai = tmp_path / "genome.fa.fai"
    fai.write_text("chr1\t5000000\t0\t0\t0\nchr2\t1000000\t0\t0\t0\n", encoding="utf-8")
    catalog = tmp_path / "catalog.tsv"
    catalog.write_text(
        "catalog_id\tclassification\tsseqid\tstart_1\tend_1\tstrand\texon_count\tor_clade\tor_class\n"
        "OR1\tfunctional\tchr1\t100\t999\t+\t1\talpha\tclass_I\n"
        "OR2\tpseudogene\tchr1\t50000\t50799\t-\t1\talpha\tclass_I\n"
        "OR3\tfunctional\tchr1\t2000000\t2000899\t+\t1\tgamma\tclass_II\n"
        "OR4\ttruncated\tchr2\t10\t400\t+\t1\teta\tnot_applicable\n",
        encoding="utf-8",
    )

    summary = analyze_or_locus_architecture(
        "test", catalog, fai, tmp_path / "test", cluster_gap_bp=100000
    )

    assert summary["counts"]["catalog_loci"] == 4
    assert summary["counts"]["operational_clusters_including_singletons"] == 3
    assert summary["counts"]["multi_locus_operational_clusters"] == 1
    assert summary["counts"]["largest_operational_cluster"] == 2
    assert summary["fractions"]["pseudogene_of_catalog"] == 0.25
