from pathlib import Path

from or_mining.or_clade_classification import classify_or_clades, infer_reference_clade


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def _blast(query: str, subject: str, pident: float = 70.0, qlen: int = 300) -> str:
    return (
        f"{query}\t{subject}\t{pident}\t250\t0\t0\t1\t250\t1\t250\t1e-50\t500\t"
        f"{qlen}\t300\n"
    )


def test_infer_reference_clade_handles_nc_identifiers() -> None:
    assert infer_reference_clade("GAMMA---XetrORs657.36") == "gamma"
    assert infer_reference_clade("Danio_rerio---ETA-1-CM002890") == "eta"
    assert infer_reference_clade("Cami-kappa.1") == "kappa"
    assert infer_reference_clade("NP_000670.1") is None


def test_classification_uses_tree_for_complete_and_blastx_for_incomplete(tmp_path: Path) -> None:
    template = _write(
        tmp_path / "template.fa",
        ">NP_OUT\nAAAA\n>ALPHA---r1\nAAAA\n>ALPHA---r2\nAAAA\n"
        ">ETA---r1\nAAAA\n>ETA---r2\nAAAA\n",
    )
    type1 = _write(tmp_path / "type1.txt", "ETA---r1\nETA---r2\n")
    type2 = _write(tmp_path / "type2.txt", "ALPHA---r1\nALPHA---r2\n")
    tree = _write(
        tmp_path / "tree.nwk",
        "(NP_OUT:1,(((ALPHA---r1:0.1,F1:0.1)0.95:0.1,ALPHA---r2:0.1)0.90:0.2,"
        "((ETA---r1:0.1,F2:0.1)0.96:0.1,ETA---r2:0.1)0.91:0.2):0.5);\n",
    )
    catalog = _write(
        tmp_path / "catalog.tsv",
        "catalog_id\tclassification\tsseqid\tstart_1\tend_1\tstrand\texon_count\n"
        "F1\tfunctional\tchr1\t1\t900\t+\t1\n"
        "F2\tfunctional\tchr1\t2000\t2900\t-\t1\n"
        "P1\tpseudogene\tchr2\t10\t809\t+\t1\n",
    )
    functional = _write(tmp_path / "functional.fa", ">F1\nMKK\n>F2\nMKK\n")
    functional_blast = _write(
        tmp_path / "functional.blastp.tsv",
        _blast("F1", "ALPHA---r1") + _blast("F2", "ETA---r1"),
    )
    pseudo_blast = _write(
        tmp_path / "pseudo.blastx.tsv", _blast("P1", "ALPHA---r2", qlen=800)
    )
    empty = _write(tmp_path / "empty.tsv", "")
    prefix = tmp_path / "test"

    summary = classify_or_clades(
        "test", tree, template, type1, type2, catalog, functional,
        functional_blast, pseudo_blast, empty, empty, prefix,
    )

    assert summary["status"] == "pass"
    assert summary["counts"]["by_clade"] == {"alpha": 2, "eta": 1}
    lines = (tmp_path / "test.OR_clade_assignments.tsv").read_text().splitlines()
    assert any(line.startswith("F1\tfunctional\talpha\tclass_I") for line in lines)
    assert any(line.startswith("F2\tfunctional\teta\tnot_applicable") for line in lines)
    assert any("P1\tpseudogene\talpha\tclass_I\tbest_blastx_propagation" in line for line in lines)
