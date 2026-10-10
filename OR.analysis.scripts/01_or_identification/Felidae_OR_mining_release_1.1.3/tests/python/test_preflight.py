from pathlib import Path

import yaml

from or_mining.preflight import load_params, run_preflight


def _write_fasta(path: Path, records: dict[str, str]) -> None:
    path.write_text(
        "".join(f">{identifier}\n{sequence}\n" for identifier, sequence in records.items()),
        encoding="utf-8",
    )


def test_preflight_accepts_cross_file_contract(tmp_path: Path) -> None:
    genome = tmp_path / "genome.fa"
    _write_fasta(genome, {"chrTest": "A" * 1200})
    Path(str(genome) + ".fai").write_text("chrTest\t1200\t9\t1200\t1201\n", encoding="utf-8")

    or_db = tmp_path / "or.prot.fa"
    gpcr_db = tmp_path / "gpcr.prot.fa"
    template = tmp_path / "template.fa"
    clade_db = tmp_path / "clade.prot.fa"
    _write_fasta(or_db, {"query1": "M" * 300})
    _write_fasta(gpcr_db, {"OR-Receptor-test": "M" * 300})
    _write_fasta(template, {"known1": "M" * 300, "known2": "M" * 300})
    _write_fasta(
        clade_db,
        {
                f"{clade.upper()}---test": "M" * 300
            for clade in (
                "alpha",
                "beta",
                "gamma",
                "delta",
                "epsilon",
                "eta",
                "zeta",
                "kappa",
                "lambda",
                "theta",
            )
        },
    )
    type1 = tmp_path / "type1.ids"
    type2 = tmp_path / "type2.ids"
    type1.write_text("known1\n", encoding="utf-8")
    type2.write_text("known2\n", encoding="utf-8")
    ipr = tmp_path / "interproscan"
    ipr.mkdir()
    (ipr / "interproscan.sh").write_text("#!/usr/bin/env bash\n", encoding="utf-8")

    merged = tmp_path / "merged.tsv"
    merged.write_text(
        "query1\tchrTest\t90\t100\t10\t0\t1\t100\t100\t399\t1e-20\t200\t0\t1\t300\t1200\n",
        encoding="utf-8",
    )
    params = {
        "taxon_genus": "Canis",
        "taxon_species": "lupus",
        "taxon_subspecies": "familiaris",
        "species": "Clu_fa",
        "assembly_id": "test",
        "genome": str(genome),
        "input_merged": str(merged),
        "or_db": str(or_db),
        "gpcr_db": str(gpcr_db),
        "template_alignment": str(template),
        "known_type1_ids": str(type1),
        "known_type2_ids": str(type2),
        "or_clade_db": str(clade_db),
        "interproscan_dir": str(ipr),
        "outdir": str(tmp_path / "Clu_fa"),
        "slurm_partition": "debug",
        "slurm_account": "test-account",
        "slurm_qos": "test-qos",
        "max_cpus": 16,
        "max_memory_gb": 128,
        "max_time_h": 24,
        "candidate_min_width": 100,
        "candidate_flank": 1000,
        "min_orf_nt": 750,
        "maximum_intron_length": 10000,
        "exonerate_top_queries_per_region": 3,
        "exonerate_pilot_pairs": 6,
        "exonerate_soft_timeout_minutes": 5,
        "strict_refinement_timeout_minutes": 5,
        "or_operational_cluster_gap_bp": 1000000,
    }
    params_path = tmp_path / "params.yml"
    params_path.write_text(yaml.safe_dump(params), encoding="utf-8")

    report = run_preflight(params_path)
    assert report["status"] == "pass"
    assert report["counts"].get("fail", 0) == 0


def test_load_params_rejects_duplicate_yaml_keys(tmp_path: Path) -> None:
    params_path = tmp_path / "duplicate.yml"
    params_path.write_text(
        "species: Fca\ninterproscan_version: 5.78-109.0\n"
        "interproscan_version: 5.79-110.0\n",
        encoding="utf-8",
    )

    try:
        load_params(params_path)
    except ValueError as exc:
        assert "duplicate YAML key 'interproscan_version'" in str(exc)
    else:
        raise AssertionError("duplicate YAML key was silently accepted")


def test_preflight_allows_generated_merged_input_to_be_pending(tmp_path: Path) -> None:
    # Reuse the fully valid fixture, then switch only the initial input contract.
    test_preflight_accepts_cross_file_contract(tmp_path)
    params_path = tmp_path / "params.yml"
    params = yaml.safe_load(params_path.read_text(encoding="utf-8"))
    params["initial_search_mode"] = "generate"
    params["input_merged"] = None
    params_path.write_text(yaml.safe_dump(params), encoding="utf-8")

    report = run_preflight(params_path)
    assert report["status"] == "pass_with_warnings"
    input_check = next(item for item in report["checks"] if item["name"] == "path:input_merged")
    assert input_check["status"] == "warning"
    assert "will be generated" in input_check["detail"]
