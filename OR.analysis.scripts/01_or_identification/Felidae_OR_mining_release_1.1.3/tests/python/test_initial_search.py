import json
from pathlib import Path

import yaml

from or_mining.initial_search import (
    merge_initial_tblastn,
    resolve_merged_input,
    split_initial_queries,
    write_resolved_params,
)


def test_initial_query_split_is_complete_balanced_and_deterministic(tmp_path: Path) -> None:
    query = tmp_path / "queries.fa"
    query.write_text(
        ">q1 first\n" + "M" * 10 + "\n"
        ">q2\n" + "M" * 8 + "\n"
        ">q3\n" + "M" * 6 + "\n"
        ">q4\n" + "M" * 4 + "\n",
        encoding="utf-8",
    )
    first = tmp_path / "first"
    second = tmp_path / "second"
    payload = split_initial_queries(query, first, str(tmp_path / "Fca"), 2)
    split_initial_queries(query, second, str(tmp_path / "Fca.second"), 2)

    assert payload["status"] == "pass"
    assert payload["counts"]["query_records"] == 4
    assert payload["counts"]["query_shards"] == 2
    assert payload["counts"]["maximum_shard_residues"] == 14
    assert payload["counts"]["minimum_shard_residues"] == 14
    assert [p.read_text() for p in sorted(first.glob("*.fa"))] == [
        p.read_text() for p in sorted(second.glob("*.fa"))
    ]


def test_initial_tblastn_merge_is_ordered_and_validated(tmp_path: Path) -> None:
    second = tmp_path / "Fca.shard_002.initial.tblastn.tsv"
    first = tmp_path / "Fca.shard_001.initial.tblastn.tsv"
    first.write_text(
        "q1\tchr1\t90\t2\t0\t0\t1\t2\t1\t6\t1e-20\t50\t0\t1\t2\t100\n",
        encoding="utf-8",
    )
    second.write_text(
        "q2\tchr1\t80\t2\t0\t0\t1\t2\t20\t15\t1e-10\t40\t0\t-1\t2\t100\n",
        encoding="utf-8",
    )
    merged = tmp_path / "Fca.initial_merged.tblastn.tsv"
    validation = tmp_path / "Fca.initial_merged.validation.json"
    summary = tmp_path / "Fca.initial_search.summary.json"
    manifest = tmp_path / "Fca.initial_query_manifest.tsv"
    manifest.write_text(
        "shard_id\tfasta\trecords\tresidues\tfirst_query\tlast_query\n"
        "shard_001\tshard_001.prot.fa\t1\t2\tq1\tq1\n"
        "shard_002\tshard_002.prot.fa\t1\t2\tq2\tq2\n",
        encoding="utf-8",
    )
    payload = merge_initial_tblastn(
        [second, first],
        manifest,
        merged,
        validation,
        summary,
        species="Fca",
        evalue="1e-5",
    )

    assert merged.read_text(encoding="utf-8").splitlines()[0].startswith("q1\t")
    assert payload["search_completeness"] == "complete"
    assert payload["counts"]["completed_query_shards"] == 2
    assert payload["counts"]["merged_rows"] == 2
    assert payload["checks"]["invalid_rows"] == 0
    assert json.loads(validation.read_text())["status"] == "pass"


def test_resolve_generated_input_and_write_snapshot(tmp_path: Path) -> None:
    outdir = tmp_path / "Fca"
    merged = outdir / "00_initial_search/merged/Fca.initial_merged.tblastn.tsv"
    merged.parent.mkdir(parents=True)
    merged.write_text("placeholder\n", encoding="utf-8")
    params = {
        "species": "Fca",
        "outdir": str(outdir),
        "initial_search_mode": "generate",
        "input_merged": None,
    }
    resolved, source = resolve_merged_input(params)
    assert resolved == merged
    assert source == "generated"

    params_path = tmp_path / "params.yml"
    params_path.write_text(yaml.safe_dump(params), encoding="utf-8")
    snapshot = tmp_path / "resolved.yml"
    written = write_resolved_params(params_path, snapshot)
    assert written["input_merged"] == str(merged.resolve())
    assert written["resolved_initial_input_source"] == "generated"
