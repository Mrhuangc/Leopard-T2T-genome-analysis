"""Validate the primary audit summary after each production stage."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from or_mining.preflight import load_params

SUMMARY_PATTERNS = {
    0: "00_initial_search/merged/{s}.initial_search.summary.json",
    1: "03_OR_Finder_results/qc/{s}.merged.validation.json",
    2: "03_OR_Finder_results/candidate_regions/{s}.candidate_regions.summary.json",
    3: "03_OR_Finder_results/orfs/{s}.orfs.summary.json",
    4: "03_OR_Finder_results/gpcr_classification/{s}.gpcr_classification.summary.json",
    5: "03_OR_Finder_results/phylogeny/filtered/{s}.tree_filter.summary.json",
    6: "03_OR_Finder_results/phylogeny/nc_r_parity/{s}.nc_r_parity.summary.json",
    7: "04_species_OR_expansion/merged/{s}.species_expansion.summary.json",
    8: "05_multiexon_OR/candidate_regions/{s}.multiexon_regions.summary.json",
    9: "05_multiexon_OR/candidate_sequences/{s}.multiexon_extraction.summary.json",
    10: "05_multiexon_OR/exonerate/merged_targeted/{s}.exonerate.raw.summary.json",
    11: "05_multiexon_OR/exonerate/model_selection/{s}.multiexon_model_selection.summary.json",
    12: "05_multiexon_OR/exonerate/strict_refinement/merged/{s}.strict_refinement.summary.json",
    13: "05_multiexon_OR/classification/{s}.multiexon_classification.summary.json",
    14: "06_OR_catalog/pre7tm/{s}.pre7tm_OR_catalog.summary.json",
    15: "08_final_OR_catalog/{s}.final_OR_catalog.summary.json",
    16: "09_NC_final_homology/catalog/{s}.NC_final_OR_catalog.summary.json",
    17: "11_complete_OR_catalog/{s}.complete_OR_catalog.summary.json",
    18: "12_OR_clade_classification/catalog/{s}.OR_clade_classification.summary.json",
    19: "13_single_species_OR_analysis/architecture/{s}.OR_locus_architecture.summary.json",
}


def check_stage(params_file: Path, stage: int) -> tuple[Path, dict[str, Any]]:
    if stage not in SUMMARY_PATTERNS:
        raise ValueError(f"stage must be one of 0..19, got {stage}")
    params = load_params(params_file)
    species = str(params.get("species") or "")
    outdir = Path(str(params.get("outdir") or "results"))
    summary = outdir / SUMMARY_PATTERNS[stage].format(s=species)
    if not summary.is_file():
        raise FileNotFoundError(f"stage {stage} summary is missing: {summary}")
    payload = json.loads(summary.read_text(encoding="utf-8"))
    allowed = {"pass", "pass_with_warnings"}
    if payload.get("status") not in allowed:
        raise ValueError(
            f"stage {stage} summary status is {payload.get('status')!r}, expected {sorted(allowed)}"
        )
    if stage == 0:
        if payload.get("search_completeness") != "complete":
            raise ValueError("stage 0 initial TBLASTN search is not complete")
        counts = payload.get("counts", {})
        if int(counts.get("completed_query_shards", 0)) != int(
            counts.get("input_query_shards", -1)
        ):
            raise ValueError("stage 0 did not complete every initial TBLASTN query shard")
        checks = payload.get("checks", {})
        if int(checks.get("invalid_rows", -1)) != 0:
            raise ValueError("stage 0 merged TBLASTN contains invalid rows")
        required_checks = (
            "all_query_shards_completed",
            "query_manifest_matches_parts",
            "merged_output_nonempty",
            "canonical_16_column_contract",
        )
        failed = [key for key in required_checks if checks.get(key) is not True]
        if failed:
            raise ValueError(f"stage 0 checks failed: {failed}")
    if stage == 10:
        if payload.get("search_completeness") != "complete":
            raise ValueError("stage 10 Exonerate search is not complete")
        counts = payload.get("counts", {})
        if int(counts.get("soft_timed_out_chunks", 0)) != 0:
            raise ValueError("stage 10 contains soft-timed-out Exonerate cells")
    if stage in {15, 16, 17} and payload.get("reportable_final_catalog") is False:
        raise ValueError(f"stage {stage} catalog is explicitly marked non-reportable")
    return summary, payload
