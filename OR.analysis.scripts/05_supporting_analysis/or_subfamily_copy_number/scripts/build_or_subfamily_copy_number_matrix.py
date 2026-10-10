#!/usr/bin/env python3
"""Expand representative OR subfamily labels to loci and build copy-number matrices."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path


UNIT_ORDER = ["Fca", "Lge", "Pbe", "Ple", "Ppa", "Ppa_or", "Pti"]
NAMED_CONFIDENCE = {"high", "medium", "low"}
PRIMARY_CONFIDENCE = {"high", "medium"}
UNNAMED = {"", "na", "unknown", "unresolved", "none"}


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Missing header: {path}")
        return list(reader.fieldnames), list(reader)


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require_columns(label: str, fields: list[str], required: set[str]) -> None:
    missing = sorted(required - set(fields))
    if missing:
        raise ValueError(f"{label} missing required columns: {', '.join(missing)}")


def is_named(row: dict[str, str]) -> bool:
    subfamily = row["final_subfamily"].strip()
    return subfamily.lower() not in UNNAMED and row["assignment_confidence"].strip().lower() in NAMED_CONFIDENCE


def matrix_rows(
    loci: list[dict[str, str]],
    subfamily_meta: dict[str, tuple[str, str]],
    confidence: set[str],
) -> list[dict[str, object]]:
    counts: Counter[tuple[str, str]] = Counter()
    for row in loci:
        if is_named(row) and row["assignment_confidence"].lower() in confidence:
            counts[(row["final_subfamily"], row["analysis_unit"])] += 1

    output: list[dict[str, object]] = []
    for subfamily in sorted(subfamily_meta):
        family, broad_clade = subfamily_meta[subfamily]
        values = [counts[(subfamily, unit)] for unit in UNIT_ORDER]
        record: dict[str, object] = {
            "or_family": family,
            "or_subfamily": subfamily,
            "broad_clade": broad_clade,
        }
        record.update(dict(zip(UNIT_ORDER, values)))
        record.update(
            {
                "total": sum(values),
                "minimum": min(values),
                "maximum": max(values),
                "count_range": max(values) - min(values),
            }
        )
        output.append(record)
    return output


def main() -> int:
    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--representatives",
        type=Path,
        default=project_root / "analysis/OR_subfamily_comparison_v1_20260903/input_snapshot/representative_OR_subfamily_assignments.tsv",
    )
    parser.add_argument(
        "--locus-map",
        type=Path,
        default=project_root / "analysis/OR_multispecies_20260825/joint_phylogeny_input/functional_OR_identical_sequence_map.tsv",
    )
    parser.add_argument(
        "--species-metadata",
        type=Path,
        default=project_root / "analysis/OR_subfamily_comparison_v1_20260903/input_snapshot/species_metadata.tsv",
    )
    parser.add_argument(
        "--previous-locus-table",
        type=Path,
        default=project_root / "analysis/OR_subfamily_comparison_v1_20260903/input_snapshot/locus_OR_subfamily_assignments.tsv",
    )
    parser.add_argument("--output-dir", type=Path, default=script_dir / "results")
    args = parser.parse_args()

    inputs = [args.representatives, args.locus_map, args.species_metadata]
    for path in inputs:
        if not path.is_file():
            raise FileNotFoundError(path)

    rep_fields, representatives = read_tsv(args.representatives)
    map_fields, locus_map = read_tsv(args.locus_map)
    metadata_fields, metadata = read_tsv(args.species_metadata)
    require_columns(
        "representative table",
        rep_fields,
        {
            "representative_id", "joint_broad_clade", "final_family", "final_subfamily",
            "assignment_confidence", "assignment_decision", "tree_nearest_reference_id",
        },
    )
    require_columns(
        "locus map",
        map_fields,
        {
            "representative_id", "global_id", "analysis_unit", "catalog_id", "protein_sha256",
            "protein_length_aa", "identical_group_size", "analysis_units_in_group",
        },
    )
    require_columns("species metadata", metadata_fields, {"analysis_unit", "scientific_name"})

    rep_ids = [row["representative_id"].strip() for row in representatives]
    map_rep_ids = [row["representative_id"].strip() for row in locus_map]
    global_ids = [row["global_id"].strip() for row in locus_map]
    catalog_ids = [row["catalog_id"].strip() for row in locus_map]
    rep_duplicates = sorted(k for k, v in Counter(rep_ids).items() if v != 1)
    global_duplicates = sorted(k for k, v in Counter(global_ids).items() if v != 1)
    catalog_duplicates = sorted(k for k, v in Counter(catalog_ids).items() if v != 1)
    missing_assignments = sorted(set(map_rep_ids) - set(rep_ids))
    unused_assignments = sorted(set(rep_ids) - set(map_rep_ids))
    units = sorted(set(row["analysis_unit"] for row in locus_map))

    hard_checks = {
        "representative_rows_are_5338": len(representatives) == 5338,
        "representative_ids_are_unique": not rep_duplicates and len(set(rep_ids)) == 5338,
        "locus_map_rows_are_6419": len(locus_map) == 6419,
        "locus_global_ids_are_unique": not global_duplicates and len(set(global_ids)) == 6419,
        "locus_catalog_ids_are_unique": not catalog_duplicates and len(set(catalog_ids)) == 6419,
        "representative_sets_match_exactly": not missing_assignments and not unused_assignments,
        "analysis_units_match_expected_seven": units == sorted(UNIT_ORDER),
        "metadata_matches_expected_seven": sorted(row["analysis_unit"] for row in metadata) == sorted(UNIT_ORDER),
    }
    if not all(hard_checks.values()):
        raise ValueError("Input integrity checks failed: " + ", ".join(k for k, v in hard_checks.items() if not v))

    assignment_by_rep = {row["representative_id"]: row for row in representatives}
    assignment_fields = [
        "joint_broad_clade", "final_family", "final_subfamily", "assignment_confidence",
        "assignment_decision", "tree_nearest_reference_id",
    ]
    expanded: list[dict[str, str]] = []
    for locus in locus_map:
        assignment = assignment_by_rep[locus["representative_id"]]
        expanded.append({**locus, **{field: assignment[field] for field in assignment_fields}})

    expanded_duplicates = sorted(k for k, v in Counter(row["global_id"] for row in expanded).items() if v != 1)
    hard_checks.update(
        {
            "expanded_rows_are_6419": len(expanded) == 6419,
            "each_locus_occurs_exactly_once": not expanded_duplicates and len({r["global_id"] for r in expanded}) == 6419,
            "all_expanded_rows_have_assignment_record": all(r["assignment_confidence"].strip() for r in expanded),
        }
    )

    subfamily_families: defaultdict[str, set[str]] = defaultdict(set)
    subfamily_clades: defaultdict[str, set[str]] = defaultdict(set)
    for row in representatives:
        if is_named(row):
            subfamily_families[row["final_subfamily"]].add(row["final_family"])
            subfamily_clades[row["final_subfamily"]].add(row["joint_broad_clade"])
    family_conflicts = {k: sorted(v) for k, v in subfamily_families.items() if len(v) != 1}
    clade_conflicts = {k: sorted(v) for k, v in subfamily_clades.items() if len(v) != 1}
    hard_checks["one_family_per_subfamily"] = not family_conflicts
    hard_checks["one_broad_clade_per_subfamily"] = not clade_conflicts
    if not all(hard_checks.values()):
        raise ValueError("Expansion integrity checks failed: " + ", ".join(k for k, v in hard_checks.items() if not v))

    subfamily_meta = {
        subfamily: (next(iter(subfamily_families[subfamily])), next(iter(subfamily_clades[subfamily])))
        for subfamily in subfamily_families
    }
    primary_matrix = matrix_rows(expanded, subfamily_meta, PRIMARY_CONFIDENCE)
    all_named_matrix = matrix_rows(expanded, subfamily_meta, NAMED_CONFIDENCE)
    matrix_fields = [
        "or_family", "or_subfamily", "broad_clade", *UNIT_ORDER,
        "total", "minimum", "maximum", "count_range",
    ]

    confidence_counts: Counter[tuple[str, str, str, str, str]] = Counter()
    for row in expanded:
        if is_named(row):
            confidence_counts[(
                row["analysis_unit"], row["final_family"], row["final_subfamily"],
                row["joint_broad_clade"], row["assignment_confidence"].lower(),
            )] += 1
    counts_long = [
        {
            "analysis_unit": key[0], "or_family": key[1], "or_subfamily": key[2],
            "broad_clade": key[3], "assignment_confidence": key[4], "functional_loci": value,
        }
        for key, value in sorted(confidence_counts.items())
    ]

    completeness_rows: list[dict[str, object]] = []
    for unit in UNIT_ORDER:
        unit_rows = [row for row in expanded if row["analysis_unit"] == unit]
        by_confidence = Counter(row["assignment_confidence"].lower() for row in unit_rows)
        primary = sum(1 for row in unit_rows if is_named(row) and row["assignment_confidence"].lower() in PRIMARY_CONFIDENCE)
        named = sum(1 for row in unit_rows if is_named(row))
        unresolved = len(unit_rows) - named
        completeness_rows.append(
            {
                "analysis_unit": unit,
                "functional_loci": len(unit_rows),
                "high_confidence_named": by_confidence["high"],
                "medium_confidence_named": by_confidence["medium"],
                "low_confidence_named": by_confidence["low"],
                "primary_high_medium_named": primary,
                "all_named": named,
                "unresolved_or_unnamed": unresolved,
                "primary_completeness_fraction": f"{primary / len(unit_rows):.6f}",
                "all_named_completeness_fraction": f"{named / len(unit_rows):.6f}",
            }
        )

    primary_sums = {unit: sum(int(row[unit]) for row in primary_matrix) for unit in UNIT_ORDER}
    all_named_sums = {unit: sum(int(row[unit]) for row in all_named_matrix) for unit in UNIT_ORDER}
    completeness_by_unit = {row["analysis_unit"]: row for row in completeness_rows}
    hard_checks.update(
        {
            "primary_matrix_column_sums_match_loci": all(
                primary_sums[u] == completeness_by_unit[u]["primary_high_medium_named"] for u in UNIT_ORDER
            ),
            "all_named_matrix_column_sums_match_loci": all(
                all_named_sums[u] == completeness_by_unit[u]["all_named"] for u in UNIT_ORDER
            ),
            "completeness_rows_sum_to_6419": sum(int(r["functional_loci"]) for r in completeness_rows) == 6419,
        }
    )

    previous_crosscheck = {"performed": False, "status": "not_available", "mismatched_cells": None}
    if args.previous_locus_table.is_file():
        old_fields, old_rows = read_tsv(args.previous_locus_table)
        common = [
            "representative_id", "global_id", "analysis_unit", "catalog_id", "protein_sha256",
            "protein_length_aa", "identical_group_size", "analysis_units_in_group", "final_family",
            "final_subfamily", "assignment_confidence", "assignment_decision", "tree_nearest_reference_id",
        ]
        require_columns("previous locus table", old_fields, set(common))
        old_by_global = {row["global_id"]: row for row in old_rows}
        mismatch_count = sum(
            1
            for row in expanded
            for field in common
            if old_by_global.get(row["global_id"], {}).get(field) != row.get(field)
        )
        previous_crosscheck = {
            "performed": True,
            "status": "pass" if len(old_rows) == len(expanded) and mismatch_count == 0 else "fail",
            "previous_rows": len(old_rows),
            "mismatched_cells": mismatch_count,
            "previous_sha256": sha256(args.previous_locus_table),
        }
        hard_checks["previous_locus_table_crosscheck_passes"] = previous_crosscheck["status"] == "pass"

    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    expanded_fields = map_fields + assignment_fields
    write_tsv(output_dir / "locus_OR_subfamily_assignments.expanded.tsv", expanded_fields, expanded)
    write_tsv(output_dir / "OR_subfamily_copy_number_matrix.high_medium.tsv", matrix_fields, primary_matrix)
    write_tsv(output_dir / "OR_subfamily_copy_number_matrix.all_named.tsv", matrix_fields, all_named_matrix)
    write_tsv(
        output_dir / "OR_subfamily_counts_by_confidence.tsv",
        ["analysis_unit", "or_family", "or_subfamily", "broad_clade", "assignment_confidence", "functional_loci"],
        counts_long,
    )
    completeness_fields = list(completeness_rows[0])
    write_tsv(output_dir / "OR_subfamily_assignment_completeness.tsv", completeness_fields, completeness_rows)
    write_tsv(output_dir / "duplicate_locus_ids.tsv", ["global_id"], [{"global_id": x} for x in expanded_duplicates])
    write_tsv(output_dir / "unmatched_representative_ids.tsv", ["representative_id", "issue"], [
        *({"representative_id": x, "issue": "map_without_assignment"} for x in missing_assignments),
        *({"representative_id": x, "issue": "assignment_without_map"} for x in unused_assignments),
    ])

    ranked = sorted(
        (row for row in primary_matrix if int(row["total"]) > 0),
        key=lambda row: (-int(row["count_range"]), -int(row["maximum"]), -int(row["total"]), str(row["or_subfamily"])),
    )[:30]
    heatmap_rows: list[dict[str, object]] = []
    selection_rows: list[dict[str, object]] = []
    for rank, row in enumerate(ranked, start=1):
        selection_rows.append({"display_rank": rank, **row})
        for unit in UNIT_ORDER:
            heatmap_rows.append(
                {
                    "display_rank": rank,
                    "or_family": row["or_family"],
                    "or_subfamily": row["or_subfamily"],
                    "broad_clade": row["broad_clade"],
                    "analysis_unit": unit,
                    "absolute_copy_number": row[unit],
                }
            )
    write_tsv(output_dir / "absolute_heatmap_top30.source_data.tsv", list(heatmap_rows[0]), heatmap_rows)
    write_tsv(output_dir / "absolute_heatmap_top30.selection.tsv", ["display_rank", *matrix_fields], selection_rows)
    write_tsv(output_dir / "species_metadata.tsv", metadata_fields, metadata)

    total_by_confidence = Counter(row["assignment_confidence"].lower() for row in expanded)
    qc = {
        "schema_version": "1.0",
        "status": "pass" if all(hard_checks.values()) else "fail",
        "observational_unit": "functional OR locus in a genome-assembly analysis unit",
        "interpretation_boundary": (
            "Descriptive copy-number differences only; expansion, contraction, gain, and loss require gene-tree and synteny evidence."
        ),
        "input_files": {
            "representative_assignments": {"path": str(args.representatives.resolve()), "sha256": sha256(args.representatives)},
            "locus_map": {"path": str(args.locus_map.resolve()), "sha256": sha256(args.locus_map)},
            "species_metadata": {"path": str(args.species_metadata.resolve()), "sha256": sha256(args.species_metadata)},
        },
        "counts": {
            "representative_rows": len(representatives),
            "functional_loci": len(expanded),
            "analysis_units": len(UNIT_ORDER),
            "standard_named_subfamilies": len(subfamily_meta),
            "high_confidence_loci": total_by_confidence["high"],
            "medium_confidence_loci": total_by_confidence["medium"],
            "low_confidence_loci": total_by_confidence["low"],
            "unresolved_loci": len(expanded) - sum(total_by_confidence[c] for c in NAMED_CONFIDENCE),
            "primary_high_medium_named_loci": sum(primary_sums.values()),
            "all_named_loci": sum(all_named_sums.values()),
            "primary_displayed_subfamilies": len(ranked),
        },
        "checks": hard_checks,
        "previous_locus_table_crosscheck": previous_crosscheck,
        "family_conflicts": family_conflicts,
        "broad_clade_conflicts": clade_conflicts,
    }
    with (output_dir / "OR_subfamily_copy_number.QC.json").open("w", encoding="utf-8") as handle:
        json.dump(qc, handle, ensure_ascii=False, indent=2)
        handle.write("\n")

    summary_lines = [
        f"status={qc['status']}",
        f"representative_rows={len(representatives)}",
        f"functional_loci={len(expanded)}",
        f"unique_functional_loci={len(set(global_ids))}",
        f"analysis_units={len(UNIT_ORDER)}",
        f"standard_named_subfamilies={len(subfamily_meta)}",
        f"high_confidence_loci={total_by_confidence['high']}",
        f"medium_confidence_loci={total_by_confidence['medium']}",
        f"low_confidence_loci={total_by_confidence['low']}",
        f"unresolved_loci={qc['counts']['unresolved_loci']}",
        f"primary_high_medium_named_loci={sum(primary_sums.values())}",
        f"all_named_loci={sum(all_named_sums.values())}",
        f"duplicate_locus_ids={len(expanded_duplicates)}",
        f"unmatched_representative_ids={len(missing_assignments) + len(unused_assignments)}",
        f"previous_locus_table_crosscheck={previous_crosscheck['status']}",
        "interpretation=descriptive_copy_number_differences_not_expansion_or_contraction",
    ]
    (output_dir / "OR_subfamily_copy_number.QC.txt").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")

    manifest_rows = []
    for path in sorted(output_dir.iterdir()):
        if path.is_file() and path.name != "output_manifest.tsv":
            manifest_rows.append({"file": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)})
    write_tsv(output_dir / "output_manifest.tsv", ["file", "bytes", "sha256"], manifest_rows)

    print("OR_SUBFAMILY_COPY_NUMBER_BUILD_COMPLETE")
    print("\n".join(summary_lines))
    print(f"output_dir={output_dir}")
    return 0 if qc["status"] == "pass" else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
