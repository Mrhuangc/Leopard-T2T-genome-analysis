#!/usr/bin/env python3
"""Regenerate the frozen four-track PAR/Xdeg input from Y orders and X loci.

The bundled ChrX catalog is an extracted, curated subset of historical GFF3s;
provide both --fca-x-gff and --ppa-x-gff to reselect loci from full GFF3s.
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path


TRACK_ORDER = ("FCA_ChrX", "FCA_ChrY", "PPA_ChrY", "PPA_ChrX")
REGION_ORDER = ("PAR", "X-degenerate")

# The ChrX GFF3 files contain additional paralogous PLCXD1/GTPBP6/STAG2
# models outside the contiguous PAR gene block. Restricting the PAR facet to
# these curated blocks prevents unrelated same-name X loci from being linked.
PAR_X_INTERVALS = {
    "FCA": (3_000_000, 11_000_000),
    "PPA": (120_000_000, 128_000_000),
}

# PPA ChrX stores its contiguous PAR block in the opposite coordinate
# orientation to both ChrY assemblies. Reverse only this schematic facet so
# conserved PAR order is readable; source coordinates and source strands are
# retained in separate columns.
DISPLAY_ORIENTATION = {
    ("PPA_ChrX", "PAR"): "reverse",
}

# PAR names are already shared X/Y symbols except this curated Y-specific label.
PAR_Y_TO_X = {
    "SHROOM2Y": "SHROOM2",
}

# Explicit ancestral X/Y homology aliases. These are the only aliases used.
XDG_Y_TO_X = {
    "CLDN34Y": "CLDN34",
    "FELIDY_ORF1": "FELIDY_ORF1",
    "UTY": "KDM6A",
    "DDX3Y": "DDX3X",
    "USP9Y": "USP9X",
    "AP1S2Y": "AP1S2",
    "EEF1A1Y": "EEF1A1",
    "UBE1Y": "UBA1",
    "KDM5D": "KDM5C",
    "EIF1AY": "EIF1AX",
    "MBTPS2Y": "MBTPS2",
    "ZFY": "ZFX",
    "EIF2S3Y": "EIF2S3",
    "AMELY": "AMELX",
    "TXLNGY": "TXLNG",
    "OFD1Y": "OFD1",
    "TMSB4Y": "TMSB4X",
    "RAB9AY": "RAB9A",
    "RBMY": "RBMX",
}

OUTPUT_COLUMNS = (
    "track",
    "region",
    "order",
    "slot",
    "gene",
    "link_id",
    "homolog_group",
    "strand",
    "source_strand",
    "display_orientation",
    "status",
    "label",
    "label_show",
    "chrom",
    "start0",
    "end0",
    "source_id",
    "annotation_class",
    "evidence",
)


def parse_attributes(text: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for field in text.strip().strip(";").split(";"):
        if not field:
            continue
        if "=" in field:
            key, value = field.split("=", 1)
            attrs[key] = value
    return attrs


def base_gene_name(name: str) -> str:
    """Collapse annotation copy suffixes such as CSF2RA.c1 for family matching."""
    return re.sub(r"\.c\d+$", "", name)


def read_x_genes(path: Path) -> list[dict[str, object]]:
    genes: list[dict[str, object]] = []
    with path.open() as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            attrs = parse_attributes(fields[8])
            source_id = attrs.get("ID") or attrs.get("Name")
            name = attrs.get("Name") or attrs.get("ID")
            if not source_id or not name:
                continue
            genes.append(
                {
                    "chrom": fields[0],
                    "start0": int(fields[3]) - 1,
                    "end0": int(fields[4]),
                    "strand": fields[6] if fields[6] in {"+", "-"} else "+",
                    "source_id": source_id,
                    "gene": base_gene_name(name),
                }
            )
    return genes


def read_y_order(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def group_for_y(region: str, family: str) -> str:
    if region == "PAR":
        return PAR_Y_TO_X.get(family, family)
    try:
        return XDG_Y_TO_X[family]
    except KeyError as exc:
        raise ValueError(f"No explicit X-degenerate X/Y alias for {family}") from exc


def make_link_id(region: str, group: str, copy_ordinal: int) -> str:
    prefix = "PAR" if region == "PAR" else "XDG"
    return f"{prefix}_{group}_{copy_ordinal}"


def build_y_rows(
    order_rows: list[dict[str, str]],
    rab9_review: Path,
) -> list[dict[str, object]]:
    selected: list[dict[str, object]] = []
    for source in order_rows:
        region = source["region_class"]
        if region not in REGION_ORDER:
            continue
        # The two unnamed provisional PAR ORFs are not comparable genes.
        if region == "PAR" and source["annotation_class"] == "provisional":
            continue
        species = source["species"]
        family = source["family"]
        selected.append(
            {
                "track": f"{species}_ChrY",
                "region": region,
                "gene": family,
                "homolog_group": group_for_y(region, family),
                "source_strand": source["y_strand"],
                "status": (
                    "candidate"
                    if source["annotation_class"] == "provisional"
                    else "intact"
                ),
                "chrom": "ChrY",
                "start0": int(source["y_start0"]),
                "end0": int(source["y_end0"]),
                "source_id": source["gene_id"],
                "annotation_class": source["annotation_class"],
                "evidence": "curated_ChrY_annotation",
            }
        )

    with rab9_review.open(newline="") as handle:
        review_rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(review_rows) != 1:
        raise ValueError("RAB9AY review must contain exactly one reviewed FCA locus")
    review = review_rows[0]
    if review["classification"] != "subthreshold_RAB9A_like_fragment":
        raise ValueError("Unexpected FCA RAB9A/Y review classification")
    selected.append(
        {
            "track": "FCA_ChrY",
            "region": "X-degenerate",
            "gene": "RAB9A/Y-remnant",
            "homolog_group": "RAB9A",
            "source_strand": review["fca_strand"],
            "status": "remnant",
            "chrom": review["fca_chrom"],
            "start0": int(review["fca_start0"]),
            "end0": int(review["fca_end0"]),
            "source_id": "FCA_RAB9A_Y_like_fragment",
            "annotation_class": review["classification"],
            "evidence": "manual_RAB9AY_synteny_review",
        }
    )
    return selected


def build_x_rows(
    y_rows: list[dict[str, object]],
    x_gffs: dict[str, Path] | None,
    x_catalog: Path,
) -> list[dict[str, object]]:
    target_groups = {
        region: {
            str(row["homolog_group"])
            for row in y_rows
            if row["region"] == region
        }
        for region in REGION_ORDER
    }
    rows: list[dict[str, object]] = []
    if x_gffs is None:
        with x_catalog.open(newline="") as handle:
            catalog = list(csv.DictReader(handle, delimiter="\t"))
        for entry in catalog:
            species = entry["track"].split("_", 1)[0]
            region = entry["region"]
            group = entry["gene"]
            if species not in {"FCA", "PPA"} or region not in REGION_ORDER or group not in target_groups[region]:
                raise ValueError(f"Unexpected ChrX catalog entry: {entry}")
            rows.append({
                "track": entry["track"], "region": region, "gene": group,
                "homolog_group": group, "source_strand": entry["source_strand"],
                "status": "intact", "chrom": entry["chrom"],
                "start0": int(entry["start0"]), "end0": int(entry["end0"]),
                "source_id": entry["source_id"], "annotation_class": "ChrX_gene",
                "evidence": "ChrX_GFF3_exact_alias",
            })
        return rows
    for species, path in x_gffs.items():
        genes = read_x_genes(path)
        for region in REGION_ORDER:
            for gene in genes:
                group = str(gene["gene"])
                if group not in target_groups[region]:
                    continue
                if region == "PAR":
                    block_start, block_end = PAR_X_INTERVALS[species]
                    if not (
                        int(gene["start0"]) < block_end
                        and int(gene["end0"]) > block_start
                    ):
                        continue
                rows.append(
                    {
                        "track": f"{species}_ChrX",
                        "region": region,
                        "gene": group,
                        "homolog_group": group,
                        "source_strand": gene["strand"],
                        "status": "intact",
                        "chrom": gene["chrom"],
                        "start0": gene["start0"],
                        "end0": gene["end0"],
                        "source_id": gene["source_id"],
                        "annotation_class": "ChrX_gene",
                        "evidence": "ChrX_GFF3_exact_alias",
                    }
                )
    return rows


def assign_order_links_and_labels(
    rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[(str(row["track"]), str(row["region"]))].append(row)

    finished: list[dict[str, object]] = []
    for track in TRACK_ORDER:
        for region in REGION_ORDER:
            orientation = DISPLAY_ORIENTATION.get((track, region), "forward")
            group_rows = sorted(
                grouped[(track, region)],
                key=lambda r: (int(r["start0"]), int(r["end0"]), str(r["source_id"])),
                reverse=orientation == "reverse",
            )
            if track == "PPA_ChrX" and region == "PAR":
                # Later diagram review ordered the two OBP copies by source
                # coordinate inside the otherwise reversed ChrX display.
                positions = [i for i, r in enumerate(group_rows) if r["gene"] == "OBP"]
                copies = sorted((group_rows[i] for i in positions), key=lambda r: int(r["start0"]))
                for i, copy in zip(positions, copies):
                    group_rows[i] = copy
            copy_count = Counter(str(row["homolog_group"]) for row in group_rows)
            copy_seen: Counter[str] = Counter()
            for order, row in enumerate(group_rows, start=1):
                group = str(row["homolog_group"])
                copy_seen[group] += 1
                copy_ordinal = copy_seen[group]
                row["order"] = order
                row["slot"] = order
                row["link_id"] = make_link_id(region, group, copy_ordinal)
                row["display_orientation"] = orientation
                source_strand = str(row["source_strand"])
                row["strand"] = (
                    source_strand
                    if orientation == "forward"
                    else ("-" if source_strand == "+" else "+")
                )
                label = str(row["gene"])
                if copy_count[group] > 1:
                    label = f"{label}-{copy_ordinal}"
                row["label"] = label
                # Frozen plot-table decisions from the later 208-row review:
                # X OBP labels are left unsuffixed.
                if track.endswith("_ChrX") and region == "PAR" and group == "OBP":
                    row["label"] = "OBP"
                # Kept verbatim for the frozen base input. The subsequent
                # reclassification table normalizes this link to CLDN34.
                if track == "PPA_ChrY" and region == "X-degenerate" and row["source_id"] == "CLDN34.c2":
                    row["link_id"] = "XDG_CLDN34Y_1"
                row["label_show"] = int(
                    region == "X-degenerate"
                    or track in {"FCA_ChrX", "PPA_ChrX"}
                )
                finished.append(row)
    return finished


def validate_rows(rows: list[dict[str, object]]) -> None:
    duplicates = Counter(
        (str(r["track"]), str(r["region"]), str(r["link_id"])) for r in rows
    )
    repeated = [key for key, count in duplicates.items() if count > 1]
    if repeated:
        raise ValueError(f"Duplicate copy-specific link_id within track/region: {repeated}")

    link_groups: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        link_groups[str(row["link_id"])].add(str(row["homolog_group"]))
    inconsistent = {key: values for key, values in link_groups.items() if len(values) != 1}
    if inconsistent:
        raise ValueError(f"Link IDs map to multiple homolog groups: {inconsistent}")

    expected = {
        ("FCA_ChrX", "PAR"): 34,
        ("FCA_ChrY", "PAR"): 36,
        ("PPA_ChrY", "PAR"): 36,
        ("PPA_ChrX", "PAR"): 36,
        ("FCA_ChrX", "X-degenerate"): 15,
        ("FCA_ChrY", "X-degenerate"): 16,
        ("PPA_ChrY", "X-degenerate"): 18,
        ("PPA_ChrX", "X-degenerate"): 17,
    }
    observed = Counter((str(r["track"]), str(r["region"])) for r in rows)
    if observed != expected:
        raise ValueError(f"Unexpected track/region counts: {observed}; expected {expected}")


def write_rows(rows: list[dict[str, object]], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=OUTPUT_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row[column] for column in OUTPUT_COLUMNS})


def build(
    y_order: Path,
    fca_x_gff: Path | None,
    ppa_x_gff: Path | None,
    x_catalog: Path,
    rab9_review: Path,
    output: Path,
) -> list[dict[str, object]]:
    y_rows = build_y_rows(read_y_order(y_order), rab9_review)
    if (fca_x_gff is None) != (ppa_x_gff is None):
        raise ValueError("Pass both --fca-x-gff and --ppa-x-gff, or neither")
    x_gffs = {"FCA": fca_x_gff, "PPA": ppa_x_gff} if fca_x_gff else None
    x_rows = build_x_rows(y_rows, x_gffs, x_catalog)
    rows = assign_order_links_and_labels(y_rows + x_rows)
    validate_rows(rows)
    write_rows(rows, output)
    return rows


def main() -> None:
    project = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--y-order",
        type=Path,
        default=project / "par_xdeg_4tracks/PAR_XDG_gene_order.tsv",
    )
    parser.add_argument(
        "--fca-x-gff",
        type=Path,
        default=None,
    )
    parser.add_argument(
        "--ppa-x-gff",
        type=Path,
        default=None,
    )
    parser.add_argument("--x-catalog", type=Path,
                        default=project / "par_xdeg_4tracks/ChrX_gene_catalog.tsv")
    parser.add_argument(
        "--rab9-review",
        type=Path,
        default=project / "results/02_fca_pseudogenes/RAB9AY.review.tsv",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=project / "par_xdeg_4tracks/base_input.tsv",
    )
    args = parser.parse_args()

    rows = build(
        args.y_order,
        args.fca_x_gff,
        args.ppa_x_gff,
        args.x_catalog,
        args.rab9_review,
        args.output,
    )
    print(f"Wrote {len(rows)} rows to {args.output}")


if __name__ == "__main__":
    main()
