#!/usr/bin/env python3
"""Audit whether high pairwise SAT113 divergence follows from copy-to-consensus divergence.

This quick check uses only existing v2.1 tables; raw monomer sequences are not
required.  It is a diagnostic relationship, not a test of haplotype structure.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Dict, List

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-sat113-audit")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        required=True,
        type=Path,
        help="TSV with species,copy_table,array_summary columns.",
    )
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(args.manifest, sep="\t", dtype=str)
    required = {"species", "copy_table", "array_summary"}
    missing = sorted(required - set(manifest.columns))
    if missing:
        raise ValueError(f"Manifest is missing columns: {', '.join(missing)}")
    base = args.manifest.resolve().parent

    rows: List[pd.DataFrame] = []
    for dataset in manifest.itertuples(index=False):
        copy_path = Path(dataset.copy_table).expanduser()
        array_path = Path(dataset.array_summary).expanduser()
        if not copy_path.is_absolute():
            copy_path = base / copy_path
        if not array_path.is_absolute():
            array_path = base / array_path
        copies = pd.read_csv(copy_path, sep="\t", dtype={"copy_id": str})
        arrays = pd.read_csv(array_path, sep="\t")
        copies = copies.loc[
            copies["phase_status"].eq("PASS") & copies["quality_tier"].eq("full_length")
        ].copy()
        copies = copies.drop_duplicates(["chrom", "start", "end", "strand"], keep="first")
        p = (
            copies.groupby("array_key")["substitution_divergence_pct"].mean().div(100.0).rename("p")
        )
        frame = arrays.set_index("array_key").join(p, how="inner").reset_index()
        frame["species"] = dataset.species
        frame["observed_pairwise_divergence_pct"] = frame[
            "intra_array_pairwise_sub_div_mean_full_length_pct"
        ]
        # If the non-consensus mass p is equally split across three alternatives:
        frame["predicted_three_alternative_pct"] = 100.0 * (
            2.0 * frame["p"] - (4.0 / 3.0) * frame["p"] ** 2
        )
        # If every variable site has only one alternative allele:
        frame["predicted_one_alternative_pct"] = 100.0 * (
            2.0 * frame["p"] * (1.0 - frame["p"])
        )
        frame["three_alternative_residual_pct"] = (
            frame["observed_pairwise_divergence_pct"]
            - frame["predicted_three_alternative_pct"]
        )
        rows.append(frame)

    combined = pd.concat(rows, ignore_index=True)
    combined.to_csv(args.outdir / "divergence_relationship.by_array.tsv", sep="\t", index=False)
    summary_rows: List[Dict[str, object]] = []
    for species, subset in combined.groupby("species", sort=False):
        valid = subset[["observed_pairwise_divergence_pct", "predicted_three_alternative_pct"]].dropna()
        correlation = np.corrcoef(
            valid["observed_pairwise_divergence_pct"],
            valid["predicted_three_alternative_pct"],
        )[0, 1]
        summary_rows.append(
            {
                "species": species,
                "n_arrays": len(valid),
                "mean_copy_to_consensus_divergence_pct": 100.0 * subset["p"].mean(),
                "mean_observed_pairwise_divergence_pct": valid[
                    "observed_pairwise_divergence_pct"
                ].mean(),
                "mean_predicted_three_alternative_pct": valid[
                    "predicted_three_alternative_pct"
                ].mean(),
                "median_observed_minus_predicted_pct": subset[
                    "three_alternative_residual_pct"
                ].median(),
                "pearson_r_observed_vs_predicted": correlation,
            }
        )
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(args.outdir / "divergence_relationship.summary.tsv", sep="\t", index=False)

    sns.set_theme(style="whitegrid", context="paper", font_scale=1.1)
    species_order = list(summary["species"])
    fig, axes = plt.subplots(1, len(species_order), figsize=(5.1 * len(species_order), 4.6), squeeze=False)
    all_values = combined[
        ["observed_pairwise_divergence_pct", "predicted_three_alternative_pct"]
    ].to_numpy()
    lower = max(0.0, float(np.nanmin(all_values)) - 2.0)
    upper = float(np.nanmax(all_values)) + 2.0
    for ax, species in zip(axes[0], species_order):
        subset = combined.loc[combined["species"] == species]
        sns.scatterplot(
            data=subset,
            x="predicted_three_alternative_pct",
            y="observed_pairwise_divergence_pct",
            s=30,
            alpha=0.65,
            linewidth=0,
            ax=ax,
        )
        ax.plot([lower, upper], [lower, upper], color="0.25", ls="--", lw=1)
        r = summary.loc[summary["species"] == species, "pearson_r_observed_vs_predicted"].iloc[0]
        ax.set(
            xlim=(lower, upper),
            ylim=(lower, upper),
            xlabel="Predicted pairwise divergence (%)",
            ylabel="Observed pairwise divergence (%)",
            title=f"{species}  r={r:.3f}",
        )
    fig.suptitle("SAT113: pairwise divergence expected from copy-to-consensus divergence", y=1.02)
    fig.tight_layout()
    fig.savefig(args.outdir / "divergence_relationship.pdf", bbox_inches="tight")
    fig.savefig(args.outdir / "divergence_relationship.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
