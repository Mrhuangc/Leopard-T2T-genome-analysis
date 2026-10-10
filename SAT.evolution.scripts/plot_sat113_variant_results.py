#!/usr/bin/env python3
"""Create publication-oriented diagnostic plots for SAT113 variant analysis."""

from __future__ import annotations

import argparse
import math
import os
from pathlib import Path
from typing import Dict, List, Sequence

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-sat113")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D


def save_figure(fig: plt.Figure, prefix: Path, dpi: int) -> None:
    fig.savefig(prefix.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(prefix.with_suffix(".png"), dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def species_palette(species: Sequence[str]) -> Dict[str, str]:
    preferred = ["#D55E00", "#0072B2", "#009E73", "#CC79A7", "#E69F00"]
    return {name: preferred[i % len(preferred)] for i, name in enumerate(species)}


def finite_series(values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    return numeric[np.isfinite(numeric)]


def panel_label(axis: plt.Axes, label: str) -> None:
    axis.text(-0.12, 1.06, label, transform=axis.transAxes, fontsize=14, fontweight="bold", va="top")


def diagnostics_figure(result_dir: Path, outdir: Path, main_threshold: float, dpi: int) -> None:
    threshold = pd.read_csv(result_dir / "threshold_sensitivity.tsv", sep="\t")
    arrays = pd.read_csv(result_dir / "array_variant_summary.tsv", sep="\t")
    pairs = pd.read_csv(result_dir / "array_pair_similarity.tsv", sep="\t")
    ld = pd.read_csv(result_dir / "site_pair_linkage.tsv", sep="\t")
    metrics = pd.read_csv(result_dir / "hypothesis_diagnostic_metrics.tsv", sep="\t")
    metric = dict(zip(metrics["metric"], metrics["value"]))

    sns.set_theme(style="whitegrid", context="paper", font_scale=1.05)
    fig, axes = plt.subplots(2, 3, figsize=(15.5, 9.2), constrained_layout=True)

    ax = axes[0, 0]
    lines = [
        ("fraction_copies_cross_array_cluster", "Multi-array cluster", "#0072B2"),
        ("fraction_copies_cross_species_cluster", "Cross-species cluster", "#D55E00"),
        ("fraction_copies_species_specific_cluster", "Species-specific cluster", "#009E73"),
    ]
    for column, label, color in lines:
        ax.plot(threshold["identity_threshold"] * 100, threshold[column] * 100, marker="o", lw=2, label=label, color=color)
    ax.axvline(main_threshold * 100, color="0.35", ls="--", lw=1)
    ax.set(xlabel="Variant clustering identity (%)", ylabel="Copies represented (%)", title="Threshold sensitivity")
    ax.legend(frameon=False, fontsize=8)
    panel_label(ax, "A")

    ax = axes[0, 1]
    exact_array = 100 * float(metric.get("fraction_copies_in_cross_array_exact_haplotypes", np.nan))
    exact_species = 100 * float(metric.get("fraction_copies_in_cross_species_exact_haplotypes", np.nan))
    cluster_array = 100 * float(metric.get("fraction_copies_in_cross_array_variant_clusters", np.nan))
    cluster_species = 100 * float(metric.get("fraction_copies_in_cross_species_variant_clusters", np.nan))
    x = np.arange(2)
    width = 0.36
    ax.bar(x - width / 2, [exact_array, exact_species], width, label="Exact haplotype", color="#999999")
    ax.bar(x + width / 2, [cluster_array, cluster_species], width, label=f"Variant cluster ({main_threshold:.0%})", color="#56B4E9")
    ax.set_xticks(x, ["Across arrays", "Across species"])
    ax.set_ylabel("Copies represented (%)")
    ax.set_title("Exact versus similarity-based sharing")
    ax.legend(frameon=False, fontsize=8)
    panel_label(ax, "B")

    ax = axes[0, 2]
    same = pairs.loc[pairs["same_species"].astype(str).str.lower().isin({"true", "1"})].copy()
    if same.empty:
        same = pairs.loc[pairs["same_species"] == True].copy()  # noqa: E712
    sns.scatterplot(
        data=same,
        x="consensus_divergence_pct",
        y="allele_frequency_mean_jsd",
        hue="species_1",
        s=26,
        alpha=0.60,
        linewidth=0,
        ax=ax,
    )
    ax.set(
        xlabel="Array-consensus divergence (%)",
        ylabel="Minor-allele frequency JSD",
        title="Same consensus can hide different copy pools",
    )
    if ax.legend_:
        ax.legend_.set_title(None)
        ax.legend_.set_frame_on(False)
    panel_label(ax, "C")

    ax = axes[1, 0]
    if not ld.empty:
        plot_ld = ld.copy()
        plot_ld["NMI"] = pd.to_numeric(
            plot_ld["categorical_normalized_mutual_information"], errors="coerce"
        )
        plot_ld = plot_ld.loc[np.isfinite(plot_ld["NMI"])]
        sns.violinplot(data=plot_ld, x="group", y="NMI", inner="quart", cut=0, color="#CC79A7", ax=ax)
        ax.tick_params(axis="x", rotation=30)
    ax.axhline(0.2, color="0.35", ls="--", lw=1)
    ax.set(xlabel=None, ylabel="Site-pair normalized mutual information", title="Mutation co-occurrence within monomers")
    panel_label(ax, "D")

    ax = axes[1, 1]
    spatial = arrays.dropna(subset=["adjacent_same_expected_random_order", "adjacent_same_observed"]).copy()
    palette = species_palette(list(dict.fromkeys(spatial["species"].astype(str))))
    for species, subset in spatial.groupby("species"):
        ax.scatter(
            subset["adjacent_same_expected_random_order"],
            subset["adjacent_same_observed"],
            s=np.clip(np.sqrt(subset["n_copies"]) * 4, 12, 110),
            alpha=0.65,
            label=species,
            color=palette.get(str(species)),
            edgecolor="white",
            linewidth=0.3,
        )
    limit = max(
        0.05,
        float(spatial[["adjacent_same_expected_random_order", "adjacent_same_observed"]].max().max()) * 1.04,
    ) if not spatial.empty else 1.0
    ax.plot([0, limit], [0, limit], color="0.25", ls="--", lw=1)
    ax.set_xlim(0, limit)
    ax.set_ylim(0, limit)
    ax.set(
        xlabel="Expected adjacent-same fraction after shuffling",
        ylabel="Observed adjacent-same fraction",
        title="Spatial blocks versus interspersion",
    )
    ax.legend(frameon=False, fontsize=8)
    panel_label(ax, "E")

    ax = axes[1, 2]
    values = arrays[
        ["species", "within_array_pi_pct", "fraction_copies_in_cross_species_clusters"]
    ].dropna()
    sns.scatterplot(
        data=values,
        x="within_array_pi_pct",
        y="fraction_copies_in_cross_species_clusters",
        hue="species",
        size=None,
        alpha=0.7,
        s=35,
        ax=ax,
    )
    ax.set(
        xlabel="Within-array copy diversity, π (%)",
        ylabel="Copies in cross-species clusters",
        title="Array diversity and shared variant pool",
    )
    ax.yaxis.set_major_formatter(lambda x, pos: f"{100*x:.0f}%")
    if ax.legend_:
        ax.legend_.set_title(None)
        ax.legend_.set_frame_on(False)
    panel_label(ax, "F")

    save_figure(fig, outdir / "SAT113_variant_diagnostics", dpi)


def origin_figure(result_dir: Path, outdir: Path, main_threshold: float, dpi: int) -> None:
    threshold_label = f"{main_threshold:.3f}".rstrip("0").rstrip(".")
    clusters = pd.read_csv(
        result_dir / f"variant_cluster_summary.identity_{threshold_label}.tsv", sep="\t"
    )
    pca = pd.read_csv(result_dir / "haplotype_pca.tsv", sep="\t")
    species_pairs = pd.read_csv(result_dir / "species_pair_origin_metrics.tsv", sep="\t")
    rarefaction = pd.read_csv(result_dir / "rarefaction_species_sharing.tsv", sep="\t")
    composition = pd.read_csv(result_dir / "array_variant_composition.long.tsv", sep="\t")
    species = [column[2:] for column in clusters.columns if column.startswith("n_") and column not in {
        "n_copies", "n_haplotypes", "n_arrays", "n_chromosomes", "n_species_observed", "n_species_qualified"
    }]
    species = [name for name in species if name in set(pca["dominant_species"].astype(str))]
    palette = species_palette(species)

    sns.set_theme(style="whitegrid", context="paper", font_scale=1.05)
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 10.2), constrained_layout=True)

    ax = axes[0, 0]
    plot_pca = pca.loc[pca["n_copies"] >= 2].copy() if "n_copies" in pca else pca.copy()
    if len(plot_pca) > 6000:
        plot_pca = plot_pca.nlargest(6000, "n_copies")
    for name, subset in plot_pca.groupby("dominant_species"):
        shared = subset["sharing_class"].eq("shared_exact_haplotype")
        sizes = np.clip(np.log2(subset["n_copies"] + 1) * 8, 8, 90)
        ax.scatter(
            subset["PC1"],
            subset["PC2"],
            s=sizes,
            c=palette.get(str(name), "0.5"),
            alpha=np.where(shared, 0.95, 0.35),
            marker="o",
            linewidth=np.where(shared, 0.6, 0),
            edgecolor=np.where(shared, "black", "none"),
            label=name,
        )
    pc1_ratio = pca.get("PC1_explained_variance_ratio", pd.Series([np.nan])).iloc[0]
    pc2_ratio = pca.get("PC2_explained_variance_ratio", pd.Series([np.nan])).iloc[0]
    ax.set(
        xlabel=f"PC1 ({100*pc1_ratio:.1f}%)" if np.isfinite(pc1_ratio) else "PC1",
        ylabel=f"PC2 ({100*pc2_ratio:.1f}%)" if np.isfinite(pc2_ratio) else "PC2",
        title="Exact monomer haplotypes",
    )
    ax.legend(frameon=False, fontsize=8)
    panel_label(ax, "A")

    ax = axes[0, 1]
    class_counts = clusters.groupby("classification")["n_copies"].sum().sort_values()
    colors = ["#56B4E9" if str(name).startswith("shared") else "#E69F00" if str(name).endswith("specific") else "#999999" for name in class_counts.index]
    ax.barh(class_counts.index, class_counts.values, color=colors)
    ax.set(xlabel="Copies", ylabel=None, title=f"Variant-cluster origin classes ({main_threshold:.0%})")
    ax.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
    panel_label(ax, "B")

    ax = axes[1, 0]
    long_rows: List[Dict[str, object]] = []
    for row in species_pairs.itertuples(index=False):
        pair = f"{row.species_1}–{row.species_2}"
        long_rows.extend(
            [
                {"pair": pair, "metric": "Consensus divergence", "value": row.consensus_divergence_pct},
                {"pair": pair, "metric": "Within-species π", "value": row.mean_within_species_pi_pct},
                {"pair": pair, "metric": "Between-species Dxy", "value": row.copy_pool_dxy_pct},
                {"pair": pair, "metric": "Net divergence Da", "value": row.net_divergence_da_pct},
            ]
        )
    long = pd.DataFrame(long_rows)
    sns.barplot(data=long, x="pair", y="value", hue="metric", ax=ax, palette="colorblind")
    ax.set(xlabel=None, ylabel="Substitution divergence (%)", title="Consensus versus full copy pools")
    ax.tick_params(axis="x", rotation=20)
    ax.legend(frameon=False, fontsize=8)
    panel_label(ax, "C")

    ax = axes[1, 1]
    rare = rarefaction.loc[np.isclose(rarefaction["identity_threshold"], main_threshold)].copy()
    jaccard_columns = [column for column in rare.columns if column.startswith("jaccard_")]
    rare_long = rare.melt(
        id_vars=["replicate"],
        value_vars=jaccard_columns,
        var_name="species_pair",
        value_name="Jaccard",
    )
    rare_long["species_pair"] = rare_long["species_pair"].str.replace("jaccard_", "", regex=False).str.replace("_", "–", regex=False)
    sns.boxplot(data=rare_long, x="species_pair", y="Jaccard", color="#009E73", width=0.55, fliersize=0, ax=ax)
    sns.stripplot(data=rare_long, x="species_pair", y="Jaccard", color="black", alpha=0.18, size=2, ax=ax)
    ax.set(xlabel=None, ylabel="Cluster-presence Jaccard", title="Rarefied cross-species sharing")
    ax.tick_params(axis="x", rotation=20)
    panel_label(ax, "D")

    save_figure(fig, outdir / "SAT113_species_origin", dpi)

    # A separate wide heatmap keeps the main figure readable.
    top_variants = clusters.nlargest(min(40, len(clusters)), "n_copies")["variant_cluster"]
    top_arrays = (
        composition.groupby(["array_uid", "species"])["n_copies"]
        .sum()
        .reset_index()
        .sort_values(["species", "n_copies"], ascending=[True, False])
        .groupby("species", sort=False)
        .head(25)["array_uid"]
    )
    heat = composition.loc[
        composition["variant_cluster"].isin(top_variants)
        & composition["array_uid"].isin(top_arrays)
    ].pivot_table(index="variant_cluster", columns="array_uid", values="fraction", fill_value=0)
    heat = heat.reindex(index=[variant for variant in top_variants if variant in heat.index])
    heat = heat.reindex(columns=[array for array in top_arrays if array in heat.columns])
    if not heat.empty:
        fig_width = max(12, 0.18 * heat.shape[1] + 5)
        fig_height = max(8, 0.20 * heat.shape[0] + 2)
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        sns.heatmap(np.log10(heat + 1e-4), cmap="mako", ax=ax, cbar_kws={"label": "log10(array fraction + 1e-4)"})
        ax.set(xlabel="Array", ylabel=f"Top {heat.shape[0]} variant clusters", title="SAT113 variant composition across arrays")
        ax.tick_params(axis="x", labelrotation=90, labelsize=6)
        ax.tick_params(axis="y", labelsize=7)
        fig.tight_layout()
        save_figure(fig, outdir / "SAT113_top_variant_array_heatmap", dpi)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result-dir", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--main-threshold", type=float, default=0.85)
    parser.add_argument("--dpi", type=int, default=300)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    diagnostics_figure(args.result_dir, args.outdir, args.main_threshold, args.dpi)
    origin_figure(args.result_dir, args.outdir, args.main_threshold, args.dpi)


if __name__ == "__main__":
    main()
