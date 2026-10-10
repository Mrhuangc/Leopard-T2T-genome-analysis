#!/usr/bin/env python3
"""Model-based validation of discrete SAT113 monomer subtypes.

A one-component categorical model already allows every monomer position to
have its own A/C/G/T frequencies.  Therefore K>1 is useful only when mutations
co-occur in reproducible multi-site combinations.  Models are evaluated by
holding out complete array/region units, so an apparent subtype must generalize
to other arrays rather than merely describe one large local expansion.
"""

from __future__ import annotations

import argparse
import math
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-sat113-subtypes")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.special import logsumexp
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score


BASES = np.array(list("ACGT"))
BASE_TO_CODE = {base: i for i, base in enumerate(BASES)}


@dataclass
class MixtureResult:
    k: int
    log_likelihood: float
    pi: np.ndarray
    theta: np.ndarray
    responsibilities: np.ndarray
    iterations: int
    converged: bool
    stability_mean_ari: float = float("nan")
    stability_min_ari: float = float("nan")


def eprint(*args: object) -> None:
    print(*args, file=sys.stderr)


def normalize_fasta_id(name: str) -> str:
    return name.split()[0].split(";")[0]


def read_fasta(path: Path) -> Dict[str, str]:
    records: Dict[str, str] = {}
    name: Optional[str] = None
    chunks: List[str] = []
    with open(path, encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records[normalize_fasta_id(name)] = "".join(chunks).upper()
                name = line[1:].strip()
                chunks = []
            else:
                if name is None:
                    raise ValueError(f"Sequence before FASTA header in {path}")
                chunks.append(line)
    if name is not None:
        records[normalize_fasta_id(name)] = "".join(chunks).upper()
    return records


def write_fasta(records: Iterable[Tuple[str, str]], path: Path, width: int = 80) -> None:
    with open(path, "w", encoding="utf-8") as out:
        for name, sequence in records:
            out.write(f">{name}\n")
            for i in range(0, len(sequence), width):
                out.write(sequence[i : i + width] + "\n")


def encode_alignment(alignment: Mapping[str, str], hap_order: Sequence[str]) -> np.ndarray:
    lengths = {len(alignment[hap]) for hap in hap_order}
    if len(lengths) != 1:
        raise ValueError("Alignment sequences do not all have the same length")
    length = lengths.pop()
    encoded = np.full((len(hap_order), length), -1, dtype=np.int8)
    for i, hap in enumerate(hap_order):
        chars = np.frombuffer(alignment[hap].encode("ascii"), dtype="S1")
        for base, code in BASE_TO_CODE.items():
            encoded[i, chars == base.encode("ascii")] = code
    return encoded


def build_one_hot(encoded: np.ndarray, sites: np.ndarray) -> np.ndarray:
    selected = encoded[:, sites]
    one_hot = np.zeros((encoded.shape[0], len(sites) * 4), dtype=np.float32)
    for base in range(4):
        one_hot[:, base::4] = selected == base
    return one_hot


def unit_capped_haplotype_weights(
    copies: pd.DataFrame,
    hap_to_index: Mapping[str, int],
    units: Optional[set[str]],
    cap_per_unit: int,
) -> np.ndarray:
    weights = np.zeros(len(hap_to_index), dtype=float)
    frame = copies if units is None else copies.loc[copies["unit_uid"].isin(units)]
    for _unit, subset in frame.groupby("unit_uid", sort=False):
        scale = 1.0
        if cap_per_unit > 0 and len(subset) > cap_per_unit:
            scale = cap_per_unit / len(subset)
        counts = subset["haplotype_id"].value_counts()
        for hap, count in counts.items():
            weights[hap_to_index[hap]] += count * scale
    return weights


def parameters_from_responsibilities(
    one_hot: np.ndarray,
    weights: np.ndarray,
    responsibilities: np.ndarray,
    alpha: float,
) -> Tuple[np.ndarray, np.ndarray]:
    weighted_resp = responsibilities * weights[:, None]
    component_weight = weighted_resp.sum(axis=0)
    pi = (component_weight + alpha) / (weights.sum() + alpha * responsibilities.shape[1])
    flat_counts = weighted_resp.T @ one_hot
    theta = flat_counts.reshape(responsibilities.shape[1], -1, 4) + alpha
    theta /= theta.sum(axis=2, keepdims=True)
    return pi, theta


def sequence_log_probabilities(one_hot: np.ndarray, pi: np.ndarray, theta: np.ndarray) -> np.ndarray:
    log_theta = np.log(np.clip(theta, 1e-12, 1.0)).reshape(theta.shape[0], -1)
    return one_hot @ log_theta.T + np.log(np.clip(pi, 1e-12, 1.0))[None, :]


def fit_em_once(
    one_hot: np.ndarray,
    weights: np.ndarray,
    k: int,
    seed: int,
    alpha: float,
    max_iter: int,
    tolerance: float,
) -> MixtureResult:
    positive = np.where(weights > 0)[0]
    if len(positive) < k:
        raise ValueError(f"Only {len(positive)} positive-weight haplotypes for K={k}")
    if k == 1:
        labels = np.zeros(len(weights), dtype=int)
    else:
        kmeans = KMeans(n_clusters=k, n_init=1, max_iter=100, random_state=seed)
        kmeans.fit(one_hot[positive], sample_weight=weights[positive])
        labels = kmeans.predict(one_hot)
    responsibilities = np.full((len(weights), k), 1e-3 / max(k - 1, 1), dtype=float)
    responsibilities[np.arange(len(weights)), labels] = 1.0 - (1e-3 if k > 1 else 0.0)
    if k == 1:
        responsibilities[:, 0] = 1.0

    previous = -np.inf
    converged = False
    log_likelihood = -np.inf
    for iteration in range(1, max_iter + 1):
        pi, theta = parameters_from_responsibilities(one_hot, weights, responsibilities, alpha)
        log_probability = sequence_log_probabilities(one_hot, pi, theta)
        normalization = logsumexp(log_probability, axis=1)
        log_likelihood = float(np.sum(weights * normalization))
        responsibilities = np.exp(log_probability - normalization[:, None])
        if np.isfinite(previous):
            improvement = log_likelihood - previous
            if improvement >= -1e-7 and improvement <= tolerance * (1.0 + abs(previous)):
                converged = True
                break
        previous = log_likelihood
    pi, theta = parameters_from_responsibilities(one_hot, weights, responsibilities, alpha)
    log_probability = sequence_log_probabilities(one_hot, pi, theta)
    normalization = logsumexp(log_probability, axis=1)
    log_likelihood = float(np.sum(weights * normalization))
    responsibilities = np.exp(log_probability - normalization[:, None])

    # Canonical labels: most abundant component first.
    expected = (responsibilities * weights[:, None]).sum(axis=0)
    order = np.argsort(-expected)
    pi = pi[order]
    theta = theta[order]
    responsibilities = responsibilities[:, order]
    return MixtureResult(k, log_likelihood, pi, theta, responsibilities, iteration, converged)


def fit_multistart(
    one_hot: np.ndarray,
    weights: np.ndarray,
    k: int,
    starts: int,
    seed: int,
    alpha: float,
    max_iter: int,
    tolerance: float,
    stability_sample: np.ndarray,
) -> MixtureResult:
    fits: List[MixtureResult] = []
    for start in range(starts):
        fits.append(
            fit_em_once(
                one_hot,
                weights,
                k,
                seed + 1009 * start + 7919 * k,
                alpha,
                max_iter,
                tolerance,
            )
        )
    fits.sort(key=lambda result: result.log_likelihood, reverse=True)
    best = fits[0]
    if k == 1 or len(fits) == 1:
        best.stability_mean_ari = 1.0
        best.stability_min_ari = 1.0
    else:
        best_labels = np.argmax(best.responsibilities, axis=1)[stability_sample]
        ari = [
            adjusted_rand_score(
                best_labels,
                np.argmax(other.responsibilities, axis=1)[stability_sample],
            )
            for other in fits[1:]
        ]
        best.stability_mean_ari = float(np.mean(ari))
        best.stability_min_ari = float(np.min(ari))
    return best


def evaluate_log_likelihood(
    one_hot: np.ndarray,
    weights: np.ndarray,
    model: MixtureResult,
) -> float:
    log_probability = sequence_log_probabilities(one_hot, model.pi, model.theta)
    return float(np.sum(weights * logsumexp(log_probability, axis=1)))


def stratified_unit_folds(copies: pd.DataFrame, folds: int, rng: np.random.Generator) -> Dict[str, int]:
    assignment: Dict[str, int] = {}
    units = copies[["unit_uid", "species"]].drop_duplicates()
    for _species, subset in units.groupby("species", sort=False):
        values = subset["unit_uid"].to_numpy().copy()
        rng.shuffle(values)
        for index, unit in enumerate(values):
            assignment[str(unit)] = index % folds
    return assignment


def component_diagnostics(
    model: MixtureResult,
    selected_sites: np.ndarray,
    reference_sequence: Optional[str] = None,
) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    for local_site, alignment_site in enumerate(selected_sites):
        profile = model.theta[:, local_site, :]
        max_difference = max(
            float(profile[:, base].max() - profile[:, base].min()) for base in range(4)
        )
        consensus_codes = np.argmax(profile, axis=1)
        row: Dict[str, object] = {
            "alignment_position_1based": int(alignment_site + 1),
            "reference_position_1based": int(alignment_site + 1),
            "reference_base": (
                reference_sequence[alignment_site]
                if reference_sequence is not None
                else ""
            ),
            "max_component_base_frequency_difference": max_difference,
            "n_distinct_component_consensus_bases": len(set(consensus_codes.tolist())),
        }
        for component in range(model.k):
            row[f"subtype_{component + 1}_consensus_base"] = BASES[consensus_codes[component]]
            row[f"subtype_{component + 1}_consensus_frequency"] = profile[
                component, consensus_codes[component]
            ]
            for base, code in BASE_TO_CODE.items():
                row[f"subtype_{component + 1}_{base}_frequency"] = profile[component, code]
        rows.append(row)
    return pd.DataFrame(rows).sort_values(
        ["max_component_base_frequency_difference", "n_distinct_component_consensus_bases"],
        ascending=[False, False],
    )


def subtype_consensuses(
    encoded: np.ndarray,
    hard_labels: np.ndarray,
    raw_weights: np.ndarray,
    k: int,
) -> Tuple[List[str], List[str]]:
    aligned: List[str] = []
    ungapped: List[str] = []
    for component in range(k):
        mask = hard_labels == component
        weights = raw_weights[mask]
        matrix = encoded[mask]
        chars: List[str] = []
        for site in range(encoded.shape[1]):
            counts = np.array(
                [weights[matrix[:, site] == base].sum() for base in range(4)], dtype=float
            )
            chars.append(BASES[int(np.argmax(counts))] if counts.sum() > 0 else "-")
        sequence = "".join(chars)
        aligned.append(sequence)
        # Preserve the fixed reference coordinate in the ordinary FASTA too.
        # N marks a reference position without a canonical base in this
        # component; deleting it would shorten and shift the subtype consensus.
        ungapped.append(sequence.replace("-", "N"))
    return aligned, ungapped


def spatial_subtype_table(copies: pd.DataFrame) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    for unit, subset in copies.groupby("unit_uid", sort=False):
        subset = subset.sort_values(["start", "end", "copy_id"])
        labels = subset["subtype"].to_numpy()
        if len(labels) < 2:
            continue
        counts = np.array(list(pd.Series(labels).value_counts().values), dtype=float)
        expected = np.sum(counts * (counts - 1)) / (len(labels) * (len(labels) - 1))
        observed = np.mean(labels[:-1] == labels[1:])
        first = subset.iloc[0]
        rows.append(
            {
                "unit_uid": unit,
                "species": first["species"],
                "unit_key": first["unit_key"],
                "chrom": first["chrom"],
                "n_copies": len(labels),
                "n_subtypes": len(counts),
                "adjacent_same_observed": observed,
                "adjacent_same_expected_random": expected,
                "adjacent_same_excess": observed - expected,
                "run_count": int(1 + np.sum(labels[:-1] != labels[1:])),
            }
        )
    return pd.DataFrame(rows)


def plot_validation(
    model_selection: pd.DataFrame,
    hap_assignments: pd.DataFrame,
    subtype_summary: pd.DataFrame,
    diagnostics: pd.DataFrame,
    array_composition: pd.DataFrame,
    recommended_k: int,
    output_prefix: Path,
) -> None:
    sns.set_theme(style="whitegrid", context="paper", font_scale=1.05)
    fig, axes = plt.subplots(2, 2, figsize=(13.2, 9.4), constrained_layout=True)

    ax = axes[0, 0]
    ax.plot(model_selection["k"], model_selection["delta_bic_vs_k1"], marker="o", color="#0072B2", label="ΔBIC vs K=1")
    ax.axhline(0, color="0.35", ls="--", lw=1)
    ax.axvline(recommended_k, color="#D55E00", ls=":", lw=2, label=f"Diagnostic K={recommended_k}")
    ax.set(xlabel="Number of categorical components (K)", ylabel="ΔBIC (lower is better)", title="Model complexity")
    ax.legend(frameon=False, fontsize=8)
    ax.text(-0.12, 1.05, "A", transform=ax.transAxes, fontsize=14, fontweight="bold")

    ax = axes[0, 1]
    ax.errorbar(
        model_selection["k"],
        model_selection["cv_loglik_per_monomer_mean"],
        yerr=model_selection["cv_loglik_per_monomer_se"],
        marker="o",
        capsize=3,
        color="#009E73",
    )
    ax.axvline(recommended_k, color="#D55E00", ls=":", lw=2)
    ax.set(xlabel="K", ylabel="Held-out log likelihood per monomer", title="Generalization to held-out arrays")
    ax.text(-0.12, 1.05, "B", transform=ax.transAxes, fontsize=14, fontweight="bold")

    ax = axes[1, 0]
    species_columns = [column for column in subtype_summary.columns if column.startswith("n_") and column not in {"n_copies", "n_haplotypes", "n_units"}]
    bottom = np.zeros(len(subtype_summary))
    colors = sns.color_palette("colorblind", len(species_columns))
    for color, column in zip(colors, species_columns):
        values = subtype_summary[column].to_numpy(dtype=float)
        ax.bar(subtype_summary["subtype"], values, bottom=bottom, label=column[2:], color=color)
        bottom += values
    ax.set(xlabel=None, ylabel="Hard-assigned copies", title="Species composition of inferred subtypes")
    ax.legend(frameon=False, fontsize=8)
    ax.tick_params(axis="x", rotation=25)
    ax.text(-0.12, 1.05, "C", transform=ax.transAxes, fontsize=14, fontweight="bold")

    ax = axes[1, 1]
    top = diagnostics.head(min(25, len(diagnostics))).copy()
    if not top.empty:
        matrix = np.zeros((recommended_k, len(top)))
        annotations = np.empty((recommended_k, len(top)), dtype=object)
        for component in range(1, recommended_k + 1):
            matrix[component - 1] = top[f"subtype_{component}_consensus_frequency"]
            annotations[component - 1] = top[f"subtype_{component}_consensus_base"]
        sns.heatmap(
            matrix,
            vmin=0.25,
            vmax=1.0,
            cmap="mako",
            annot=annotations,
            fmt="",
            xticklabels=top["alignment_position_1based"],
            yticklabels=[f"Subtype {i}" for i in range(1, recommended_k + 1)],
            cbar_kws={"label": "Consensus-base frequency"},
            ax=ax,
        )
        ax.tick_params(axis="x", rotation=90, labelsize=7)
    ax.set(xlabel="Top diagnostic alignment positions", ylabel=None, title="Subtype-defining multi-site pattern")
    ax.text(-0.12, 1.05, "D", transform=ax.transAxes, fontsize=14, fontweight="bold")

    if (
        "k_search_status" in model_selection.columns
        and str(model_selection["k_search_status"].iloc[0]).startswith("UPPER_BOUNDARY")
    ):
        fig.suptitle(
            "K search reached the tested upper boundary; subtype count is unresolved",
            color="#B22222",
            fontsize=12,
            fontweight="bold",
        )

    fig.savefig(output_prefix.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(output_prefix.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)

    if not array_composition.empty and recommended_k > 1:
        heat = array_composition.pivot_table(
            index="unit_uid", columns="subtype", values="fraction", fill_value=0
        )
        order = (
            array_composition[["unit_uid", "species", "chrom"]]
            .drop_duplicates()
            .sort_values(["species", "chrom", "unit_uid"])["unit_uid"]
        )
        heat = heat.reindex(index=[unit for unit in order if unit in heat.index])
        fig, ax = plt.subplots(figsize=(7.5, max(5, 0.12 * len(heat) + 2)))
        sns.heatmap(heat, cmap="rocket_r", vmin=0, vmax=1, ax=ax, cbar_kws={"label": "Array fraction"})
        ax.set(xlabel=None, ylabel="Array/merged region", title="Model-based SAT113 subtype composition")
        ax.tick_params(axis="y", labelsize=6)
        fig.tight_layout()
        fig.savefig(output_prefix.parent / "subtype_array_composition.pdf", bbox_inches="tight")
        fig.savefig(output_prefix.parent / "subtype_array_composition.png", dpi=300, bbox_inches="tight")
        plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--haplotype-summary", required=True, type=Path)
    parser.add_argument("--alignment", required=True, type=Path)
    parser.add_argument(
        "--reference-fasta",
        type=Path,
        help="One ungapped fixed-coordinate reference used to build --alignment.",
    )
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--group-by", choices=["array_key", "region_key"], default="region_key")
    parser.add_argument("--max-k", type=int, default=12)
    parser.add_argument("--starts", type=int, default=8)
    parser.add_argument("--cv-folds", type=int, default=5)
    parser.add_argument("--cv-starts", type=int, default=2)
    parser.add_argument("--max-sites", type=int, default=80)
    parser.add_argument("--min-maf", type=float, default=0.03)
    parser.add_argument("--min-occupancy", type=float, default=0.90)
    parser.add_argument("--exclude-core-edge-sites", type=int, default=0)
    parser.add_argument("--cap-copies-per-unit", type=int, default=100)
    parser.add_argument("--min-component-fraction", type=float, default=0.01)
    parser.add_argument("--min-component-copies", type=float, default=20.0)
    parser.add_argument(
        "--boundary-buffer",
        type=int,
        default=2,
        help="Treat a CV optimum within this many K values of Kmax as unresolved.",
    )
    parser.add_argument(
        "--min-cv-increment-per-site",
        type=float,
        default=0.001,
        help="Smallest held-out K-to-K+1 improvement considered practically meaningful.",
    )
    parser.add_argument("--min-stability-ari", type=float, default=0.80)
    parser.add_argument("--dirichlet-alpha", type=float, default=0.5)
    parser.add_argument("--max-iterations", type=int, default=150)
    parser.add_argument("--tolerance", type=float, default=1e-6)
    parser.add_argument("--seed", type=int, default=20260721)
    args = parser.parse_args()
    if args.max_k < 1:
        raise ValueError("--max-k must be >= 1")
    if args.starts < 2:
        raise ValueError("--starts must be >= 2 so stability can be assessed")
    if args.cv_folds < 2:
        raise ValueError("--cv-folds must be >= 2")
    if args.cv_starts < 1:
        raise ValueError("--cv-starts must be >= 1")
    if args.boundary_buffer < 1:
        raise ValueError("--boundary-buffer must be >= 1")
    if args.min_cv_increment_per_site < 0:
        raise ValueError("--min-cv-increment-per-site must be >= 0")
    args.outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    metadata = pd.read_csv(args.metadata, sep="\t", dtype={"copy_id": str})
    hap_summary = pd.read_csv(args.haplotype_summary, sep="\t")
    if "haplotype_id" not in metadata.columns:
        raise ValueError("copy_metadata.tsv must contain haplotype_id")
    if args.group_by not in metadata.columns:
        raise ValueError(f"Metadata lacks {args.group_by}")
    metadata["unit_key"] = metadata[args.group_by].astype(str)
    metadata["unit_uid"] = metadata["species"].astype(str) + "|" + metadata["unit_key"]
    species_order = list(dict.fromkeys(metadata["species"].astype(str)))

    alignment = read_fasta(args.alignment)
    hap_order = list(hap_summary["haplotype_id"].astype(str))
    missing = sorted(set(hap_order) - set(alignment))
    if missing:
        raise ValueError(f"Alignment is missing {len(missing)} haplotypes")
    hap_to_index = {hap: i for i, hap in enumerate(hap_order)}
    encoded = encode_alignment(alignment, hap_order)
    reference_sequence: Optional[str] = None
    if args.reference_fasta is not None:
        reference_records = read_fasta(args.reference_fasta)
        if len(reference_records) != 1:
            raise ValueError("--reference-fasta must contain exactly one sequence")
        reference_sequence = next(iter(reference_records.values())).replace("-", "")
        if len(reference_sequence) != encoded.shape[1]:
            raise ValueError(
                "Fixed reference and alignment lengths differ: "
                f"reference={len(reference_sequence)}, alignment={encoded.shape[1]}"
            )
    raw_counts_series = metadata["haplotype_id"].value_counts()
    raw_weights = np.array([raw_counts_series.get(hap, 0) for hap in hap_order], dtype=float)
    full_weights = unit_capped_haplotype_weights(
        metadata, hap_to_index, None, args.cap_copies_per_unit
    )

    canonical = encoded >= 0
    occupancy = (canonical * full_weights[:, None]).sum(axis=0) / full_weights.sum()
    counts = np.zeros((encoded.shape[1], 4), dtype=float)
    for base in range(4):
        counts[:, base] = ((encoded == base) * full_weights[:, None]).sum(axis=0)
    frequencies = counts / np.maximum(counts.sum(axis=1, keepdims=True), 1e-12)
    maf = 1.0 - frequencies.max(axis=1)
    candidate_sites = np.where((occupancy >= args.min_occupancy) & (maf >= args.min_maf))[0]
    if args.exclude_core_edge_sites > 0:
        core = np.where(occupancy >= args.min_occupancy)[0]
        excluded = set(core[: args.exclude_core_edge_sites]) | set(core[-args.exclude_core_edge_sites :])
        candidate_sites = np.array([site for site in candidate_sites if site not in excluded], dtype=int)
    if len(candidate_sites) > args.max_sites:
        entropy = -np.sum(
            np.where(frequencies > 0, frequencies * np.log(np.clip(frequencies, 1e-12, 1.0)), 0),
            axis=1,
        )
        candidate_sites = candidate_sites[np.argsort(-entropy[candidate_sites])[: args.max_sites]]
        candidate_sites.sort()
    if len(candidate_sites) < 2:
        raise ValueError("Fewer than two informative sites remain; discrete subtype fitting is not meaningful")
    one_hot = build_one_hot(encoded, candidate_sites)
    eprint(
        f"Using {len(candidate_sites)} informative sites, {len(hap_order)} exact haplotypes, "
        f"effective capped weight={full_weights.sum():.1f}"
    )

    stability_size = min(50000, max(5000, int(round(full_weights.sum() * 3))))
    stability_sample = rng.choice(
        len(hap_order), size=stability_size, replace=True, p=full_weights / full_weights.sum()
    )
    max_k = min(args.max_k, int(np.sum(full_weights > 0)))
    models: Dict[int, MixtureResult] = {}
    selection_rows: List[Dict[str, object]] = []
    for k in range(1, max_k + 1):
        eprint(f"Fitting full model K={k}")
        model = fit_multistart(
            one_hot,
            full_weights,
            k,
            args.starts,
            args.seed,
            args.dirichlet_alpha,
            args.max_iterations,
            args.tolerance,
            stability_sample,
        )
        models[k] = model
        hard = np.argmax(model.responsibilities, axis=1)
        hard_weight = np.array([full_weights[hard == component].sum() for component in range(k)])
        expected_weight = (model.responsibilities * full_weights[:, None]).sum(axis=0)
        parameters = (k - 1) + k * len(candidate_sites) * 3
        bic = -2.0 * model.log_likelihood + parameters * math.log(full_weights.sum())
        aic = -2.0 * model.log_likelihood + 2.0 * parameters
        classification_entropy = -float(
            np.sum(
                full_weights[:, None]
                * model.responsibilities
                * np.log(np.clip(model.responsibilities, 1e-12, 1.0))
            )
        )
        icl = bic + 2.0 * classification_entropy
        max_posterior = np.max(model.responsibilities, axis=1)
        row = {
            "k": k,
            "log_likelihood": model.log_likelihood,
            "n_parameters": parameters,
            "bic": bic,
            "aic": aic,
            "icl": icl,
            "classification_entropy": classification_entropy,
            "mean_max_posterior": np.average(max_posterior, weights=full_weights),
            "fraction_weight_posterior_ge_0_90": np.sum(
                full_weights[max_posterior >= 0.90]
            )
            / full_weights.sum(),
            "min_hard_component_weight": hard_weight.min(),
            "min_expected_component_weight": expected_weight.min(),
            "min_component_fraction": hard_weight.min() / full_weights.sum(),
            "stability_mean_ari": model.stability_mean_ari,
            "stability_min_ari": model.stability_min_ari,
            "iterations": model.iterations,
            "converged": model.converged,
        }
        selection_rows.append(row)

    selection = pd.DataFrame(selection_rows)

    # Cross-validation holds out entire region/array units.
    fold_assignment = stratified_unit_folds(metadata, args.cv_folds, rng)
    all_units = set(fold_assignment)
    cv_rows: List[Dict[str, object]] = []
    for fold in range(args.cv_folds):
        test_units = {unit for unit, assigned in fold_assignment.items() if assigned == fold}
        train_units = all_units - test_units
        train_weights = unit_capped_haplotype_weights(
            metadata, hap_to_index, train_units, args.cap_copies_per_unit
        )
        test_weights = unit_capped_haplotype_weights(
            metadata, hap_to_index, test_units, args.cap_copies_per_unit
        )
        if train_weights.sum() <= 0 or test_weights.sum() <= 0:
            continue
        sample = rng.choice(
            len(hap_order),
            size=min(10000, max(1000, int(train_weights.sum()))),
            replace=True,
            p=train_weights / train_weights.sum(),
        )
        for k in range(1, max_k + 1):
            if np.sum(train_weights > 0) < k:
                continue
            model = fit_multistart(
                one_hot,
                train_weights,
                k,
                args.cv_starts,
                args.seed + 100003 * (fold + 1),
                args.dirichlet_alpha,
                args.max_iterations,
                args.tolerance,
                sample,
            )
            test_ll = evaluate_log_likelihood(one_hot, test_weights, model)
            cv_rows.append(
                {
                    "fold": fold + 1,
                    "k": k,
                    "n_train_units": len(train_units),
                    "n_test_units": len(test_units),
                    "test_effective_weight": test_weights.sum(),
                    "test_log_likelihood": test_ll,
                    "test_loglik_per_monomer": test_ll / test_weights.sum(),
                    "test_loglik_per_site": test_ll / (test_weights.sum() * len(candidate_sites)),
                }
            )
    cv = pd.DataFrame(cv_rows)
    cv.to_csv(args.outdir / "subtype_cross_validation.tsv", sep="\t", index=False)
    cv_summary = (
        cv.groupby("k")
        .agg(
            cv_loglik_per_monomer_mean=("test_loglik_per_monomer", "mean"),
            cv_loglik_per_monomer_sd=("test_loglik_per_monomer", "std"),
            cv_folds=("test_loglik_per_monomer", "size"),
        )
        .reset_index()
    )
    cv_summary["cv_loglik_per_monomer_se"] = cv_summary[
        "cv_loglik_per_monomer_sd"
    ] / np.sqrt(cv_summary["cv_folds"])

    # Paired fold-wise increments remove most of the irrelevant variation in
    # absolute log likelihood between easy and difficult held-out folds.  The
    # previous implementation applied a one-SE rule to the unpaired fold SD;
    # with only three folds that could select K=4 even when every fold clearly
    # improved from K=4 through the tested upper boundary K=6.
    cv_pivot = cv.pivot(index="fold", columns="k", values="test_loglik_per_monomer")
    increment_rows: List[Dict[str, object]] = []
    for k in sorted(cv_summary["k"].astype(int)):
        row: Dict[str, object] = {"k": k}
        if k > 1 and k in cv_pivot.columns and (k - 1) in cv_pivot.columns:
            differences = (cv_pivot[k] - cv_pivot[k - 1]).dropna()
            mean_increment = float(differences.mean()) if len(differences) else float("nan")
            sd_increment = float(differences.std(ddof=1)) if len(differences) > 1 else 0.0
            se_increment = sd_increment / math.sqrt(len(differences)) if len(differences) else float("nan")
            row.update(
                {
                    "cv_increment_vs_previous_per_monomer_mean": mean_increment,
                    "cv_increment_vs_previous_per_monomer_se": se_increment,
                    "cv_increment_vs_previous_per_site_mean": mean_increment / len(candidate_sites),
                    "cv_increment_vs_previous_per_site_se": se_increment / len(candidate_sites),
                }
            )
        increment_rows.append(row)
    increment_summary = pd.DataFrame(increment_rows)
    cv_summary = cv_summary.merge(increment_summary, on="k", how="left")
    selection = selection.merge(cv_summary, on="k", how="left")
    selection["delta_bic_vs_k1"] = selection["bic"] - selection.loc[
        selection["k"] == 1, "bic"
    ].iloc[0]
    baseline_cv = selection.loc[selection["k"] == 1, "cv_loglik_per_monomer_mean"].iloc[0]
    selection["cv_gain_vs_k1_per_monomer"] = (
        selection["cv_loglik_per_monomer_mean"] - baseline_cv
    )
    selection["cv_gain_vs_k1_per_site"] = selection["cv_gain_vs_k1_per_monomer"] / len(
        candidate_sites
    )
    min_weight = max(args.min_component_copies, args.min_component_fraction * full_weights.sum())
    selection["passes_min_component_size"] = selection["min_hard_component_weight"] >= min_weight
    expected_cv_folds = int(cv["fold"].nunique())
    selection["passes_complete_cv"] = selection["cv_folds"].fillna(0).astype(int) == expected_cv_folds
    selection["passes_model_validity"] = (
        selection["passes_min_component_size"] & selection["passes_complete_cv"]
    )

    valid = selection.loc[selection["passes_model_validity"]].copy()
    if valid.empty:
        raise RuntimeError(
            "No fitted K passed both the component-size and complete-CV rules"
        )
    if valid["cv_loglik_per_monomer_mean"].notna().any():
        cv_best_k = int(valid.loc[valid["cv_loglik_per_monomer_mean"].idxmax(), "k"])
        paired_rows: List[Dict[str, object]] = []
        for k in selection["k"].astype(int):
            if cv_best_k not in cv_pivot.columns or k not in cv_pivot.columns:
                continue
            losses = (cv_pivot[cv_best_k] - cv_pivot[k]).dropna()
            mean_loss = float(losses.mean()) if len(losses) else float("nan")
            sd_loss = float(losses.std(ddof=1)) if len(losses) > 1 else 0.0
            se_loss = sd_loss / math.sqrt(len(losses)) if len(losses) else float("nan")
            paired_rows.append(
                {
                    "k": k,
                    "paired_cv_loss_to_best_per_monomer_mean": mean_loss,
                    "paired_cv_loss_to_best_per_monomer_se": se_loss,
                    "paired_cv_loss_to_best_per_site_mean": mean_loss / len(candidate_sites),
                }
            )
        selection = selection.merge(pd.DataFrame(paired_rows), on="k", how="left")
        valid = selection.loc[selection["passes_model_validity"]].copy()
        practical_loss = args.min_cv_increment_per_site * len(candidate_sites)
        tolerance_to_best = np.maximum(
            valid["paired_cv_loss_to_best_per_monomer_se"].fillna(0.0),
            practical_loss,
        )
        eligible = valid.loc[
            valid["paired_cv_loss_to_best_per_monomer_mean"] <= tolerance_to_best
        ]
        if eligible.empty:
            eligible = valid.loc[valid["k"] == cv_best_k]
        recommended_k = int(eligible["k"].min())
        selection_rule = "paired array-held-out CV one-SE/practical-loss rule"
    else:
        recommended_k = int(valid.loc[valid["bic"].idxmin(), "k"])
        cv_best_k = recommended_k
        selection_rule = "minimum BIC"

    max_valid_k = int(valid["k"].max())
    bic_best_k = int(valid.loc[valid["bic"].idxmin(), "k"])
    boundary_start = max(int(valid["k"].min()), max_valid_k - max(args.boundary_buffer, 1) + 1)
    boundary_probe = selection.loc[selection["k"] == cv_best_k]
    boundary_increment = (
        float(boundary_probe["cv_increment_vs_previous_per_site_mean"].iloc[0])
        if not boundary_probe.empty
        and pd.notna(boundary_probe["cv_increment_vs_previous_per_site_mean"].iloc[0])
        else 0.0
    )
    # A small component becoming invalid is a biological/statistical stopping
    # rule, not evidence that the user-specified K ceiling was too low.
    search_boundary_reached = bool(
        max_valid_k == max_k
        and max_k == args.max_k
        and cv_best_k >= boundary_start
        and boundary_increment >= args.min_cv_increment_per_site
    )
    suggested_next_max_k = (
        max_k + max(4, int(math.ceil(max_k / 2.0)))
        if search_boundary_reached
        else max_k
    )
    k_search_status = (
        "UPPER_BOUNDARY_REACHED_EXTEND_K"
        if search_boundary_reached
        else "INTERNAL_CV_PLATEAU_OR_OPTIMUM"
    )
    selection["near_search_boundary"] = selection["k"] >= boundary_start
    selection["cv_best_k"] = cv_best_k
    selection["k_search_status"] = k_search_status
    selected_model = models[recommended_k]
    selection["selected_for_diagnostics"] = selection["k"] == recommended_k
    # Backward-compatible alias.  When k_search_status reports a boundary,
    # this K is a diagnostic simplification rather than a resolved subtype count.
    selection["recommended"] = selection["selected_for_diagnostics"]
    selection.to_csv(args.outdir / "subtype_model_selection.tsv", sep="\t", index=False)

    responsibilities = selected_model.responsibilities
    hard_labels = np.argmax(responsibilities, axis=1)
    posterior = np.max(responsibilities, axis=1)
    hap_assignments = pd.DataFrame(
        {
            "haplotype_id": hap_order,
            "subtype": [f"Subtype_{label + 1}" for label in hard_labels],
            "posterior_probability": posterior,
            "raw_copy_count": raw_weights,
            "capped_analysis_weight": full_weights,
        }
    )
    for component in range(recommended_k):
        hap_assignments[f"posterior_Subtype_{component + 1}"] = responsibilities[:, component]
    hap_assignments.to_csv(args.outdir / "subtype_haplotype_assignments.tsv", sep="\t", index=False)

    copy_assign = metadata.merge(
        hap_assignments[["haplotype_id", "subtype", "posterior_probability"]],
        on="haplotype_id",
        how="left",
        validate="many_to_one",
    )
    copy_assign.to_csv(args.outdir / "subtype_copy_assignments.tsv", sep="\t", index=False)

    subtype_rows: List[Dict[str, object]] = []
    for component in range(recommended_k):
        name = f"Subtype_{component + 1}"
        subset = copy_assign.loc[copy_assign["subtype"] == name]
        species_counts = subset.groupby("species").size().to_dict()
        row: Dict[str, object] = {
            "subtype": name,
            "n_copies": len(subset),
            "n_haplotypes": subset["haplotype_id"].nunique(),
            "n_units": subset["unit_uid"].nunique(),
            "mean_posterior_probability": subset["posterior_probability"].mean(),
            "mixture_proportion": selected_model.pi[component],
        }
        for species in species_order:
            row[f"n_{species}"] = species_counts.get(species, 0)
        subtype_rows.append(row)
    subtype_summary = pd.DataFrame(subtype_rows)
    subtype_summary.to_csv(args.outdir / "subtype_summary.tsv", sep="\t", index=False)

    composition = (
        copy_assign.groupby(["unit_uid", "species", "unit_key", "chrom", "subtype"])
        .size()
        .rename("n_copies")
        .reset_index()
    )
    composition["fraction"] = composition["n_copies"] / composition.groupby("unit_uid")[
        "n_copies"
    ].transform("sum")
    composition.to_csv(args.outdir / "subtype_array_composition.tsv", sep="\t", index=False)
    spatial = spatial_subtype_table(copy_assign)
    spatial.to_csv(args.outdir / "subtype_spatial_structure.tsv", sep="\t", index=False)

    diagnostics = component_diagnostics(
        selected_model, candidate_sites, reference_sequence=reference_sequence
    )
    diagnostics.to_csv(args.outdir / "subtype_diagnostic_sites.tsv", sep="\t", index=False)
    n_strong_diagnostic = int(
        (
            (diagnostics["max_component_base_frequency_difference"] >= 0.50)
            & (diagnostics["n_distinct_component_consensus_bases"] >= 2)
        ).sum()
    )

    aligned_consensus, ungapped_consensus = subtype_consensuses(
        encoded, hard_labels, raw_weights, recommended_k
    )
    write_fasta(
        ((f"Subtype_{i + 1}", sequence) for i, sequence in enumerate(aligned_consensus)),
        args.outdir / "subtype_consensus.aln.fa",
    )
    write_fasta(
        ((f"Subtype_{i + 1}", sequence) for i, sequence in enumerate(ungapped_consensus)),
        args.outdir / "subtype_consensus.fa",
    )

    selected_row = selection.loc[selection["k"] == recommended_k].iloc[0]
    if recommended_k == 1:
        evidence = "no_discrete_subtype_support"
    else:
        strong = (
            selected_row["cv_gain_vs_k1_per_site"] >= 0.002
            and selected_row["delta_bic_vs_k1"] <= -10.0
            and selected_row["stability_mean_ari"] >= args.min_stability_ari
            and selected_row["mean_max_posterior"] >= 0.80
            and n_strong_diagnostic >= 3
            and not search_boundary_reached
        )
        if search_boundary_reached:
            evidence = "provisional_subtype_support_k_unresolved"
        else:
            evidence = "strong_discrete_subtype_support" if strong else "provisional_subtype_support"
    last_valid_row = selection.loc[selection["k"] == max_valid_k].iloc[0]
    conclusion = pd.DataFrame(
        [
            {"metric": "recommended_k", "value": recommended_k},
            {
                "metric": "recommended_k_interpretation",
                "value": (
                    "diagnostic_selection_only_k_unresolved"
                    if search_boundary_reached
                    else "resolved_within_tested_range"
                ),
            },
            {"metric": "selection_rule", "value": selection_rule},
            {"metric": "evidence_class", "value": evidence},
            {"metric": "k_search_status", "value": k_search_status},
            {"metric": "max_k_tested", "value": max_k},
            {"metric": "max_k_passing_component_size", "value": max_valid_k},
            {"metric": "cv_best_k", "value": cv_best_k},
            {"metric": "bic_best_k", "value": bic_best_k},
            {"metric": "suggested_next_max_k", "value": suggested_next_max_k},
            {
                "metric": "last_k_cv_increment_per_site_mean",
                "value": last_valid_row.get("cv_increment_vs_previous_per_site_mean", float("nan")),
            },
            {
                "metric": "last_k_cv_increment_per_site_se",
                "value": last_valid_row.get("cv_increment_vs_previous_per_site_se", float("nan")),
            },
            {"metric": "informative_sites", "value": len(candidate_sites)},
            {"metric": "fixed_reference_length", "value": encoded.shape[1]},
            {"metric": "strong_diagnostic_sites", "value": n_strong_diagnostic},
            {"metric": "cv_gain_vs_k1_per_monomer", "value": selected_row["cv_gain_vs_k1_per_monomer"]},
            {"metric": "cv_gain_vs_k1_per_site", "value": selected_row["cv_gain_vs_k1_per_site"]},
            {"metric": "delta_bic_vs_k1", "value": selected_row["delta_bic_vs_k1"]},
            {"metric": "mean_max_posterior", "value": selected_row["mean_max_posterior"]},
            {"metric": "stability_mean_ari", "value": selected_row["stability_mean_ari"]},
            {"metric": "minimum_component_weight_threshold", "value": min_weight},
            {"metric": "group_by", "value": args.group_by},
        ]
    )
    conclusion.to_csv(args.outdir / "subtype_validation_conclusion.tsv", sep="\t", index=False)

    plot_validation(
        selection,
        hap_assignments,
        subtype_summary,
        diagnostics,
        composition,
        recommended_k,
        args.outdir / "SAT113_subtype_validation",
    )
    eprint(
        f"Diagnostic K={recommended_k}; evidence={evidence}; "
        f"K search={k_search_status}; strong diagnostic sites={n_strong_diagnostic}"
    )


if __name__ == "__main__":
    main()
