#!/usr/bin/env python3
"""Test why SAT113 has high within-array but low between-consensus divergence.

The analysis separates four logically different properties that a consensus
sequence alone collapses:

* exact monomer haplotype sharing;
* similarity-cluster (variant-family) sharing at several identity thresholds;
* per-site minor-allele frequency similarity among arrays;
* linkage/co-occurrence of mutations within monomers and spatial arrangement of
  variant families along each genomic array.

It also summarizes joint three-species clusters to distinguish shared/ancestral
variant candidates from species-enriched or species-specific candidates.  The
statistics are descriptive for repeat copies and are not population-genetic FST
estimates from independently sampled individuals.
"""

from __future__ import annotations

import argparse
import gzip
import itertools
import math
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, MutableMapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from sklearn.utils.extmath import randomized_svd


CANONICAL = np.array(list("ACGT"))
BASE_TO_CODE = {"A": 0, "C": 1, "G": 2, "T": 3}


def eprint(*args: object) -> None:
    print(*args, file=sys.stderr)


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r", encoding="utf-8")


def normalize_fasta_id(name: str) -> str:
    return name.split()[0].split(";")[0]


def read_fasta_dict(path: Path) -> Dict[str, str]:
    records: Dict[str, str] = {}
    name: Optional[str] = None
    chunks: List[str] = []
    with open_text(path) as handle:
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


def parse_uc(path: Path) -> Dict[str, str]:
    """Return haplotype -> centroid haplotype from a VSEARCH .uc file."""
    assignment: Dict[str, str] = {}
    with open_text(path) as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 10 or fields[0] not in {"S", "H"}:
                continue
            query = normalize_fasta_id(fields[8])
            target = query if fields[0] == "S" or fields[9] == "*" else normalize_fasta_id(fields[9])
            assignment[query] = target
    if not assignment:
        raise ValueError(f"No S/H assignments found in {path}")
    return assignment


def parse_threshold_path(value: str) -> Tuple[float, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("--uc must be IDENTITY=PATH, for example 0.85=clusters.0.85.uc")
    identity, path = value.split("=", 1)
    try:
        threshold = float(identity)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"Invalid identity threshold: {identity}") from exc
    if threshold <= 0 or threshold > 1:
        raise argparse.ArgumentTypeError("Identity threshold must be in (0, 1]")
    return threshold, Path(path)


def threshold_label(threshold: float) -> str:
    return f"{threshold:.3f}".rstrip("0").rstrip(".")


def encode_alignment(alignment: Mapping[str, str], hap_order: Sequence[str]) -> np.ndarray:
    lengths = {len(alignment[hap]) for hap in hap_order}
    if len(lengths) != 1:
        raise ValueError(f"Aligned FASTA contains unequal sequence lengths: {sorted(lengths)}")
    length = lengths.pop()
    encoded = np.full((len(hap_order), length), -1, dtype=np.int8)
    for i, hap in enumerate(hap_order):
        sequence = alignment[hap]
        for base, code in BASE_TO_CODE.items():
            encoded[i, np.fromiter((char == base for char in sequence), dtype=bool, count=length)] = code
    return encoded


def entropy_from_counts(counts: np.ndarray) -> float:
    total = float(np.sum(counts))
    if total <= 0:
        return float("nan")
    probabilities = counts[counts > 0] / total
    return float(-np.sum(probabilities * np.log2(probabilities)))


def normalized_mutual_information(contingency: np.ndarray) -> float:
    contingency = np.asarray(contingency, dtype=float)
    total = contingency.sum()
    if total <= 0:
        return float("nan")
    row = contingency.sum(axis=1)
    col = contingency.sum(axis=0)
    expected = np.outer(row, col) / total
    mask = contingency > 0
    mi = float(np.sum((contingency[mask] / total) * np.log2(contingency[mask] / expected[mask])))
    h_row = entropy_from_counts(row)
    h_col = entropy_from_counts(col)
    if not np.isfinite(h_row) or not np.isfinite(h_col) or h_row <= 0 or h_col <= 0:
        return 0.0
    return mi / math.sqrt(h_row * h_col)


def js_divergence(p: np.ndarray, q: np.ndarray) -> float:
    """Jensen-Shannon divergence in bits for two probability vectors."""
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)
    p_sum, q_sum = p.sum(), q.sum()
    if p_sum <= 0 or q_sum <= 0:
        return float("nan")
    p = p / p_sum
    q = q / q_sum
    middle = 0.5 * (p + q)
    p_mask = p > 0
    q_mask = q > 0
    return float(
        0.5 * np.sum(p[p_mask] * np.log2(p[p_mask] / middle[p_mask]))
        + 0.5 * np.sum(q[q_mask] * np.log2(q[q_mask] / middle[q_mask]))
    )


def mean_site_jsd(p: np.ndarray, q: np.ndarray, valid: np.ndarray) -> float:
    values: List[float] = []
    for site in np.where(valid)[0]:
        values.append(js_divergence(p[site], q[site]))
    return float(np.nanmean(values)) if values else float("nan")


def base_frequencies(
    encoded: np.ndarray,
    hap_indices: np.ndarray,
    weights: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    length = encoded.shape[1]
    counts = np.zeros((length, 4), dtype=float)
    selected = encoded[hap_indices]
    for base in range(4):
        counts[:, base] = ((selected == base) * weights[:, None]).sum(axis=0)
    totals = counts.sum(axis=1)
    frequencies = np.divide(
        counts,
        totals[:, None],
        out=np.zeros_like(counts),
        where=totals[:, None] > 0,
    )
    return frequencies, totals


def diversity_from_frequencies(freq: np.ndarray, valid: np.ndarray, n: int) -> float:
    if not np.any(valid):
        return float("nan")
    diversity = 1.0 - np.sum(freq[valid] ** 2, axis=1)
    if n > 1:
        diversity *= n / (n - 1)
    return float(np.mean(diversity))


def dxy_from_frequencies(p: np.ndarray, q: np.ndarray, valid: np.ndarray) -> float:
    if not np.any(valid):
        return float("nan")
    return float(np.mean(1.0 - np.sum(p[valid] * q[valid], axis=1)))


def consensus_divergence(p: np.ndarray, q: np.ndarray, valid: np.ndarray) -> float:
    if not np.any(valid):
        return float("nan")
    return float(np.mean(np.argmax(p[valid], axis=1) != np.argmax(q[valid], axis=1)))


def cluster_classification(
    counts_by_species: Mapping[str, int],
    species_totals: Mapping[str, int],
    species_order: Sequence[str],
    min_copies: int,
    min_fraction: float,
) -> Tuple[str, List[str]]:
    present = [
        species
        for species in species_order
        if counts_by_species.get(species, 0) >= min_copies
        and counts_by_species.get(species, 0) / species_totals[species] >= min_fraction
    ]
    if len(present) == len(species_order):
        return f"shared_{len(species_order)}_species", present
    if len(present) >= 2:
        return "shared_" + "_".join(present), present
    if len(present) == 1:
        return present[0] + "_specific", present
    return "rare_unclassified", present


def build_cluster_summary(
    copy_frame: pd.DataFrame,
    variant_column: str,
    centroid_by_variant: Mapping[str, str],
    species_order: Sequence[str],
    min_species_copies: int,
    min_species_fraction: float,
) -> pd.DataFrame:
    species_totals = copy_frame.groupby("species").size().to_dict()
    rows: List[Dict[str, object]] = []
    grouped = copy_frame.groupby(variant_column, sort=False)
    for variant, subset in grouped:
        species_counts = subset.groupby("species").size().to_dict()
        classification, present = cluster_classification(
            species_counts,
            species_totals,
            species_order,
            min_species_copies,
            min_species_fraction,
        )
        counts = np.array([species_counts.get(species, 0) for species in species_order], dtype=float)
        dominant_index = int(np.argmax(counts))
        row: Dict[str, object] = {
            "variant_cluster": variant,
            "centroid_haplotype": centroid_by_variant.get(str(variant), ""),
            "n_copies": len(subset),
            "n_haplotypes": subset["haplotype_id"].nunique(),
            "n_arrays": subset["array_uid"].nunique(),
            "n_chromosomes": subset[["species", "chrom"]].drop_duplicates().shape[0],
            "n_species_observed": int(np.sum(counts > 0)),
            "n_species_qualified": len(present),
            "classification": classification,
            "qualified_species": ",".join(present),
            "dominant_species": species_order[dominant_index],
            "dominant_species_purity": counts[dominant_index] / counts.sum(),
        }
        for species in species_order:
            row[f"n_{species}"] = int(species_counts.get(species, 0))
            row[f"fraction_of_{species}"] = species_counts.get(species, 0) / species_totals[species]
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["n_copies", "variant_cluster"], ascending=[False, True])


def association_nmi(copy_frame: pd.DataFrame, variant_column: str, label_column: str) -> float:
    contingency = pd.crosstab(copy_frame[variant_column], copy_frame[label_column]).to_numpy()
    return normalized_mutual_information(contingency)


def calculate_threshold_metrics(
    copy_frame: pd.DataFrame,
    cluster_summary: pd.DataFrame,
    variant_column: str,
    threshold: float,
) -> Dict[str, object]:
    total = len(copy_frame)
    cross_array = set(cluster_summary.loc[cluster_summary["n_arrays"] >= 2, "variant_cluster"])
    cross_species = set(
        cluster_summary.loc[cluster_summary["n_species_qualified"] >= 2, "variant_cluster"]
    )
    all_species = set(
        cluster_summary.loc[
            cluster_summary["classification"].str.match(r"shared_\d+_species"),
            "variant_cluster",
        ]
    )
    specific = set(
        cluster_summary.loc[cluster_summary["classification"].str.endswith("_specific"), "variant_cluster"]
    )
    return {
        "identity_threshold": threshold,
        "n_clusters": len(cluster_summary),
        "median_cluster_copies": float(cluster_summary["n_copies"].median()),
        "fraction_copies_cross_array_cluster": copy_frame[variant_column].isin(cross_array).mean(),
        "fraction_copies_cross_species_cluster": copy_frame[variant_column].isin(cross_species).mean(),
        "fraction_copies_all_species_cluster": copy_frame[variant_column].isin(all_species).mean(),
        "fraction_copies_species_specific_cluster": copy_frame[variant_column].isin(specific).mean(),
        "cluster_species_nmi": association_nmi(copy_frame, variant_column, "species"),
        "cluster_array_nmi": association_nmi(copy_frame, variant_column, "array_uid"),
        "n_copies": total,
    }


def make_array_profiles(
    copy_frame: pd.DataFrame,
    variant_column: str,
    variant_order: Sequence[str],
    encoded: np.ndarray,
    hap_to_index: Mapping[str, int],
    global_core: np.ndarray,
    min_array_occupancy: float,
) -> Tuple[pd.DataFrame, np.ndarray, np.ndarray, np.ndarray, Dict[str, int]]:
    arrays = (
        copy_frame.groupby("array_uid", sort=False)
        .agg(
            species=("species", "first"),
            unit_key=("unit_key", "first"),
            chrom=("chrom", "first"),
            n_source_arrays=("source_array_key", "nunique"),
            constituent_array_keys=(
                "source_array_key",
                lambda values: ",".join(sorted(set(map(str, values)))),
            ),
        )
        .reset_index()
        .sort_values(["species", "chrom", "unit_key"])
        .reset_index(drop=True)
    )
    array_index = {uid: i for i, uid in enumerate(arrays["array_uid"])}
    variant_index = {variant: i for i, variant in enumerate(variant_order)}
    composition = np.zeros((len(arrays), len(variant_order)), dtype=float)
    frequencies = np.zeros((len(arrays), encoded.shape[1], 4), dtype=float)
    occupancies = np.zeros((len(arrays), encoded.shape[1]), dtype=float)
    n_copies = np.zeros(len(arrays), dtype=int)

    for uid, subset in copy_frame.groupby("array_uid", sort=False):
        i = array_index[uid]
        n_copies[i] = len(subset)
        variant_counts = subset[variant_column].value_counts()
        for variant, count in variant_counts.items():
            composition[i, variant_index[variant]] = count
        hap_counts = subset["haplotype_id"].value_counts()
        hap_indices = np.array([hap_to_index[hap] for hap in hap_counts.index], dtype=int)
        weights = hap_counts.to_numpy(dtype=float)
        freq, totals = base_frequencies(encoded, hap_indices, weights)
        frequencies[i] = freq
        occupancies[i] = totals / len(subset)

    composition = np.divide(
        composition,
        composition.sum(axis=1, keepdims=True),
        out=np.zeros_like(composition),
        where=composition.sum(axis=1, keepdims=True) > 0,
    )
    valid = occupancies >= min_array_occupancy
    valid &= global_core[None, :]
    arrays["n_copies"] = n_copies
    arrays["n_variant_clusters"] = (composition > 0).sum(axis=1)
    arrays["variant_hill_simpson"] = np.divide(
        1.0,
        np.sum(composition**2, axis=1),
        out=np.full(len(arrays), np.nan),
        where=np.sum(composition**2, axis=1) > 0,
    )
    arrays["top_variant_fraction"] = composition.max(axis=1)
    arrays["within_array_pi_pct"] = [
        100.0 * diversity_from_frequencies(frequencies[i], valid[i], int(n_copies[i]))
        for i in range(len(arrays))
    ]
    arrays["n_core_sites"] = valid.sum(axis=1)
    return arrays, composition, frequencies, valid, array_index


def array_pair_table(
    arrays: pd.DataFrame,
    composition: np.ndarray,
    frequencies: np.ndarray,
    valid: np.ndarray,
) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    for i, j in itertools.combinations(range(len(arrays)), 2):
        jointly_valid = valid[i] & valid[j]
        pi_i = arrays.iloc[i]["within_array_pi_pct"] / 100.0
        pi_j = arrays.iloc[j]["within_array_pi_pct"] / 100.0
        dxy = dxy_from_frequencies(frequencies[i], frequencies[j], jointly_valid)
        fst_like = float("nan")
        if np.isfinite(dxy) and dxy > 0:
            fst_like = 1.0 - (pi_i + pi_j) / (2.0 * dxy)
        comp_jsd = js_divergence(composition[i], composition[j])
        row = {
            "array_uid_1": arrays.iloc[i]["array_uid"],
            "array_uid_2": arrays.iloc[j]["array_uid"],
            "species_1": arrays.iloc[i]["species"],
            "species_2": arrays.iloc[j]["species"],
            "chrom_1": arrays.iloc[i]["chrom"],
            "chrom_2": arrays.iloc[j]["chrom"],
            "same_species": arrays.iloc[i]["species"] == arrays.iloc[j]["species"],
            "same_chromosome_label": (
                arrays.iloc[i]["species"] == arrays.iloc[j]["species"]
                and arrays.iloc[i]["chrom"] == arrays.iloc[j]["chrom"]
            ),
            "n_joint_core_sites": int(jointly_valid.sum()),
            "consensus_divergence_pct": 100.0
            * consensus_divergence(frequencies[i], frequencies[j], jointly_valid),
            "copy_pool_dxy_pct": 100.0 * dxy,
            "fst_like_copy_pool": fst_like,
            "allele_frequency_mean_jsd": mean_site_jsd(
                frequencies[i], frequencies[j], jointly_valid
            ),
            "variant_composition_jsd": comp_jsd,
            "variant_composition_js_distance": math.sqrt(comp_jsd)
            if np.isfinite(comp_jsd)
            else float("nan"),
            "variant_bray_curtis_similarity": 1.0
            - 0.5 * float(np.sum(np.abs(composition[i] - composition[j]))),
        }
        rows.append(row)
    return pd.DataFrame(rows)


def spatial_statistics(
    copy_frame: pd.DataFrame,
    variant_column: str,
    permutations: int,
    rng: np.random.Generator,
) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    for uid, subset in copy_frame.groupby("array_uid", sort=False):
        subset = subset.sort_values(["start", "end", "copy_id"])
        labels = subset[variant_column].to_numpy()
        n = len(labels)
        if n < 2:
            continue
        observed = float(np.mean(labels[:-1] == labels[1:]))
        counts = np.array(list(Counter(labels).values()), dtype=float)
        expected = float(np.sum(counts * (counts - 1)) / (n * (n - 1)))
        run_count = int(1 + np.sum(labels[:-1] != labels[1:]))
        permutation_values = np.empty(permutations, dtype=float)
        for replicate in range(permutations):
            shuffled = rng.permutation(labels)
            permutation_values[replicate] = np.mean(shuffled[:-1] == shuffled[1:])
        p_upper = (1 + np.sum(permutation_values >= observed)) / (permutations + 1)
        p_lower = (1 + np.sum(permutation_values <= observed)) / (permutations + 1)
        first = subset.iloc[0]
        rows.append(
            {
                "array_uid": uid,
                "species": first["species"],
                "unit_key": first["unit_key"],
                "chrom": first["chrom"],
                "n_source_arrays": subset["source_array_key"].nunique(),
                "n_copies": n,
                "n_variant_clusters": len(counts),
                "adjacent_same_observed": observed,
                "adjacent_same_expected_random_order": expected,
                "adjacent_same_excess": observed - expected,
                "permutation_p_spatial_clustering": p_upper,
                "permutation_p_spatial_alternation": p_lower,
                "run_count": run_count,
                "run_fraction_of_max": run_count / n,
            }
        )
    return pd.DataFrame(rows)


def ld_for_group(
    group_name: str,
    encoded: np.ndarray,
    weights: np.ndarray,
    core: np.ndarray,
    min_maf: float,
    max_sites: int,
) -> pd.DataFrame:
    hap_indices = np.where(weights > 0)[0]
    if len(hap_indices) == 0:
        return pd.DataFrame()
    matrix = encoded[hap_indices]
    w = weights[hap_indices].astype(float)
    freq, totals = base_frequencies(encoded, hap_indices, w)
    major = np.argmax(freq, axis=1)
    major_freq = np.max(freq, axis=1)
    maf = 1.0 - major_freq
    occupancy = totals / w.sum()
    sites = np.where(core & (occupancy >= 0.90) & (maf >= min_maf))[0]
    if len(sites) > max_sites:
        # Keep the most informative sites, closest to 50:50.
        sites = sites[np.argsort(np.abs(maf[sites] - 0.5))[:max_sites]]
        sites.sort()
    rows: List[Dict[str, object]] = []
    for site1, site2 in itertools.combinations(sites, 2):
        x = matrix[:, site1]
        y = matrix[:, site2]
        valid = (x >= 0) & (y >= 0)
        if valid.sum() < 2 or w[valid].sum() < 5:
            continue
        xv, yv, wv = x[valid], y[valid], w[valid]
        contingency = np.bincount(xv * 4 + yv, weights=wv, minlength=16).reshape(4, 4)
        nmi = normalized_mutual_information(contingency)
        xb = (xv != major[site1]).astype(float)
        yb = (yv != major[site2]).astype(float)
        total_w = wv.sum()
        mean_x = float(np.sum(wv * xb) / total_w)
        mean_y = float(np.sum(wv * yb) / total_w)
        var_x = float(np.sum(wv * (xb - mean_x) ** 2) / total_w)
        var_y = float(np.sum(wv * (yb - mean_y) ** 2) / total_w)
        r2 = float("nan")
        if var_x > 0 and var_y > 0:
            covariance = float(np.sum(wv * (xb - mean_x) * (yb - mean_y)) / total_w)
            r2 = covariance * covariance / (var_x * var_y)
        rows.append(
            {
                "group": group_name,
                "site_1_alignment_1based": int(site1 + 1),
                "site_2_alignment_1based": int(site2 + 1),
                "maf_site_1": float(maf[site1]),
                "maf_site_2": float(maf[site2]),
                "categorical_normalized_mutual_information": nmi,
                "derived_vs_major_r_squared": r2,
                "weighted_comparable_copies": int(total_w),
            }
        )
    return pd.DataFrame(rows)


def species_profiles(
    copy_frame: pd.DataFrame,
    species_order: Sequence[str],
    variant_column: str,
    variant_order: Sequence[str],
    encoded: np.ndarray,
    hap_to_index: Mapping[str, int],
    core: np.ndarray,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    profiles: Dict[str, Dict[str, object]] = {}
    variant_index = {variant: i for i, variant in enumerate(variant_order)}
    for species in species_order:
        subset = copy_frame.loc[copy_frame["species"] == species]
        hap_counts = subset["haplotype_id"].value_counts()
        indices = np.array([hap_to_index[hap] for hap in hap_counts.index], dtype=int)
        weights = hap_counts.to_numpy(dtype=float)
        freq, totals = base_frequencies(encoded, indices, weights)
        valid = core & (totals / len(subset) >= 0.90)
        composition = np.zeros(len(variant_order), dtype=float)
        for variant, count in subset[variant_column].value_counts().items():
            composition[variant_index[variant]] = count
        composition /= composition.sum()
        pi = diversity_from_frequencies(freq, valid, len(subset))
        profiles[species] = {
            "freq": freq,
            "valid": valid,
            "composition": composition,
            "pi": pi,
            "n": len(subset),
        }

    summary_rows = [
        {
            "species": species,
            "n_copies": profiles[species]["n"],
            "within_species_copy_pi_pct": 100.0 * profiles[species]["pi"],
            "n_core_sites": int(np.sum(profiles[species]["valid"])),
        }
        for species in species_order
    ]
    pair_rows: List[Dict[str, object]] = []
    for first, second in itertools.combinations(species_order, 2):
        a, b = profiles[first], profiles[second]
        valid = a["valid"] & b["valid"]
        dxy = dxy_from_frequencies(a["freq"], b["freq"], valid)
        mean_pi = (a["pi"] + b["pi"]) / 2.0
        comp_a, comp_b = a["composition"], b["composition"]
        occupied_a, occupied_b = comp_a > 0, comp_b > 0
        union = int(np.sum(occupied_a | occupied_b))
        intersection = int(np.sum(occupied_a & occupied_b))
        pair_rows.append(
            {
                "species_1": first,
                "species_2": second,
                "n_joint_core_sites": int(valid.sum()),
                "consensus_divergence_pct": 100.0
                * consensus_divergence(a["freq"], b["freq"], valid),
                "copy_pool_dxy_pct": 100.0 * dxy,
                "mean_within_species_pi_pct": 100.0 * mean_pi,
                "net_divergence_da_pct": 100.0 * (dxy - mean_pi),
                "fst_like_copy_pool": 1.0 - mean_pi / dxy if dxy > 0 else float("nan"),
                "allele_frequency_mean_jsd": mean_site_jsd(a["freq"], b["freq"], valid),
                "variant_composition_jsd": js_divergence(comp_a, comp_b),
                "variant_bray_curtis_similarity": 1.0
                - 0.5 * float(np.sum(np.abs(comp_a - comp_b))),
                "variant_cluster_presence_jaccard": intersection / union if union else float("nan"),
                "n_variant_clusters_shared": intersection,
                "n_variant_clusters_union": union,
            }
        )
    return pd.DataFrame(summary_rows), pd.DataFrame(pair_rows)


def rarefaction_sharing(
    copy_frame: pd.DataFrame,
    species_order: Sequence[str],
    variant_column: str,
    threshold: float,
    replicates: int,
    min_cluster_copies: int,
    rng: np.random.Generator,
) -> pd.DataFrame:
    labels = {
        species: copy_frame.loc[copy_frame["species"] == species, variant_column].to_numpy()
        for species in species_order
    }
    sample_size = min(len(values) for values in labels.values())
    rows: List[Dict[str, object]] = []
    for replicate in range(1, replicates + 1):
        present: Dict[str, set[str]] = {}
        for species, values in labels.items():
            chosen = rng.choice(values, size=sample_size, replace=False)
            counts = Counter(chosen)
            present[species] = {
                cluster for cluster, count in counts.items() if count >= min_cluster_copies
            }
        all_shared = set.intersection(*(present[species] for species in species_order))
        row: Dict[str, object] = {
            "identity_threshold": threshold,
            "replicate": replicate,
            "copies_rarefied_per_species": sample_size,
            "n_clusters_shared_all_species": len(all_shared),
        }
        for species in species_order:
            row[f"n_clusters_{species}"] = len(present[species])
        for first, second in itertools.combinations(species_order, 2):
            union = present[first] | present[second]
            intersection = present[first] & present[second]
            row[f"shared_{first}_{second}"] = len(intersection)
            row[f"jaccard_{first}_{second}"] = len(intersection) / len(union) if union else float("nan")
        rows.append(row)
    return pd.DataFrame(rows)


def ortholog_locus_test(
    pair_frame: pd.DataFrame,
    locus_map_path: Path,
    permutations: int,
    rng: np.random.Generator,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Compare cross-species array pairs assigned to the same syntenic locus.

    The locus map must contain species, unit_key, and ortholog_group.  Labels
    are permuted within species, preserving how many arrays each species assigns
    to every group.
    """
    locus = pd.read_csv(locus_map_path, sep="\t", dtype=str)
    required = {"species", "unit_key", "ortholog_group"}
    missing = sorted(required - set(locus.columns))
    if missing:
        raise ValueError(f"{locus_map_path}: missing columns: {', '.join(missing)}")
    locus["array_uid"] = locus["species"] + "|" + locus["unit_key"]
    if locus["array_uid"].duplicated().any():
        raise ValueError("Each species/unit_key may occur only once in the locus map")
    uid_to_group = dict(zip(locus["array_uid"], locus["ortholog_group"]))
    annotated = pair_frame.copy()
    annotated["ortholog_group_1"] = annotated["array_uid_1"].map(uid_to_group)
    annotated["ortholog_group_2"] = annotated["array_uid_2"].map(uid_to_group)
    annotated["both_loci_mapped"] = annotated["ortholog_group_1"].notna() & annotated[
        "ortholog_group_2"
    ].notna()
    annotated["same_ortholog_group"] = (
        annotated["both_loci_mapped"]
        & (annotated["ortholog_group_1"] == annotated["ortholog_group_2"])
        & (~annotated["same_species"])
    )
    cross = annotated.loc[annotated["both_loci_mapped"] & (~annotated["same_species"])].copy()
    metrics = [
        "copy_pool_dxy_pct",
        "allele_frequency_mean_jsd",
        "variant_composition_jsd",
        "variant_bray_curtis_similarity",
    ]
    rows: List[Dict[str, object]] = []
    for metric in metrics:
        same = pd.to_numeric(cross.loc[cross["same_ortholog_group"], metric], errors="coerce").dropna()
        different = pd.to_numeric(cross.loc[~cross["same_ortholog_group"], metric], errors="coerce").dropna()
        observed = float(same.median() - different.median()) if len(same) and len(different) else float("nan")
        permutation_deltas: List[float] = []
        if np.isfinite(observed):
            base = locus[["array_uid", "species", "ortholog_group"]].copy()
            for _ in range(permutations):
                shuffled_map: Dict[str, str] = {}
                for species, subset in base.groupby("species", sort=False):
                    labels = rng.permutation(subset["ortholog_group"].to_numpy())
                    shuffled_map.update(zip(subset["array_uid"], labels))
                group1 = cross["array_uid_1"].map(shuffled_map)
                group2 = cross["array_uid_2"].map(shuffled_map)
                same_mask = group1 == group2
                perm_same = pd.to_numeric(cross.loc[same_mask, metric], errors="coerce").dropna()
                perm_diff = pd.to_numeric(cross.loc[~same_mask, metric], errors="coerce").dropna()
                if len(perm_same) and len(perm_diff):
                    permutation_deltas.append(float(perm_same.median() - perm_diff.median()))
        p_lower = float("nan")
        p_upper = float("nan")
        if permutation_deltas:
            values = np.asarray(permutation_deltas)
            p_lower = (1 + np.sum(values <= observed)) / (len(values) + 1)
            p_upper = (1 + np.sum(values >= observed)) / (len(values) + 1)
        rows.append(
            {
                "metric": metric,
                "n_same_ortholog_pairs": len(same),
                "n_different_locus_pairs": len(different),
                "median_same_ortholog": same.median() if len(same) else float("nan"),
                "median_different_locus": different.median() if len(different) else float("nan"),
                "same_minus_different": observed,
                "permutation_p_lower": p_lower,
                "permutation_p_upper": p_upper,
                "n_valid_permutations": len(permutation_deltas),
            }
        )
    return annotated, pd.DataFrame(rows)


def pca_haplotypes(
    encoded: np.ndarray,
    core: np.ndarray,
    hap_order: Sequence[str],
    hap_copy_counts: np.ndarray,
    hap_summary: pd.DataFrame,
    n_components: int,
    seed: int,
) -> pd.DataFrame:
    sites = np.where(core)[0]
    one_hot = np.zeros((len(hap_order), len(sites) * 4), dtype=np.float32)
    for base in range(4):
        one_hot[:, base::4] = encoded[:, sites] == base
    mean = np.average(one_hot, axis=0, weights=hap_copy_counts)
    centered = one_hot - mean
    weighted = centered * np.sqrt(hap_copy_counts[:, None])
    components = min(n_components, weighted.shape[0] - 1, weighted.shape[1])
    if components < 1:
        return hap_summary.copy()
    _u, singular_values, components_matrix = randomized_svd(
        weighted, n_components=components, random_state=seed
    )
    scores = centered @ components_matrix.T
    total_weighted_ss = float(np.sum(weighted**2))
    explained_ratios = (
        singular_values**2 / total_weighted_ss
        if total_weighted_ss > 0
        else np.full(components, np.nan)
    )
    result = hap_summary.copy()
    for i in range(components):
        result[f"PC{i + 1}"] = scores[:, i]
        result[f"PC{i + 1}_explained_variance_ratio"] = explained_ratios[i]
    return result


def exact_haplotype_species_summary(
    copy_frame: pd.DataFrame,
    species_order: Sequence[str],
) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    for hap, subset in copy_frame.groupby("haplotype_id", sort=False):
        counts = subset.groupby("species").size().to_dict()
        row: Dict[str, object] = {
            "haplotype_id": hap,
            "n_copies": len(subset),
            "n_arrays": subset["array_uid"].nunique(),
            "n_species": subset["species"].nunique(),
            "species_list": ",".join(sorted(subset["species"].unique())),
        }
        for species in species_order:
            row[f"n_{species}"] = counts.get(species, 0)
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["n_copies", "haplotype_id"], ascending=[False, True])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--copy-haplotypes", required=True, type=Path)
    parser.add_argument("--haplotype-summary", required=True, type=Path)
    parser.add_argument("--haplotype-fasta", required=True, type=Path)
    parser.add_argument("--alignment", required=True, type=Path)
    parser.add_argument(
        "--uc",
        action="append",
        required=True,
        type=parse_threshold_path,
        metavar="IDENTITY=PATH",
    )
    parser.add_argument("--main-threshold", type=float, default=0.85)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--min-global-occupancy", type=float, default=0.90)
    parser.add_argument("--min-array-occupancy", type=float, default=0.80)
    parser.add_argument("--min-species-copies", type=int, default=3)
    parser.add_argument("--min-species-fraction", type=float, default=0.0005)
    parser.add_argument("--min-ld-maf", type=float, default=0.05)
    parser.add_argument("--max-ld-sites", type=int, default=160)
    parser.add_argument("--spatial-permutations", type=int, default=200)
    parser.add_argument("--rarefaction-replicates", type=int, default=100)
    parser.add_argument("--rarefaction-min-cluster-copies", type=int, default=3)
    parser.add_argument("--max-tree-clusters", type=int, default=240)
    parser.add_argument("--min-tree-cluster-copies", type=int, default=10)
    parser.add_argument(
        "--group-by",
        choices=["array_key", "region_key"],
        default="region_key",
        help="Biological comparison unit; region_key avoids counting split pieces of one locus as arrays.",
    )
    parser.add_argument(
        "--locus-map",
        type=Path,
        help="Optional TSV with species,unit_key,ortholog_group for a syntenic-locus test.",
    )
    parser.add_argument("--locus-permutations", type=int, default=1000)
    parser.add_argument(
        "--exclude-core-edge-sites",
        type=int,
        default=0,
        help="Exclude this many occupied sites from each phased monomer edge (sensitivity test).",
    )
    parser.add_argument("--seed", type=int, default=20260721)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    metadata = pd.read_csv(args.metadata, sep="\t", dtype={"copy_id": str})
    copy_hap = pd.read_csv(args.copy_haplotypes, sep="\t")
    hap_summary_input = pd.read_csv(args.haplotype_summary, sep="\t")
    if "haplotype_id" in metadata.columns:
        check = metadata[["seq_id", "haplotype_id"]].merge(
            copy_hap,
            on="seq_id",
            how="left",
            validate="one_to_one",
            suffixes=("_metadata", "_mapping"),
        )
        if check["haplotype_id_mapping"].isna().any() or not (
            check["haplotype_id_metadata"] == check["haplotype_id_mapping"]
        ).all():
            raise ValueError("haplotype assignments disagree between metadata and copy_haplotype.tsv")
        copy_frame = metadata.copy()
    else:
        copy_frame = metadata.merge(copy_hap, on="seq_id", how="inner", validate="one_to_one")
        if len(copy_frame) != len(metadata):
            raise ValueError("copy_haplotype.tsv does not cover every metadata row")
    if args.group_by not in copy_frame.columns:
        raise ValueError(f"Metadata lacks requested grouping column: {args.group_by}")
    if copy_frame[args.group_by].isna().any() or copy_frame[args.group_by].astype(str).eq("").any():
        raise ValueError(f"Metadata contains empty {args.group_by} values")
    copy_frame["source_array_key"] = copy_frame["array_key"].astype(str)
    copy_frame["unit_key"] = copy_frame[args.group_by].astype(str)
    copy_frame["array_uid"] = copy_frame["species"].astype(str) + "|" + copy_frame["unit_key"]
    species_order = list(dict.fromkeys(copy_frame["species"].astype(str)))
    if len(species_order) < 2:
        raise ValueError("At least two species are required for origin analysis")

    hap_sequences = read_fasta_dict(args.haplotype_fasta)
    alignment = read_fasta_dict(args.alignment)
    hap_order = list(hap_summary_input["haplotype_id"].astype(str))
    missing_sequence = sorted(set(hap_order) - set(hap_sequences))
    missing_alignment = sorted(set(hap_order) - set(alignment))
    if missing_sequence or missing_alignment:
        raise ValueError(
            f"Missing haplotypes: raw_fasta={missing_sequence[:5]}, alignment={missing_alignment[:5]}"
        )
    hap_to_index = {hap: i for i, hap in enumerate(hap_order)}
    encoded = encode_alignment(alignment, hap_order)
    hap_copy_counts_series = copy_frame["haplotype_id"].value_counts()
    hap_copy_counts = np.array([hap_copy_counts_series.get(hap, 0) for hap in hap_order], dtype=float)
    global_occupancy = ((encoded >= 0) * hap_copy_counts[:, None]).sum(axis=0) / hap_copy_counts.sum()
    global_core = global_occupancy >= args.min_global_occupancy
    if args.exclude_core_edge_sites > 0:
        occupied_positions = np.where(global_core)[0]
        edge = args.exclude_core_edge_sites
        if len(occupied_positions) <= 2 * edge + 20:
            raise ValueError("Too few core sites remain after --exclude-core-edge-sites")
        global_core[occupied_positions[:edge]] = False
        global_core[occupied_positions[-edge:]] = False
    if global_core.sum() < 20:
        raise ValueError(
            f"Only {global_core.sum()} alignment columns pass global occupancy; inspect phasing/alignment"
        )

    uc_items = sorted(args.uc, key=lambda item: item[0])
    thresholds = [item[0] for item in uc_items]
    if not any(abs(threshold - args.main_threshold) < 1e-9 for threshold in thresholds):
        raise ValueError(f"--main-threshold {args.main_threshold} has no matching --uc input")

    threshold_metrics: List[Dict[str, object]] = []
    cluster_summaries: Dict[float, pd.DataFrame] = {}
    variant_columns: Dict[float, str] = {}
    rarefaction_tables: List[pd.DataFrame] = []
    assignment_columns = [
        "seq_id",
        "species",
        "array_uid",
        "unit_key",
        "source_array_key",
        "haplotype_id",
    ]

    for threshold, uc_path in uc_items:
        eprint(f"Analyzing identity threshold {threshold_label(threshold)}")
        assignment = parse_uc(uc_path)
        missing = sorted(set(hap_order) - set(assignment))
        if missing:
            raise ValueError(f"{uc_path}: {len(missing)} haplotypes lack cluster assignment")
        centroids = sorted(set(assignment.values()))
        centroid_to_variant = {
            centroid: f"VAR{int(round(threshold * 100)):02d}_{i:05d}"
            for i, centroid in enumerate(centroids, 1)
        }
        hap_to_variant = {hap: centroid_to_variant[centroid] for hap, centroid in assignment.items()}
        variant_to_centroid = {variant: centroid for centroid, variant in centroid_to_variant.items()}
        column = f"variant_id_{threshold_label(threshold).replace('.', '_')}"
        variant_columns[threshold] = column
        copy_frame[column] = copy_frame["haplotype_id"].map(hap_to_variant)
        assignment_columns.append(column)
        summary = build_cluster_summary(
            copy_frame,
            column,
            variant_to_centroid,
            species_order,
            args.min_species_copies,
            args.min_species_fraction,
        )
        cluster_summaries[threshold] = summary
        summary.to_csv(
            args.outdir / f"variant_cluster_summary.identity_{threshold_label(threshold)}.tsv",
            sep="\t",
            index=False,
        )
        threshold_metrics.append(calculate_threshold_metrics(copy_frame, summary, column, threshold))
        rarefaction_tables.append(
            rarefaction_sharing(
                copy_frame,
                species_order,
                column,
                threshold,
                args.rarefaction_replicates,
                args.rarefaction_min_cluster_copies,
                rng,
            )
        )

    copy_frame[assignment_columns].to_csv(
        args.outdir / "copy_variant_assignments.tsv", sep="\t", index=False
    )
    threshold_metrics_frame = pd.DataFrame(threshold_metrics)
    threshold_metrics_frame.to_csv(
        args.outdir / "threshold_sensitivity.tsv", sep="\t", index=False
    )
    pd.concat(rarefaction_tables, ignore_index=True).to_csv(
        args.outdir / "rarefaction_species_sharing.tsv", sep="\t", index=False
    )

    main_threshold = min(thresholds, key=lambda value: abs(value - args.main_threshold))
    main_column = variant_columns[main_threshold]
    main_summary = cluster_summaries[main_threshold]
    variant_order = list(main_summary["variant_cluster"])
    arrays, composition, frequencies, valid, array_index = make_array_profiles(
        copy_frame,
        main_column,
        variant_order,
        encoded,
        hap_to_index,
        global_core,
        args.min_array_occupancy,
    )

    cluster_flags = main_summary.set_index("variant_cluster")[["n_arrays", "n_species_qualified"]]
    cross_array_variants = set(cluster_flags.index[cluster_flags["n_arrays"] >= 2])
    cross_species_variants = set(cluster_flags.index[cluster_flags["n_species_qualified"] >= 2])
    array_extra: List[Dict[str, object]] = []
    for uid, subset in copy_frame.groupby("array_uid", sort=False):
        counts = subset[main_column].value_counts(normalize=True)
        array_extra.append(
            {
                "array_uid": uid,
                "fraction_copies_in_cross_array_clusters": counts[
                    counts.index.isin(cross_array_variants)
                ].sum(),
                "fraction_copies_in_cross_species_clusters": counts[
                    counts.index.isin(cross_species_variants)
                ].sum(),
                "n_exact_haplotypes": subset["haplotype_id"].nunique(),
                "exact_haplotype_hill_simpson": 1.0
                / np.sum(subset["haplotype_id"].value_counts(normalize=True).to_numpy() ** 2),
            }
        )
    arrays = arrays.merge(pd.DataFrame(array_extra), on="array_uid", how="left")

    spatial = spatial_statistics(
        copy_frame, main_column, args.spatial_permutations, rng
    )
    arrays = arrays.merge(
        spatial,
        on=["array_uid", "species", "unit_key", "chrom", "n_copies", "n_source_arrays"],
        how="left",
    )
    arrays.to_csv(args.outdir / "array_variant_summary.tsv", sep="\t", index=False)
    spatial.to_csv(args.outdir / "array_spatial_variant_structure.tsv", sep="\t", index=False)

    composition_rows: List[Dict[str, object]] = []
    for i, array_row in arrays.iterrows():
        for j, variant in enumerate(variant_order):
            if composition[i, j] > 0:
                composition_rows.append(
                    {
                        "array_uid": array_row["array_uid"],
                        "species": array_row["species"],
                        "unit_key": array_row["unit_key"],
                        "chrom": array_row["chrom"],
                        "variant_cluster": variant,
                        "fraction": composition[i, j],
                        "n_copies": int(round(composition[i, j] * array_row["n_copies"])),
                    }
                )
    pd.DataFrame(composition_rows).to_csv(
        args.outdir / "array_variant_composition.long.tsv", sep="\t", index=False
    )
    pair_frame = array_pair_table(arrays, composition, frequencies, valid)
    pair_frame.to_csv(args.outdir / "array_pair_similarity.tsv", sep="\t", index=False)
    if args.locus_map is not None:
        annotated_pairs, locus_test = ortholog_locus_test(
            pair_frame, args.locus_map, args.locus_permutations, rng
        )
        annotated_pairs.to_csv(
            args.outdir / "array_pair_similarity.with_ortholog_loci.tsv", sep="\t", index=False
        )
        locus_test.to_csv(
            args.outdir / "ortholog_locus_origin_test.tsv", sep="\t", index=False
        )

    species_summary, species_pairs = species_profiles(
        copy_frame,
        species_order,
        main_column,
        variant_order,
        encoded,
        hap_to_index,
        global_core,
    )
    species_summary.to_csv(args.outdir / "species_copy_pool_summary.tsv", sep="\t", index=False)
    species_pairs.to_csv(args.outdir / "species_pair_origin_metrics.tsv", sep="\t", index=False)

    hap_species = exact_haplotype_species_summary(copy_frame, species_order)
    hap_species.to_csv(args.outdir / "exact_haplotype_species_sharing.tsv", sep="\t", index=False)

    # Mutation linkage: stable divergent variants should create correlated
    # multi-site substitutions, whereas independent old mutations should not.
    ld_tables: List[pd.DataFrame] = []
    all_weights = hap_copy_counts.copy()
    ld_tables.append(
        ld_for_group(
            "ALL",
            encoded,
            all_weights,
            global_core,
            args.min_ld_maf,
            args.max_ld_sites,
        )
    )
    for species in species_order:
        species_counts = copy_frame.loc[copy_frame["species"] == species, "haplotype_id"].value_counts()
        weights = np.array([species_counts.get(hap, 0) for hap in hap_order], dtype=float)
        ld_tables.append(
            ld_for_group(
                species,
                encoded,
                weights,
                global_core,
                args.min_ld_maf,
                args.max_ld_sites,
            )
        )
    nonempty_ld = [table for table in ld_tables if not table.empty]
    if nonempty_ld:
        ld_frame = pd.concat(nonempty_ld, ignore_index=True)
    else:
        ld_frame = pd.DataFrame(
            columns=[
                "group",
                "site_1_alignment_1based",
                "site_2_alignment_1based",
                "maf_site_1",
                "maf_site_2",
                "categorical_normalized_mutual_information",
                "derived_vs_major_r_squared",
                "weighted_comparable_copies",
            ]
        )
    ld_frame.to_csv(args.outdir / "site_pair_linkage.tsv", sep="\t", index=False)
    if not ld_frame.empty:
        ld_summary = (
            ld_frame.groupby("group")
            .agg(
                n_site_pairs=("categorical_normalized_mutual_information", "size"),
                median_nmi=("categorical_normalized_mutual_information", "median"),
                q95_nmi=("categorical_normalized_mutual_information", lambda x: x.quantile(0.95)),
                fraction_nmi_ge_0_2=("categorical_normalized_mutual_information", lambda x: np.mean(x >= 0.2)),
                median_r_squared=("derived_vs_major_r_squared", "median"),
                q95_r_squared=("derived_vs_major_r_squared", lambda x: x.quantile(0.95)),
            )
            .reset_index()
        )
    else:
        ld_summary = pd.DataFrame()
    ld_summary.to_csv(args.outdir / "site_pair_linkage.summary.tsv", sep="\t", index=False)

    # Haplotype PCA is a visualization aid; it is not used to classify origin.
    counts_by_hap_species = pd.crosstab(copy_frame["haplotype_id"], copy_frame["species"])
    hap_annotation = hap_summary_input.set_index("haplotype_id").reindex(hap_order).reset_index()
    for species in species_order:
        hap_annotation[f"n_{species}"] = [
            counts_by_hap_species.loc[hap, species]
            if hap in counts_by_hap_species.index and species in counts_by_hap_species.columns
            else 0
            for hap in hap_order
        ]
    species_count_matrix = hap_annotation[[f"n_{species}" for species in species_order]].to_numpy()
    hap_annotation["dominant_species"] = [
        species_order[int(np.argmax(row))] for row in species_count_matrix
    ]
    hap_annotation["dominant_species_purity"] = np.max(species_count_matrix, axis=1) / np.maximum(
        species_count_matrix.sum(axis=1), 1
    )
    hap_annotation["sharing_class"] = np.where(
        np.sum(species_count_matrix > 0, axis=1) >= 2, "shared_exact_haplotype", "one_species_only"
    )
    pca_frame = pca_haplotypes(
        encoded,
        global_core,
        hap_order,
        hap_copy_counts,
        hap_annotation,
        n_components=5,
        seed=args.seed,
    )
    pca_frame.to_csv(args.outdir / "haplotype_pca.tsv", sep="\t", index=False)

    # Select a balanced set of variant centroids for an optional short-sequence
    # tree/network.  Tree topology must be interpreted cautiously at 113 bp.
    eligible = main_summary.loc[
        main_summary["n_copies"] >= args.min_tree_cluster_copies
    ].copy()
    selected: List[str] = []
    classes = list(eligible["classification"].drop_duplicates())
    per_class = max(1, args.max_tree_clusters // max(1, len(classes)))
    for classification in classes:
        selected.extend(
            eligible.loc[eligible["classification"] == classification, "variant_cluster"]
            .head(per_class)
            .tolist()
        )
    for variant in eligible["variant_cluster"]:
        if len(selected) >= args.max_tree_clusters:
            break
        if variant not in selected:
            selected.append(variant)
    selected = selected[: args.max_tree_clusters]
    tree_meta = main_summary.loc[main_summary["variant_cluster"].isin(selected)].copy()
    tree_meta.to_csv(args.outdir / "tree_tip_metadata.tsv", sep="\t", index=False)
    write_fasta(
        (
            (row.variant_cluster, hap_sequences[row.centroid_haplotype])
            for row in tree_meta.itertuples(index=False)
        ),
        args.outdir / "tree_representatives.fa",
    )

    exact_cross_array_fraction = hap_species.loc[
        hap_species["n_arrays"] >= 2, "n_copies"
    ].sum() / len(copy_frame)
    exact_cross_species_fraction = hap_species.loc[
        hap_species["n_species"] >= 2, "n_copies"
    ].sum() / len(copy_frame)
    main_metrics = threshold_metrics_frame.loc[
        np.isclose(threshold_metrics_frame["identity_threshold"], main_threshold)
    ].iloc[0]
    same_species_pairs = pair_frame.loc[pair_frame["same_species"]]
    different_species_pairs = pair_frame.loc[~pair_frame["same_species"]]
    all_ld = ld_summary.loc[ld_summary["group"] == "ALL"] if not ld_summary.empty else pd.DataFrame()
    diagnostic_rows = [
        {"metric": "main_identity_threshold", "value": main_threshold},
        {"metric": "n_copies", "value": len(copy_frame)},
        {"metric": "n_exact_haplotypes", "value": copy_frame["haplotype_id"].nunique()},
        {"metric": "fraction_copies_in_cross_array_exact_haplotypes", "value": exact_cross_array_fraction},
        {"metric": "fraction_copies_in_cross_species_exact_haplotypes", "value": exact_cross_species_fraction},
        {
            "metric": "fraction_copies_in_cross_array_variant_clusters",
            "value": main_metrics["fraction_copies_cross_array_cluster"],
        },
        {
            "metric": "fraction_copies_in_cross_species_variant_clusters",
            "value": main_metrics["fraction_copies_cross_species_cluster"],
        },
        {
            "metric": "median_same_species_array_consensus_divergence_pct",
            "value": same_species_pairs["consensus_divergence_pct"].median(),
        },
        {
            "metric": "median_same_species_array_copy_pool_dxy_pct",
            "value": same_species_pairs["copy_pool_dxy_pct"].median(),
        },
        {
            "metric": "median_same_species_array_allele_frequency_jsd",
            "value": same_species_pairs["allele_frequency_mean_jsd"].median(),
        },
        {
            "metric": "median_same_species_array_variant_composition_jsd",
            "value": same_species_pairs["variant_composition_jsd"].median(),
        },
        {
            "metric": "median_different_species_array_variant_composition_jsd",
            "value": different_species_pairs["variant_composition_jsd"].median(),
        },
        {
            "metric": "median_spatial_adjacent_same_excess",
            "value": spatial["adjacent_same_excess"].median(),
        },
    ]
    if not all_ld.empty:
        diagnostic_rows.extend(
            [
                {"metric": "global_median_site_pair_nmi", "value": all_ld.iloc[0]["median_nmi"]},
                {"metric": "global_q95_site_pair_nmi", "value": all_ld.iloc[0]["q95_nmi"]},
                {
                    "metric": "global_fraction_site_pairs_nmi_ge_0_2",
                    "value": all_ld.iloc[0]["fraction_nmi_ge_0_2"],
                },
            ]
        )
    pd.DataFrame(diagnostic_rows).to_csv(
        args.outdir / "hypothesis_diagnostic_metrics.tsv", sep="\t", index=False
    )

    core_rows = pd.DataFrame(
        {
            "alignment_position_1based": np.arange(1, encoded.shape[1] + 1),
            "weighted_occupancy": global_occupancy,
            "is_core_column": global_core,
        }
    )
    core_rows.to_csv(args.outdir / "alignment_core_columns.tsv", sep="\t", index=False)
    eprint(
        f"Done: {len(copy_frame):,} copies, {len(hap_order):,} exact haplotypes, "
        f"{int(global_core.sum())} core alignment columns"
    )


if __name__ == "__main__":
    main()
