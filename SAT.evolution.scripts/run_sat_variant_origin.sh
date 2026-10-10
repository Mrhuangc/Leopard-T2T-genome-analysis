#!/usr/bin/env bash
set -euo pipefail

usage() {
    printf '%s\n' \
        "Usage: bash run_sat_variant_origin.sh MANIFEST.tsv OUTDIR [THREADS]" \
        "" \
        "Environment variables:" \
        "  IDENTITIES=0.80,0.85,0.90,0.95   VSEARCH clustering thresholds" \
        "  MAIN_IDENTITY=0.85                threshold used for detailed plots" \
        "  RUN_TREE=0                        run optional MAFFT + IQ-TREE2 tree" \
        "  MIN_LENGTH_RATIO=0.80              complete-copy lower bound" \
        "  MAX_LENGTH_RATIO=1.20              complete-copy upper bound" \
        "  EXCLUDE_CORE_EDGE_SITES=0          set to 5 for a boundary sensitivity run" \
        "  LOCUS_MAP=/path/loci.tsv            optional syntenic-array grouping" \
        "  GROUP_BY=region_key                 merge split array_id pieces by region (recommended)" \
        "  REFERENCE_FASTA=/path/SAT113.fa     optional validated one-monomer reference" \
        "  REFERENCE_REFINEMENT_ROUNDS=2       copy-weighted reference refinement rounds" \
        "  REFERENCE_ALIGNMENT_STRATEGY=multipair  accurate MAFFT mode for divergent copies" \
        "  MIN_ALIGNMENT_COVERAGE=0.80         fixed-reference alignment QC threshold" \
        "  ALLOW_POOR_ALIGNMENT=0              diagnostic override; keep 0 for formal runs" \
        "  CROSS_SPECIES_ORIENTATION_MODE=auto  allow species-level reverse-complement normalization" \
        "  PHASE_SEARCH_RADIUS=12              per-haplotype search around the species phase anchor" \
        "  PHASE_KMER_SIZE=7                   circular k-mer size used to rank phase candidates" \
        "  PHASE_CANDIDATE_COUNT=6             rotations evaluated by banded edit distance" \
        "  PHASE_ALIGNMENT_BAND=16             edit-distance alignment half-band" \
        "  MIN_PHASE_ANCHOR_SIMILARITY=0.70    minimum species-consensus/reference similarity" \
        "  MAX_PHASE_TIED_WEIGHT_FRACTION=0.25 maximum copy-capped weight with tied local phase" \
        "  ALLOW_POOR_CROSS_SPECIES_PHASE=0    diagnostic override; keep 0 for formal runs" \
        "  MAX_SUBTYPES=18                     categorical mixture K=1..Kmax" \
        "  MODEL_STARTS=24                     independent initializations per full model" \
        "  CV_FOLDS=5                          held-out region folds" \
        "  CV_STARTS=4                         initializations per CV fit" \
        "  MIN_CV_INCREMENT_PER_SITE=0.001     practical K-to-K+1 gain threshold" \
        "  SEED=20260721                     random seed" \
        "" \
        "The manifest format is documented in sat113_manifest.example.tsv."
}

if [[ $# -lt 2 || $# -gt 3 ]]; then
    usage >&2
    exit 2
fi

manifest=$1
outdir=$2
threads=${3:-8}
identities=${IDENTITIES:-0.80,0.85,0.90,0.95}
main_identity=${MAIN_IDENTITY:-0.85}
run_tree=${RUN_TREE:-0}
seed=${SEED:-20260721}
min_length_ratio=${MIN_LENGTH_RATIO:-0.80}
max_length_ratio=${MAX_LENGTH_RATIO:-1.20}
exclude_core_edge_sites=${EXCLUDE_CORE_EDGE_SITES:-0}
locus_map=${LOCUS_MAP:-}
group_by=${GROUP_BY:-region_key}
reference_fasta=${REFERENCE_FASTA:-}
reference_refinement_rounds=${REFERENCE_REFINEMENT_ROUNDS:-2}
reference_alignment_strategy=${REFERENCE_ALIGNMENT_STRATEGY:-multipair}
min_alignment_coverage=${MIN_ALIGNMENT_COVERAGE:-0.80}
allow_poor_alignment=${ALLOW_POOR_ALIGNMENT:-0}
cross_species_orientation_mode=${CROSS_SPECIES_ORIENTATION_MODE:-auto}
phase_search_radius=${PHASE_SEARCH_RADIUS:-12}
phase_kmer_size=${PHASE_KMER_SIZE:-7}
phase_candidate_count=${PHASE_CANDIDATE_COUNT:-6}
phase_alignment_band=${PHASE_ALIGNMENT_BAND:-16}
min_phase_anchor_similarity=${MIN_PHASE_ANCHOR_SIMILARITY:-0.70}
max_phase_tied_weight_fraction=${MAX_PHASE_TIED_WEIGHT_FRACTION:-0.25}
allow_poor_cross_species_phase=${ALLOW_POOR_CROSS_SPECIES_PHASE:-0}
max_subtypes=${MAX_SUBTYPES:-18}
model_starts=${MODEL_STARTS:-24}
cv_folds=${CV_FOLDS:-5}
cv_starts=${CV_STARTS:-4}
min_cv_increment_per_site=${MIN_CV_INCREMENT_PER_SITE:-0.001}

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
python_bin=${PYTHON_BIN:-python3}
vsearch_bin=${VSEARCH_BIN:-vsearch}
mafft_bin=${MAFFT_BIN:-mafft}
iqtree_bin=${IQTREE_BIN:-iqtree2}

for executable in "$python_bin" "$vsearch_bin" "$mafft_bin"; do
    if ! command -v "$executable" >/dev/null 2>&1; then
        printf 'Required executable not found: %s\n' "$executable" >&2
        exit 1
    fi
done

mkdir -p "$outdir" "$outdir/prepared" "$outdir/clusters" "$outdir/results" "$outdir/results/subtypes" "$outdir/figures" "$outdir/logs"
export MPLCONFIGDIR=${MPLCONFIGDIR:-"$outdir/.matplotlib"}
mkdir -p "$MPLCONFIGDIR"

printf '[1/8] Preparing copies in each species-specific global phase\n'
"$python_bin" "$script_dir/prepare_sat113_variant_sequences.py" \
    --manifest "$manifest" \
    --outdir "$outdir/prepared" \
    --min-length-ratio "$min_length_ratio" \
    --max-length-ratio "$max_length_ratio" \
    2>"$outdir/logs/01_prepare.log"

printf '[2/8] Inferring a bootstrap SAT113 reference coordinate\n'
reference_args=()
if [[ -n "$reference_fasta" ]]; then
    reference_args+=(--reference-fasta "$reference_fasta")
fi
alignment_qc_args=()
if [[ "$allow_poor_alignment" == "1" ]]; then
    alignment_qc_args+=(--allow-poor-alignment)
fi
"$python_bin" "$script_dir/align_sat113_to_reference.py" \
    --haplotype-fasta "$outdir/prepared/exact_haplotypes.fa" \
    --haplotype-summary "$outdir/prepared/exact_haplotype_summary.tsv" \
    --metadata "$outdir/prepared/copy_metadata.tsv" \
    --output-alignment "$outdir/prepared/exact_haplotypes.bootstrap_reference.aln.fa" \
    --output-reference "$outdir/prepared/SAT113.bootstrap_reference.fa" \
    --output-qc "$outdir/prepared/bootstrap_reference_alignment_qc.tsv" \
    --output-summary "$outdir/prepared/bootstrap_reference_alignment_summary.tsv" \
    --output-coordinate-map "$outdir/prepared/bootstrap_reference_coordinate_map.tsv" \
    --log "$outdir/logs/02_mafft_bootstrap_reference.log" \
    --mafft "$mafft_bin" \
    --mafft-strategy "$reference_alignment_strategy" \
    --threads "$threads" \
    --group-by "$group_by" \
    --refinement-rounds "$reference_refinement_rounds" \
    --min-alignment-coverage "$min_alignment_coverage" \
    --allow-poor-alignment \
    ${reference_args[@]+"${reference_args[@]}"} \
    2>"$outdir/logs/02_bootstrap_reference.wrapper.log"

printf '[3/8] Normalizing all species to the same circular phase\n'
phase_qc_args=()
if [[ "$allow_poor_cross_species_phase" == "1" ]]; then
    phase_qc_args+=(--allow-poor-phase)
fi
"$python_bin" "$script_dir/normalize_sat113_cross_species_phase.py" \
    --metadata "$outdir/prepared/copy_metadata.tsv" \
    --copy-haplotypes "$outdir/prepared/copy_haplotype.tsv" \
    --haplotype-summary "$outdir/prepared/exact_haplotype_summary.tsv" \
    --haplotype-fasta "$outdir/prepared/exact_haplotypes.fa" \
    --reference-fasta "$outdir/prepared/SAT113.bootstrap_reference.fa" \
    --output-metadata "$outdir/prepared/copy_metadata.reference_phase.tsv" \
    --output-copy-haplotypes "$outdir/prepared/copy_haplotype.reference_phase.tsv" \
    --output-haplotype-summary "$outdir/prepared/exact_haplotype_summary.reference_phase.tsv" \
    --output-haplotype-fasta "$outdir/prepared/exact_haplotypes.reference_phase.fa" \
    --output-all-copies-fasta "$outdir/prepared/all_copies.reference_phase.fa" \
    --output-phase-map "$outdir/prepared/cross_species_phase_map.tsv" \
    --output-species-summary "$outdir/prepared/cross_species_phase_species_summary.tsv" \
    --output-summary "$outdir/prepared/cross_species_phase_summary.tsv" \
    --group-by "$group_by" \
    --orientation-mode "$cross_species_orientation_mode" \
    --phase-search-radius "$phase_search_radius" \
    --kmer-size "$phase_kmer_size" \
    --candidate-count "$phase_candidate_count" \
    --alignment-band "$phase_alignment_band" \
    --min-anchor-similarity "$min_phase_anchor_similarity" \
    --max-ambiguous-weight-fraction "$max_phase_tied_weight_fraction" \
    ${phase_qc_args[@]+"${phase_qc_args[@]}"} \
    2>"$outdir/logs/03_cross_species_phase.log"

printf '[4/8] Re-aligning normalized haplotypes to the fixed reference\n'
"$python_bin" "$script_dir/align_sat113_to_reference.py" \
    --haplotype-fasta "$outdir/prepared/exact_haplotypes.reference_phase.fa" \
    --haplotype-summary "$outdir/prepared/exact_haplotype_summary.reference_phase.tsv" \
    --metadata "$outdir/prepared/copy_metadata.reference_phase.tsv" \
    --reference-fasta "$outdir/prepared/SAT113.bootstrap_reference.fa" \
    --output-alignment "$outdir/prepared/exact_haplotypes.reference.aln.fa" \
    --output-reference "$outdir/prepared/SAT113.fixed_reference.fa" \
    --output-qc "$outdir/prepared/reference_alignment_qc.tsv" \
    --output-summary "$outdir/prepared/reference_alignment_summary.tsv" \
    --output-coordinate-map "$outdir/prepared/reference_coordinate_map.tsv" \
    --log "$outdir/logs/04_mafft_reference_alignment.log" \
    --mafft "$mafft_bin" \
    --mafft-strategy "$reference_alignment_strategy" \
    --threads "$threads" \
    --group-by "$group_by" \
    --refinement-rounds 0 \
    --min-alignment-coverage "$min_alignment_coverage" \
    ${alignment_qc_args[@]+"${alignment_qc_args[@]}"} \
    2>"$outdir/logs/04_reference_alignment.wrapper.log"

printf '[5/8] Validating whether discrete multi-site subtypes exist\n'
"$python_bin" "$script_dir/validate_sat113_subtypes.py" \
    --metadata "$outdir/prepared/copy_metadata.reference_phase.tsv" \
    --haplotype-summary "$outdir/prepared/exact_haplotype_summary.reference_phase.tsv" \
    --alignment "$outdir/prepared/exact_haplotypes.reference.aln.fa" \
    --reference-fasta "$outdir/prepared/SAT113.fixed_reference.fa" \
    --group-by "$group_by" \
    --max-k "$max_subtypes" \
    --starts "$model_starts" \
    --cv-folds "$cv_folds" \
    --cv-starts "$cv_starts" \
    --min-cv-increment-per-site "$min_cv_increment_per_site" \
    --exclude-core-edge-sites "$exclude_core_edge_sites" \
    --seed "$seed" \
    --outdir "$outdir/results/subtypes" \
    2>"$outdir/logs/05_subtype_validation.log"

printf '[6/8] Clustering reference-phased exact haplotypes at multiple identities\n'
IFS=',' read -r -a identity_array <<< "$identities"
uc_args=()
for identity in "${identity_array[@]}"; do
    clean_identity=$(printf '%s' "$identity" | sed 's/0*$//; s/\.$//')
    uc_file="$outdir/clusters/identity_${clean_identity}.uc"
    centroid_file="$outdir/clusters/identity_${clean_identity}.centroids.fa"
    "$vsearch_bin" \
        --cluster_size "$outdir/prepared/exact_haplotypes.reference_phase.fa" \
        --id "$identity" \
        --iddef 1 \
        --strand plus \
        --qmask none \
        --query_cov 0.80 \
        --target_cov 0.80 \
        --sizein \
        --sizeout \
        --threads "$threads" \
        --centroids "$centroid_file" \
        --uc "$uc_file" \
        2>"$outdir/logs/06_vsearch_identity_${clean_identity}.log"
    uc_args+=(--uc "${identity}=${uc_file}")
done

printf '[7/8] Testing sharing, linkage, array composition, and spatial order\n'
analysis_extra=()
if [[ -n "$locus_map" ]]; then
    analysis_extra+=(--locus-map "$locus_map")
fi
"$python_bin" "$script_dir/analyze_sat113_variants.py" \
    --metadata "$outdir/prepared/copy_metadata.reference_phase.tsv" \
    --copy-haplotypes "$outdir/prepared/copy_haplotype.reference_phase.tsv" \
    --haplotype-summary "$outdir/prepared/exact_haplotype_summary.reference_phase.tsv" \
    --haplotype-fasta "$outdir/prepared/exact_haplotypes.reference_phase.fa" \
    --alignment "$outdir/prepared/exact_haplotypes.reference.aln.fa" \
    ${uc_args[@]+"${uc_args[@]}"} \
    --main-threshold "$main_identity" \
    --seed "$seed" \
    --exclude-core-edge-sites "$exclude_core_edge_sites" \
    --group-by "$group_by" \
    ${analysis_extra[@]+"${analysis_extra[@]}"} \
    --outdir "$outdir/results" \
    2>"$outdir/logs/07_analysis.log"

printf '[8/8] Plotting diagnostic and origin figures\n'
"$python_bin" "$script_dir/plot_sat113_variant_results.py" \
    --result-dir "$outdir/results" \
    --outdir "$outdir/figures" \
    --main-threshold "$main_identity" \
    2>"$outdir/logs/08_plot.log"

if [[ "$run_tree" == "1" ]]; then
    if ! command -v "$iqtree_bin" >/dev/null 2>&1; then
        printf 'RUN_TREE=1, but IQ-TREE executable was not found: %s\n' "$iqtree_bin" >&2
        exit 1
    fi
    printf '[optional] Building a short-sequence centroid tree\n'
    "$mafft_bin" --thread "$threads" --auto \
        "$outdir/results/tree_representatives.fa" \
        >"$outdir/results/tree_representatives.aln.fa" \
        2>"$outdir/logs/09_mafft_tree.log"
    "$iqtree_bin" \
        -s "$outdir/results/tree_representatives.aln.fa" \
        -m GTR+F+G4 \
        -B 1000 \
        -alrt 1000 \
        -T AUTO \
        --prefix "$outdir/results/SAT113_variant_centroids" \
        >"$outdir/logs/10_iqtree.stdout.log" \
        2>"$outdir/logs/10_iqtree.stderr.log"
fi

printf 'Completed. Key tables: %s\n' "$outdir/results/hypothesis_diagnostic_metrics.tsv"
printf 'Cross-species phase QC: %s\n' "$outdir/prepared/cross_species_phase_summary.tsv"
printf 'Reference-alignment QC: %s\n' "$outdir/prepared/reference_alignment_summary.tsv"
printf 'Subtype validation: %s\n' "$outdir/results/subtypes/subtype_validation_conclusion.tsv"
printf 'Key figures: %s\n' "$outdir/figures/SAT113_variant_diagnostics.pdf"
