process SELECT_EXONERATE_MODELS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(raw_exonerate), path(regions_tsv), path(target_fasta), path(manifest_tsv), path(merged_hits)
    path package_src
    val maximum_exons
    val minimum_intron_bp
    val maximum_intron_bp
    val refinement_flank_bp
    val blast_clearance_bp
    val minimum_blast_segment_bp

    output:
    tuple val(meta), path("${meta.id}.multiexon_models.all.tsv"), emit: all_models
    tuple val(meta), path("${meta.id}.multiexon_refinement_candidates.tsv"), emit: candidates
    tuple val(meta), path("${meta.id}.multiexon_screened_models.gff3"), emit: gff3
    tuple val(meta), path("${meta.id}.multiexon_model_selection.summary.json"), emit: summary
    path 'versions.yml', emit: versions

    script:
    def prefix = meta.id
    """
    PYTHONPATH="\$PWD" python -m or_mining.exonerate_model_selection \
        "$raw_exonerate" \
        "$regions_tsv" \
        "$target_fasta" \
        "$manifest_tsv" \
        "$merged_hits" \
        --output-prefix "$prefix" \
        --maximum-exons "$maximum_exons" \
        --minimum-intron-bp "$minimum_intron_bp" \
        --maximum-intron-bp "$maximum_intron_bp" \
        --refinement-flank-bp "$refinement_flank_bp" \
        --blast-clearance-bp "$blast_clearance_bp" \
        --minimum-blast-segment-bp "$minimum_blast_segment_bp"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    echo -e 'candidate_id\tfilter_status' > "${prefix}.multiexon_models.all.tsv"
    echo -e 'candidate_id\tfilter_status' > "${prefix}.multiexon_refinement_candidates.tsv"
    echo '##gff-version 3' > "${prefix}.multiexon_screened_models.gff3"
    echo '{"status":"stub"}' > "${prefix}.multiexon_model_selection.summary.json"
    echo '"${task.process}": {python: "3.11"}' > versions.yml
    """
}
