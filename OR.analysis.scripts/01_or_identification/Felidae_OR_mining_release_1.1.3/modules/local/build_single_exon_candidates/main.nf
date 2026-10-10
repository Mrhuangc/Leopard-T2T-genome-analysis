process BUILD_SINGLE_EXON_CANDIDATES {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(merged_hits), path(catalog), path(reference_or), path(species_or)
    path package_src
    val clearance
    val flank

    output:
    tuple val(meta), path("${meta.id}.single_exon_residual_candidates.tsv"), path("${meta.id}.single_exon_residual_queries.prot.fa"), path("${meta.id}.single_exon_residual_candidates.faidx.txt"), emit: candidates
    tuple val(meta), path("${meta.id}.single_exon_residual_candidates.summary.json"), emit: summary
    path 'versions.yml', emit: versions

    script:
    def prefix = meta.id
    """
    PYTHONPATH="\$PWD" python -m or_mining.run_single_exon_recovery prepare \
        "$merged_hits" \
        "$catalog" \
        "$reference_or" \
        "$species_or" \
        --species "$prefix" \
        --output-prefix "$prefix" \
        --clearance "$clearance" \
        --flank "$flank"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    echo -e 'candidate_id\tsseqid\tstrand\tquery_id\ttrue_start_1\ttrue_end_1\texpanded_start_1\texpanded_end_1\tslen\talignment_aa\tpident\tevalue\tbitscore\tfaidx_region\tsource_hit_count' > "${prefix}.single_exon_residual_candidates.tsv"
    touch "${prefix}.single_exon_residual_queries.prot.fa"
    touch "${prefix}.single_exon_residual_candidates.faidx.txt"
    echo '{"status":"stub"}' > "${prefix}.single_exon_residual_candidates.summary.json"
    echo '"${task.process}": {python: "3.11"}' > versions.yml
    """
}
