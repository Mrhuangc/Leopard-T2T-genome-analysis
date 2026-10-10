process ANALYZE_OR_LOCUS_ARCHITECTURE {
    tag "$meta.id"
    label 'process_single'

    // The stage-13 parser only needs the already validated Python runtime.
    conda "${moduleDir}/../filter_or_tree/environment.yml"

    input:
    tuple val(meta), path(annotated_catalog)
    path genome_fai
    path package_src
    val cluster_gap_bp

    output:
    tuple val(meta), path("${meta.id}.OR_locus_architecture.tsv"), path("${meta.id}.OR_operational_clusters.tsv"), path("${meta.id}.OR_contig_counts.tsv"), path("${meta.id}.OR_clade_composition.tsv"), emit: tables
    tuple val(meta), path("${meta.id}.OR_locus_architecture.summary.json"), emit: summary
    path 'versions.yml', emit: versions

    script:
    def prefix = meta.id
    """
    PYTHONPATH="\$PWD" python -m or_mining.run_or_locus_architecture \
        "$annotated_catalog" \
        "$genome_fai" \
        --species "${meta.id}" \
        --output-prefix "$prefix" \
        --cluster-gap-bp "$cluster_gap_bp"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    touch "${prefix}.OR_locus_architecture.tsv" "${prefix}.OR_operational_clusters.tsv"
    touch "${prefix}.OR_contig_counts.tsv" "${prefix}.OR_clade_composition.tsv"
    echo '{"status":"stub"}' > "${prefix}.OR_locus_architecture.summary.json"
    echo '"${task.process}": {python: "3.11"}' > versions.yml
    """
}
