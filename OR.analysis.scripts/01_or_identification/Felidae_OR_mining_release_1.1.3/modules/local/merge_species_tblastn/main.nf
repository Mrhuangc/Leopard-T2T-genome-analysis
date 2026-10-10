process MERGE_SPECIES_TBLASTN {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/../classify_gpcr/environment.yml"

    input:
    tuple val(meta), path(species_hits)
    path initial_hits
    path package_src
    val max_hits_per_query

    output:
    tuple val(meta), path("${meta.id}.expanded.merged.tblastn.tsv"),        emit: merged
    tuple val(meta), path("${meta.id}.species.filtered.tblastn.tsv"),       emit: filtered
    tuple val(meta), path("${meta.id}.species.kept_query_ids.txt"),         emit: kept_queries
    tuple val(meta), path("${meta.id}.species.outlier_query_ids.txt"),      emit: outliers
    tuple val(meta), path("${meta.id}.species_expansion.summary.json"),     emit: summary
    path 'versions.yml',                                                     emit: versions

    script:
    def prefix = meta.id
    """
    PYTHONPATH="\$PWD" python -m or_mining.merge_species_hits \
        "$initial_hits" \
        "$species_hits" \
        --output-prefix "$prefix" \
        --max-hits-per-query "$max_hits_per_query"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    touch "${prefix}.expanded.merged.tblastn.tsv"
    touch "${prefix}.species.filtered.tblastn.tsv"
    touch "${prefix}.species.kept_query_ids.txt"
    touch "${prefix}.species.outlier_query_ids.txt"
    echo '{"status":"stub"}' > "${prefix}.species_expansion.summary.json"
    echo '"${task.process}": {python: "3.11"}' > versions.yml
    """
}
