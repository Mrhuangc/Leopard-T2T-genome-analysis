process PREPARE_EXONERATE_QUERIES {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(reference_or), path(species_or)
    path package_src
    val requested_chunks

    output:
    tuple val(meta), path('query_chunks/chunk_*.prot.fa'), emit: chunks
    tuple val(meta), path("${meta.id}.exonerate_queries.summary.json"), emit: summary
    path 'versions.yml', emit: versions

    script:
    def prefix = meta.id
    """
    PYTHONPATH="\$PWD" python -m or_mining.prepare_exonerate_queries \
        "$reference_or" \
        "$species_or" \
        --output-dir query_chunks \
        --summary "${prefix}.exonerate_queries.summary.json" \
        --chunks "$requested_chunks"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    mkdir -p query_chunks
    touch query_chunks/chunk_001.prot.fa
    echo '{"status":"stub"}' > "${prefix}.exonerate_queries.summary.json"
    echo '"${task.process}": {python: "3.11"}' > versions.yml
    """
}
