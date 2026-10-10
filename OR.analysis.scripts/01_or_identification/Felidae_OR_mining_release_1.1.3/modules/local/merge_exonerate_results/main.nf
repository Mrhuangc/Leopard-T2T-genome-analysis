process MERGE_EXONERATE_RESULTS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/../prepare_exonerate_queries/environment.yml"

    input:
    tuple val(meta), path(exonerate_results), path(chunk_statuses)
    path package_src

    output:
    tuple val(meta), path("${meta.id}.exonerate.raw.txt"), emit: merged
    tuple val(meta), path("${meta.id}.exonerate.raw.summary.json"), emit: summary
    path 'versions.yml', emit: versions

    script:
    def prefix = meta.id
    def result_args = exonerate_results.collect { path -> "\"${path}\"" }.join(' ')
    def status_args = chunk_statuses.collect { path -> "\"${path}\"" }.join(' ')
    """
    PYTHONPATH="\$PWD" python -m or_mining.merge_exonerate \
        --inputs $result_args \
        --statuses $status_args \
        --output "${prefix}.exonerate.raw.txt" \
        --summary "${prefix}.exonerate.raw.summary.json"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    touch "${prefix}.exonerate.raw.txt"
    echo '{"status":"stub"}' > "${prefix}.exonerate.raw.summary.json"
    echo '"${task.process}": {python: "3.11"}' > versions.yml
    """
}
