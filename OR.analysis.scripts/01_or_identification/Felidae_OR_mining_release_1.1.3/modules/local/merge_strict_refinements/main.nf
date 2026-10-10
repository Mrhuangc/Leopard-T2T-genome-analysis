process MERGE_STRICT_REFINEMENTS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(results), path(statuses), path(manifest)
    path package_src

    output:
    tuple val(meta), path("${meta.id}.strict_refinement.raw.txt"), emit: merged
    tuple val(meta), path("${meta.id}.strict_refinement.summary.json"), emit: summary
    path 'versions.yml', emit: versions

    script:
    def result_args = results.collect { result -> "\"${result}\"" }.join(' ')
    def status_args = statuses.collect { status -> "\"${status}\"" }.join(' ')
    """
    PYTHONPATH="\$PWD" python -m or_mining.merge_strict_refinements \
        --results $result_args \
        --statuses $status_args \
        --manifest "$manifest" \
        --output "${meta.id}.strict_refinement.raw.txt" \
        --summary "${meta.id}.strict_refinement.summary.json"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch "${meta.id}.strict_refinement.raw.txt"
    echo '{"status":"stub"}' > "${meta.id}.strict_refinement.summary.json"
    echo '"${task.process}": {python: "3.11"}' > versions.yml
    """
}
