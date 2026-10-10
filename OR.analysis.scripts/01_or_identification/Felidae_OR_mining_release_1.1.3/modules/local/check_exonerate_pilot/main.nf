process CHECK_EXONERATE_PILOT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(results), path(statuses)
    path package_src
    val expected_pairs
    val minimum_positive_pairs

    output:
    path "${meta.id}.exonerate_pilot_gate.ok", emit: gate
    tuple val(meta), path("${meta.id}.exonerate_pilot.summary.json"), emit: summary
    path 'versions.yml', emit: versions

    script:
    def prefix = meta.id
    def result_args = results.collect { result -> "\"${result}\"" }.join(' ')
    def status_args = statuses.collect { status -> "\"${status}\"" }.join(' ')
    """
    PYTHONPATH="\$PWD" python -m or_mining.check_exonerate_pilot \
        --results $result_args \
        --statuses $status_args \
        --output "${prefix}.exonerate_pilot.summary.json" \
        --gate "${prefix}.exonerate_pilot_gate.ok" \
        --expected-pairs "$expected_pairs" \
        --minimum-positive-pairs "$minimum_positive_pairs"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    echo 'status=pass' > "${prefix}.exonerate_pilot_gate.ok"
    echo '{"status":"stub"}' > "${prefix}.exonerate_pilot.summary.json"
    echo '"${task.process}": {python: "3.11"}' > versions.yml
    """
}
