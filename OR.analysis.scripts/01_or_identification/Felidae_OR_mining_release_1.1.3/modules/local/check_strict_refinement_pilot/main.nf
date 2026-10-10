process CHECK_STRICT_REFINEMENT_PILOT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(results), path(statuses)
    path package_src
    val expected_pairs

    output:
    path "${meta.id}.strict_refinement_pilot_gate.ok", emit: gate
    tuple val(meta), path("${meta.id}.strict_refinement_pilot.summary.json"), emit: summary
    path 'versions.yml', emit: versions

    script:
    def result_args = results.collect { result -> "\"${result}\"" }.join(' ')
    def status_args = statuses.collect { status -> "\"${status}\"" }.join(' ')
    """
    PYTHONPATH="\$PWD" python -m or_mining.check_strict_refinement_pilot \
        --results $result_args \
        --statuses $status_args \
        --output "${meta.id}.strict_refinement_pilot.summary.json" \
        --gate "${meta.id}.strict_refinement_pilot_gate.ok" \
        --expected-pairs "$expected_pairs"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    echo 'status=pass' > "${meta.id}.strict_refinement_pilot_gate.ok"
    echo '{"status":"stub"}' > "${meta.id}.strict_refinement_pilot.summary.json"
    echo '"${task.process}": {python: "3.11"}' > versions.yml
    """
}
