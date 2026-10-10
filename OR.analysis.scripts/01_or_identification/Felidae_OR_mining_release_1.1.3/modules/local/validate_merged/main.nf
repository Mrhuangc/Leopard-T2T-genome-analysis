process VALIDATE_MERGED {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(input_merged)
    path package_src

    output:
    tuple val(meta), path("${meta.id}.merged.validation.json"), emit: report
    path 'versions.yml',                                      emit: versions

    script:
    def prefix = meta.id
    """
    PYTHONPATH="\$PWD" python -m or_mining.validate_merged \
        "$input_merged" \
        --report "${prefix}.merged.validation.json"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    echo '{"status":"stub"}' > "${prefix}.merged.validation.json"
    echo '"${task.process}": {python: "3.11"}' > versions.yml
    """
}

