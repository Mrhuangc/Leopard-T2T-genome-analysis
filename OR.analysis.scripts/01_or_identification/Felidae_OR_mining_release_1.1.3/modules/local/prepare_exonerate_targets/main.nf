process PREPARE_EXONERATE_TARGETS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(target_fasta)
    path package_src
    val requested_chunks

    output:
    tuple val(meta), path('target_chunks/target_*.nucl.fa'), emit: chunks
    tuple val(meta), path("${meta.id}.exonerate_targets.summary.json"), emit: summary
    path 'versions.yml', emit: versions

    script:
    def prefix = meta.id
    """
    PYTHONPATH="\$PWD" python -m or_mining.prepare_exonerate_targets \
        "$target_fasta" \
        --output-dir target_chunks \
        --summary "${prefix}.exonerate_targets.summary.json" \
        --chunks "$requested_chunks"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    mkdir -p target_chunks
    touch target_chunks/target_001.nucl.fa
    echo '{"status":"stub"}' > "${prefix}.exonerate_targets.summary.json"
    echo '"${task.process}": {python: "3.11"}' > versions.yml
    """
}
