process EXONERATE_SMOKE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/../exonerate_protein2genome/environment.yml"

    input:
    val meta

    output:
    tuple val(meta), path("${meta.id}.exonerate_toolchain.txt"), emit: report
    path 'versions.yml', emit: versions

    script:
    def prefix = meta.id
    """
    command -v exonerate
    # Exonerate 2.4.0 prints a valid version banner but returns exit code 1.
    exonerate --version > "${prefix}.exonerate_toolchain.txt" 2>&1 || true
    grep -qi 'exonerate' "${prefix}.exonerate_toolchain.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        exonerate: \$(exonerate --version 2>&1 | sed -n '1s/.*version //p' || true)
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    echo 'exonerate version 2.4.0' > "${prefix}.exonerate_toolchain.txt"
    echo '"${task.process}": {exonerate: "2.4.0"}' > versions.yml
    """
}
