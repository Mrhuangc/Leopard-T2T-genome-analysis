process MAFFT_ADD {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(candidate_proteins)
    path template_alignment

    output:
    tuple val(meta), path("${meta.id}.OR_plus_reference.aln.fa"), path(candidate_proteins), emit: alignment
    tuple val(meta), path("${meta.id}.mafft.log"),                                      emit: log
    path 'versions.yml',                                                               emit: versions

    script:
    def prefix = meta.id
    """
    mafft \
        --thread "$task.cpus" \
        --add "$candidate_proteins" \
        --keeplength \
        "$template_alignment" \
        > "${prefix}.OR_plus_reference.aln.fa" \
        2> "${prefix}.mafft.log"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mafft: \$(mafft --version 2>&1 | sed 's/^v//')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    cp "$template_alignment" "${prefix}.OR_plus_reference.aln.fa"
    touch "${prefix}.mafft.log"
    echo '"${task.process}": {mafft: "7.526"}' > versions.yml
    """
}
