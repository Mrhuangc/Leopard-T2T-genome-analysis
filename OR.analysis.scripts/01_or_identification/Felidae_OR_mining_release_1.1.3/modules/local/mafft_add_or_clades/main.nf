process MAFFT_ADD_OR_CLADES {
    tag "$meta.id"
    label 'process_medium'

    // Reuse the already validated/cached MAFFT toolchain from stage 03.
    conda "${moduleDir}/../mafft_add/environment.yml"

    input:
    tuple val(meta), path(functional_protein)
    path template_alignment

    output:
    tuple val(meta), path("${meta.id}.final_OR_plus_clade_reference.aln.fa"), emit: alignment
    tuple val(meta), path("${meta.id}.OR_clades.mafft.log"), emit: log
    path 'versions.yml', emit: versions

    script:
    def prefix = meta.id
    """
    mafft \
        --thread "$task.cpus" \
        --add "$functional_protein" \
        --keeplength \
        "$template_alignment" \
        > "${prefix}.final_OR_plus_clade_reference.aln.fa" \
        2> "${prefix}.OR_clades.mafft.log"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mafft: \$(mafft --version 2>&1 | sed 's/^v//')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    cp "$template_alignment" "${prefix}.final_OR_plus_clade_reference.aln.fa"
    touch "${prefix}.OR_clades.mafft.log"
    echo '"${task.process}": {mafft: "7.526"}' > versions.yml
    """
}
