process GETORF {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(candidate_fasta)
    val min_orf_nt

    output:
    tuple val(meta), path("${meta.id}.getorf.raw.fa"), emit: orfs
    path 'versions.yml',                              emit: versions

    script:
    def prefix = meta.id
    """
    getorf \
        -sequence "$candidate_fasta" \
        -outseq "${prefix}.getorf.raw.fa" \
        -minsize "$min_orf_nt" \
        -find 3 \
        -auto

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        emboss: \$(getorf -version 2>&1 | sed 's/^EMBOSS://')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    touch "${prefix}.getorf.raw.fa"
    echo '"${task.process}": {emboss: "6.6.0.0"}' > versions.yml
    """
}

