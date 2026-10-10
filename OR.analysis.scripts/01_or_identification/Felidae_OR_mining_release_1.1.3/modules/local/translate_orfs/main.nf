process TRANSLATE_ORFS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(nucleotide_orfs)
    path package_src

    output:
    tuple val(meta), path("${meta.id}.orfs.prot.fa"),                emit: proteins
    tuple val(meta), path("${meta.id}.translation.summary.json"),   emit: summary
    path 'versions.yml',                                            emit: versions

    script:
    def prefix = meta.id
    """
    transeq \
        -sequence "$nucleotide_orfs" \
        -outseq raw_transeq.prot.fa \
        -frame 1 \
        -table 0 \
        -auto

    PYTHONPATH="\$PWD" python -m or_mining.normalize_translation \
        "$nucleotide_orfs" \
        raw_transeq.prot.fa \
        --output "${prefix}.orfs.prot.fa" \
        --summary "${prefix}.translation.summary.json"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        emboss: \$(transeq -version 2>&1 | sed 's/^EMBOSS://')
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    touch "${prefix}.orfs.prot.fa"
    echo '{"status":"stub"}' > "${prefix}.translation.summary.json"
    echo '"${task.process}": {emboss: "6.6.0.0", python: "3.11"}' > versions.yml
    """
}

