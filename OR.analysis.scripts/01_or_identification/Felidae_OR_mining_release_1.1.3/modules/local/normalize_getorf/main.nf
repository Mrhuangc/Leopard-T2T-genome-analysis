process NORMALIZE_GETORF {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(raw_orfs)
    path regions_tsv
    path package_src

    output:
    tuple val(meta), path("${meta.id}.orfs.nt.fa"),       emit: fasta
    tuple val(meta), path("${meta.id}.orfs.tsv"),         emit: table
    tuple val(meta), path("${meta.id}.orfs.summary.json"), emit: summary
    path 'versions.yml',                                  emit: versions

    script:
    def prefix = meta.id
    """
    PYTHONPATH="\$PWD" python -m or_mining.normalize_orfs \
        "$raw_orfs" \
        "$regions_tsv" \
        --species "$meta.id" \
        --output-fasta "${prefix}.orfs.nt.fa" \
        --output-tsv "${prefix}.orfs.tsv" \
        --summary "${prefix}.orfs.summary.json"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    touch "${prefix}.orfs.nt.fa"
    touch "${prefix}.orfs.tsv"
    echo '{"status":"stub"}' > "${prefix}.orfs.summary.json"
    echo '"${task.process}": {python: "3.11"}' > versions.yml
    """
}

