process CLASSIFY_GPCR {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(blastp_tsv), path(protein_orfs)
    path nucleotide_orfs
    path orf_table
    path package_src

    output:
    tuple val(meta), path("${meta.id}.gpcr_classification.tsv"),          emit: classification
    tuple val(meta), path("${meta.id}.OR_candidates.prot.fa"),           emit: proteins
    tuple val(meta), path("${meta.id}.OR_candidates.nt.fa"),             emit: nucleotides
    tuple val(meta), path("${meta.id}.OR_candidates.tsv"),               emit: coordinates
    tuple val(meta), path("${meta.id}.OR_candidates.ids.txt"),           emit: ids
    tuple val(meta), path("${meta.id}.gpcr_classification.summary.json"), emit: summary
    path 'versions.yml',                                                   emit: versions

    script:
    def prefix = meta.id
    """
    PYTHONPATH="\$PWD" python -m or_mining.classify_gpcr \
        "$blastp_tsv" \
        "$protein_orfs" \
        "$nucleotide_orfs" \
        "$orf_table" \
        --output-prefix "$prefix"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    touch "${prefix}.gpcr_classification.tsv"
    touch "${prefix}.OR_candidates.prot.fa"
    touch "${prefix}.OR_candidates.nt.fa"
    touch "${prefix}.OR_candidates.tsv"
    touch "${prefix}.OR_candidates.ids.txt"
    echo '{"status":"stub"}' > "${prefix}.gpcr_classification.summary.json"
    echo '"${task.process}": {python: "3.11"}' > versions.yml
    """
}

