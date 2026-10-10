process FILTER_OR_TREE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(tree_file), path(candidate_proteins)
    path candidate_nucleotides
    path candidate_coordinates
    path template_alignment
    path type1_ids
    path type2_ids
    path package_src

    output:
    tuple val(meta), path("${meta.id}.tree_filtered.prot.fa"),                emit: proteins
    tuple val(meta), path("${meta.id}.tree_filtered.nt.fa"),                  emit: nucleotides
    tuple val(meta), path("${meta.id}.tree_filtered.tsv"),                    emit: coordinates
    tuple val(meta), path("${meta.id}.tree_filter.decisions.tsv"),            emit: decisions
    tuple val(meta), path("${meta.id}.tree_filtered.ids.txt"),                emit: ids
    tuple val(meta), path("${meta.id}.nc_legacy.rooted.treefile"),           emit: rooted_tree
    tuple val(meta), path("${meta.id}.tree_filter.summary.json"),             emit: summary
    path 'versions.yml',                                                       emit: versions

    script:
    def prefix = meta.id
    """
    PYTHONPATH="\$PWD" python -m or_mining.filter_tree \
        "$tree_file" \
        "$template_alignment" \
        "$type1_ids" \
        "$type2_ids" \
        "$candidate_proteins" \
        "$candidate_nucleotides" \
        "$candidate_coordinates" \
        --output-prefix "$prefix"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biopython: \$(python -c 'import Bio; print(Bio.__version__)')
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    touch "${prefix}.tree_filtered.prot.fa"
    touch "${prefix}.tree_filtered.nt.fa"
    touch "${prefix}.tree_filtered.tsv"
    touch "${prefix}.tree_filter.decisions.tsv"
    touch "${prefix}.tree_filtered.ids.txt"
    touch "${prefix}.nc_legacy.rooted.treefile"
    echo '{"status":"stub"}' > "${prefix}.tree_filter.summary.json"
    echo '"${task.process}": {biopython: "1.85", python: "3.11"}' > versions.yml
    """
}
