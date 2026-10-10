process CLASSIFY_OR_CLADES {
    tag "$meta.id"
    label 'process_single'

    // Reuse the validated Biopython tree environment from stage 03.
    conda "${moduleDir}/../filter_or_tree/environment.yml"

    input:
    tuple val(meta), path(tree_file)
    tuple val(meta2), path(functional_blastp), path(pseudogene_blastx), path(truncated_blastx), path(ambiguous_blastx)
    path template_alignment
    path type1_ids
    path type2_ids
    path catalog
    path functional_protein
    path package_src
    val minimum_tree_support

    output:
    tuple val(meta), path("${meta.id}.OR_clade_assignments.tsv"), path("${meta.id}.OR_catalog_with_clades.tsv"), path("${meta.id}.OR_clade_counts.tsv"), emit: tables
    tuple val(meta), path("${meta.id}.OR_clades.audited.unrooted.treefile"), emit: tree
    tuple val(meta), path("${meta.id}.OR_clade_classification.summary.json"), emit: summary
    path 'versions.yml', emit: versions

    script:
    def prefix = meta.id
    if (meta.id != meta2.id) {
        error "tree and BLAST species IDs differ: ${meta.id}/${meta2.id}"
    }
    """
    PYTHONPATH="\$PWD" python -m or_mining.run_or_clade_classification \
        "$tree_file" \
        "$template_alignment" \
        "$type1_ids" \
        "$type2_ids" \
        "$catalog" \
        "$functional_protein" \
        "$functional_blastp" \
        "$pseudogene_blastx" \
        "$truncated_blastx" \
        "$ambiguous_blastx" \
        --species "${meta.id}" \
        --output-prefix "$prefix" \
        --minimum-tree-support "$minimum_tree_support"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biopython: \$(python -c 'import Bio; print(Bio.__version__)')
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    touch "${prefix}.OR_clade_assignments.tsv" "${prefix}.OR_catalog_with_clades.tsv"
    touch "${prefix}.OR_clade_counts.tsv" "${prefix}.OR_clades.audited.unrooted.treefile"
    echo '{"status":"stub"}' > "${prefix}.OR_clade_classification.summary.json"
    echo '"${task.process}": {python: "3.11", biopython: "1.85"}' > versions.yml
    """
}
