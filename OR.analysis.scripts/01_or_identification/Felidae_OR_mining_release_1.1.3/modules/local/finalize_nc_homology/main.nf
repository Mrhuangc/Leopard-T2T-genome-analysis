process FINALIZE_NC_HOMOLOGY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(catalog), path(gff3), path(functional_nucleotide), path(functional_protein), path(pseudogene_nucleotide), path(truncated_nucleotide), path(ambiguous_nucleotide)
    tuple val(meta2), path(functional_blastp), path(pseudogene_blastx), path(truncated_blastx), path(ambiguous_blastx)
    val blast_evalue
    path package_src

    output:
    tuple val(meta), path("${meta.id}.NC_final_OR_catalog.tsv"), path("${meta.id}.NC_final_OR_catalog.gff3"), emit: catalog
    tuple val(meta), path("${meta.id}.NC_final_OR_homology_audit.tsv"), emit: audit
    tuple val(meta), path("${meta.id}.NC_excluded_nonOR_catalog.tsv"), path("${meta.id}.NC_excluded_nonOR.nucl.fa"), path("${meta.id}.NC_excluded_nonOR.prot.fa"), emit: rejected
    tuple val(meta), path("${meta.id}.NC_final_functional_OR.nucl.fa"), path("${meta.id}.NC_final_functional_OR.prot.fa"), emit: functional
    tuple val(meta), path("${meta.id}.NC_final_pseudogene_OR.nucl.fa"), emit: pseudogene
    tuple val(meta), path("${meta.id}.NC_final_truncated_OR.nucl.fa"), emit: truncated
    tuple val(meta), path("${meta.id}.NC_final_ambiguous_OR.nucl.fa"), emit: ambiguous
    tuple val(meta), path("${meta.id}.NC_final_OR_catalog.summary.json"), emit: summary
    path 'versions.yml', emit: versions

    script:
    def prefix = meta.id
    if (meta.id != meta2.id) {
        error "catalog and BLAST species IDs differ: ${meta.id}/${meta2.id}"
    }
    """
    PYTHONPATH="\$PWD" python -m or_mining.run_final_homology \
        "$catalog" \
        "$gff3" \
        "$functional_nucleotide" \
        "$functional_protein" \
        "$pseudogene_nucleotide" \
        "$truncated_nucleotide" \
        "$ambiguous_nucleotide" \
        "$functional_blastp" \
        "$pseudogene_blastx" \
        "$truncated_blastx" \
        "$ambiguous_blastx" \
        --species "${meta.id}" \
        --blast-evalue "$blast_evalue" \
        --output-prefix "$prefix"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    echo -e 'catalog_id\tclassification\tfinal_homology_decision' > "${prefix}.NC_final_OR_catalog.tsv"
    echo '##gff-version 3' > "${prefix}.NC_final_OR_catalog.gff3"
    echo -e 'catalog_id\tdecision\treason' > "${prefix}.NC_final_OR_homology_audit.tsv"
    echo -e 'catalog_id\tclassification\tfinal_homology_decision' > "${prefix}.NC_excluded_nonOR_catalog.tsv"
    touch "${prefix}.NC_excluded_nonOR.nucl.fa" "${prefix}.NC_excluded_nonOR.prot.fa"
    touch "${prefix}.NC_final_functional_OR.nucl.fa" "${prefix}.NC_final_functional_OR.prot.fa"
    touch "${prefix}.NC_final_pseudogene_OR.nucl.fa" "${prefix}.NC_final_truncated_OR.nucl.fa"
    touch "${prefix}.NC_final_ambiguous_OR.nucl.fa"
    echo '{"status":"stub","reportable_final_catalog":false}' > "${prefix}.NC_final_OR_catalog.summary.json"
    echo '"${task.process}": {python: "stub"}' > versions.yml
    """
}
