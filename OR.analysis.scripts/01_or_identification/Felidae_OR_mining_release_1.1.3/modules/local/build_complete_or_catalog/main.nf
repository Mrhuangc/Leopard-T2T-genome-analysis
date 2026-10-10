process BUILD_COMPLETE_OR_CATALOG {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(existing_catalog), path(existing_gff3), path(existing_functional_nucleotide), path(existing_functional_protein), path(existing_pseudogene), path(existing_truncated), path(existing_ambiguous)
    tuple val(meta2), path(recovery_calls), path(recovered_pseudogene), path(recovered_truncated), path(recovered_ambiguous), path(recovery_blastx)
    path package_src

    output:
    tuple val(meta), path("${meta.id}.complete_OR_catalog.tsv"), path("${meta.id}.complete_OR_catalog.gff3"), emit: catalog
    tuple val(meta), path("${meta.id}.complete_functional_OR.nucl.fa"), path("${meta.id}.complete_functional_OR.prot.fa"), emit: functional
    tuple val(meta), path("${meta.id}.complete_pseudogene_OR.nucl.fa"), emit: pseudogene
    tuple val(meta), path("${meta.id}.complete_truncated_OR.nucl.fa"), emit: truncated
    tuple val(meta), path("${meta.id}.complete_ambiguous_OR.nucl.fa"), emit: ambiguous
    tuple val(meta), path("${meta.id}.single_exon_recovery_audit.tsv"), path("${meta.id}.single_exon_nonOR_rejected.nucl.fa"), emit: audit
    tuple val(meta), path("${meta.id}.complete_OR_catalog.summary.json"), emit: summary
    path 'versions.yml', emit: versions

    script:
    def prefix = meta.id
    if (meta.id != meta2.id) {
        error "existing and recovered catalog species IDs differ: ${meta.id}/${meta2.id}"
    }
    """
    PYTHONPATH="\$PWD" python -m or_mining.run_single_exon_recovery complete \
        "$existing_catalog" \
        "$existing_gff3" \
        "$existing_functional_nucleotide" \
        "$existing_functional_protein" \
        "$existing_pseudogene" \
        "$existing_truncated" \
        "$existing_ambiguous" \
        "$recovery_calls" \
        "$recovered_pseudogene" \
        "$recovered_truncated" \
        "$recovered_ambiguous" \
        "$recovery_blastx" \
        --species "$prefix" \
        --output-prefix "$prefix"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    echo -e 'catalog_id\tclassification' > "${prefix}.complete_OR_catalog.tsv"
    echo '##gff-version 3' > "${prefix}.complete_OR_catalog.gff3"
    touch "${prefix}.complete_functional_OR.nucl.fa" "${prefix}.complete_functional_OR.prot.fa"
    touch "${prefix}.complete_pseudogene_OR.nucl.fa" "${prefix}.complete_truncated_OR.nucl.fa"
    touch "${prefix}.complete_ambiguous_OR.nucl.fa"
    echo -e 'candidate_id\tdecision' > "${prefix}.single_exon_recovery_audit.tsv"
    touch "${prefix}.single_exon_nonOR_rejected.nucl.fa"
    echo '{"status":"stub","reportable_final_catalog":false}' > "${prefix}.complete_OR_catalog.summary.json"
    echo '"${task.process}": {python: "3.11"}' > versions.yml
    """
}
