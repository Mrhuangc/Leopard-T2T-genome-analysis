process FINALIZE_OR_CATALOG {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(pre7tm_catalog), path(pre7tm_gff3), path(functional_nucleotide), path(functional_protein), path(pseudogene_nucleotide), path(truncated_nucleotide), path(ambiguous_nucleotide)
    tuple val(meta2), path(interpro_tsv), path(version_report)
    val expected_version
    val expected_ipr
    path package_src

    output:
    tuple val(meta), path("${meta.id}.final_OR_catalog.tsv"), emit: catalog
    tuple val(meta), path("${meta.id}.final_OR_catalog.gff3"), emit: gff3
    tuple val(meta), path("${meta.id}.interpro_7tm_audit.tsv"), emit: audit
    tuple val(meta), path("${meta.id}.final_functional_OR.nucl.fa"), path("${meta.id}.final_functional_OR.prot.fa"), emit: functional
    tuple val(meta), path("${meta.id}.final_pseudogene_OR.nucl.fa"), emit: pseudogene
    tuple val(meta), path("${meta.id}.final_truncated_OR.nucl.fa"), emit: truncated
    tuple val(meta), path("${meta.id}.final_ambiguous_OR.nucl.fa"), emit: ambiguous
    tuple val(meta), path("${meta.id}.non7tm_rejected.nucl.fa"), path("${meta.id}.non7tm_rejected.prot.fa"), emit: rejected
    tuple val(meta), path("${meta.id}.final_OR_catalog.summary.json"), emit: summary
    path 'versions.yml', emit: versions

    script:
    def prefix = meta.id
    if (meta.id != meta2.id) {
        error "pre-7TM and InterProScan species IDs differ: ${meta.id}/${meta2.id}"
    }
    """
    PYTHONPATH="\$PWD" python -m or_mining.run_final_catalog \
        "$pre7tm_catalog" \
        "$pre7tm_gff3" \
        "$functional_nucleotide" \
        "$functional_protein" \
        "$pseudogene_nucleotide" \
        "$truncated_nucleotide" \
        "$ambiguous_nucleotide" \
        "$interpro_tsv" \
        "$version_report" \
        --species "${meta.id}" \
        --expected-version "$expected_version" \
        --expected-ipr "$expected_ipr" \
        --output-prefix "$prefix"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
        expected_interpro_accession: $expected_ipr
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    echo -e 'catalog_id\\tclassification\\tinterpro_7tm_status' > "${prefix}.final_OR_catalog.tsv"
    echo '##gff-version 3' > "${prefix}.final_OR_catalog.gff3"
    echo -e 'catalog_id\\tstatus' > "${prefix}.interpro_7tm_audit.tsv"
    touch "${prefix}.final_functional_OR.nucl.fa" "${prefix}.final_functional_OR.prot.fa"
    touch "${prefix}.final_pseudogene_OR.nucl.fa" "${prefix}.final_truncated_OR.nucl.fa"
    touch "${prefix}.final_ambiguous_OR.nucl.fa"
    touch "${prefix}.non7tm_rejected.nucl.fa" "${prefix}.non7tm_rejected.prot.fa"
    echo '{"status":"stub","reportable_final_catalog":false}' > "${prefix}.final_OR_catalog.summary.json"
    echo '"${task.process}": {python: "stub"}' > versions.yml
    """
}
