process BUILD_PRE7TM_CATALOG {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(genome_fai), path(single_nucl), path(single_prot), path(single_tsv), path(multiexon_calls), path(multiexon_gff3), path(multiexon_functional_nucl), path(multiexon_functional_prot), path(multiexon_pseudogene_nucl), path(multiexon_truncated_nucl), path(multiexon_ambiguous_nucl)
    path package_src

    output:
    tuple val(meta), path("${meta.id}.pre7tm_OR_catalog.tsv"), emit: catalog
    tuple val(meta), path("${meta.id}.pre7tm_OR_catalog.gff3"), emit: gff3
    tuple val(meta), path("${meta.id}.pre7tm_OR_id_map.tsv"), emit: id_map
    tuple val(meta), path("${meta.id}.pre7tm_overlap_audit.tsv"), emit: overlap_audit
    tuple val(meta), path("${meta.id}.pre7tm_functional_candidates.nucl.fa"), emit: functional_nucl
    tuple val(meta), path("${meta.id}.pre7tm_functional_candidates.prot.fa"), emit: functional_prot
    tuple val(meta), path("${meta.id}.pre7tm_pseudogene.nucl.fa"), emit: pseudogene_nucl
    tuple val(meta), path("${meta.id}.pre7tm_truncated.nucl.fa"), emit: truncated_nucl
    tuple val(meta), path("${meta.id}.pre7tm_ambiguous.nucl.fa"), emit: ambiguous_nucl
    tuple val(meta), path("${meta.id}.pre7tm_OR_catalog.summary.json"), emit: summary
    path 'versions.yml', emit: versions

    script:
    def prefix = meta.id
    """
    PYTHONPATH="\$PWD" python -m or_mining.run_preliminary_catalog \
        "$genome_fai" \
        "$single_nucl" \
        "$single_prot" \
        "$single_tsv" \
        "$multiexon_calls" \
        "$multiexon_gff3" \
        "$multiexon_functional_nucl" \
        "$multiexon_functional_prot" \
        "$multiexon_pseudogene_nucl" \
        "$multiexon_truncated_nucl" \
        "$multiexon_ambiguous_nucl" \
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
    echo -e 'catalog_id\tclassification' > "${prefix}.pre7tm_OR_catalog.tsv"
    echo '##gff-version 3' > "${prefix}.pre7tm_OR_catalog.gff3"
    echo -e 'catalog_id\tsource_id\tsource_stage' > "${prefix}.pre7tm_OR_id_map.tsv"
    echo -e 'left_id\tright_id\tsseqid\taction\tselected_id\tremoved_ids' > "${prefix}.pre7tm_overlap_audit.tsv"
    touch "${prefix}.pre7tm_functional_candidates.nucl.fa"
    touch "${prefix}.pre7tm_functional_candidates.prot.fa"
    touch "${prefix}.pre7tm_pseudogene.nucl.fa"
    touch "${prefix}.pre7tm_truncated.nucl.fa"
    touch "${prefix}.pre7tm_ambiguous.nucl.fa"
    echo '{"status":"stub","reportable_final_catalog":false}' > "${prefix}.pre7tm_OR_catalog.summary.json"
    echo '"${task.process}": {python: "3.11"}' > versions.yml
    """
}
