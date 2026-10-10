process CLASSIFY_SINGLE_EXON_RECOVERY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(candidates_tsv), path(targets), path(standard_tsv), path(gapless_tsv)
    path package_src

    output:
    tuple val(meta), path("${meta.id}.single_exon_recovery_calls.tsv"), path("${meta.id}.single_exon_recovery_candidates.nucl.fa"), path("${meta.id}.single_exon_pseudogene.nucl.fa"), path("${meta.id}.single_exon_truncated.nucl.fa"), path("${meta.id}.single_exon_ambiguous.nucl.fa"), emit: classified
    tuple val(meta), path("${meta.id}.single_exon_no_targeted_hit.nucl.fa"), emit: rejected
    tuple val(meta), path("${meta.id}.single_exon_recovery_classification.summary.json"), emit: summary
    path 'versions.yml', emit: versions

    script:
    def prefix = meta.id
    """
    PYTHONPATH="\$PWD" python -m or_mining.run_single_exon_recovery classify \
        "$candidates_tsv" \
        "$targets" \
        "$standard_tsv" \
        "$gapless_tsv" \
        --output-prefix "$prefix"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    echo -e 'candidate_id\tclassification\tevidence' > "${prefix}.single_exon_recovery_calls.tsv"
    touch "${prefix}.single_exon_recovery_candidates.nucl.fa"
    touch "${prefix}.single_exon_pseudogene.nucl.fa" "${prefix}.single_exon_truncated.nucl.fa"
    touch "${prefix}.single_exon_ambiguous.nucl.fa" "${prefix}.single_exon_no_targeted_hit.nucl.fa"
    echo '{"status":"stub"}' > "${prefix}.single_exon_recovery_classification.summary.json"
    echo '"${task.process}": {python: "3.11"}' > versions.yml
    """
}
