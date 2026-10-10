process EXTRACT_SINGLE_EXON_TARGETS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(candidates_tsv), path(queries), path(faidx_regions)
    path genome
    path genome_fai
    path package_src

    output:
    tuple val(meta), path(candidates_tsv), path(queries), path("${meta.id}.single_exon_residual_targets.nucl.fa"), emit: targets
    tuple val(meta), path("${meta.id}.single_exon_residual_targets.summary.json"), emit: summary
    path 'versions.yml', emit: versions

    script:
    def prefix = meta.id
    """
    test -s "$genome_fai"
    if [[ -s "$faidx_regions" ]]; then
        samtools faidx -r "$faidx_regions" "$genome" > raw_targets.fa
    else
        : > raw_targets.fa
    fi

    PYTHONPATH="\$PWD" python -m or_mining.run_single_exon_recovery normalize-targets \
        raw_targets.fa \
        "$candidates_tsv" \
        --output "${prefix}.single_exon_residual_targets.nucl.fa" \
        --summary "${prefix}.single_exon_residual_targets.summary.json"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | sed '1!d; s/samtools //')
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    touch "${prefix}.single_exon_residual_targets.nucl.fa"
    echo '{"status":"stub"}' > "${prefix}.single_exon_residual_targets.summary.json"
    echo '"${task.process}": {samtools: "1.21", python: "3.11"}' > versions.yml
    """
}
