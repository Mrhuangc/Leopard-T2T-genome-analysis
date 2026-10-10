process EXTRACT_MULTIEXON_SEQUENCES {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/../extract_candidate_sequences/environment.yml"

    input:
    tuple val(meta), path(genome), path(regions_tsv), path(faidx_regions)
    path package_src

    output:
    tuple val(meta), path("${meta.id}.multiexon_regions.fa"), emit: sequences
    tuple val(meta), path("${meta.id}.multiexon_extraction.summary.json"), emit: summary
    path 'versions.yml', emit: versions

    script:
    def prefix = meta.id
    """
    samtools faidx "$genome"
    samtools faidx -r "$faidx_regions" "$genome" > raw_multiexon_regions.fa

    PYTHONPATH="\$PWD" python -m or_mining.extract_candidates \
        raw_multiexon_regions.fa \
        "$regions_tsv" \
        --output "${prefix}.multiexon_regions.fa" \
        --summary "${prefix}.multiexon_extraction.summary.json"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | sed '1!d; s/samtools //')
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    touch "${prefix}.multiexon_regions.fa"
    echo '{"status":"stub"}' > "${prefix}.multiexon_extraction.summary.json"
    echo '"${task.process}": {samtools: "1.21", python: "3.11"}' > versions.yml
    """
}
