process EXTRACT_CANDIDATE_SEQUENCES {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(genome), path(regions_tsv), path(faidx_regions)
    path package_src

    output:
    tuple val(meta), path("${meta.id}.candidate_regions.fa"),           emit: sequences
    tuple val(meta), path("${meta.id}.candidate_extraction.summary.json"), emit: summary
    path 'versions.yml',                                                emit: versions

    script:
    def prefix = meta.id
    """
    samtools faidx "$genome"
    samtools faidx -r "$faidx_regions" "$genome" > raw_candidate_regions.fa

    PYTHONPATH="\$PWD" python -m or_mining.extract_candidates \
        raw_candidate_regions.fa \
        "$regions_tsv" \
        --output "${prefix}.candidate_regions.fa" \
        --summary "${prefix}.candidate_extraction.summary.json"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | sed '1!d; s/samtools //')
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    touch "${prefix}.candidate_regions.fa"
    echo '{"status":"stub"}' > "${prefix}.candidate_extraction.summary.json"
    echo '"${task.process}": {samtools: "1.21", python: "3.11"}' > versions.yml
    """
}

