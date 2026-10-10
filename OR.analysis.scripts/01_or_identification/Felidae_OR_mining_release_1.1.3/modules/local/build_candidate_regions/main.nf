process BUILD_CANDIDATE_REGIONS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(input_merged)
    path package_src
    val min_width
    val flank

    output:
    tuple val(meta), path("${meta.id}.candidate_regions.tsv"),          emit: regions
    tuple val(meta), path("${meta.id}.candidate_regions.bed"),          emit: bed
    tuple val(meta), path("${meta.id}.candidate_regions.faidx.txt"),    emit: faidx
    tuple val(meta), path("${meta.id}.candidate_regions.summary.json"), emit: summary
    path 'versions.yml',                                                 emit: versions

    script:
    def prefix = meta.id
    """
    PYTHONPATH="\$PWD" python -m or_mining.build_regions \
        "$input_merged" \
        --species "$meta.id" \
        --output-prefix "$prefix" \
        --min-width "$min_width" \
        --flank "$flank"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    touch "${prefix}.candidate_regions.tsv"
    touch "${prefix}.candidate_regions.bed"
    touch "${prefix}.candidate_regions.faidx.txt"
    echo '{"status":"stub"}' > "${prefix}.candidate_regions.summary.json"
    echo '"${task.process}": {python: "3.11"}' > versions.yml
    """
}

