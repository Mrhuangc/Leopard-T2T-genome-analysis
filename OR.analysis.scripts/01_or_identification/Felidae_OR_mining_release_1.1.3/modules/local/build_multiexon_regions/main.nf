process BUILD_MULTIEXON_REGIONS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/../classify_gpcr/environment.yml"

    input:
    tuple val(meta), path(merged_hits)
    path known_coordinates
    path package_src
    val maximum_intron_length
    val gene_clearance

    output:
    tuple val(meta), path("${meta.id}.multiexon_regions.tsv"),       emit: regions
    tuple val(meta), path("${meta.id}.multiexon_regions.bed"),       emit: bed
    tuple val(meta), path("${meta.id}.multiexon_regions.faidx.txt"), emit: faidx
    tuple val(meta), path("${meta.id}.multiexon_regions.summary.json"), emit: summary
    path 'versions.yml', emit: versions

    script:
    def prefix = meta.id
    """
    PYTHONPATH="\$PWD" python -m or_mining.build_multiexon_regions \
        "$merged_hits" \
        "$known_coordinates" \
        --output-prefix "$prefix" \
        --species "$prefix" \
        --maximum-intron-length "$maximum_intron_length" \
        --gene-clearance "$gene_clearance"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    touch "${prefix}.multiexon_regions.tsv"
    touch "${prefix}.multiexon_regions.bed"
    touch "${prefix}.multiexon_regions.faidx.txt"
    echo '{"status":"stub"}' > "${prefix}.multiexon_regions.summary.json"
    echo '"${task.process}": {python: "3.11"}' > versions.yml
    """
}
