process PREPARE_TARGETED_EXONERATE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(reference_or), path(species_or), path(target_fasta), path(regions_tsv), path(merged_hits)
    path package_src
    val top_queries_per_region

    output:
    tuple val(meta), path('pairs/*.query.prot.fa'), path('pairs/*.target.nucl.fa'), emit: pairs
    tuple val(meta), path("${meta.id}.targeted_exonerate_pairs.tsv"), emit: manifest
    tuple val(meta), path("${meta.id}.targeted_exonerate_pairs.summary.json"), emit: summary
    path 'versions.yml', emit: versions

    script:
    def prefix = meta.id
    """
    PYTHONPATH="\$PWD" python -m or_mining.prepare_targeted_exonerate \
        "$reference_or" \
        "$species_or" \
        "$target_fasta" \
        "$regions_tsv" \
        "$merged_hits" \
        --output-dir pairs \
        --manifest "${prefix}.targeted_exonerate_pairs.tsv" \
        --summary "${prefix}.targeted_exonerate_pairs.summary.json" \
        --top-queries-per-region "$top_queries_per_region"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    mkdir -p pairs
    echo -e '>q\nM' > pairs/pair_000001.query.prot.fa
    echo -e '>t\nATG' > pairs/pair_000001.target.nucl.fa
    echo -e 'pair_id\tregion_id\n' > "${prefix}.targeted_exonerate_pairs.tsv"
    echo '{"status":"stub"}' > "${prefix}.targeted_exonerate_pairs.summary.json"
    echo '"${task.process}": {python: "3.11"}' > versions.yml
    """
}
