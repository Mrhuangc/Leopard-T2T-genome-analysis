process PREPARE_STRICT_REFINEMENT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(reference_or), path(species_or), path(region_fasta), path(regions_tsv), path(candidates_tsv)
    path package_src

    output:
    tuple val(meta), path('pairs/*.query.prot.fa'), path('pairs/*.target.nucl.fa'), emit: pairs
    tuple val(meta), path("${meta.id}.strict_refinement_pairs.tsv"), emit: manifest
    tuple val(meta), path("${meta.id}.strict_refinement_pairs.summary.json"), emit: summary
    path 'versions.yml', emit: versions

    script:
    def prefix = meta.id
    """
    PYTHONPATH="\$PWD" python -m or_mining.prepare_strict_refinement \
        "$reference_or" \
        "$species_or" \
        "$region_fasta" \
        "$regions_tsv" \
        "$candidates_tsv" \
        --output-dir pairs \
        --manifest "${prefix}.strict_refinement_pairs.tsv" \
        --summary "${prefix}.strict_refinement_pairs.summary.json"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    mkdir -p pairs
    echo -e '>q\nM' > pairs/strict_000001__ORMEC000001.query.prot.fa
    echo -e '>ORMEC000001\nATG' > pairs/strict_000001__ORMEC000001.target.nucl.fa
    echo -e 'strict_pair_id\tcandidate_id\nstrict_000001__ORMEC000001\tORMEC000001' > "${prefix}.strict_refinement_pairs.tsv"
    echo '{"status":"stub"}' > "${prefix}.strict_refinement_pairs.summary.json"
    echo '"${task.process}": {python: "3.11"}' > versions.yml
    """
}
