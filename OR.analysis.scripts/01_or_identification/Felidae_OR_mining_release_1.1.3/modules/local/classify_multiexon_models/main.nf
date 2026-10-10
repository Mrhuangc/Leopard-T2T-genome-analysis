process CLASSIFY_MULTIEXON_MODELS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(strict_raw), path(strict_manifest), path(regions_tsv), path(regions_fasta), path(reference_or), path(species_or)
    path package_src
    val minimum_orf_nt
    val minimum_model_overlap
    val maximum_exons
    val edge_distance_bp
    val edge_n_window_bp
    val edge_n_run_bp

    output:
    tuple val(meta), path("${meta.id}.multiexon_calls.tsv"), emit: calls
    tuple val(meta), path("${meta.id}.multiexon_calls.gff3"), emit: gff3
    tuple val(meta), path("${meta.id}.multiexon_functional.nucl.fa"), emit: functional_nucl
    tuple val(meta), path("${meta.id}.multiexon_functional.prot.fa"), emit: functional_prot
    tuple val(meta), path("${meta.id}.multiexon_pseudogene.nucl.fa"), emit: pseudogene_nucl
    tuple val(meta), path("${meta.id}.multiexon_truncated.nucl.fa"), emit: truncated_nucl
    tuple val(meta), path("${meta.id}.multiexon_ambiguous.nucl.fa"), emit: ambiguous_nucl
    tuple val(meta), path("${meta.id}.multiexon_classification.summary.json"), emit: summary
    path 'versions.yml', emit: versions

    script:
    def prefix = meta.id
    """
    PYTHONPATH="\$PWD" python -m or_mining.run_multiexon_classification \
        "$strict_raw" \
        "$strict_manifest" \
        "$regions_tsv" \
        "$regions_fasta" \
        "$reference_or" \
        "$species_or" \
        --output-prefix "$prefix" \
        --minimum-orf-nt "$minimum_orf_nt" \
        --minimum-model-overlap "$minimum_model_overlap" \
        --maximum-exons "$maximum_exons" \
        --edge-distance-bp "$edge_distance_bp" \
        --edge-n-window-bp "$edge_n_window_bp" \
        --edge-n-run-bp "$edge_n_run_bp"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    echo -e 'candidate_id\tclassification' > "${prefix}.multiexon_calls.tsv"
    echo '##gff-version 3' > "${prefix}.multiexon_calls.gff3"
    touch "${prefix}.multiexon_functional.nucl.fa"
    touch "${prefix}.multiexon_functional.prot.fa"
    touch "${prefix}.multiexon_pseudogene.nucl.fa"
    touch "${prefix}.multiexon_truncated.nucl.fa"
    touch "${prefix}.multiexon_ambiguous.nucl.fa"
    echo '{"status":"stub"}' > "${prefix}.multiexon_classification.summary.json"
    echo '"${task.process}": {python: "3.11"}' > versions.yml
    """
}
