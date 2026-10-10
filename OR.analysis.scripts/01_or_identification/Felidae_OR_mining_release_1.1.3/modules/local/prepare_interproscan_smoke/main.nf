process PREPARE_INTERPROSCAN_SMOKE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(reference_or)
    val record_count
    path package_src

    output:
    tuple val(meta), path("${meta.id}.interproscan_smoke.prot.fa"), emit: fasta
    tuple val(meta), path("${meta.id}.interproscan_smoke.input.summary.json"), emit: summary

    script:
    def prefix = meta.id
    """
    PYTHONPATH="\$PWD" python -m or_mining.run_interproscan_smoke prepare \
        "$reference_or" \
        --records "$record_count" \
        --output "${prefix}.interproscan_smoke.prot.fa" \
        --summary "${prefix}.interproscan_smoke.input.summary.json"
    """

    stub:
    def prefix = meta.id
    """
    echo -e '>known_OR_smoke\nMTEST' > "${prefix}.interproscan_smoke.prot.fa"
    echo '{"status":"stub"}' > "${prefix}.interproscan_smoke.input.summary.json"
    """
}
