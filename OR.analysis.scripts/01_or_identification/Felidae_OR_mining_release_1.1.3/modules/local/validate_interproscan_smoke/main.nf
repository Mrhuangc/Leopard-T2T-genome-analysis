process VALIDATE_INTERPROSCAN_SMOKE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/../prepare_interproscan_smoke/environment.yml"

    input:
    tuple val(meta), path(proteins), path(result_tsv), path(version_report)
    val expected_version
    path package_src

    output:
    tuple val(meta), path("${meta.id}.interproscan_functional_smoke.summary.json"), emit: summary
    tuple val(meta), path("${meta.id}.interproscan_functional_smoke.ok"), emit: gate

    script:
    def prefix = meta.id
    """
    PYTHONPATH="\$PWD" python -m or_mining.run_interproscan_smoke validate \
        "$proteins" \
        "$result_tsv" \
        "$version_report" \
        --expected-version "$expected_version" \
        --output "${prefix}.interproscan_functional_smoke.summary.json" \
        --gate "${prefix}.interproscan_functional_smoke.ok"
    """

    stub:
    def prefix = meta.id
    """
    echo '{"status":"stub"}' > "${prefix}.interproscan_functional_smoke.summary.json"
    echo 'status=stub' > "${prefix}.interproscan_functional_smoke.ok"
    """
}
