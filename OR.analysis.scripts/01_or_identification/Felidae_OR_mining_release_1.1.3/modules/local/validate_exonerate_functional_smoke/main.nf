process VALIDATE_EXONERATE_FUNCTIONAL_SMOKE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/../prepare_exonerate_functional_smoke/environment.yml"

    input:
    tuple val(meta), path(result), path(chunk_summary)
    path package_src

    output:
    tuple val(meta), path("${meta.id}.exonerate_functional_smoke.summary.json"), emit: summary

    script:
    """
    PYTHONPATH="\$PWD" python -m or_mining.run_exonerate_smoke validate \
        "$result" \
        "$chunk_summary" \
        --output "${meta.id}.exonerate_functional_smoke.summary.json"
    """

    stub:
    """
    echo '{"status":"stub"}' > "${meta.id}.exonerate_functional_smoke.summary.json"
    """
}
