process CHECK_EXONERATE_SMOKE_GATE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/../prepare_exonerate_functional_smoke/environment.yml"

    input:
    tuple val(meta), path(smoke_summary)
    path package_src

    output:
    path "${meta.id}.exonerate_smoke_gate.ok", emit: gate

    script:
    """
    PYTHONPATH="\$PWD" python -m or_mining.run_exonerate_smoke gate \
        "$smoke_summary" \
        --output "${meta.id}.exonerate_smoke_gate.ok"
    """

    stub:
    """
    echo 'status=stub' > "${meta.id}.exonerate_smoke_gate.ok"
    """
}
