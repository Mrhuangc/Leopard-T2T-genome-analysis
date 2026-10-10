process PREPARE_EXONERATE_FUNCTIONAL_SMOKE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    val meta
    path package_src

    output:
    tuple val(meta), path('synthetic_spliced_or.prot.fa'), path('synthetic_target.nucl.fa'), emit: inputs
    path "${meta.id}.exonerate_smoke_inputs.summary.json", emit: summary

    script:
    """
    PYTHONPATH="\$PWD" python -m or_mining.run_exonerate_smoke generate \
        --query synthetic_spliced_or.prot.fa \
        --target synthetic_target.nucl.fa \
        --summary "${meta.id}.exonerate_smoke_inputs.summary.json"
    """

    stub:
    """
    touch synthetic_spliced_or.prot.fa synthetic_target.nucl.fa
    echo '{"status":"stub"}' > "${meta.id}.exonerate_smoke_inputs.summary.json"
    """
}
