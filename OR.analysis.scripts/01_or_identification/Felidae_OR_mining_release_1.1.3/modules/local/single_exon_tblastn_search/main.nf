process SINGLE_EXON_TBLASTN_SEARCH {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(candidates_tsv), path(queries), path(targets)
    tuple val(meta2), path(pilot_gate)
    path package_src
    val evalue
    val pair_timeout_seconds

    output:
    tuple val(meta), path(candidates_tsv), path(targets), path("${meta.id}.single_exon_tblastn.standard.tsv"), path("${meta.id}.single_exon_tblastn.gapless.tsv"), emit: results
    tuple val(meta), path("${meta.id}.single_exon_tblastn.summary.json"), emit: summary
    path 'versions.yml', emit: versions

    script:
    def prefix = meta.id
    if (meta.id != meta2.id) {
        error "target and pilot species IDs differ: ${meta.id}/${meta2.id}"
    }
    """
    grep -Fxq 'status=pass' "$pilot_gate"

    PYTHONPATH="\$PWD" python -m or_mining.run_single_exon_recovery search \
        "$candidates_tsv" \
        "$queries" \
        "$targets" \
        --standard-output "${prefix}.single_exon_tblastn.standard.tsv" \
        --gapless-output "${prefix}.single_exon_tblastn.gapless.tsv" \
        --summary "${prefix}.single_exon_tblastn.summary.json" \
        --evalue "$evalue" \
        --workers "$task.cpus" \
        --timeout-seconds "$pair_timeout_seconds"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(tblastn -version | sed -n '1s/^tblastn: //p')
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    touch "${prefix}.single_exon_tblastn.standard.tsv" "${prefix}.single_exon_tblastn.gapless.tsv"
    echo '{"status":"stub"}' > "${prefix}.single_exon_tblastn.summary.json"
    echo '"${task.process}": {blast: "2.16.0+", python: "3.11"}' > versions.yml
    """
}
