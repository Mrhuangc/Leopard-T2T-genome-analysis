process SINGLE_EXON_TBLASTN_PILOT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(candidates_tsv), path(queries), path(targets)
    path package_src
    val evalue
    val pilot_candidates
    val pair_timeout_seconds

    output:
    tuple val(meta), path("${meta.id}.single_exon_tblastn_pilot.ok"), emit: gate
    tuple val(meta), path("${meta.id}.single_exon_tblastn_pilot.summary.json"), path("${meta.id}.single_exon_tblastn_pilot.standard.tsv"), path("${meta.id}.single_exon_tblastn_pilot.gapless.tsv"), emit: results
    path 'versions.yml', emit: versions

    script:
    def prefix = meta.id
    """
    PYTHONPATH="\$PWD" python -m or_mining.run_single_exon_recovery search \
        "$candidates_tsv" \
        "$queries" \
        "$targets" \
        --standard-output "${prefix}.single_exon_tblastn_pilot.standard.tsv" \
        --gapless-output "${prefix}.single_exon_tblastn_pilot.gapless.tsv" \
        --summary "${prefix}.single_exon_tblastn_pilot.summary.json" \
        --gate "${prefix}.single_exon_tblastn_pilot.ok" \
        --evalue "$evalue" \
        --workers "$task.cpus" \
        --timeout-seconds "$pair_timeout_seconds" \
        --limit "$pilot_candidates"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(tblastn -version | sed -n '1s/^tblastn: //p')
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    echo 'status=pass' > "${prefix}.single_exon_tblastn_pilot.ok"
    echo '{"status":"stub"}' > "${prefix}.single_exon_tblastn_pilot.summary.json"
    touch "${prefix}.single_exon_tblastn_pilot.standard.tsv" "${prefix}.single_exon_tblastn_pilot.gapless.tsv"
    echo '"${task.process}": {blast: "2.16.0+", python: "3.11"}' > versions.yml
    """
}
