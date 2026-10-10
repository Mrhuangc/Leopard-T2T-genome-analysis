process EXONERATE_STRICT_REFINEMENT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(query), path(target)
    val min_intron
    val max_intron
    val timeout_seconds
    path package_src

    output:
    tuple val(meta), path("${query.simpleName}__${target.simpleName}.strict.exonerate.txt"), path("${query.simpleName}__${target.simpleName}.strict.exonerate.summary.json"), emit: artifacts
    path "${query.simpleName}__${target.simpleName}.strict.exonerate.stderr.txt", emit: logs
    path "${query.simpleName}__${target.simpleName}.strict.exonerate.versions.yml", emit: versions

    script:
    def cell_name = "${query.simpleName}__${target.simpleName}"
    def output_name = "${cell_name}.strict.exonerate.txt"
    def raw_name = "${cell_name}.strict.exonerate.raw.txt"
    def summary_name = "${cell_name}.strict.exonerate.summary.json"
    def stderr_name = "${cell_name}.strict.exonerate.stderr.txt"
    def versions_name = "${cell_name}.strict.exonerate.versions.yml"
    """
    set +e
    started_epoch=\$(date +%s)
    timeout --signal=TERM --kill-after=30 "${timeout_seconds}s" exonerate \
        -E FALSE \
        -S no \
        --bestn 1 \
        --showtargetgff TRUE \
        --showalignment TRUE \
        --showvulgar TRUE \
        --model protein2genome \
        --ryo '%tcs' \
        --minintron "$min_intron" \
        --maxintron "$max_intron" \
        "$query" \
        "$target" \
        > "$raw_name" \
        2> "$stderr_name"
    command_exit=\$?
    finished_epoch=\$(date +%s)
    elapsed_seconds=\$((finished_epoch - started_epoch))
    set -e

    if [[ "\$command_exit" -ne 0 && "\$command_exit" -ne 124 ]]; then
        exit "\$command_exit"
    fi

    PYTHONPATH="\$PWD" python -m or_mining.finalize_exonerate_chunk \
        "$raw_name" \
        --output "$output_name" \
        --summary "$summary_name" \
        --exit-code "\$command_exit" \
        --timeout-seconds "$timeout_seconds" \
        --exhaustive FALSE \
        --suboptimal no \
        --timeout-is-fatal true \
        --elapsed-seconds "\$elapsed_seconds" \
        --bestn 1

    if [[ "\$command_exit" -eq 124 ]]; then
        echo 'FATAL: strict Exonerate refinement exceeded its bounded timeout' >&2
        exit 42
    fi

    cat <<-END_VERSIONS > "$versions_name"
    "${task.process}":
        exonerate: \$(exonerate --version 2>&1 | sed -n '1s/.*version //p' || true)
    END_VERSIONS
    """

    stub:
    def cell_name = "${query.simpleName}__${target.simpleName}"
    def output_name = "${cell_name}.strict.exonerate.txt"
    def summary_name = "${cell_name}.strict.exonerate.summary.json"
    def stderr_name = "${cell_name}.strict.exonerate.stderr.txt"
    def versions_name = "${cell_name}.strict.exonerate.versions.yml"
    """
    echo 'C4 Alignment:' > "$output_name"
    echo 'vulgar: q 0 1 . t 0 3 + 1 M 1 3' >> "$output_name"
    echo '# --- START OF GFF DUMP ---' >> "$output_name"
    echo '# --- END OF GFF DUMP ---' >> "$output_name"
    echo '{"status":"pass","termination":"completed","elapsed_seconds":0,"output_file":"$output_name","search_parameters":{"model":"protein2genome","exhaustive":false,"suboptimal":false,"bestn":1,"timeout_is_fatal":true},"counts":{"retained_complete_alignment_blocks":1,"retained_vulgar_alignments":1,"retained_gff_dumps":1,"dropped_incomplete_alignment_blocks":0}}' > "$summary_name"
    touch "$stderr_name"
    echo '"${task.process}": {exonerate: "2.4.0"}' > "$versions_name"
    """
}
