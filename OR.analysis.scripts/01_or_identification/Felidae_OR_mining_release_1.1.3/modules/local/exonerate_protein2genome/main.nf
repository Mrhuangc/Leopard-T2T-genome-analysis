process EXONERATE_PROTEIN2GENOME {
    tag "$meta.id:${query_chunk.simpleName}:${target_chunk.simpleName}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(query_chunk), path(target_chunk)
    val exhaustive
    val suboptimal
    val min_intron
    val max_intron
    val timeout_seconds
    val fail_on_timeout
    path package_src

    output:
    tuple val(meta), path("${query_chunk.simpleName}__${target_chunk.simpleName}.exonerate.txt"), path("${query_chunk.simpleName}__${target_chunk.simpleName}.exonerate.summary.json"), emit: artifacts
    path "${query_chunk.simpleName}__${target_chunk.simpleName}.exonerate.stderr.txt", emit: logs
    path "${query_chunk.simpleName}__${target_chunk.simpleName}.exonerate.versions.yml", emit: versions

    script:
    def exhaustive_arg = exhaustive ? 'TRUE' : 'FALSE'
    def suboptimal_arg = suboptimal ? 'yes' : 'no'
    def cell_name = "${query_chunk.simpleName}__${target_chunk.simpleName}"
    def output_name = "${cell_name}.exonerate.txt"
    def raw_name = "${cell_name}.exonerate.raw.txt"
    def summary_name = "${cell_name}.exonerate.summary.json"
    def stderr_name = "${cell_name}.exonerate.stderr.txt"
    def versions_name = "${cell_name}.exonerate.versions.yml"
    """
    set +e
    started_epoch=\$(date +%s)
    timeout --signal=TERM --kill-after=30 "${timeout_seconds}s" exonerate \
        -E "$exhaustive_arg" \
        -S "$suboptimal_arg" \
        --showtargetgff TRUE \
        --showalignment TRUE \
        --showvulgar TRUE \
        --model protein2genome \
        --ryo '%tcs' \
        --minintron "$min_intron" \
        --maxintron "$max_intron" \
        "$query_chunk" \
        "$target_chunk" \
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
        --exhaustive "$exhaustive_arg" \
        --suboptimal "$suboptimal_arg" \
        --timeout-is-fatal "$fail_on_timeout" \
        --elapsed-seconds "\$elapsed_seconds"

    if [[ "\$command_exit" -eq 124 && "$fail_on_timeout" == "true" ]]; then
        echo "FATAL: Exonerate exceeded the bounded timeout; remaining tasks will not be released" >&2
        exit 42
    fi

    cat <<-END_VERSIONS > "$versions_name"
    "${task.process}":
        exonerate: \$(exonerate --version 2>&1 | sed -n '1s/.*version //p' || true)
    END_VERSIONS
    """

    stub:
    def exhaustive_arg = exhaustive ? 'true' : 'false'
    def suboptimal_arg = suboptimal ? 'true' : 'false'
    def cell_name = "${query_chunk.simpleName}__${target_chunk.simpleName}"
    def output_name = "${cell_name}.exonerate.txt"
    def summary_name = "${cell_name}.exonerate.summary.json"
    def stderr_name = "${cell_name}.exonerate.stderr.txt"
    def versions_name = "${cell_name}.exonerate.versions.yml"
    """
    echo '# FELIDAE_OR_CHUNK_STATUS completed' > "$output_name"
    echo '{"schema_version":"1.2","status":"pass","termination":"completed","output_file":"$output_name","search_parameters":{"model":"protein2genome","exhaustive":$exhaustive_arg,"suboptimal":$suboptimal_arg,"bestn":null,"timeout_is_fatal":true},"counts":{"retained_complete_alignment_blocks":0,"dropped_incomplete_alignment_blocks":0}}' > "$summary_name"
    touch "$stderr_name"
    echo '"${task.process}": {exonerate: "2.4.0"}' > "$versions_name"
    """
}
