process MERGE_INITIAL_TBLASTN {
    tag "$species"
    label 'process_low'

    conda "${moduleDir}/../classify_gpcr/environment.yml"

    input:
    path part_files
    path query_manifest
    path package_src
    val species
    val evalue

    output:
    path "${species}.initial_merged.tblastn.tsv",        emit: merged
    path "${species}.initial_merged.validation.json",    emit: validation
    path "${species}.initial_search.summary.json",       emit: summary
    path 'versions.yml',                                  emit: versions

    script:
    def part_args = part_files.collect { path -> "\"${path}\"" }.join(' ')
    """
    PYTHONPATH="\$PWD" python -m or_mining.run_initial_search merge \
        $part_args \
        --manifest "$query_manifest" \
        --output-merged "${species}.initial_merged.tblastn.tsv" \
        --validation "${species}.initial_merged.validation.json" \
        --summary "${species}.initial_search.summary.json" \
        --species "$species" \
        --evalue "$evalue"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def part_args = part_files.collect { path -> "\"${path}\"" }.join(' ')
    """
    cat $part_args > "${species}.initial_merged.tblastn.tsv"
    echo '{"status":"pass","counts":{"invalid_rows":0}}' \
        > "${species}.initial_merged.validation.json"
    echo '{"status":"pass","search_completeness":"complete"}' \
        > "${species}.initial_search.summary.json"
    echo '"${task.process}": {python: "3.11"}' > versions.yml
    """
}
