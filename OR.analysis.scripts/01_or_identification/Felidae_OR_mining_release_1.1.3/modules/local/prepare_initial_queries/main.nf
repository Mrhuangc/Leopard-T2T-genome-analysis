process PREPARE_INITIAL_QUERIES {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/../classify_gpcr/environment.yml"

    input:
    tuple val(meta), path(query_fasta)
    path package_src
    val requested_shards

    output:
    tuple val(meta), path('query_shards/*.prot.fa'),                    emit: shards
    tuple val(meta), path("${meta.id}.initial_query_manifest.tsv"),     emit: manifest
    tuple val(meta), path("${meta.id}.initial_query_shards.summary.json"), emit: summary
    path 'versions.yml',                                                 emit: versions

    script:
    def prefix = meta.id
    """
    PYTHONPATH="\$PWD" python -m or_mining.run_initial_search split \
        "$query_fasta" \
        --output-dir query_shards \
        --output-prefix "$prefix" \
        --shards "$requested_shards"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    mkdir -p query_shards
    printf '>stub_OR\nM\n' > query_shards/shard_001.prot.fa
    printf 'shard_id\tfasta\trecords\tresidues\tfirst_query\tlast_query\n' \
        > "${prefix}.initial_query_manifest.tsv"
    echo '{"status":"pass","counts":{"query_shards":1}}' \
        > "${prefix}.initial_query_shards.summary.json"
    echo '"${task.process}": {python: "3.11"}' > versions.yml
    """
}
