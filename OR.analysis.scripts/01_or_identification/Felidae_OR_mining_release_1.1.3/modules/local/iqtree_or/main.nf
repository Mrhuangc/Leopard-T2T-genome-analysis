process IQTREE_OR {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(alignment), path(candidate_proteins)
    val model
    val iterations
    val seed

    output:
    tuple val(meta), path("${meta.id}.or_phylogeny.treefile"), path(candidate_proteins), emit: tree
    tuple val(meta), path("${meta.id}.or_phylogeny.iqtree"),                        emit: report
    tuple val(meta), path("${meta.id}.or_phylogeny.log"),                           emit: log
    path 'versions.yml',                                                             emit: versions

    script:
    def prefix = meta.id
    """
    if command -v iqtree2 >/dev/null 2>&1; then
        IQTREE_BIN=iqtree2
    else
        IQTREE_BIN=iqtree
    fi

    "\$IQTREE_BIN" \
        -s "$alignment" \
        -st AA \
        -m "$model" \
        -nt "$task.cpus" \
        -n "$iterations" \
        -seed "$seed" \
        -keep-ident \
        -redo \
        --prefix "${prefix}.or_phylogeny"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iqtree: \$("\$IQTREE_BIN" --version | sed -n '1s/.*version //p' | awk '{print \$1}')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    touch "${prefix}.or_phylogeny.treefile"
    touch "${prefix}.or_phylogeny.iqtree"
    touch "${prefix}.or_phylogeny.log"
    echo '"${task.process}": {iqtree: "2.3.6"}' > versions.yml
    """
}

