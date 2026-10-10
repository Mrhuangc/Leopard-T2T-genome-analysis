process FASTTREE_OR_CLADES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(alignment)

    output:
    tuple val(meta), path("${meta.id}.OR_clades.fasttree.treefile"), emit: tree
    tuple val(meta), path("${meta.id}.OR_clades.fasttree.log"),      emit: log
    path 'versions.yml',                                             emit: versions

    script:
    def prefix = meta.id
    """
    if command -v FastTree >/dev/null 2>&1; then
        FASTTREE_BIN=FastTree
    else
        FASTTREE_BIN=fasttree
    fi

    "\$FASTTREE_BIN" "$alignment" \
        > "${prefix}.OR_clades.fasttree.treefile" \
        2> "${prefix}.OR_clades.fasttree.log"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fasttree: \$("\$FASTTREE_BIN" 2>&1 | sed -n '1s/.*Version //p' | awk '{print \$1}')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    touch "${prefix}.OR_clades.fasttree.treefile" "${prefix}.OR_clades.fasttree.log"
    echo '"${task.process}": {fasttree: "2.1.11"}' > versions.yml
    """
}
