process NC_R_TREE_PARSER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(tree_file)
    path type1_ids
    path type2_ids
    path python_ids
    path parser_script

    output:
    tuple val(meta), path("${meta.id}.nc_r.ids.txt"),                 emit: ids
    tuple val(meta), path("${meta.id}.nc_r_only.ids.txt"),            emit: r_only
    tuple val(meta), path("${meta.id}.python_only.ids.txt"),          emit: python_only
    tuple val(meta), path("${meta.id}.nc_r.rooted.treefile"),         emit: rooted_tree
    tuple val(meta), path("${meta.id}.nc_r_parity.summary.json"),     emit: summary
    path 'versions.yml',                                               emit: versions

    script:
    def prefix = meta.id
    """
    Rscript "$parser_script" \
        "$tree_file" \
        "$type1_ids" \
        "$type2_ids" \
        "$python_ids" \
        "$prefix"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r: \$(Rscript -e 'cat(paste(R.version\$major, R.version\$minor, sep="."))')
        ape: \$(Rscript -e 'cat(as.character(packageVersion("ape")))')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    touch "${prefix}.nc_r.ids.txt"
    touch "${prefix}.nc_r_only.ids.txt"
    touch "${prefix}.python_only.ids.txt"
    touch "${prefix}.nc_r.rooted.treefile"
    echo '{"status":"stub"}' > "${prefix}.nc_r_parity.summary.json"
    echo '"${task.process}": {r: "4.4", ape: "5.8-1"}' > versions.yml
    """
}
