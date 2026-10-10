process MAKEBLASTDB_NUCL {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/../blastp_gpcr/environment.yml"

    input:
    tuple val(meta), path(genome)

    output:
    tuple val(meta), path("${meta.id}.genome_db.*"), emit: db
    path 'versions.yml',                              emit: versions

    script:
    def prefix = meta.id
    """
    makeblastdb \
        -in "$genome" \
        -dbtype nucl \
        -out "${prefix}.genome_db"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(makeblastdb -version | sed -n '1s/^makeblastdb: //p')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    touch "${prefix}.genome_db.ndb"
    echo '"${task.process}": {blast: "2.16.0+"}' > versions.yml
    """
}
