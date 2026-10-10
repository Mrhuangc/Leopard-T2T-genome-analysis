process TBLASTN_SPECIES_OR {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/../blastp_gpcr/environment.yml"

    input:
    tuple val(meta), path(query_fasta)
    tuple val(db_meta), path(genome_db_files)
    val evalue

    output:
    tuple val(meta), path("${meta.id}.species_OR_vs_genome.tblastn.tsv"), emit: hits
    path 'versions.yml',                                                       emit: versions

    script:
    if (meta.id != db_meta.id) {
        error "Query and genome database metadata IDs differ: ${meta.id} versus ${db_meta.id}"
    }
    def prefix = meta.id
    """
    tblastn \
        -query "$query_fasta" \
        -db "${prefix}.genome_db" \
        -evalue "$evalue" \
        -num_threads "$task.cpus" \
        -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe sframe qlen slen' \
        -out "${prefix}.species_OR_vs_genome.tblastn.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(tblastn -version | sed -n '1s/^tblastn: //p')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    touch "${prefix}.species_OR_vs_genome.tblastn.tsv"
    echo '"${task.process}": {blast: "2.16.0+"}' > versions.yml
    """
}
