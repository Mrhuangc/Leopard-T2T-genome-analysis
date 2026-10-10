process TBLASTN_INITIAL {
    tag "$meta.id:$query_fasta.simpleName"
    label 'process_high'

    conda "${moduleDir}/../blastp_gpcr/environment.yml"

    input:
    tuple val(meta), path(query_fasta)
    tuple val(db_meta), path(genome_db_files)
    val evalue

    output:
    tuple val(meta), path("${meta.id}.*.initial.tblastn.tsv"), emit: hits
    path 'versions.yml',                                      emit: versions

    script:
    if (meta.id != db_meta.id) {
        error "Query and genome database metadata IDs differ: ${meta.id} versus ${db_meta.id}"
    }
    def shard = query_fasta.name.replace('.prot.fa', '')
    def output_name = "${meta.id}.${shard}.initial.tblastn.tsv"
    """
    tblastn \
        -query "$query_fasta" \
        -db "${meta.id}.genome_db" \
        -evalue "$evalue" \
        -num_threads "$task.cpus" \
        -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe sframe qlen slen' \
        -out "$output_name"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(tblastn -version | sed -n '1s/^tblastn: //p')
    END_VERSIONS
    """

    stub:
    def shard = query_fasta.name.replace('.prot.fa', '')
    def output_name = "${meta.id}.${shard}.initial.tblastn.tsv"
    """
    printf 'stub_OR\tstub_contig\t100\t1\t0\t0\t1\t1\t1\t3\t1e-10\t20\t0\t1\t1\t3\n' \
        > "$output_name"
    echo '"${task.process}": {blast: "2.16.0+"}' > versions.yml
    """
}
