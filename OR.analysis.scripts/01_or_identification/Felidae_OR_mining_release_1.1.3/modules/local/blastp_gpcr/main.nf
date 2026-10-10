process BLASTP_GPCR {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(protein_orfs)
    path gpcr_db
    val evalue

    output:
    tuple val(meta), path("${meta.id}.gpcr.blastp.tsv"), path(protein_orfs), emit: results
    path 'versions.yml',                                                  emit: versions

    script:
    def prefix = meta.id
    """
    makeblastdb \
        -in "$gpcr_db" \
        -dbtype prot \
        -out gpcr_database

    blastp \
        -query "$protein_orfs" \
        -db gpcr_database \
        -evalue "$evalue" \
        -max_target_seqs 1 \
        -num_threads "$task.cpus" \
        -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle' \
        -out "${prefix}.gpcr.blastp.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastp -version | sed -n '1s/^blastp: //p')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    touch "${prefix}.gpcr.blastp.tsv"
    echo '"${task.process}": {blast: "2.16.0+"}' > versions.yml
    """
}

