process BLASTX_SINGLE_EXON_RECOVERY {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(calls), path(candidates), path(pseudogene), path(truncated), path(ambiguous)
    path gpcr_db
    val evalue

    output:
    tuple val(meta), path(calls), path(pseudogene), path(truncated), path(ambiguous), path("${meta.id}.single_exon_recovery.blastx.tsv"), emit: results
    path 'versions.yml', emit: versions

    script:
    def prefix = meta.id
    """
    makeblastdb -in "$gpcr_db" -dbtype prot -out gpcr_database

    if grep -q '^>' "$candidates"; then
        blastx \
            -query "$candidates" \
            -db gpcr_database \
            -evalue "$evalue" \
            -max_target_seqs 1 \
            -num_threads "$task.cpus" \
            -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle' \
            -out "${prefix}.single_exon_recovery.blastx.tsv"
    else
        : > "${prefix}.single_exon_recovery.blastx.tsv"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastx -version | sed -n '1s/^blastx: //p')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    touch "${prefix}.single_exon_recovery.blastx.tsv"
    echo '"${task.process}": {blast: "2.16.0+"}' > versions.yml
    """
}
