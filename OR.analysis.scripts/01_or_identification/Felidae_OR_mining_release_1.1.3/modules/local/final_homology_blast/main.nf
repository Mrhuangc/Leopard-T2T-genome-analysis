process FINAL_HOMOLOGY_BLAST {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(functional_protein), path(pseudogene_nucleotide), path(truncated_nucleotide), path(ambiguous_nucleotide)
    path gpcr_db
    val evalue

    output:
    tuple val(meta), path("${meta.id}.final_functional.blastp.tsv"), path("${meta.id}.final_pseudogene.blastx.tsv"), path("${meta.id}.final_truncated.blastx.tsv"), path("${meta.id}.final_ambiguous.blastx.tsv"), emit: results
    path 'versions.yml', emit: versions

    script:
    def prefix = meta.id
    """
    makeblastdb \
        -in "$gpcr_db" \
        -dbtype prot \
        -out gpcr_database

    run_blastp() {
        local query="\$1"
        local output="\$2"
        if grep -q '^>' "\$query"; then
            blastp \
                -query "\$query" \
                -db gpcr_database \
                -evalue "$evalue" \
                -max_target_seqs 1 \
                -num_threads "$task.cpus" \
                -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle' \
                -out "\$output"
        else
            : > "\$output"
        fi
    }

    run_blastx() {
        local query="\$1"
        local output="\$2"
        if grep -q '^>' "\$query"; then
            blastx \
                -query "\$query" \
                -db gpcr_database \
                -evalue "$evalue" \
                -max_target_seqs 1 \
                -num_threads "$task.cpus" \
                -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle' \
                -out "\$output"
        else
            : > "\$output"
        fi
    }

    run_blastp "$functional_protein" "${prefix}.final_functional.blastp.tsv"
    run_blastx "$pseudogene_nucleotide" "${prefix}.final_pseudogene.blastx.tsv"
    run_blastx "$truncated_nucleotide" "${prefix}.final_truncated.blastx.tsv"
    run_blastx "$ambiguous_nucleotide" "${prefix}.final_ambiguous.blastx.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blastp: \$(blastp -version | sed -n '1s/^blastp: //p')
        blastx: \$(blastx -version | sed -n '1s/^blastx: //p')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    touch "${prefix}.final_functional.blastp.tsv"
    touch "${prefix}.final_pseudogene.blastx.tsv"
    touch "${prefix}.final_truncated.blastx.tsv"
    touch "${prefix}.final_ambiguous.blastx.tsv"
    echo '"${task.process}": {blastp: "2.16.0+", blastx: "2.16.0+"}' > versions.yml
    """
}
