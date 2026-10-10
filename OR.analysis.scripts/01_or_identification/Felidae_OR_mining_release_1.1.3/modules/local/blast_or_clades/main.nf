process BLAST_OR_CLADES {
    tag "$meta.id"
    label 'process_medium'

    // Reuse the stage-09 BLAST toolchain instead of solving a duplicate env.
    conda "${moduleDir}/../final_homology_blast/environment.yml"

    input:
    tuple val(meta), path(functional_protein), path(pseudogene_nucleotide), path(truncated_nucleotide), path(ambiguous_nucleotide)
    path clade_reference_protein
    val blast_evalue

    output:
    tuple val(meta), path("${meta.id}.functional.clade_blastp.tsv"), path("${meta.id}.pseudogene.clade_blastx.tsv"), path("${meta.id}.truncated.clade_blastx.tsv"), path("${meta.id}.ambiguous.clade_blastx.tsv"), emit: hits
    path 'versions.yml', emit: versions

    script:
    def prefix = meta.id
    def fields = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen'
    """
    makeblastdb -in "$clade_reference_protein" -dbtype prot -out clade_reference >/dev/null

    blastp -query "$functional_protein" -db clade_reference \
        -evalue "$blast_evalue" -max_target_seqs 5 -num_threads "$task.cpus" \
        -outfmt '6 $fields' -out "${prefix}.functional.clade_blastp.tsv"

    for CLASS in pseudogene truncated ambiguous; do
        case "\$CLASS" in
            pseudogene) QUERY="$pseudogene_nucleotide" ;;
            truncated) QUERY="$truncated_nucleotide" ;;
            ambiguous) QUERY="$ambiguous_nucleotide" ;;
        esac
        OUTPUT="${prefix}.\${CLASS}.clade_blastx.tsv"
        if grep -q '^>' "\$QUERY"; then
            blastx -query "\$QUERY" -db clade_reference \
                -evalue "$blast_evalue" -max_target_seqs 5 -num_threads "$task.cpus" \
                -outfmt '6 $fields' -out "\$OUTPUT"
        else
            touch "\$OUTPUT"
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastp -version | sed -n '1s/.*: //p')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    touch "${prefix}.functional.clade_blastp.tsv"
    touch "${prefix}.pseudogene.clade_blastx.tsv"
    touch "${prefix}.truncated.clade_blastx.tsv"
    touch "${prefix}.ambiguous.clade_blastx.tsv"
    echo '"${task.process}": {blast: "2.16.0"}' > versions.yml
    """
}
