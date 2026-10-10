process PREPARE_SPECIES_OR_QUERIES {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(protein_fasta)
    path nucleotide_fasta
    path package_src
    val cluster_threshold
    val cluster_identity

    output:
    tuple val(meta), path("${meta.id}.species_OR_queries.prot.fa"), emit: queries
    tuple val(meta), path("${meta.id}.clear_OR.prot.fa"),          emit: clear_proteins
    tuple val(meta), path("${meta.id}.clear_OR.ids.txt"),          emit: clear_ids
    tuple val(meta), path("${meta.id}.ambiguous_OR.ids.txt"),      emit: ambiguous_ids
    tuple val(meta), path("${meta.id}.species_OR_queries.prot.fa.clstr"), emit: clusters
    tuple val(meta), path("${meta.id}.species_OR_queries.summary.json"),  emit: summary
    path 'versions.yml',                                                   emit: versions

    script:
    def prefix = meta.id
    """
    PYTHONPATH="\$PWD" python -m or_mining.prepare_species_queries \
        "$protein_fasta" \
        "$nucleotide_fasta" \
        --output-prefix "$prefix" \
        --cluster-threshold "$cluster_threshold" \
        --cluster-identity "$cluster_identity" \
        --cpus "$task.cpus"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cd-hit: \$(cd-hit 2>&1 | sed -n 's/.*CD-HIT version \([^ ]*\).*/\1/p' | head -n 1)
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    touch "${prefix}.species_OR_queries.prot.fa"
    touch "${prefix}.clear_OR.prot.fa"
    touch "${prefix}.clear_OR.ids.txt"
    touch "${prefix}.ambiguous_OR.ids.txt"
    touch "${prefix}.species_OR_queries.prot.fa.clstr"
    echo '{"status":"stub"}' > "${prefix}.species_OR_queries.summary.json"
    echo '"${task.process}": {cd-hit: "4.8.1", python: "3.11"}' > versions.yml
    """
}
