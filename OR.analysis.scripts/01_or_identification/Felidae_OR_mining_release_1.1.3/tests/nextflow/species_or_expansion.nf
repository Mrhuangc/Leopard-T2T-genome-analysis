nextflow.enable.dsl = 2

include { PREPARE_SPECIES_OR_QUERIES } from '../../modules/local/prepare_species_or_queries/main'
include { MAKEBLASTDB_NUCL }           from '../../modules/local/makeblastdb_nucl/main'
include { TBLASTN_SPECIES_OR }         from '../../modules/local/tblastn_species_or/main'
include { MERGE_SPECIES_TBLASTN }      from '../../modules/local/merge_species_tblastn/main'

workflow {
    if (!params.species) {
        error "Missing required parameter: species"
    }
    if (!params.genome) {
        error "Missing required parameter: genome"
    }
    if (!params.input_merged) {
        error "Missing required parameter: input_merged"
    }

    def meta = [id: params.species]
    def filtered_dir = file(
        "${params.outdir}/03_OR_Finder_results/phylogeny/filtered",
        checkIfExists: true,
    )
    def protein_fasta = file(
        "${filtered_dir}/${params.species}.tree_filtered.prot.fa",
        checkIfExists: true,
    )
    def nucleotide_fasta = file(
        "${filtered_dir}/${params.species}.tree_filtered.nt.fa",
        checkIfExists: true,
    )
    def genome = file(params.genome, checkIfExists: true)
    def initial_hits = file(params.input_merged, checkIfExists: true)
    def package_src = file("${launchDir}/src/or_mining", checkIfExists: true)

    PREPARE_SPECIES_OR_QUERIES(
        channel.of([meta, protein_fasta]),
        channel.value(nucleotide_fasta),
        channel.value(package_src),
        channel.value(params.species_or_cluster_threshold),
        channel.value(params.species_or_cluster_identity),
    )
    MAKEBLASTDB_NUCL(channel.of([meta, genome]))
    TBLASTN_SPECIES_OR(
        PREPARE_SPECIES_OR_QUERIES.out.queries,
        MAKEBLASTDB_NUCL.out.db,
        channel.value(params.species_tblastn_evalue),
    )
    MERGE_SPECIES_TBLASTN(
        TBLASTN_SPECIES_OR.out.hits,
        channel.value(initial_hits),
        channel.value(package_src),
        channel.value(params.max_hits_per_species_query),
    )
}
