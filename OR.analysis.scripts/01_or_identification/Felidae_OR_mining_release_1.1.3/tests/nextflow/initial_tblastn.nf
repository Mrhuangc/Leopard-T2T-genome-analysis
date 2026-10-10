nextflow.enable.dsl = 2

include { PREPARE_INITIAL_QUERIES } from '../../modules/local/prepare_initial_queries/main'
include { MAKEBLASTDB_NUCL }       from '../../modules/local/makeblastdb_nucl/main'
include { TBLASTN_INITIAL }         from '../../modules/local/tblastn_initial/main'
include { MERGE_INITIAL_TBLASTN }   from '../../modules/local/merge_initial_tblastn/main'

workflow {
    if (!params.species) {
        error "Missing required parameter: species"
    }
    if (!params.genome) {
        error "Missing required parameter: genome"
    }
    if (!params.or_db) {
        error "Missing required parameter: or_db"
    }

    def meta = [id: params.species]
    def genome = file(params.genome, checkIfExists: true)
    def query_fasta = file(params.or_db, checkIfExists: true)
    def package_src = file("${launchDir}/src/or_mining", checkIfExists: true)

    PREPARE_INITIAL_QUERIES(
        channel.of([meta, query_fasta]),
        channel.value(package_src),
        channel.value(params.initial_query_shards),
    )
    MAKEBLASTDB_NUCL(channel.of([meta, genome]))

    def query_shards = PREPARE_INITIAL_QUERIES.out.shards.flatMap { meta_map, shard_files ->
        def files = shard_files instanceof List ? shard_files : [shard_files]
        files.collect { shard_file -> tuple(meta_map, shard_file) }
    }
    def reusable_genome_db = MAKEBLASTDB_NUCL.out.db.first()
    TBLASTN_INITIAL(
        query_shards,
        reusable_genome_db,
        channel.value(params.initial_tblastn_evalue),
    )

    def completed_parts = TBLASTN_INITIAL.out.hits
        .map { meta_map, part_file -> part_file }
        .collect()
    def query_manifest = PREPARE_INITIAL_QUERIES.out.manifest
        .map { meta_map, manifest_file -> manifest_file }
        .first()
    MERGE_INITIAL_TBLASTN(
        completed_parts,
        query_manifest,
        channel.value(package_src),
        channel.value(params.species),
        channel.value(params.initial_tblastn_evalue),
    )
}
