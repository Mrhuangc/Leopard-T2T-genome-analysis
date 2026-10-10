nextflow.enable.dsl = 2

include { BUILD_SINGLE_EXON_CANDIDATES } from '../../modules/local/build_single_exon_candidates/main'
include { EXTRACT_SINGLE_EXON_TARGETS } from '../../modules/local/extract_single_exon_targets/main'
include { SINGLE_EXON_TBLASTN_PILOT } from '../../modules/local/single_exon_tblastn_pilot/main'

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
    def prefix = params.species
    def stage04 = "${params.outdir}/04_species_OR_expansion"
    def stage09 = "${params.outdir}/09_NC_final_homology/catalog"
    def package_src = file("${launchDir}/src/or_mining", checkIfExists: true)

    BUILD_SINGLE_EXON_CANDIDATES(
        channel.value(tuple(
            meta,
            file("${stage04}/merged/${prefix}.expanded.merged.tblastn.tsv", checkIfExists: true),
            file("${stage09}/${prefix}.NC_final_OR_catalog.tsv", checkIfExists: true),
            file(params.or_db, checkIfExists: true),
            file("${stage04}/queries/${prefix}.clear_OR.prot.fa", checkIfExists: true),
        )),
        channel.value(package_src),
        channel.value(params.single_exon_catalog_clearance),
        channel.value(params.single_exon_recovery_flank),
    )
    EXTRACT_SINGLE_EXON_TARGETS(
        BUILD_SINGLE_EXON_CANDIDATES.out.candidates,
        channel.value(file(params.genome, checkIfExists: true)),
        channel.value(file("${params.genome}.fai", checkIfExists: true)),
        channel.value(package_src),
    )
    SINGLE_EXON_TBLASTN_PILOT(
        EXTRACT_SINGLE_EXON_TARGETS.out.targets,
        channel.value(package_src),
        channel.value(params.single_exon_targeted_evalue),
        channel.value(params.single_exon_pilot_candidates),
        channel.value(params.single_exon_pair_timeout_seconds),
    )
}
