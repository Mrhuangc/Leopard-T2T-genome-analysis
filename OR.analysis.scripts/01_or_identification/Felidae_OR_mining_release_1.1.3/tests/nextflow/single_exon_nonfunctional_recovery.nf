nextflow.enable.dsl = 2

include { BUILD_SINGLE_EXON_CANDIDATES } from '../../modules/local/build_single_exon_candidates/main'
include { EXTRACT_SINGLE_EXON_TARGETS } from '../../modules/local/extract_single_exon_targets/main'
include { SINGLE_EXON_TBLASTN_PILOT } from '../../modules/local/single_exon_tblastn_pilot/main'
include { SINGLE_EXON_TBLASTN_SEARCH } from '../../modules/local/single_exon_tblastn_search/main'
include { CLASSIFY_SINGLE_EXON_RECOVERY } from '../../modules/local/classify_single_exon_recovery/main'
include { BLASTX_SINGLE_EXON_RECOVERY } from '../../modules/local/blastx_single_exon_recovery/main'
include { BUILD_COMPLETE_OR_CATALOG } from '../../modules/local/build_complete_or_catalog/main'

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
    if (!params.gpcr_db) {
        error "Missing required parameter: gpcr_db"
    }

    def meta = [id: params.species]
    def prefix = params.species
    def stage04 = "${params.outdir}/04_species_OR_expansion"
    def stage09 = "${params.outdir}/09_NC_final_homology/catalog"
    def package_src = file("${launchDir}/src/or_mining", checkIfExists: true)
    def genome = file(params.genome, checkIfExists: true)
    def genome_fai = file("${params.genome}.fai", checkIfExists: true)
    def existing_catalog = channel.value(tuple(
        meta,
        file("${stage09}/${prefix}.NC_final_OR_catalog.tsv", checkIfExists: true),
        file("${stage09}/${prefix}.NC_final_OR_catalog.gff3", checkIfExists: true),
        file("${stage09}/${prefix}.NC_final_functional_OR.nucl.fa", checkIfExists: true),
        file("${stage09}/${prefix}.NC_final_functional_OR.prot.fa", checkIfExists: true),
        file("${stage09}/${prefix}.NC_final_pseudogene_OR.nucl.fa", checkIfExists: true),
        file("${stage09}/${prefix}.NC_final_truncated_OR.nucl.fa", checkIfExists: true),
        file("${stage09}/${prefix}.NC_final_ambiguous_OR.nucl.fa", checkIfExists: true),
    ))

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
        channel.value(genome),
        channel.value(genome_fai),
        channel.value(package_src),
    )
    SINGLE_EXON_TBLASTN_PILOT(
        EXTRACT_SINGLE_EXON_TARGETS.out.targets,
        channel.value(package_src),
        channel.value(params.single_exon_targeted_evalue),
        channel.value(params.single_exon_pilot_candidates),
        channel.value(params.single_exon_pair_timeout_seconds),
    )
    SINGLE_EXON_TBLASTN_SEARCH(
        EXTRACT_SINGLE_EXON_TARGETS.out.targets,
        SINGLE_EXON_TBLASTN_PILOT.out.gate,
        channel.value(package_src),
        channel.value(params.single_exon_targeted_evalue),
        channel.value(params.single_exon_pair_timeout_seconds),
    )
    CLASSIFY_SINGLE_EXON_RECOVERY(
        SINGLE_EXON_TBLASTN_SEARCH.out.results,
        channel.value(package_src),
    )
    BLASTX_SINGLE_EXON_RECOVERY(
        CLASSIFY_SINGLE_EXON_RECOVERY.out.classified,
        file(params.gpcr_db, checkIfExists: true),
        channel.value(params.gpcr_blast_evalue),
    )
    BUILD_COMPLETE_OR_CATALOG(
        existing_catalog,
        BLASTX_SINGLE_EXON_RECOVERY.out.results,
        channel.value(package_src),
    )
}
