nextflow.enable.dsl = 2

include { PREPARE_TARGETED_EXONERATE } from '../../modules/local/prepare_targeted_exonerate/main'

workflow {
    if (!params.species) {
        error "Missing required parameter: species"
    }
    if (!params.or_db) {
        error "Missing required parameter: or_db"
    }

    def meta = [id: params.species]
    def reference_or = file(params.or_db, checkIfExists: true)
    def species_or = file(
        "${params.outdir}/04_species_OR_expansion/queries/${params.species}.clear_OR.prot.fa",
        checkIfExists: true,
    )
    def target_fasta = file(
        "${params.outdir}/05_multiexon_OR/candidate_sequences/${params.species}.multiexon_regions.fa",
        checkIfExists: true,
    )
    def regions_tsv = file(
        "${params.outdir}/05_multiexon_OR/candidate_regions/${params.species}.multiexon_regions.tsv",
        checkIfExists: true,
    )
    def merged_hits = file(
        "${params.outdir}/04_species_OR_expansion/merged/${params.species}.expanded.merged.tblastn.tsv",
        checkIfExists: true,
    )
    def package_src = file("${launchDir}/src/or_mining", checkIfExists: true)

    PREPARE_TARGETED_EXONERATE(
        channel.of([meta, reference_or, species_or, target_fasta, regions_tsv, merged_hits]),
        channel.value(package_src),
        channel.value(params.exonerate_top_queries_per_region),
    )
}
