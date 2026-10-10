nextflow.enable.dsl = 2

include { PREPARE_EXONERATE_QUERIES } from '../../modules/local/prepare_exonerate_queries/main'

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
    def package_src = file("${launchDir}/src/or_mining", checkIfExists: true)

    PREPARE_EXONERATE_QUERIES(
        channel.of([meta, reference_or, species_or]),
        channel.value(package_src),
        channel.value(params.exonerate_query_chunks),
    )
}
