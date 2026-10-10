nextflow.enable.dsl = 2

include { EXTRACT_MULTIEXON_SEQUENCES } from '../../modules/local/extract_multiexon_sequences/main'

workflow {
    if (!params.species) {
        error "Missing required parameter: species"
    }
    if (!params.genome) {
        error "Missing required parameter: genome"
    }

    def meta = [id: params.species]
    def candidate_dir = file(
        "${params.outdir}/05_multiexon_OR/candidate_regions",
        checkIfExists: true,
    )
    def genome = file(params.genome, checkIfExists: true)
    def regions_tsv = file(
        "${candidate_dir}/${params.species}.multiexon_regions.tsv",
        checkIfExists: true,
    )
    def faidx_regions = file(
        "${candidate_dir}/${params.species}.multiexon_regions.faidx.txt",
        checkIfExists: true,
    )
    def package_src = file("${launchDir}/src/or_mining", checkIfExists: true)

    EXTRACT_MULTIEXON_SEQUENCES(
        channel.of([meta, genome, regions_tsv, faidx_regions]),
        channel.value(package_src),
    )
}
