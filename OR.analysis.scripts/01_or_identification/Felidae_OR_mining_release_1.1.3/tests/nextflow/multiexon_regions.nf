nextflow.enable.dsl = 2

include { BUILD_MULTIEXON_REGIONS } from '../../modules/local/build_multiexon_regions/main'

workflow {
    if (!params.species) {
        error "Missing required parameter: species"
    }
    def meta = [id: params.species]
    def merged = file(
        "${params.outdir}/04_species_OR_expansion/merged/${params.species}.expanded.merged.tblastn.tsv",
        checkIfExists: true,
    )
    def known_coordinates = file(
        "${params.outdir}/03_OR_Finder_results/phylogeny/filtered/${params.species}.tree_filtered.tsv",
        checkIfExists: true,
    )
    def package_src = file("${launchDir}/src/or_mining", checkIfExists: true)

    BUILD_MULTIEXON_REGIONS(
        channel.of([meta, merged]),
        channel.value(known_coordinates),
        channel.value(package_src),
        channel.value(params.maximum_intron_length),
        channel.value(params.confirmed_or_clearance),
    )
}
