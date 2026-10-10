nextflow.enable.dsl = 2

include { BUILD_CANDIDATE_REGIONS } from '../../modules/local/build_candidate_regions/main'

workflow {
    if (!params.input_merged) {
        error "Missing required parameter: input_merged"
    }
    if (!params.species) {
        error "Missing required parameter: species"
    }

    def meta = [id: params.species]
    def merged = file(params.input_merged, checkIfExists: true)
    def package_src = file("${launchDir}/src/or_mining", checkIfExists: true)

    BUILD_CANDIDATE_REGIONS(
        channel.of([meta, merged]),
        channel.value(package_src),
        channel.value(params.candidate_min_width),
        channel.value(params.candidate_flank),
    )
}

