nextflow.enable.dsl = 2

include { VALIDATE_MERGED } from '../../modules/local/validate_merged/main'

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

    VALIDATE_MERGED(channel.of([meta, merged]), channel.value(package_src))
}

