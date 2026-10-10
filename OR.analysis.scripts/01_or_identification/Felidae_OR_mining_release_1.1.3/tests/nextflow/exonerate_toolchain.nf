nextflow.enable.dsl = 2

include { EXONERATE_SMOKE } from '../../modules/local/exonerate_smoke/main'

workflow {
    if (!params.species) {
        error "Missing required parameter: species"
    }
    EXONERATE_SMOKE(channel.of([id: params.species]))
}
