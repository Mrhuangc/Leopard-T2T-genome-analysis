nextflow.enable.dsl = 2

include { PREPARE_EXONERATE_FUNCTIONAL_SMOKE } from '../../modules/local/prepare_exonerate_functional_smoke/main'
include { EXONERATE_PROTEIN2GENOME as EXONERATE_FUNCTIONAL_SMOKE } from '../../modules/local/exonerate_protein2genome/main'
include { VALIDATE_EXONERATE_FUNCTIONAL_SMOKE } from '../../modules/local/validate_exonerate_functional_smoke/main'

workflow {
    if (!params.species) {
        error "Missing required parameter: species"
    }

    def meta = [id: params.species]
    def package_src = file("${launchDir}/src/or_mining", checkIfExists: true)
    def timeout_seconds = (params.exonerate_smoke_timeout_minutes as int) * 60

    PREPARE_EXONERATE_FUNCTIONAL_SMOKE(
        channel.value(meta),
        channel.value(package_src),
    )
    EXONERATE_FUNCTIONAL_SMOKE(
        PREPARE_EXONERATE_FUNCTIONAL_SMOKE.out.inputs,
        channel.value(params.exonerate_smoke_exhaustive),
        channel.value(false),
        channel.value(params.exonerate_min_intron),
        channel.value(params.maximum_intron_length),
        channel.value(timeout_seconds),
        channel.value(true),
        channel.value(package_src),
    )
    VALIDATE_EXONERATE_FUNCTIONAL_SMOKE(
        EXONERATE_FUNCTIONAL_SMOKE.out.artifacts,
        channel.value(package_src),
    )
}
