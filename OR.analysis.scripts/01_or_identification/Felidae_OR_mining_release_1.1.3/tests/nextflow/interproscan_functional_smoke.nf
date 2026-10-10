nextflow.enable.dsl = 2

include { PREPARE_INTERPROSCAN_SMOKE } from '../../modules/local/prepare_interproscan_smoke/main'
include { INTERPROSCAN_FUNCTIONAL_SMOKE } from '../../modules/local/interproscan_functional_smoke/main'
include { VALIDATE_INTERPROSCAN_SMOKE } from '../../modules/local/validate_interproscan_smoke/main'

workflow {
    if (!params.species) {
        error "Missing required parameter: species"
    }
    if (!params.or_db) {
        error "Missing required parameter: or_db"
    }
    if (!params.interproscan_dir) {
        error "Missing required parameter: interproscan_dir"
    }

    def meta = [id: params.species]
    def reference_or = file(params.or_db, checkIfExists: true)
    def package_src = file("${launchDir}/src/or_mining", checkIfExists: true)

    PREPARE_INTERPROSCAN_SMOKE(
        channel.of(tuple(meta, reference_or)),
        channel.value(params.interproscan_smoke_records),
        channel.value(package_src),
    )
    INTERPROSCAN_FUNCTIONAL_SMOKE(
        PREPARE_INTERPROSCAN_SMOKE.out.fasta,
        channel.value(params.interproscan_dir),
        channel.value(params.interproscan_version),
        channel.value(params.interproscan_applications),
        channel.value(params.interproscan_smoke_disable_precalc),
    )

    def validation_inputs = PREPARE_INTERPROSCAN_SMOKE.out.fasta
        .join(INTERPROSCAN_FUNCTIONAL_SMOKE.out.result)
        .map { meta_map, proteins, result_tsv, version_report ->
            tuple(meta_map, proteins, result_tsv, version_report)
        }
    VALIDATE_INTERPROSCAN_SMOKE(
        validation_inputs,
        channel.value(params.interproscan_version),
        channel.value(package_src),
    )
}
