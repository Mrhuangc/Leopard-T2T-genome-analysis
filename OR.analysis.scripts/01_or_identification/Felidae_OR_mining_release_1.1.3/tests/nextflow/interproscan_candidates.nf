nextflow.enable.dsl = 2

include { PREPARE_INTERPROSCAN_SMOKE } from '../../modules/local/prepare_interproscan_smoke/main'
include { INTERPROSCAN_FUNCTIONAL_SMOKE } from '../../modules/local/interproscan_functional_smoke/main'
include { VALIDATE_INTERPROSCAN_SMOKE } from '../../modules/local/validate_interproscan_smoke/main'
include { INTERPROSCAN_CANDIDATE_SCAN } from '../../modules/local/interproscan_candidate_scan/main'
include { FINALIZE_OR_CATALOG } from '../../modules/local/finalize_or_catalog/main'

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
    def pre7tm_dir = "${params.outdir}/06_OR_catalog/pre7tm"
    def prefix = params.species
    def reference_or = file(params.or_db, checkIfExists: true)
    def package_src = file("${launchDir}/src/or_mining", checkIfExists: true)
    def candidate_proteins = file(
        "${pre7tm_dir}/${prefix}.pre7tm_functional_candidates.prot.fa",
        checkIfExists: true,
    )
    def pre7tm_inputs = channel.value(tuple(
        meta,
        file("${pre7tm_dir}/${prefix}.pre7tm_OR_catalog.tsv", checkIfExists: true),
        file("${pre7tm_dir}/${prefix}.pre7tm_OR_catalog.gff3", checkIfExists: true),
        file("${pre7tm_dir}/${prefix}.pre7tm_functional_candidates.nucl.fa", checkIfExists: true),
        candidate_proteins,
        file("${pre7tm_dir}/${prefix}.pre7tm_pseudogene.nucl.fa", checkIfExists: true),
        file("${pre7tm_dir}/${prefix}.pre7tm_truncated.nucl.fa", checkIfExists: true),
        file("${pre7tm_dir}/${prefix}.pre7tm_ambiguous.nucl.fa", checkIfExists: true),
    ))

    PREPARE_INTERPROSCAN_SMOKE(
        channel.value(tuple(meta, reference_or)),
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
    def smoke_validation_inputs = PREPARE_INTERPROSCAN_SMOKE.out.fasta
        .join(INTERPROSCAN_FUNCTIONAL_SMOKE.out.result)
        .map { meta_map, proteins, result_tsv, version_report ->
            tuple(meta_map, proteins, result_tsv, version_report)
        }
    VALIDATE_INTERPROSCAN_SMOKE(
        smoke_validation_inputs,
        channel.value(params.interproscan_version),
        channel.value(package_src),
    )

    INTERPROSCAN_CANDIDATE_SCAN(
        channel.value(tuple(meta, candidate_proteins)),
        VALIDATE_INTERPROSCAN_SMOKE.out.gate.map { _meta_map, gate -> gate },
        channel.value(params.interproscan_dir),
        channel.value(params.interproscan_version),
        channel.value(params.interproscan_expected_ipr),
        channel.value(params.interproscan_applications),
        channel.value(params.interproscan_candidate_disable_precalc),
    )

    FINALIZE_OR_CATALOG(
        pre7tm_inputs,
        INTERPROSCAN_CANDIDATE_SCAN.out.result,
        channel.value(params.interproscan_version),
        channel.value(params.interproscan_expected_ipr),
        channel.value(package_src),
    )
}
