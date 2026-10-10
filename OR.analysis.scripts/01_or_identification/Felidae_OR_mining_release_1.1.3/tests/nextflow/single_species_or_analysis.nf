nextflow.enable.dsl = 2

include { ANALYZE_OR_LOCUS_ARCHITECTURE } from '../../modules/local/analyze_or_locus_architecture/main'

workflow {
    if (!params.species) {
        error "Missing required parameter: species"
    }
    if (!params.genome) {
        error "Missing required parameter: genome"
    }

    def meta = [id: params.species]
    def prefix = params.species
    def stage12 = "${params.outdir}/12_OR_clade_classification/catalog"
    def package_src = file("${launchDir}/src/or_mining", checkIfExists: true)
    ANALYZE_OR_LOCUS_ARCHITECTURE(
        channel.of(tuple(
            meta,
            file("${stage12}/${prefix}.OR_catalog_with_clades.tsv", checkIfExists: true),
        )),
        channel.value(file("${params.genome}.fai", checkIfExists: true)),
        channel.value(package_src),
        channel.value(params.or_operational_cluster_gap_bp),
    )
}
