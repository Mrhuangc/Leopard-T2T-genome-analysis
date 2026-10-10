nextflow.enable.dsl = 2

include { FINAL_HOMOLOGY_BLAST } from '../../modules/local/final_homology_blast/main'
include { FINALIZE_NC_HOMOLOGY } from '../../modules/local/finalize_nc_homology/main'

workflow {
    if (!params.species) {
        error "Missing required parameter: species"
    }
    if (!params.gpcr_db) {
        error "Missing required parameter: gpcr_db"
    }

    def meta = [id: params.species]
    def prefix = params.species
    def stage08 = "${params.outdir}/08_final_OR_catalog"
    def catalog_inputs = channel.value(tuple(
        meta,
        file("${stage08}/${prefix}.final_OR_catalog.tsv", checkIfExists: true),
        file("${stage08}/${prefix}.final_OR_catalog.gff3", checkIfExists: true),
        file("${stage08}/${prefix}.final_functional_OR.nucl.fa", checkIfExists: true),
        file("${stage08}/${prefix}.final_functional_OR.prot.fa", checkIfExists: true),
        file("${stage08}/${prefix}.final_pseudogene_OR.nucl.fa", checkIfExists: true),
        file("${stage08}/${prefix}.final_truncated_OR.nucl.fa", checkIfExists: true),
        file("${stage08}/${prefix}.final_ambiguous_OR.nucl.fa", checkIfExists: true),
    ))
    def blast_inputs = channel.value(tuple(
        meta,
        file("${stage08}/${prefix}.final_functional_OR.prot.fa", checkIfExists: true),
        file("${stage08}/${prefix}.final_pseudogene_OR.nucl.fa", checkIfExists: true),
        file("${stage08}/${prefix}.final_truncated_OR.nucl.fa", checkIfExists: true),
        file("${stage08}/${prefix}.final_ambiguous_OR.nucl.fa", checkIfExists: true),
    ))
    def package_src = file("${launchDir}/src/or_mining", checkIfExists: true)

    FINAL_HOMOLOGY_BLAST(
        blast_inputs,
        file(params.gpcr_db, checkIfExists: true),
        channel.value(params.gpcr_blast_evalue),
    )
    FINALIZE_NC_HOMOLOGY(
        catalog_inputs,
        FINAL_HOMOLOGY_BLAST.out.results,
        channel.value(params.gpcr_blast_evalue),
        channel.value(package_src),
    )
}
