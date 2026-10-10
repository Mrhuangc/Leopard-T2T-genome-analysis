nextflow.enable.dsl = 2

include { MAFFT_ADD_OR_CLADES } from '../../modules/local/mafft_add_or_clades/main'
include { FASTTREE_OR_CLADES } from '../../modules/local/fasttree_or_clades/main'
include { BLAST_OR_CLADES } from '../../modules/local/blast_or_clades/main'
include { CLASSIFY_OR_CLADES } from '../../modules/local/classify_or_clades/main'

workflow {
    if (!params.species) {
        error "Missing required parameter: species"
    }

    def meta = [id: params.species]
    def prefix = params.species
    def stage11 = "${params.outdir}/11_complete_OR_catalog"
    def package_src = file("${launchDir}/src/or_mining", checkIfExists: true)
    def template = file(params.template_alignment, checkIfExists: true)
    def type1 = file(params.known_type1_ids, checkIfExists: true)
    def type2 = file(params.known_type2_ids, checkIfExists: true)
    def reference = file(params.or_clade_db ?: params.or_db, checkIfExists: true)
    def catalog = file("${stage11}/${prefix}.complete_OR_catalog.tsv", checkIfExists: true)
    def functional = file("${stage11}/${prefix}.complete_functional_OR.prot.fa", checkIfExists: true)
    def pseudogene = file("${stage11}/${prefix}.complete_pseudogene_OR.nucl.fa", checkIfExists: true)
    def truncated = file("${stage11}/${prefix}.complete_truncated_OR.nucl.fa", checkIfExists: true)
    def ambiguous = file("${stage11}/${prefix}.complete_ambiguous_OR.nucl.fa", checkIfExists: true)

    MAFFT_ADD_OR_CLADES(
        channel.of(tuple(meta, functional)),
        channel.value(template),
    )
    FASTTREE_OR_CLADES(MAFFT_ADD_OR_CLADES.out.alignment)
    BLAST_OR_CLADES(
        channel.of(tuple(meta, functional, pseudogene, truncated, ambiguous)),
        channel.value(reference),
        channel.value(params.clade_blast_evalue),
    )
    CLASSIFY_OR_CLADES(
        FASTTREE_OR_CLADES.out.tree,
        BLAST_OR_CLADES.out.hits,
        channel.value(template),
        channel.value(type1),
        channel.value(type2),
        channel.value(catalog),
        channel.value(functional),
        channel.value(package_src),
        channel.value(params.clade_minimum_tree_support),
    )
}
