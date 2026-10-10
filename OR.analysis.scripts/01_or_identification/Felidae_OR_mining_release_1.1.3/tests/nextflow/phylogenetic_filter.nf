nextflow.enable.dsl = 2

include { MAFFT_ADD }      from '../../modules/local/mafft_add/main'
include { IQTREE_OR }      from '../../modules/local/iqtree_or/main'
include { FILTER_OR_TREE } from '../../modules/local/filter_or_tree/main'

workflow {
    if (!params.species) {
        error "Missing required parameter: species"
    }

    def meta = [id: params.species]
    def classification_dir = file(
        "${params.outdir}/03_OR_Finder_results/gpcr_classification",
        checkIfExists: true,
    )
    def candidate_proteins = file(
        "${classification_dir}/${params.species}.OR_candidates.prot.fa",
        checkIfExists: true,
    )
    def candidate_nucleotides = file(
        "${classification_dir}/${params.species}.OR_candidates.nt.fa",
        checkIfExists: true,
    )
    def candidate_coordinates = file(
        "${classification_dir}/${params.species}.OR_candidates.tsv",
        checkIfExists: true,
    )
    def template_alignment = file(params.template_alignment, checkIfExists: true)
    def type1_ids = file(params.known_type1_ids, checkIfExists: true)
    def type2_ids = file(params.known_type2_ids, checkIfExists: true)
    def package_src = file("${launchDir}/src/or_mining", checkIfExists: true)

    MAFFT_ADD(
        channel.of([meta, candidate_proteins]),
        channel.value(template_alignment),
    )
    IQTREE_OR(
        MAFFT_ADD.out.alignment,
        channel.value(params.iqtree_model),
        channel.value(params.iqtree_iterations),
        channel.value(params.iqtree_seed),
    )
    FILTER_OR_TREE(
        IQTREE_OR.out.tree,
        channel.value(candidate_nucleotides),
        channel.value(candidate_coordinates),
        channel.value(template_alignment),
        channel.value(type1_ids),
        channel.value(type2_ids),
        channel.value(package_src),
    )
}

