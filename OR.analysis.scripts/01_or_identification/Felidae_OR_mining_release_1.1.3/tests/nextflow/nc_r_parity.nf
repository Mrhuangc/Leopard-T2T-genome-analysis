nextflow.enable.dsl = 2

include { NC_R_TREE_PARSER } from '../../modules/local/nc_r_tree_parser/main'

workflow {
    if (!params.species) {
        error "Missing required parameter: species"
    }

    def meta = [id: params.species]
    def phylogeny_dir = file(
        "${params.outdir}/03_OR_Finder_results/phylogeny",
        checkIfExists: true,
    )
    def tree_file = file(
        "${phylogeny_dir}/tree/${params.species}.or_phylogeny.treefile",
        checkIfExists: true,
    )
    def python_ids = file(
        "${phylogeny_dir}/filtered/${params.species}.tree_filtered.ids.txt",
        checkIfExists: true,
    )
    def type1_ids = file(params.known_type1_ids, checkIfExists: true)
    def type2_ids = file(params.known_type2_ids, checkIfExists: true)
    def parser_script = file(
        "${launchDir}/modules/local/nc_r_tree_parser/nc_tree_parser.R",
        checkIfExists: true,
    )

    NC_R_TREE_PARSER(
        channel.of([meta, tree_file]),
        channel.value(type1_ids),
        channel.value(type2_ids),
        channel.value(python_ids),
        channel.value(parser_script),
    )
}
