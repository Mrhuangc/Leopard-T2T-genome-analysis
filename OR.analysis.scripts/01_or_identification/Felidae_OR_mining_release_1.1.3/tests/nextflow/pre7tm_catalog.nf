nextflow.enable.dsl = 2

include { BUILD_PRE7TM_CATALOG } from '../../modules/local/build_pre7tm_catalog/main'

workflow {
    if (!params.species) {
        error "Missing required parameter: species"
    }
    if (!params.genome) {
        error "Missing required parameter: genome"
    }

    def meta = [id: params.species]
    def genome_fai = file("${params.genome}.fai", checkIfExists: true)
    def filtered_dir = "${params.outdir}/03_OR_Finder_results/phylogeny/filtered"
    def classified_dir = "${params.outdir}/05_multiexon_OR/classification"
    def package_src = file("${launchDir}/src/or_mining", checkIfExists: true)

    BUILD_PRE7TM_CATALOG(
        channel.of(tuple(
            meta,
            genome_fai,
            file("${filtered_dir}/${params.species}.tree_filtered.nt.fa", checkIfExists: true),
            file("${filtered_dir}/${params.species}.tree_filtered.prot.fa", checkIfExists: true),
            file("${filtered_dir}/${params.species}.tree_filtered.tsv", checkIfExists: true),
            file("${classified_dir}/${params.species}.multiexon_calls.tsv", checkIfExists: true),
            file("${classified_dir}/${params.species}.multiexon_calls.gff3", checkIfExists: true),
            file("${classified_dir}/${params.species}.multiexon_functional.nucl.fa", checkIfExists: true),
            file("${classified_dir}/${params.species}.multiexon_functional.prot.fa", checkIfExists: true),
            file("${classified_dir}/${params.species}.multiexon_pseudogene.nucl.fa", checkIfExists: true),
            file("${classified_dir}/${params.species}.multiexon_truncated.nucl.fa", checkIfExists: true),
            file("${classified_dir}/${params.species}.multiexon_ambiguous.nucl.fa", checkIfExists: true),
        )),
        channel.value(package_src),
    )
}
