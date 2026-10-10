nextflow.enable.dsl = 2

include { SELECT_EXONERATE_MODELS } from '../../modules/local/select_exonerate_models/main'

workflow {
    if (!params.species) {
        error "Missing required parameter: species"
    }
    if ((params.maximum_or_exons as int) < 2) {
        error "maximum_or_exons must be at least 2"
    }
    if ((params.exonerate_min_intron as int) < 1) {
        error "exonerate_min_intron must be positive"
    }
    if ((params.maximum_intron_length as int) < (params.exonerate_min_intron as int)) {
        error "maximum_intron_length must be >= exonerate_min_intron"
    }

    def meta = [id: params.species]
    def raw_exonerate = file(
        "${params.outdir}/05_multiexon_OR/exonerate/merged_targeted/${params.species}.exonerate.raw.txt",
        checkIfExists: true,
    )
    def regions_tsv = file(
        "${params.outdir}/05_multiexon_OR/candidate_regions/${params.species}.multiexon_regions.tsv",
        checkIfExists: true,
    )
    def target_fasta = file(
        "${params.outdir}/05_multiexon_OR/candidate_sequences/${params.species}.multiexon_regions.fa",
        checkIfExists: true,
    )
    def manifest_tsv = file(
        "${params.outdir}/05_multiexon_OR/exonerate/targeted_pairs/${params.species}.targeted_exonerate_pairs.tsv",
        checkIfExists: true,
    )
    def merged_hits = file(
        "${params.outdir}/04_species_OR_expansion/merged/${params.species}.expanded.merged.tblastn.tsv",
        checkIfExists: true,
    )
    def package_src = file("${launchDir}/src/or_mining", checkIfExists: true)

    SELECT_EXONERATE_MODELS(
        channel.of(tuple(meta, raw_exonerate, regions_tsv, target_fasta, manifest_tsv, merged_hits)),
        channel.value(package_src),
        channel.value(params.maximum_or_exons),
        channel.value(params.exonerate_min_intron),
        channel.value(params.maximum_intron_length),
        channel.value(params.exonerate_refinement_flank),
        channel.value(params.exonerate_model_blast_clearance),
        channel.value(params.exonerate_min_reduced_blast_segment),
    )
}
