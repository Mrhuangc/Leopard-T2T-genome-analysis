nextflow.enable.dsl = 2

include { CLASSIFY_MULTIEXON_MODELS } from '../../modules/local/classify_multiexon_models/main'

workflow {
    if (!params.species) {
        error "Missing required parameter: species"
    }
    if (!params.or_db) {
        error "Missing required parameter: or_db"
    }
    if ((params.min_orf_nt as int) < 3 || (params.min_orf_nt as int) % 3 != 0) {
        error "min_orf_nt must be a positive multiple of three"
    }
    if ((params.multiexon_min_model_overlap as double) <= 0 ||
        (params.multiexon_min_model_overlap as double) > 1) {
        error "multiexon_min_model_overlap must be in (0, 1]"
    }

    def meta = [id: params.species]
    def strict_raw = file(
        "${params.outdir}/05_multiexon_OR/exonerate/strict_refinement/merged/${params.species}.strict_refinement.raw.txt",
        checkIfExists: true,
    )
    def strict_manifest = file(
        "${params.outdir}/05_multiexon_OR/exonerate/strict_refinement/pairs/${params.species}.strict_refinement_pairs.tsv",
        checkIfExists: true,
    )
    def regions_tsv = file(
        "${params.outdir}/05_multiexon_OR/candidate_regions/${params.species}.multiexon_regions.tsv",
        checkIfExists: true,
    )
    def regions_fasta = file(
        "${params.outdir}/05_multiexon_OR/candidate_sequences/${params.species}.multiexon_regions.fa",
        checkIfExists: true,
    )
    def reference_or = file(params.or_db, checkIfExists: true)
    def species_or = file(
        "${params.outdir}/04_species_OR_expansion/queries/${params.species}.clear_OR.prot.fa",
        checkIfExists: true,
    )
    def package_src = file("${launchDir}/src/or_mining", checkIfExists: true)

    CLASSIFY_MULTIEXON_MODELS(
        channel.of(tuple(meta, strict_raw, strict_manifest, regions_tsv, regions_fasta, reference_or, species_or)),
        channel.value(package_src),
        channel.value(params.min_orf_nt),
        channel.value(params.multiexon_min_model_overlap),
        channel.value(params.maximum_or_exons),
        channel.value(params.multiexon_edge_distance),
        channel.value(params.multiexon_edge_n_window),
        channel.value(params.multiexon_edge_n_run),
    )
}
