nextflow.enable.dsl = 2

include { EXTRACT_CANDIDATE_SEQUENCES } from '../../modules/local/extract_candidate_sequences/main'
include { GETORF }                      from '../../modules/local/getorf/main'
include { NORMALIZE_GETORF }            from '../../modules/local/normalize_getorf/main'

workflow {
    if (!params.genome) {
        error "Missing required parameter: genome"
    }
    if (!params.species) {
        error "Missing required parameter: species"
    }

    def meta = [id: params.species]
    def genome = file(params.genome, checkIfExists: true)
    def regions_dir = file(
        "${params.outdir}/03_OR_Finder_results/candidate_regions",
        checkIfExists: true,
    )
    def regions_tsv = file(
        "${regions_dir}/${params.species}.candidate_regions.tsv",
        checkIfExists: true,
    )
    def faidx_regions = file(
        "${regions_dir}/${params.species}.candidate_regions.faidx.txt",
        checkIfExists: true,
    )
    def package_src = file("${launchDir}/src/or_mining", checkIfExists: true)

    EXTRACT_CANDIDATE_SEQUENCES(
        channel.of([meta, genome, regions_tsv, faidx_regions]),
        channel.value(package_src),
    )
    GETORF(
        EXTRACT_CANDIDATE_SEQUENCES.out.sequences,
        channel.value(params.min_orf_nt),
    )
    NORMALIZE_GETORF(
        GETORF.out.orfs,
        channel.value(regions_tsv),
        channel.value(package_src),
    )
}

