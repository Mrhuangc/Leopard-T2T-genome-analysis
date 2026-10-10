nextflow.enable.dsl = 2

include { PREPARE_STRICT_REFINEMENT } from '../../modules/local/prepare_strict_refinement/main'
include { EXONERATE_STRICT_REFINEMENT as EXONERATE_STRICT_PILOT } from '../../modules/local/exonerate_strict_refinement/main'
include { EXONERATE_STRICT_REFINEMENT as EXONERATE_STRICT_SEARCH } from '../../modules/local/exonerate_strict_refinement/main'
include { CHECK_STRICT_REFINEMENT_PILOT } from '../../modules/local/check_strict_refinement_pilot/main'
include { MERGE_STRICT_REFINEMENTS } from '../../modules/local/merge_strict_refinements/main'

workflow {
    if (!params.species) {
        error "Missing required parameter: species"
    }
    if (!params.or_db) {
        error "Missing required parameter: or_db"
    }
    if ((params.strict_refinement_pilot_pairs as int) < 1) {
        error "strict_refinement_pilot_pairs must be positive"
    }
    if ((params.strict_refinement_timeout_minutes as int) < 1) {
        error "strict_refinement_timeout_minutes must be positive"
    }

    def meta = [id: params.species]
    def reference_or = file(params.or_db, checkIfExists: true)
    def species_or = file(
        "${params.outdir}/04_species_OR_expansion/queries/${params.species}.clear_OR.prot.fa",
        checkIfExists: true,
    )
    def region_fasta = file(
        "${params.outdir}/05_multiexon_OR/candidate_sequences/${params.species}.multiexon_regions.fa",
        checkIfExists: true,
    )
    def regions_tsv = file(
        "${params.outdir}/05_multiexon_OR/candidate_regions/${params.species}.multiexon_regions.tsv",
        checkIfExists: true,
    )
    def candidates_tsv = file(
        "${params.outdir}/05_multiexon_OR/exonerate/model_selection/${params.species}.multiexon_refinement_candidates.tsv",
        checkIfExists: true,
    )
    def package_src = file("${launchDir}/src/or_mining", checkIfExists: true)
    def timeout_seconds = (params.strict_refinement_timeout_minutes as int) * 60
    def pilot_count = params.strict_refinement_pilot_pairs as int

    PREPARE_STRICT_REFINEMENT(
        channel.of(tuple(meta, reference_or, species_or, region_fasta, regions_tsv, candidates_tsv)),
        channel.value(package_src),
    )

    def strict_pairs = PREPARE_STRICT_REFINEMENT.out.pairs.flatMap {
        pair_meta, query_files, target_files ->
        def targets_by_pair = target_files.collectEntries { target_file ->
            def pair_id = target_file.name.replace('.target.nucl.fa', '')
            [(pair_id): target_file]
        }
        query_files.sort { query_file -> query_file.name }.withIndex().collect {
            query_file, index ->
            def pair_id = query_file.name.replace('.query.prot.fa', '')
            def target_file = targets_by_pair[pair_id]
            if (target_file == null) {
                error "Missing strict-refinement target FASTA for ${pair_id}"
            }
            def candidate_id = pair_id.tokenize('__')[-1]
            tuple(
                [id: candidate_id, species: pair_meta.id, pilot: index < pilot_count],
                query_file,
                target_file,
            )
        }
    }
    def pair_branches = strict_pairs.branch { pair_meta, query_file, target_file ->
        pilot: pair_meta.pilot
        production: true
    }

    EXONERATE_STRICT_PILOT(
        pair_branches.pilot,
        channel.value(params.exonerate_min_intron),
        channel.value(params.maximum_intron_length),
        channel.value(timeout_seconds),
        channel.value(package_src),
    )
    def collected_pilot = EXONERATE_STRICT_PILOT.out.artifacts
        .map { pair_meta, result, status ->
            tuple(pair_meta.species, [id: pair_meta.species], result, status)
        }
        .groupTuple(by: 0)
        .map { species_id, metas, results, statuses ->
            tuple(metas[0], results, statuses)
        }
    CHECK_STRICT_REFINEMENT_PILOT(
        collected_pilot,
        channel.value(package_src),
        channel.value(pilot_count),
    )

    def gated_production_pairs = pair_branches.production
        .combine(CHECK_STRICT_REFINEMENT_PILOT.out.gate)
        .map { pair_meta, query_file, target_file, pilot_gate ->
            tuple(pair_meta, query_file, target_file)
        }
    EXONERATE_STRICT_SEARCH(
        gated_production_pairs,
        channel.value(params.exonerate_min_intron),
        channel.value(params.maximum_intron_length),
        channel.value(timeout_seconds),
        channel.value(package_src),
    )

    def all_artifacts = EXONERATE_STRICT_PILOT.out.artifacts
        .mix(EXONERATE_STRICT_SEARCH.out.artifacts)
    def collected_results = all_artifacts
        .map { pair_meta, result, status ->
            tuple(pair_meta.species, [id: pair_meta.species], result, status)
        }
        .groupTuple(by: 0)
        .map { species_id, metas, results, statuses ->
            tuple(metas[0], results, statuses)
        }
    def results_with_manifest = collected_results
        .combine(PREPARE_STRICT_REFINEMENT.out.manifest)
        .map { result_meta, results, statuses, manifest_meta, manifest ->
            tuple(result_meta, results, statuses, manifest)
        }
    MERGE_STRICT_REFINEMENTS(
        results_with_manifest,
        channel.value(package_src),
    )
}
