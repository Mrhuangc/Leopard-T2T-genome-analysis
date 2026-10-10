nextflow.enable.dsl = 2

include { PREPARE_TARGETED_EXONERATE } from '../../modules/local/prepare_targeted_exonerate/main'
include { EXONERATE_PROTEIN2GENOME as EXONERATE_TARGETED_PILOT } from '../../modules/local/exonerate_protein2genome/main'
include { EXONERATE_PROTEIN2GENOME as EXONERATE_TARGETED_SEARCH } from '../../modules/local/exonerate_protein2genome/main'
include { CHECK_EXONERATE_PILOT } from '../../modules/local/check_exonerate_pilot/main'
include { MERGE_EXONERATE_RESULTS } from '../../modules/local/merge_exonerate_results/main'
include { PREPARE_EXONERATE_FUNCTIONAL_SMOKE } from '../../modules/local/prepare_exonerate_functional_smoke/main'
include { EXONERATE_PROTEIN2GENOME as EXONERATE_FUNCTIONAL_SMOKE } from '../../modules/local/exonerate_protein2genome/main'
include { VALIDATE_EXONERATE_FUNCTIONAL_SMOKE } from '../../modules/local/validate_exonerate_functional_smoke/main'
include { CHECK_EXONERATE_SMOKE_GATE } from '../../modules/local/check_exonerate_smoke_gate/main'

workflow {
    if (!params.species) {
        error "Missing required parameter: species"
    }
    if (!params.or_db) {
        error "Missing required parameter: or_db"
    }
    if (params.exonerate_exhaustive) {
        error "Targeted discovery requires exonerate_exhaustive=false; exhaustive all-region search is computationally unbounded"
    }
    if (params.exonerate_suboptimal) {
        error "Targeted discovery requires exonerate_suboptimal=false"
    }
    if ((params.exonerate_pilot_pairs as int) < 1) {
        error "exonerate_pilot_pairs must be at least 1"
    }
    if ((params.exonerate_top_queries_per_region as int) < 1) {
        error "exonerate_top_queries_per_region must be at least 1"
    }
    if ((params.exonerate_soft_timeout_minutes as int) < 1) {
        error "exonerate_soft_timeout_minutes must be at least 1"
    }
    if ((params.exonerate_pilot_min_positive_pairs as int) < 1 ||
        (params.exonerate_pilot_min_positive_pairs as int) >
        (params.exonerate_pilot_pairs as int)) {
        error "exonerate_pilot_min_positive_pairs must be within 1..exonerate_pilot_pairs"
    }

    def meta = [id: params.species]
    def reference_or = file(params.or_db, checkIfExists: true)
    def species_or = file(
        "${params.outdir}/04_species_OR_expansion/queries/${params.species}.clear_OR.prot.fa",
        checkIfExists: true,
    )
    def target_fasta = file(
        "${params.outdir}/05_multiexon_OR/candidate_sequences/${params.species}.multiexon_regions.fa",
        checkIfExists: true,
    )
    def regions_tsv = file(
        "${params.outdir}/05_multiexon_OR/candidate_regions/${params.species}.multiexon_regions.tsv",
        checkIfExists: true,
    )
    def merged_hits = file(
        "${params.outdir}/04_species_OR_expansion/merged/${params.species}.expanded.merged.tblastn.tsv",
        checkIfExists: true,
    )
    def package_src = file("${launchDir}/src/or_mining", checkIfExists: true)
    def timeout_seconds = (params.exonerate_soft_timeout_minutes as int) * 60
    def smoke_timeout_seconds = (params.exonerate_smoke_timeout_minutes as int) * 60
    def pilot_count = params.exonerate_pilot_pairs as int

    PREPARE_EXONERATE_FUNCTIONAL_SMOKE(
        channel.value(meta),
        channel.value(package_src),
    )
    EXONERATE_FUNCTIONAL_SMOKE(
        PREPARE_EXONERATE_FUNCTIONAL_SMOKE.out.inputs,
        channel.value(params.exonerate_smoke_exhaustive),
        channel.value(false),
        channel.value(params.exonerate_min_intron),
        channel.value(params.maximum_intron_length),
        channel.value(smoke_timeout_seconds),
        channel.value(true),
        channel.value(package_src),
    )
    VALIDATE_EXONERATE_FUNCTIONAL_SMOKE(
        EXONERATE_FUNCTIONAL_SMOKE.out.artifacts,
        channel.value(package_src),
    )
    CHECK_EXONERATE_SMOKE_GATE(
        VALIDATE_EXONERATE_FUNCTIONAL_SMOKE.out.summary,
        channel.value(package_src),
    )

    PREPARE_TARGETED_EXONERATE(
        channel.of([meta, reference_or, species_or, target_fasta, regions_tsv, merged_hits]),
        channel.value(package_src),
        channel.value(params.exonerate_top_queries_per_region),
    )
    def targeted_pairs = PREPARE_TARGETED_EXONERATE.out.pairs.flatMap {
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
                error "Missing target FASTA for ${pair_id}"
            }
            tuple(pair_meta + [pair_id: pair_id, pilot: index < pilot_count], query_file, target_file)
        }
    }
    def pair_branches = targeted_pairs.branch { pair_meta, query_file, target_file ->
        pilot: pair_meta.pilot
        production: true
    }
    def gated_pilot_pairs = pair_branches.pilot
        .combine(CHECK_EXONERATE_SMOKE_GATE.out.gate)
        .map { pair_meta, query_file, target_file, smoke_gate ->
            tuple(pair_meta, query_file, target_file)
        }
    EXONERATE_TARGETED_PILOT(
        gated_pilot_pairs,
        channel.value(params.exonerate_exhaustive),
        channel.value(params.exonerate_suboptimal),
        channel.value(params.exonerate_min_intron),
        channel.value(params.maximum_intron_length),
        channel.value(timeout_seconds),
        channel.value(true),
        channel.value(package_src),
    )

    def collected_pilot = EXONERATE_TARGETED_PILOT.out.artifacts
        .map { pair_meta, result, status ->
            tuple(pair_meta.id, [id: pair_meta.id], result, status)
        }
        .groupTuple(by: 0)
        .map { species_id, metas, results, statuses ->
            tuple(metas[0], results, statuses)
        }
    CHECK_EXONERATE_PILOT(
        collected_pilot,
        channel.value(package_src),
        channel.value(pilot_count),
        channel.value(params.exonerate_pilot_min_positive_pairs),
    )

    def gated_production_pairs = pair_branches.production
        .combine(CHECK_EXONERATE_PILOT.out.gate)
        .map { pair_meta, query_file, target_file, pilot_gate ->
            tuple(pair_meta, query_file, target_file)
        }
    EXONERATE_TARGETED_SEARCH(
        gated_production_pairs,
        channel.value(params.exonerate_exhaustive),
        channel.value(params.exonerate_suboptimal),
        channel.value(params.exonerate_min_intron),
        channel.value(params.maximum_intron_length),
        channel.value(timeout_seconds),
        channel.value(true),
        channel.value(package_src),
    )

    def all_artifacts = EXONERATE_TARGETED_PILOT.out.artifacts
        .mix(EXONERATE_TARGETED_SEARCH.out.artifacts)
    def collected_results = all_artifacts
        .map { result_meta, result, status ->
            tuple(result_meta.id, [id: result_meta.id], result, status)
        }
        .groupTuple(by: 0)
        .map { species_id, metas, results, statuses ->
            tuple(metas[0], results, statuses)
        }
    MERGE_EXONERATE_RESULTS(
        collected_results,
        channel.value(package_src),
    )
}
