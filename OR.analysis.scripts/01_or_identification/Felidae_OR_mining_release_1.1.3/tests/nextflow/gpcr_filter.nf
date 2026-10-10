nextflow.enable.dsl = 2

include { TRANSLATE_ORFS } from '../../modules/local/translate_orfs/main'
include { BLASTP_GPCR }    from '../../modules/local/blastp_gpcr/main'
include { CLASSIFY_GPCR }  from '../../modules/local/classify_gpcr/main'

workflow {
    if (!params.gpcr_db) {
        error "Missing required parameter: gpcr_db"
    }
    if (!params.species) {
        error "Missing required parameter: species"
    }

    def meta = [id: params.species]
    def orf_dir = file("${params.outdir}/03_OR_Finder_results/orfs", checkIfExists: true)
    def nucleotide_orfs = file("${orf_dir}/${params.species}.orfs.nt.fa", checkIfExists: true)
    def orf_table = file("${orf_dir}/${params.species}.orfs.tsv", checkIfExists: true)
    def gpcr_db = file(params.gpcr_db, checkIfExists: true)
    def package_src = file("${launchDir}/src/or_mining", checkIfExists: true)

    TRANSLATE_ORFS(
        channel.of([meta, nucleotide_orfs]),
        channel.value(package_src),
    )
    BLASTP_GPCR(
        TRANSLATE_ORFS.out.proteins,
        channel.value(gpcr_db),
        channel.value(params.gpcr_blast_evalue),
    )
    CLASSIFY_GPCR(
        BLASTP_GPCR.out.results,
        channel.value(nucleotide_orfs),
        channel.value(orf_table),
        channel.value(package_src),
    )
}

