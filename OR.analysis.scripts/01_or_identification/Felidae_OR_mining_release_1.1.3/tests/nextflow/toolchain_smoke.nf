nextflow.enable.dsl = 2

process OR_FINDER_TOOLCHAIN_SMOKE {
    tag 'blast-emboss'
    label 'process_single'
    conda "${launchDir}/envs/or_finder.yml"
    publishDir "${params.outdir}/smoke", mode: params.publish_mode, overwrite: true

    output:
    path 'or_finder_toolchain.txt'

    script:
    """
    {
        echo "status=ok"
        echo "hostname=\$(hostname)"
        echo "job_id=\${SLURM_JOB_ID:-unset}"
        echo "requested_cpus=${task.cpus}"
        tblastn -version | head -n 1
        getorf -version
        transeq -version
    } > or_finder_toolchain.txt 2>&1
    """
}

workflow {
    OR_FINDER_TOOLCHAIN_SMOKE()
}

