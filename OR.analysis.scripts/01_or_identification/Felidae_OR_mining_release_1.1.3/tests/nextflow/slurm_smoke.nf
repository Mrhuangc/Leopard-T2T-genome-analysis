nextflow.enable.dsl = 2

process SLURM_SMOKE {
    label 'process_single'
    publishDir "${params.outdir}/smoke", mode: params.publish_mode, overwrite: true

    output:
    path 'slurm_smoke.txt'

    script:
    """
    {
        echo "status=ok"
        echo "hostname=\$(hostname)"
        echo "date=\$(date --iso-8601=seconds)"
        echo "job_id=\${SLURM_JOB_ID:-unset}"
        echo "cpus=\${SLURM_CPUS_PER_TASK:-unset}"
        python --version 2>&1
    } > slurm_smoke.txt
    """
}

workflow {
    SLURM_SMOKE()
}
