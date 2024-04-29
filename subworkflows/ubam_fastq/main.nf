process ubam_to_fastq {
    publishDir "${params.outdir}"
    cpus params.t
    label "samtools"

    input:
    path ubam

    output:
    path "${ubam.baseName}.fq.gz"

    script:
    def m = params.m ? "-T '*'" : "" 
    """
    samtools fastq $m -@ $task.cpus $ubam > "${ubam.baseName}.fq.gz" 
    """
}