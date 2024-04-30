process ubam_to_fastq {
    publishDir "${params.out_dir}", mode : "copy"
    label "samtools"

    input:
    path ubam

    output:
    path "${ubam.baseName}.fq.gz"

    script:
    def mod = params.no_mod ? "" : "-T '*'" 
    """
    samtools fastq $mod -@ $task.cpus $ubam > "${ubam.baseName}.fq.gz" 
    """
}