process ubam_to_fastq {
    publishDir "${params.out_dir}/reads", mode : "copy"
    label "sam_sm"

    input:
    path ubam

    output:
    path "${ubam.baseName}.fq.gz"

    script:
    def mod = params.no_mod ? "" : "-T '*'" 
    """
    samtools fastq $mod -@ $params.threads $ubam | gzip > "${ubam.baseName}.fq.gz" 
    """
}