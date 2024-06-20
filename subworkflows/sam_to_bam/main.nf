process sam_to_bam {
    publishDir "${params.out_dir}/alignments", mode : "copy"
    label "sam_big"

    input:
    path sam

    output:
    path "${sam.baseName}.bam"

    script:
    """
    samtools view --no-PG -@ $params.threads -Sb $sam -o ${sam.baseName}.bam
    """
}