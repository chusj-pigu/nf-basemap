process qs_filter {
    publishDir "${params.out_dir}/reads", mode : "copy"
    label "sam_sm"

    input:
    path ubam

    output:
    path "${ubam.baseName}_unaligned_pass.bam", emit: ubam_pass
    path "${ubam.baseName}_unaligned_fail.bam", emit: ubam_fail

    script:
    """
    samtools view --no-PG -@ $params.threads -e '[qs] >=10' -b $ubam --output ${ubam.baseName}_unaligned_pass.bam --unoutput ${ubam.baseName}_unaligned_fail.bam
    """
}