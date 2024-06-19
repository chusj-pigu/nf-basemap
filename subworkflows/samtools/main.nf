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

process sam_sort {
    publishDir "${params.out_dir}/alignments", mode : "copy"
    label "sam_big"

    input: 
    path pass_bam

    output:
    tuple path("${pass_bam.baseName}.bam"), path("${pass_bam.baseName}.bam.bai")
    
    script:
    """
    samtools sort -@ $params.threads --write-index $pass_bam -o ${pass_bam.baseName}.bam##idx##${pass_bam.baseName}.bam.bai
    """
}

