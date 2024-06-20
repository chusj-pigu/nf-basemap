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

