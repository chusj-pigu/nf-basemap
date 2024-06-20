process sam_sort {
    publishDir "${params.out_dir}/alignments", mode : "copy"
    label "sam_big"

    input: 
    path sam

    output:
    tuple path("${sam.baseName}.bam"), path("${sam.baseName}.bam.bai")
    
    script:
    """
    samtools view --no-PG -@ $params.threads -Sb $sam \
    | samtools sort -@ $params.threads --write-index -o ${sam.baseName}.bam##idx##${sam.baseName}.bam.bai
    """
}

