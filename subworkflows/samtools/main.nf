process sam_to_bam {
    publishDir "${params.out_dir}", mode : "copy"
    label "samtools"

    input:
    path sam

    output:
    path "${sam.baseName}_all.bam"

    script:
    """
    samtools view -@ $params.threads -Sb $sam -o ${sam.baseName}_all.bam
    """
}

process sam_qs_filter {
    publishDir "${params.out_dir}", mode : "copy"
    label "samtools"

    input:
    path bam

    output:
    path "${bam.baseName}_pass.bam"

    script:
    """
    samtools view --no-PG -@ $params.threads -e '[qs] >=10' -b $bam -o ${bam.baseName}_pass.bam
    """
}

process sam_sort {
    publishDir "${params.out_dir}", mode : "copy"
    label "samtools"

    input:
    path aligned

    output:
    path "${aligned.baseName}_sorted.bam"

    script:
    """
    samtools sort -@ $params.threads $aligned -o ${aligned.baseName}_sorted.bam
    """
}

process sam_index {
    publishDir "${params.out_dir}", mode : "copy"
    label "samtools"

    input:
    path sorted

    output:
    path "${sorted.baseName}_index.bam.bai"

    script:
    """
    samtools index -@ $params.threads $sorted -o ${sorted.baseName}_index.bam.bai
    """
}

process sam_stats {
    publishDir "${params.out_dir}", mode : "copy"
    label "samtools"

    input:
    path sorted

    output:
    path "${sorted.baseName}.stats.txt"

    script:
    """
    samtools stats -@ $params.threads $sorted > ${sorted.baseName}.stats.txt
    """
}