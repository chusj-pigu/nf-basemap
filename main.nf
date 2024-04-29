process sam_to_bam {
    cpus params.t
    publishDir "${params.out_dir}"
    label "samtools"

    input:
    path sam

    output:
    path "${sam.baseName}.bam"

    script:
    """
    samtools view -@ $task.cpus -e '[qs] >=10' -Sb $sam > ${sam.baseName}.bam
    """
}

process sam_sort {
    cpus params.t
    publishDir "${params.out_dir}"
    label "samtools"

    input:
    path aligned

    output:
    path "${aligned.baseName}_sorted.bam"

    script:
    """
    samtools sort -@ $task.cpus $aligned > ${aligned.baseName}_sorted.bam
    """
}

process sam_index {
    cpus params.t
    publishDir "${params.out_dir}"
    label "samtools"

    input:
    path sorted

    output:
    path "${sorted.baseName}_index.bam"

    script:
    """
    samtools index -@ $task.cpus $sorted > ${sorted.baseName}_index.bam
    """
}

process sam_stats {
    cpus params.t
    publishDir "${params.out_dir}"
    label "samtools"

    input:
    path sorted

    output:
    path "${sorted.baseName}.stats.txt"

    script:
    """
    samtools stats -@ $task.cpus $sorted > ${sorted.baseName}.stats.txt
    """
}