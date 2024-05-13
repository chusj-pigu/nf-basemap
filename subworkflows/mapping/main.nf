process mapping {
    publishDir "${params.out_dir}", mode : "copy"
    label "minimap"

    input:
    path ref
    path fastq

    output:
    path "${params.sample_id}_aligned.sam"

    script:
    """
    minimap2 -y -ax map-ont -t $params.threads $ref $fastq > "${params.sample_id}_aligned.sam"
    """
}