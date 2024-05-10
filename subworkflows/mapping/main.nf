process mapping {
    publishDir "${params.out_dir}", mode : "copy"
    label "minimap"

    input:
    path ref
    path fastq

    output:
    path "${fastq.getBaseName(2)}.sam"

    script:
    def mod = params.no_mod ? "" : "-y"
    """
    minimap2 -ax map-ont $mod -t $params.threads $ref $fastq > "${fastq.getBaseName(2)}.sam"
    """
}