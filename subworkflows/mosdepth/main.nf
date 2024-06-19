process mosdepth {
    publishDir "${params.out_dir}/reports", mode : "copy"

    input:
    tuple path(bam), path(bai)

    output:
    path "${bam.baseName}.mosdepth.global.dist.txt"
    path "${bam.baseName}.mosdepth.summary.txt"

    script:
    """
    mosdepth -n '${bam.baseName}' $bam
    """
}