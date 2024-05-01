process multiqc {
    publishDir "${params.out_dir}", mode : "copy"
    
    input:
    path "${params.out_dir}/*"

    output:
    path "${params.out_dir}/multiqc_report.html"

    script:
    """
    multiqc . 
    """

}