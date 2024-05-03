process multiqc {
    publishDir "${params.out_dir}", mode : "copy"
    
    input:
    path "*"

    output:
    path "multiqc_report.html"
    path "multiqc_data"

    script:
    """
    multiqc . 
    """

}