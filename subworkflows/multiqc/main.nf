process multiqc {
    publishDir "${params.out_dir}", mode : "copy"
    
    input:
    path "*"

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc .
    """

}
