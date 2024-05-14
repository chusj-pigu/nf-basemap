process multiqc {
    publishDir "${params.out_dir}", mode : "copy"
    
    input:
    path multiqc_config
    path "*"

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc --config $multiqc_config .
    """

}
