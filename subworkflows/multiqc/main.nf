process multiqc {
    publishDir "${params.out_dir}", mode : "copy"
    
    input:
    path multiqc_config
    path "*"

    output:
    stdout

    script:
    """
    multiqc --config $multiqc_config .
    """

}
