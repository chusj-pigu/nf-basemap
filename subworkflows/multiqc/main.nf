process multiqc {
    publishDir "${params.outdir}" 
    
    input:
    path "${params.outdir}/*"

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc . 
    """

}