process multiqc {
    publishDir "${params.out_dir}", mode : "copy"
    
    input:
    path "*"

    output:
    path "multiqc_report.html"

    script:
    def datadir = params.git ? "--no-data-dir" : ""
    """
    multiqc $datadir -n stdout . > multiqc_report.html
    """

}
