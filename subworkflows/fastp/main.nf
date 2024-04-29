process qc_fastq {
    label "fastp" 
    publishDir "${params.outdir}", mode: "copy"

    input:
    path fastq

    output:
    path "${fastq.baseName}_log.json"

    script:
    """
    fastp -Q -i $fastq -j ${fastq.baseName}_log.json
    """
}