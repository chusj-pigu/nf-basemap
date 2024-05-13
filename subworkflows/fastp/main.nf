process cat_fastq { 

    input:
    path '*.fq.gz'

    output:
    path "all.fq.gz"

    script:
    """
    cat *.fq.gz > all.fq.gz
    """
}

process qc_fastq {
    label "fastp" 
    publishDir "${params.out_dir}", mode: "copy"

    input:
    path fq

    output:
    path "${fq.baseName}_log.json"

    script:
    """
    fastp -e 10 -i $fq -j ${fq.baseName}_log.json
    """
}