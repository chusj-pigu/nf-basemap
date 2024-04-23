// Usage help

def helpMessage() {
  log.info """
        Usage:
        The typical command for running the pipeline is as follows:
        nextflow run chusj-pigu/wf-mapping --reads SAMPLE_ID.fq.gz --ref REF.fasta

        Mandatory arguments:
         --reads                        Path to the input file in FASTQ format or uBAM if used with --bam
         --ref                          Path to the reference fasta file 

         Optional arguments:
         --out_dir                      Output directory to place mapped files and reports in [default: output]
         --t                            Number of CPUs to use [default: 4]
         --m                            Conserve base modification information
         --bam                          Take an uBAM file as input
         -profiles                      Use Docker, Singularity or Apptainer to run the workflow [default: Docker]
         --help                         This usage statement.
        """
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

process ubam_to_fastq {
    publishDir="${launchDir}/output/"
    cpus params.t
    label "samtools"

    input:
    path ubam

    output:
    path "${ubam.baseName}.fq.gz"

    script:
    def m = params.m ? "-T '*'" : "" 
    """
    samtools fastq $m -@ $task.cpus $ubam > "${ubam.baseName}.fq.gz" 
    """
}

process mapping {
    cpus params.t
    publishDir "${params.out_dir}"
    label "minimap"

    input:
    path ref
    path fastq

    output:
    path "${fastq.getBaseName(2)}.sam"

    script:
    def m = params.m ? "-y" : ""
    """
    minimap2 -ax map-ont $m -t $task.cpus $ref $fastq > "${fastq.getBaseName(2)}.sam"
    """
}

process sam_to_bam {
    cpus params.t
    publishDir "${params.out_dir}"
    label "samtools"

    input:
    path sam

    output:
    path "${sam.baseName}.bam"

    script:
    """
    samtools view -@ $task.cpus -Sb $sam > ${sam.baseName}.bam
    """
}

process sam_sort {
    cpus params.t
    publishDir "${params.out_dir}"
    label "samtools"

    input:
    path aligned

    output:
    path "${aligned.baseName}_sorted.bam"

    script:
    """
    samtools sort -@ $task.cpus $aligned > ${aligned.baseName}_sorted.bam
    """
}

process sam_index {
    cpus params.t
    publishDir "${params.out_dir}"
    label "samtools"

    input:
    path sorted

    output:
    path "${sorted.baseName}_index.bam"

    script:
    """
    samtools index -@ $task.cpus $sorted > ${sorted.baseName}_index.bam
    """
}

process sam_stats {
    cpus params.t
    publishDir "${params.out_dir}"
    label "samtools"

    input:
    path sorted

    output:
    path "${sorted.baseName}.stats.txt"

    script:
    """
    samtools stats -@ $task.cpus $sorted > ${sorted.baseName}.stats.txt
    """
}

process multiqc {
    publishDir "${params.out_dir}" 
    
    input:
    path "${params.out_dir}/*"

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc .
    """

}


workflow {
    ref_ch = Channel.fromPath(params.ref)
    reads_ch = Channel.fromPath(params.reads) 
    if ( params. bam ) {
        ubam_to_fastq(reads_ch)
        mapping(ref_ch, ubam_to_fastq.out)
    }
    else {
       mapping(ref_ch,reads_ch) 
    }
    sam_to_bam(mapping.out)
    sam_sort(sam_to_bam.out)
    sam_index(sam_sort.out)
    sam_stats(sam_sort.out)
    multi_ch = Channel.empty()
        .mix(sam_stats.out)
        .collect()
        .set { stat_files }
    multiqc(stat_files)
}