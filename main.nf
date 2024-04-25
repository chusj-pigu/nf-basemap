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

include { ubam_to_fastq } from './subworkflows/ubam_fastq'
include { qc_fastq} from './subworkflows/fastp'
include { mapping } from './subworkflows/mapping'
include { sam_to_bam } from './subworkflows/samtools'
include { sam_sort } from './subworkflows/samtools'
include { sam_index } from './subworkflows/samtools'
include { sam_stats } from './subworkflows/samtools'
include { multiqc } from './subworkflows/multiqc'

workflow {
    ref_ch = Channel.fromPath(params.ref)
    reads_ch = Channel.fromPath(params.reads) 
    if ( params. bam ) {
        ubam_to_fastq(reads_ch)
        qc_fastq(ubam_to_fastq.out)
        mapping(ref_ch, ubam_to_fastq.out)
    }
    else {
        qc_fastq(reads_ch)
        mapping(ref_ch,reads_ch)
    }
    sam_to_bam(mapping.out)
        sam_sort(sam_to_bam.out)
        sam_index(sam_sort.out)
        sam_stats(sam_sort.out)
        multi_ch = Channel.empty()
            .mix(qc_fastq.out)
            .mix(sam_stats.out)
            .collect()
        multiqc(multi_ch)
}