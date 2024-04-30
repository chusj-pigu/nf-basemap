// Usage help

def helpMessage() {
  log.info """
        Usage:
        The typical command for running the pipeline is as follows:
        nextflow run chusj-pigu/wf-mapping --reads SAMPLE_ID.fq.gz --ref REF.fasta

        Mandatory arguments:
         --pod5                         Path to the directory containing pod5 files
         --ref                          Path to the reference fasta file 

         Optional arguments:
         --outdir                       Output directory to place mapped files and reports in [default: output]
         --sample_id                    Will name output files according to sample id [default: reads]
         --m_bases                      Modified bases to be called, separated by commas if more than one is desired. Requires path to model if run with drac profile [default: 5mCG_5hmCG].
         --model                        Basecalling model to use [default: sup@v4.3.0].
         --no_mod                       Basecalling without base modification [default: false]
         --model_path                   Path for the basecalling model, required when running with drac profile [default: path to sup@v4.3.0]
         --m_bases_path                 Path for the modified basecalling model, required when running with drac profile [default: path to sup@v4.3.0_5mCG_5hmCG]
         -profile                       Use standard for running locally, or drac when running on Digital Research Alliance of Canada Narval [default: standard]
         --help                         This usage statement.
        """
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

include { basecall } from './subworkflows/dorado'
include { pod5_channel } from './subworkflows/pod5'
include { pod5_subset } from './subworkflows/pod5'
include { ubam_to_fastq } from './subworkflows/ubam_fastq'
include { qc_fastq} from './subworkflows/fastp'
include { mapping } from './subworkflows/mapping'
include { sam_to_bam } from './subworkflows/samtools'
include { sam_sort } from './subworkflows/samtools'
include { sam_index } from './subworkflows/samtools'
include { sam_stats } from './subworkflows/samtools'
include { multiqc } from './subworkflows/multiqc'

workflow {

    pod5_ch = Channel.fromPath(params.pod5)
    model_ch = params.model ? Channel.of(params.model) : Channel.fromPath(params.model)
    ref_ch = Channel.fromPath(params.ref)
    
    pod5_channel(pod5_ch)
    pod5_subset(pod5_ch,pod5_channel.out)
    
    basecall(pod5_subset.out, model_ch)
    
    ubam_to_fastq(basecall.out)
    qc_fastq(ubam_to_fastq.out)
    
    mapping(ref_ch, ubam_to_fastq.out)
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