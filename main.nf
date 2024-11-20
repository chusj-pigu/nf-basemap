// Usage help

def helpMessage() {
  log.info """
        Usage:
        The typical command for running the pipeline is as follows:
        nextflow run chusj-pigu/nf-basecall --pod5 /path/to/pod5 --ref /path/to/REF.fasta

        Mandatory arguments:
         --pod5 or --fastq              Path to the directory containing pod5 files or fastq files

         Optional arguments:
         --ref                          Path to the reference fasta file 
         --skip_basecall                Basecalling step will be skipped, input must me in fastq [default: false]
         --skip_mapping                 Mapping will be skipped [default: false]
         --duplex                       Dorado will basecall in duplex mode instead of simplex [default: false]
         --out_dir                      Output directory to place mapped files and reports in [default: output]
         --sample_id                    Will name output files according to sample id [default: reads]
         --m_bases                      Modified bases to be called, separated by commas if more than one is desired. Requires path to model if run with drac profile [default: 5mCG_5hmCG].
         --model                        Basecalling model to use [default: sup@v5.0.0].
         --no_mod                       Basecalling without base modification [default: false]
         --model_path                   Path for the basecalling model, required when running with drac profile [default: path to sup@v5.0.0]
         --m_bases_path                 Path for the modified basecalling model, required when running with drac profile [default: path to sup@v5.0.0_5mCG_5hmCG]
         -profile                       Use standard for running locally, or drac when running on Digital Research Alliance of Canada Narval [default: standard]
         --bed                          Bed file containing regions of interest to compute mapping statistics with mosdepth
         --resume                       Use when basecalling crashes and you want to resume from the existing ubam [default: false]
         --ubam                         Partial ubam that basecalling will resume from [default: NO_FILE]
         --batch                        Batchsize for basecalling, if 0 optimal batchsize will be automatically selected [default: 0]
         --help                         This usage statement.
        """
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

include { SIMPLEX } from './subworkflows/simplex'
include { DUPLEX } from './subworkflows/duplex'
include { ALIGNMENT } from './subworkflows/mapping'
include { multiqc } from './modules/multiqc'

workflow {

    // Create channel for partial ubam to make basecalling resuming possible:
    partial_ubam = Channel.fromPath(params.ubam)

    if (params.skip_basecall) {
        ref_ch = Channel.fromPath(params.ref)
        fastq_ch = Channel.fromPath(params.fastq)
        
        ALIGNMENT(fastq_ch, ref_ch)

        multi_ch = Channel.empty()
            .mix(ALIGNMENT.out.mosdepth_dist, ALIGNMENT.out.mosdepth_summary, ALIGNMENT.out.mosdepth_bed)
            .collect()
        multiqc(multi_ch)
    }

    else if (params.duplex) {
        pod5_ch = Channel.fromPath(params.pod5)
        model_ch = params.model ? Channel.of(params.model) : Channel.fromPath(params.model_path)
        ref_ch = Channel.fromPath(params.ref)
        DUPLEX(pod5_ch, model_ch, partial_ubam, ref_ch)

    } else {
        pod5_ch = Channel.fromPath(params.pod5)
        model_ch = params.model ? Channel.of(params.model) : Channel.fromPath(params.model_path)
        ref_ch = Channel.fromPath(params.ref)
        SIMPLEX(pod5_ch, model_ch, partial_ubam, ref_ch)
    }
}
