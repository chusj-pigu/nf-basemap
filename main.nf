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
         --b                            Batchsize for basecalling, if 0 optimal batchsize will be automatically selected [default: 0]
         --help                         This usage statement.
        """
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

include { basecall } from './modules/dorado'
include { pod5_channel } from './modules/pod5'
include { subset } from './modules/pod5'
include { qs_filter } from './modules/samtools'
include { ubam_to_fastq as ubam_to_fastq_p } from './modules/samtools'
include { ubam_to_fastq as ubam_to_fastq_f } from './modules/samtools'
include { nanoplot } from "./modules/nanoplot"
include { mapping } from './modules/minimap'
include { sam_sort } from './modules/samtools'
include { mosdepth } from './modules/mosdepth'
include { multiqc } from './modules/multiqc'

workflow {
    if (params.duplex) {
        pod5_in = Channel.fromPath(params.pod5)
        model_ch = params.model ? Channel.of(params.model) : Channel.fromPath(params.model_path)
        pod5_channel(pod5_in)
        pod5_ch = subset(pod5_in,pod5_channel.out)
    } else {
        pod5_ch = Channel.fromPath(params.pod5)
        model_ch = params.model ? Channel.of(params.model) : Channel.fromPath(params.model_path)
    }

    if (params.skip_basecall) {
        ref_ch = Channel.fromPath(params.ref)
        fastq_ch = Channel.fromPath(params.fastq)

        mapping(ref_ch, fastq_ch)

        sam_sort(mapping.out)
        mosdepth(sam_sort.out)

        multi_ch = Channel.empty()
            .mix(mosdepth.out)
            .collect()
        multiqc(multi_ch)
    }

    else if (params.skip_mapping) {
        
        basecall(pod5_ch, model_ch)

        qs_filter(basecall.out)
        nanoplot(basecall.out)

        fq_pass = ubam_to_fastq_p(qs_filter.out.ubam_pass)
        fq_fail = ubam_to_fastq_f(qs_filter.out.ubam_fail)

        multi_ch = Channel.empty()
            .mix(nanoplot.out)
            .collect()
        multiqc(multi_ch)

    } else {
        ref_ch = Channel.fromPath(params.ref)
        basecall(pod5_ch, model_ch)

        qs_filter(basecall.out)
        nanoplot(basecall.out)

        fq_pass = ubam_to_fastq_p(qs_filter.out.ubam_pass)
        fq_fail = ubam_to_fastq_f(qs_filter.out.ubam_fail)

        mapping(ref_ch, fq_pass)

        sam_sort(mapping.out)
        mosdepth(sam_sort.out)

        multi_ch = Channel.empty()
            .mix(nanoplot.out,mosdepth.out)
            .collect()
        multiqc(multi_ch)
    }
}
