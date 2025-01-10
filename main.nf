// Usage help

def helpMessage() {
  log.info """
        Usage:
        The typical command for running the pipeline is as follows:
        nextflow run chusj-pigu/nf-basecall --pod5 /path/to/pod5 --ref /path/to/REF.fasta

        Mandatory arguments:
         --sample_sheet             Path to the samplesheet with first column containing sample_id and second column containing the absolute path to fastq or pod5 files (MANDATORY)

         Optional arguments:
         --ref                          Path to the reference fasta file 
         --skip_basecall                Basecalling step will be skipped, input must me in fastq [default: false]
         --skip_mapping                 Mapping will be skipped [default: false]
         --duplex                       Dorado will basecall in duplex mode instead of simplex [default: false]
         --out_dir                      Output directory to place mapped files and reports in [default: output]
         --m_bases                      Modified bases to be called, separated by commas if more than one is desired. Requires path to model if run with drac profile [default: 5mCG_5hmCG].
         --model                        Basecalling model to use [default: sup@v5.0.0].
         --device                       Parameter to choose device when basecalling. Specify CPU or GPU device: 'auto', 'cpu', 'cuda:all' or 'cuda:<device_id>[,<device_id>...]'. Specifying 'auto' will choose either 'cpu', 'metal' or 'cuda:all' depending on the presence of a GPU device. [nargs=0..1] [default: "auto"]
         --no_mod                       Basecalling without base modification [default: false]
         --model_path                   Path for the basecalling model, required when running with drac profile [default: path to sup@v5.0.0]
         --m_bases_path                 Path for the modified basecalling model, required when running with drac profile [default: path to sup@v5.0.0_5mCG_5hmCG]
         -profile                       Use standard for running locally, or drac when running on Digital Research Alliance of Canada Narval [default: standard]
         --bed                          Bed file containing regions of interest to compute mapping statistics with mosdepth
         --resume                       Use when basecalling crashes and you want to resume from the existing ubam [default: false]
         --ubam                         Path to the samplesheet with first column containing sample_id and second column containing the absolute path to ubam file to append to [default: none].
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

    // Create channel for bed file 
    bed_ch = Channel.fromPath(params.bed)
    ubam_ch = Channel.fromPath("${projectDir}/assets/NO_UBAM")

    // Create channel for partial ubam to make basecalling resuming possible:
    if (params.ubam != null) {
        partial_ubam = Channel.fromPath(params.ubam)
            .splitCsv(header: true)
            .map { row -> tuple(row.sample_id, file(row.ubam)) }
    
        sheet_ch = Channel.fromPath(params.sample_sheet)
            .splitCsv(header: true)
            .map { row -> tuple(row.sample_id, file(row.path)) }
            .join(partial_ubam)
    } else {
        sheet_ch = Channel.fromPath(params.sample_sheet)
            .splitCsv(header: true)
            .map { row -> tuple(row.sample_id, file(row.path)) }
            .combine(ubam_ch)
    }

    if (params.skip_basecall) {
        ref_ch = Channel.fromPath(params.ref)
        
        ALIGNMENT(sheet_ch, bed_ch, ref_ch)

        multi_ch = Channel.empty()
            .mix(ALIGNMENT.out.mosdepth_all_out)
            .collect()
        multiqc(multi_ch)
    }

    else if (params.duplex) {
        model_ch = params.model ? Channel.of(params.model) : Channel.fromPath(params.model_path)
        ref_ch = Channel.fromPath(params.ref)
        DUPLEX(sheet_ch, model_ch, bed_ch, ref_ch)

    } else {
        model_ch = params.model ? Channel.of(params.model) : Channel.fromPath(params.model_path)
        ref_ch = Channel.fromPath(params.ref)
        SIMPLEX(sheet_ch, model_ch, bed_ch, ref_ch)
    }
}
