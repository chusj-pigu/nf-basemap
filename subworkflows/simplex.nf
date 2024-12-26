include { basecall } from '../modules/dorado'
include { demultiplex } from '../modules/dorado'
include { qs_filter } from '../modules/samtools'
include { ubam_to_fastq as ubam_to_fastq_p } from '../modules/samtools'
include { ubam_to_fastq as ubam_to_fastq_f } from '../modules/samtools'
include { nanoplot } from "../modules/nanoplot"
include { ALIGNMENT } from '../subworkflows/mapping'
include { multiqc } from '../modules/multiqc'

workflow SIMPLEX {
    take:
    sample_sheet
    model
    bed
    ref
    
    main:
    input_ch = sample_sheet
        .combine(model)

    basecall(input_ch)

    if (params.demux != null) {
        demultiplex(basecall.out)
        demultiplex_out = demultiplex.out  // Capture the output channel

        split_bams = demultiplex_out.flatMap { sample_id, bam_files ->
        bam_files.collect { bam -> tuple(sample_id, bam) }
}
        qs_filter(split_bams)

        ubam_to_fastq_p(qs_filter.out.ubam_pass)
        ubam_to_fastq_f(qs_filter.out.ubam_fail)


    } else {
        qs_filter(basecall.out)
        nanoplot(basecall.out)

        ubam_to_fastq_p(qs_filter.out.ubam_pass)
        ubam_to_fastq_f(qs_filter.out.ubam_fail)
    }

    
    if (params.skip_mapping) {
        multi_ch = Channel.empty()
            .mix(nanoplot.out)
            .collect()
        multiqc(multi_ch)

    } else {
        ALIGNMENT(ubam_to_fastq_p.out, bed, ref)

        multi_ch = Channel.empty()
            .mix(nanoplot.out,ALIGNMENT.out.mosdepth_dist,ALIGNMENT.out.mosdepth_summ)
            .collect()
        multiqc(multi_ch)
    }

}