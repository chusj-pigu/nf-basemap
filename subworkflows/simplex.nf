include { basecall } from '../modules/dorado'
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
    basecall(sample_sheet, model)

    qs_filter(basecall.out)

    nanoplot(basecall.out)

    ubam_to_fastq_p(qs_filter.out.ubam_pass)
    ubam_to_fastq_f(qs_filter.out.ubam_fail)
    
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