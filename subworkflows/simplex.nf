include { basecall } from '../modules/dorado'
include { qs_filter } from '../modules/samtools'
include { ubam_to_fastq as ubam_to_fastq_p } from '../modules/samtools'
include { ubam_to_fastq as ubam_to_fastq_f } from '../modules/samtools'
include { nanoplot } from "../modules/nanoplot"
include { ALIGNMENT } from '../subworkflows/mapping'
include { multiqc } from '../modules/multiqc'

workflow SIMPLEX {
    take:
    pod5
    model
    ubam
    ref
    
    main:
    basecall(pod5, model, ubam)

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
        ALIGNMENT(ubam_to_fastq_p.out, ref)

        multi_ch = Channel.empty()
            .mix(nanoplot.out,ALIGNMENT.out.mosdepth_all_out)
            .collect()
        multiqc(multi_ch)
    }

}