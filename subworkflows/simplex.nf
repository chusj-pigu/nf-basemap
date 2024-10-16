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
    ref
    
    main:
    basecall(pod5, model)

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
            .mix(nanoplot.out,ALIGNMENT.out.mosdepth_dist, ALIGNMENT.out.mosdepth_summary, ALIGNMENT.out.mosdepth_bed)
            .collect()
        multiqc(multi_ch)
    }

    emit:
    fq_pass = ubam_to_fastq_p.out
    fq_fail = ubam_to_fastq_f.out
    bam = ALIGNMENT.out.bam

}