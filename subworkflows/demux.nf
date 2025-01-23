include { basecall } from '../modules/dorado'
include { demultiplex } from '../modules/dorado'
include { qs_filter } from '../modules/samtools'
include { ubam_to_fastq as ubam_to_fastq_p } from '../modules/samtools'
include { ubam_to_fastq as ubam_to_fastq_f } from '../modules/samtools'
include { nanoplot } from "../modules/nanoplot"
include { ALIGNMENT } from '../subworkflows/mapping'
include { multiqc } from '../modules/multiqc'

workflow DEMULTIPLEX {
    take:
    sample_sheet
    sample_demux
    model
    bed
    ref
    
    main:
    input_ch = sample_sheet
        .combine(model)

    basecall(input_ch)

    demultiplex(basecall.out)
    demultiplex_out = demultiplex.out  // Capture the output channel

    split_bams = demultiplex_out.flatMap { sample_id, bam_files ->
    bam_files.collect { bam ->
        // Use the baseName of the BAM file for the output
        def bam_base = bam.baseName
        tuple(sample_id,bam_base,bam)  // Create a tuple with the baseName as the key and the BAM file
        }
    }
    qs_filter(split_bams)

    // We remove the higher level sample_id used for basecalling
    red_qsfilt_pass = qs_filter.out.ubam_pass
        .map { sample_id, barcode, ubam -> 
            tuple(barcode, ubam) }

    red_qsfilt_fail = qs_filter.out.ubam_fail
        .map { sample_id, barcode, ubam -> 
            tuple(barcode, ubam) }
    
    // Get the sample_ids from the sample_sheet that specify each barcode to each sample_id
    new_sample_ids_pass = sample_demux 
        .combine(red_qsfilt_pass, by:0)
        .map { barcode, sample_id, ubam ->
            tuple(sample_id, barcode, ubam) }

    new_sample_ids_fail = sample_demux 
        .combine(red_qsfilt_fail, by:0)
        .map { barcode, sample_id, ubam ->
            tuple(sample_id, barcode, ubam) }

    ubam_to_fastq_p(new_sample_ids_pass)
    ubam_to_fastq_f(new_sample_ids_fail)

    nanoplot(ubam_to_fastq_p.out)

    if (params.skip_mapping) {
            // Skip alignment and run multiqc only on nanoplot output
        multi_ch = nanoplot.out
            .collect()

        multiqc(multi_ch)
    } else {
            // Perform alignment if skip_mapping is not set
        ALIGNMENT(ubam_to_fastq_p.out, bed, ref)

        multi_ch = nanoplot.out
            .mix(ALIGNMENT.out.mosdepth_dist, ALIGNMENT.out.mosdepth_summ)
            .collect()

        multiqc(multi_ch)
    }
}