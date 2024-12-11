include { mapping } from '../modules/minimap'
include { sam_sort } from '../modules/samtools'
include { mosdepth } from '../modules/mosdepth'

workflow ALIGNMENT {
    take:
    sample_sheet
    bed
    ref

    main:
    input_ch = sample_sheet
        .combine(ref)
    mapping(input_ch)
    sam_sort(mapping.out)
    bed_ch = bed
        .map( it -> tuple(it, 1796, 0))
    mos_in = sam_sort.out
        .combine(bed_ch)
    mosdepth(mos_in)

    emit: 
    bam = sam_sort.out
    mosdepth_dist = mosdepth.out.dist
    mosdepth_summ = mosdepth.out.summary
}