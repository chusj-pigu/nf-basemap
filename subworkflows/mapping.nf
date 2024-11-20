include { mapping } from '../modules/minimap'
include { sam_sort } from '../modules/samtools'
include { mosdepth } from '../modules/mosdepth'

workflow ALIGNMENT {
    take:
    fastq
    ref

    main:
    mapping(ref, fastq)
    sam_sort(mapping.out)
    mosdepth(sam_sort.out)

    emit: 
    bam = sam_sort.out
    mosdepth_all_out = mosdepth.out
}