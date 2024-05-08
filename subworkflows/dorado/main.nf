process basecall {
    publishDir "${params.out_dir}", mode : "copy"
    label "dorado"

    input:
    path pod5
    val model

    output:
    path "${params.sample_id}.bam"

    script:
    def mod = params.no_mod ? "" : (params.m_bases_path ? "--modified-bases-models ${params.m_bases_path}" : "--modified-bases ${params.m_bases}")
    def dev = params.dorado_cpu ? '-x "cpu"' : ""
    """
    dorado duplex $dev $model $pod5 $mod > ${params.sample_id}.bam
    """
}

