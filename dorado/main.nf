process basecall {
    publishDir "${params.out_dir}"
    label "dorado"
    label "gpu"
    cpus params.t

    input:
    path pod5
    val model

    output:
    path "${params.sample_id}.bam"

    script:
    def mod = params.no_mod ? "" : (profile == "drac" ? "--modified-bases-models ${params.m_bases}" : "--modified-bases ${params.m_bases}")
    """
    dorado duplex $model $pod5 $mod > ${params.sample_id}.bam
    """
}

