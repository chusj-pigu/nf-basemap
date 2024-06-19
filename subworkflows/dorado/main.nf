process basecall {
    label "dorado"

    input:
    path pod5
    val model

    output:
    path "${params.sample_id}.bam"

    script:
    def call = params.simplex ? "basecaller" : "duplex"
    def mod = params.no_mod ? "" : (params.m_bases_path ? "--modified-bases-models ${params.m_bases_path}" : "--modified-bases ${params.m_bases}")
    def dev = params.dorado_cpu ? '-x "cpu"' : ""
    def b = params.b ? "-b $params.b" : ""
    """
    dorado $call $b $dev $model $pod5 $mod > ${params.sample_id}.bam
    """
}

