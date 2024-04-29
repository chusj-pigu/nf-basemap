process mapping {
    cpus params.t
    publishDir "${params.outdir}"
    label "minimap"

    input:
    path ref
    path fastq

    output:
    path "${fastq.getBaseName(2)}.sam"

    script:
    def m = params.m ? "-y" : ""
    """
    minimap2 -ax map-ont $m -t $task.cpus $ref $fastq > "${fastq.getBaseName(2)}.sam"
    """
}