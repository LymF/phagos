process FASTP {
    tag "$meta.id"
    publishDir { "${params.outdir}/${meta.id}/fastp" }, mode: 'copy'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}_1.clean.fq.gz"), path("${meta.id}_2.clean.fq.gz"), emit: reads
    path "fastp.{html,json}", emit: qc

    script:
    """
    fastp \\
        -i $r1 -I $r2 \\
        -o ${meta.id}_1.clean.fq.gz \\
        -O ${meta.id}_2.clean.fq.gz \\
        -h fastp.html -j fastp.json \\
        -q ${params.fastp_min_qual} \\
        --length_required ${params.fastp_min_len} \\
        -w $task.cpus
    """
}
