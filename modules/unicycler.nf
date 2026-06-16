process UNICYCLER {
    tag "$meta.id"
    publishDir "${params.outdir}/${meta.id}/assembly", mode: 'copy'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("assembly.fasta"), emit: contigs
    path "unicycler.log",                   emit: log

    script:
    """
    unicycler \\
        -1 $r1 -2 $r2 \\
        -o . \\
        -t $task.cpus \\
        2>&1 | tee unicycler.log
    """
}
