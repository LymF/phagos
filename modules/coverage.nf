process BWA_INDEX {
    tag "$meta.id"

    input:
    tuple val(meta), path(genome)

    output:
    tuple val(meta), path(genome), path("${genome}.*"), emit: indexed

    script:
    """
    bwa-mem2 index $genome
    """
}

process BWA_MEM {
    tag "$meta.id"

    input:
    tuple val(meta), path(genome), path(index)
    tuple val(meta2), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}.bam"), emit: bam

    script:
    """
    bwa-mem2 mem -t $task.cpus $genome $r1 $r2 \\
        | samtools sort -@ $task.cpus -o ${meta.id}.bam
    samtools index ${meta.id}.bam
    """
}

process SAMTOOLS_DEPTH {
    tag "$meta.id"
    publishDir { "${params.outdir}/${meta.id}/coverage" }, mode: 'copy'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("depth.tsv"), emit: depth

    script:
    """
    samtools depth -a $bam > depth.tsv
    echo "Mean coverage: \$(awk '{s+=\$3} END {printf "%.0fx", s/NR}' depth.tsv)"
    """
}
