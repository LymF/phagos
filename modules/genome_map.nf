process GENOME_MAP {
    tag "$meta.id"
    publishDir { "${params.outdir}/${meta.id}/figures" }, mode: 'copy'

    input:
    tuple val(meta), path(gbk)
    tuple val(meta2), path(genome)
    tuple val(meta3), path(depth)   // opcional — pode ser um arquivo vazio

    output:
    path "genome_map_${meta.id}.{pdf,png}", emit: figures

    script:
    def depth_flag = depth.name != 'NO_FILE' ? "--depth $depth" : ""
    """
    genome_map.py \\
        --gbk   $gbk \\
        --fasta $genome \\
        $depth_flag \\
        --name  "${params.phage_name}" \\
        --out   genome_map_${meta.id}.pdf
    """
}
