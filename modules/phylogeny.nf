process PHYLO_PREP {
    tag "$meta.id"
    publishDir "${params.outdir}/${meta.id}/phylo/markers", mode: 'copy'

    input:
    tuple val(meta), path(taxonomy)
    tuple val(meta2), path(gbk)

    output:
    tuple val(meta), path("*.faa"), emit: marker_fastas

    script:
    def api = params.ncbi_api_key ? "--api-key ${params.ncbi_api_key}" : ""
    """
    phylo_prep.py \\
        --taxmyphage $taxonomy \\
        --gbk        $gbk \\
        --msl        ${params.msl} \\
        --markers    "${params.phylo_markers}" \\
        --name       "${params.phage_name}" \\
        --email      ${params.ncbi_email} \\
        $api \\
        --outdir     .
    """
}

process MAFFT {
    tag "${meta.id}:${faa.baseName}"
    publishDir "${params.outdir}/${meta.id}/phylo", mode: 'copy'

    input:
    tuple val(meta), path(faa)

    output:
    tuple val(meta), path("${faa.baseName}_aligned.faa"), emit: aligned

    script:
    """
    mafft --localpair --maxiterate 1000 --thread $task.cpus $faa \\
        > ${faa.baseName}_aligned.faa
    """
}

process TRIMAL {
    tag "${meta.id}:${aligned.baseName}"

    input:
    tuple val(meta), path(aligned)

    output:
    tuple val(meta), path("${aligned.baseName}_trimmed.faa"), emit: trimmed

    script:
    """
    trimal -in $aligned -out ${aligned.baseName}_trimmed.faa -automated1
    """
}

process IQTREE {
    tag "${meta.id}:${trimmed.baseName}"
    publishDir "${params.outdir}/${meta.id}/phylo", mode: 'copy'

    input:
    tuple val(meta), path(trimmed)

    output:
    tuple val(meta), path("${trimmed.baseName}_tree.treefile"), emit: tree
    path "${trimmed.baseName}_tree.*",                          emit: all

    script:
    """
    iqtree2 \\
        -s $trimmed \\
        -m MFP \\
        -T $task.cpus \\
        --alrt 1000 --bb 1000 \\
        -pre ${trimmed.baseName}_tree \\
        --redo
    """
}
