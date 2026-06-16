process TAXMYPHAGE {
    tag "$meta.id"
    publishDir { "${params.outdir}/${meta.id}/taxmy" }, mode: 'copy'

    input:
    tuple val(meta), path(genome)

    output:
    tuple val(meta), path("Summary_taxonomy.tsv"), emit: taxonomy
    path "*",                                       emit: all

    script:
    """
    taxmyphage run \\
        -i $genome \\
        -o . \\
        -p ${meta.id} \\
        -t $task.cpus \\
        -db ${params.taxmyphage_db}
    """
}

process PHAROKKA {
    tag "$meta.id"
    publishDir { "${params.outdir}/${meta.id}/pharokka" }, mode: 'copy'

    input:
    tuple val(meta), path(genome)

    output:
    tuple val(meta), path("${meta.id}.gbk"),                    emit: gbk
    tuple val(meta), path("${meta.id}_cds_functions.tsv"),      emit: functions
    path "*",                                                    emit: all

    script:
    """
    pharokka.py \\
        -i $genome \\
        -o . \\
        -d ${params.pharokka_db} \\
        -t $task.cpus \\
        -p ${meta.id} \\
        --dnaapler \\
        --trna_scan_model bacterial \\
        --force
    """
}

process PHOLD {
    tag "$meta.id"
    publishDir { "${params.outdir}/${meta.id}/phold" }, mode: 'copy'

    input:
    tuple val(meta), path(gbk)

    output:
    tuple val(meta), path("phold.gbk"),                    emit: gbk
    tuple val(meta), path("phold_all_cds_functions.tsv"),  emit: functions
    path "*",                                              emit: all

    script:
    """
    phold run \\
        -i $gbk \\
        -o . \\
        -d ${params.phold_db} \\
        -t $task.cpus \\
        --force
    """
}

process VIBRANT {
    tag "$meta.id"
    publishDir { "${params.outdir}/${meta.id}/vibrant" }, mode: 'copy'

    input:
    tuple val(meta), path(genome)

    output:
    tuple val(meta), path("VIBRANT_${meta.id}/"), emit: results

    script:
    def vdb = params.vibrant_db
    """
    VIBRANT_run.py \\
        -i $genome \\
        -f nucl \\
        -t $task.cpus \\
        -no_plot \\
        -d ${vdb}/databases \\
        -m ${vdb}/files
    """
}
