process COMPARATIVE_PREP {
    tag "$meta.id"
    publishDir "${params.outdir}/${meta.id}/comparative", mode: 'copy'

    input:
    tuple val(meta), path(genome)
    tuple val(meta2), path(taxonomy)
    tuple val(meta3), path(gbk)

    output:
    tuple val(meta), path("fastani_results.tsv"),      emit: fastani,  optional: true
    tuple val(meta), path("viridic_similarity.tsv"),   emit: viridic,  optional: true
    tuple val(meta), path("clinker_genomes/"),         emit: clinker_gbks
    path "metadata_genomes.tsv",                       emit: metadata, optional: true

    script:
    def api = params.ncbi_api_key ? "--api-key ${params.ncbi_api_key}" : ""
    """
    comparative_prep.py \\
        --taxmyphage $taxonomy \\
        --genome     $genome \\
        --gbk        $gbk \\
        --msl        ${params.msl} \\
        --email      ${params.ncbi_email} \\
        $api \\
        --threads    $task.cpus \\
        --top-clinker ${params.clinker_top_n} \\
        --outdir     .
    """
}

process CLINKER {
    tag "$meta.id"
    publishDir "${params.outdir}/${meta.id}/figures", mode: 'copy'

    input:
    tuple val(meta), path(clinker_dir)

    output:
    tuple val(meta), path(clinker_dir), path("clinker_synteny.html"), emit: html
    path "clinker_synteny.svg",                                        emit: svg

    script:
    """
    clinker ${clinker_dir}/*.gb* \\
        -i ${params.clinker_identity} \\
        -o clinker_synteny.html \\
        -p clinker_synteny.svg
    """
}

process CLINKER_PHROGS {
    tag "$meta.id"
    publishDir "${params.outdir}/${meta.id}/figures", mode: 'copy'

    input:
    tuple val(meta), path(clinker_dir), path(clinker_html)

    output:
    path "clinker_phrogs.html", emit: html

    script:
    """
    clinker_phrogs_colors.py \\
        --clinker_html $clinker_html \\
        --gbk_dir      $clinker_dir \\
        --output       clinker_phrogs.html
    """
}
