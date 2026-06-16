process CHECKV {
    tag "$meta.id"
    publishDir { "${params.outdir}/${meta.id}/checkv" }, mode: 'copy'

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path("quality_summary.tsv"), emit: summary
    tuple val(meta), path("*"),                   emit: all

    script:
    """
    checkv end_to_end \\
        $contigs \\
        . \\
        -t $task.cpus \\
        -d ${params.checkv_db}
    """
}

// Seleciona o contig principal (Complete > High-quality > mais longo)
process SELECT_CONTIG {
    tag "$meta.id"
    publishDir { "${params.outdir}/${meta.id}" }, mode: 'copy'

    input:
    tuple val(meta), path(summary)
    tuple val(meta2), path(contigs)

    output:
    tuple val(meta), path("${meta.id}.fasta"), emit: genome

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    from Bio import SeqIO

    df = pd.read_csv("$summary", sep="\\t")
    complete = df[df["checkv_quality"] == "Complete"]
    if not complete.empty:
        chosen = complete.iloc[0]["contig_id"]
    else:
        hq = df[df["checkv_quality"] == "High-quality"]
        chosen = hq.iloc[0]["contig_id"] if not hq.empty else df.iloc[0]["contig_id"]

    with open("${meta.id}.fasta", "w") as out:
        for rec in SeqIO.parse("$contigs", "fasta"):
            if rec.id == chosen:
                rec.id = "${meta.id}"
                rec.description = "${params.phage_name}"
                from Bio import SeqIO as _SeqIO
                _SeqIO.write(rec, out, "fasta")
                print(f"Selected: {chosen} ({len(rec)} bp)")
                break
    """
}
