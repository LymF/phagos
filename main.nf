#!/usr/bin/env nextflow
// =============================================================================
// PHAGOS — PHAge Genome One-Stop analysis (DSL2)
// Reads -> Assembly -> Quality -> Annotation -> Comparative -> Phylogeny -> Figures
//
// Run:
//   nextflow run main.nf -params-file config.yaml -profile apptainer,local
//   nextflow run main.nf -params-file config.yaml -profile docker,local
//   nextflow run main.nf -params-file config.yaml -profile conda,local
//
// Resume after failure:
//   nextflow run main.nf -params-file config.yaml -profile apptainer,local -resume
// =============================================================================

nextflow.enable.dsl = 2

include { FASTP }                                     from './modules/fastp'
include { UNICYCLER }                                 from './modules/unicycler'
include { CHECKV; SELECT_CONTIG }                     from './modules/checkv'
include { BWA_INDEX; BWA_MEM; SAMTOOLS_DEPTH }        from './modules/coverage'
include { TAXMYPHAGE; PHAROKKA; PHOLD; VIBRANT }      from './modules/annotation'
include { COMPARATIVE_PREP; CLINKER; CLINKER_PHROGS } from './modules/comparative'
include { PHYLO_PREP; MAFFT; TRIMAL; IQTREE }        from './modules/phylogeny'
include { GENOME_MAP }                                from './modules/genome_map'

// ── Required parameter validation ─────────────────────────────────────────
def validate_params() {
    def errors = []
    if (!params.reads_r1)       errors << "--reads_r1 is required"
    if (!params.reads_r2)       errors << "--reads_r2 is required"
    if (!params.pharokka_db)    errors << "pharokka_db not set in config.yaml"
    if (!params.phold_db)       errors << "phold_db not set in config.yaml"
    if (!params.checkv_db)      errors << "checkv_db not set in config.yaml"
    if (!params.ncbi_email)     errors << "ncbi_email not set in config.yaml"
    if (errors) {
        log.error "Invalid parameters:\n  " + errors.join("\n  ")
        System.exit(1)
    }
}

// ── Main workflow ──────────────────────────────────────────────────────────
workflow {

    validate_params()

    // Startup log
    log.info """
    ╔══════════════════════════════════════════════════╗
    ║      PHAGOS — PHAge Genome One-Stop v1.0         ║
    ╠══════════════════════════════════════════════════╣
    ║  Sample   : ${params.sample}
    ║  Reads R1 : ${params.reads_r1}
    ║  Reads R2 : ${params.reads_r2}
    ║  Outdir   : ${params.outdir}
    ║  Profile  : ${workflow.profile}
    ╚══════════════════════════════════════════════════╝
    """.stripIndent()

    // Sample meta map propagated through all channels
    meta = [id: params.sample]

    // ── Input channels ────────────────────────────────────────────────────
    ch_reads = Channel.of(
        tuple(meta, file(params.reads_r1), file(params.reads_r2))
    )

    // Placeholder file for optional depth input in genome_map
    ch_no_file = Channel.value(tuple(meta, file("$projectDir/assets/NO_FILE")))

    // ── Phase 1: QC and assembly ──────────────────────────────────────────
    FASTP(ch_reads)

    UNICYCLER(FASTP.out.reads)

    CHECKV(UNICYCLER.out.contigs)

    SELECT_CONTIG(
        CHECKV.out.summary,
        UNICYCLER.out.contigs
    )
    ch_genome = SELECT_CONTIG.out.genome

    // ── Phase 2: Coverage ─────────────────────────────────────────────────
    BWA_INDEX(ch_genome)
    BWA_MEM(BWA_INDEX.out.indexed, FASTP.out.reads)
    SAMTOOLS_DEPTH(BWA_MEM.out.bam)
    ch_depth = SAMTOOLS_DEPTH.out.depth

    // ── Phase 3: Taxonomy and annotation ──────────────────────────────────
    TAXMYPHAGE(ch_genome)
    ch_taxonomy = TAXMYPHAGE.out.taxonomy

    PHAROKKA(ch_genome)

    PHOLD(PHAROKKA.out.gbk)
    ch_gbk = PHOLD.out.gbk

    if (!params.skip_vibrant) {
        VIBRANT(ch_genome)
    }

    // ── Phase 4: Comparative genomics ─────────────────────────────────────
    if (!params.skip_comparative) {
        COMPARATIVE_PREP(ch_genome, ch_taxonomy, ch_gbk)

        if (!params.skip_clinker) {
            CLINKER(COMPARATIVE_PREP.out.clinker_gbks)
            CLINKER_PHROGS(CLINKER.out.html)
        }
    }

    // ── Phase 5: Phylogenetics ────────────────────────────────────────────
    if (!params.skip_phylo) {
        PHYLO_PREP(ch_taxonomy, ch_gbk)

        // Each marker .faa becomes an independent alignment + tree job
        ch_markers = PHYLO_PREP.out.marker_fastas
            .transpose()

        MAFFT(ch_markers)
        TRIMAL(MAFFT.out.aligned)
        IQTREE(TRIMAL.out.trimmed)
    }

    // ── Phase 6: Genome map ───────────────────────────────────────────────
    GENOME_MAP(
        ch_gbk,
        ch_genome,
        ch_depth.ifEmpty(ch_no_file)
    )
}

// ── Completion summary ────────────────────────────────────────────────────
workflow.onComplete {
    log.info """
    ╔══════════════════════════════════════════════════╗
    ║  Pipeline ${workflow.success ? "completed ✓" : "FAILED ✗"}
    ║  Duration : ${workflow.duration}
    ║  Outdir   : ${params.outdir}/${params.sample}/
    ╚══════════════════════════════════════════════════╝
    """.stripIndent()
}
