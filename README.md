# PHAGOS

**PHAge Genome One-Stop analysis** — end-to-end Nextflow pipeline for bacteriophage genome characterization from Illumina paired-end reads.

```
Reads (paired-end)
  │
  ├─ [fastp]          QC + adapter trimming
  ├─ [Unicycler]      de novo assembly
  ├─ [CheckV]         genome completeness + contig selection
  │
  ├─ [bwa-mem2]       read mapping → coverage depth
  │
  ├─ [taxmyPHAGE]     ICTV taxonomy assignment
  ├─ [Pharokka]       functional annotation (PHROGS, tRNA, CRISPR)
  ├─ [Phold]          structural-homology re-annotation (improves ~15 pp)
  ├─ [VIBRANT]        lifestyle prediction (lytic / lysogenic)
  │
  ├─ [comparative_prep]   FastANI + VIRIDIC-equivalent + reference download
  ├─ [Clinker]            synteny visualization
  ├─ [clinker_phrogs]     PHROGS-colored synteny HTML
  │
  ├─ [phylo_prep]     marker protein download (ICTV MSL + NCBI)
  ├─ [MAFFT]          multiple sequence alignment
  ├─ [trimAl]         alignment trimming
  ├─ [IQ-TREE2]       maximum-likelihood phylogeny (per marker)
  │
  └─ [genome_map]     publication-quality circular genome figure
```

---

## Quick start

```bash
# 1. Copy and edit the config
cp config.yaml my_run.yaml
nano my_run.yaml          # set reads_r1/r2, phage_name, db paths

# 2. Run with Apptainer (recommended on HPC)
nextflow run main.nf \
  -params-file my_run.yaml \
  -profile apptainer,local

# 3. Resume after a failed step
nextflow run main.nf \
  -params-file my_run.yaml \
  -profile apptainer,local \
  -resume
```

---

## Configuration

All user-facing parameters live in `config.yaml`. Copy it and edit before running.

### Sample and reads

| Parameter   | Description                          | Example            |
|-------------|--------------------------------------|--------------------|
| `sample`    | Sample ID (used in all output paths) | `"B3"`             |
| `reads_r1`  | Path to R1 FASTQ (can be gzipped)   | `"B3_1.fq.gz"`     |
| `reads_r2`  | Path to R2 FASTQ                     | `"B3_2.fq.gz"`     |
| `outdir`    | Root output directory                | `"results_B3"`     |
| `phage_name`| Display name for figure titles       | `"Ralstonia phage B3"` |
| `host`      | Host organism (metadata only)        | `"Ralstonia solanacearum"` |

### Database paths

| Parameter       | Tool        | Default NAS path                                  |
|-----------------|-------------|---------------------------------------------------|
| `pharokka_db`   | Pharokka    | `/media/nas1/LITRP.DBs/metagenDBs/pharokka_db`   |
| `phold_db`      | Phold       | `/media/nas1/LITRP.DBs/metagenDBs/phold_db`      |
| `checkv_db`     | CheckV      | `/media/nas1/LITRP.DBs/metagenDBs/checkv_db`     |
| `taxmyphage_db` | taxmyPHAGE  | `/media/nas1/LITRP.DBs/metagenDBs/taxmyphage`    |
| `vibrant_db`    | VIBRANT     | `/media/nas1/LITRP.DBs/metagenDBs/vibrant-1.0.1` |
| `msl`           | ICTV MSL    | path to `ICTV_Master_Species_List_*.xlsx`         |

### NCBI Entrez

| Parameter      | Description                                        |
|----------------|----------------------------------------------------|
| `ncbi_email`   | Required for all NCBI downloads                    |
| `ncbi_api_key` | Optional — raises rate limit from 3 to 10 req/s   |

### Analysis parameters

| Parameter          | Default | Description                                           |
|--------------------|---------|-------------------------------------------------------|
| `fastp_min_qual`   | 20      | Minimum base quality for trimming                     |
| `fastp_min_len`    | 50      | Minimum read length after trimming                    |
| `phylo_markers`    | `"terminase large subunit,major capsid protein"` | Comma-separated marker proteins for the phylogeny |
| `clinker_top_n`    | 6       | Number of closest reference genomes in Clinker        |
| `clinker_identity` | 0.3     | Minimum link identity shown in Clinker (0–1)          |

### Skip flags

Set any of these to `true` to skip that phase:

| Flag               | Skips                                                |
|--------------------|------------------------------------------------------|
| `skip_vibrant`     | VIBRANT lifestyle prediction                         |
| `skip_comparative` | FastANI + VIRIDIC similarity + Clinker (entire phase 4) |
| `skip_clinker`     | Clinker + PHROGS coloring only (FastANI still runs)  |
| `skip_phylo`       | MAFFT + trimAl + IQ-TREE2 (entire phase 5)          |

---

## Output structure

```
results_B3/
├── B3/
│   ├── fastp/
│   │   ├── fastp.html              QC report
│   │   └── fastp.json
│   ├── B3.fasta                    Selected genome contig
│   ├── checkv/
│   │   └── quality_summary.tsv    Completeness and quality metrics
│   ├── taxmy/
│   │   └── Summary_taxonomy.tsv   ICTV taxonomy assignment
│   ├── pharokka/
│   │   └── B3.gbk                 Initial annotation
│   ├── phold/
│   │   ├── phold.gbk              Final annotation ← used by all downstream steps
│   │   └── phold_all_cds_functions.tsv
│   ├── vibrant/
│   │   └── VIBRANT_B3/            Lifestyle prediction results
│   ├── comparative/
│   │   ├── fastani_results.tsv    All-vs-all ANI
│   │   ├── viridic_similarity.tsv VIRIDIC-equivalent intergenomic similarity
│   │   ├── metadata_genomes.tsv   Reference genome metadata
│   │   └── clinker_genomes/       GBK files fed to Clinker
│   ├── phylo/
│   │   ├── markers/               Per-marker FASTA files
│   │   ├── *_aligned.faa          MAFFT alignments
│   │   ├── *_trimmed.faa          trimAl-trimmed alignments
│   │   └── *_tree.treefile        IQ-TREE2 best trees
│   └── figures/
│       ├── genome_map_B3.pdf      Circular genome map (vector)
│       ├── genome_map_B3.png      Circular genome map (600 DPI)
│       ├── clinker_synteny.html   Synteny browser (default Clinker colors)
│       ├── clinker_synteny.svg
│       └── clinker_phrogs.html    Synteny browser colored by PHROGS category ★
└── pipeline_info/
    ├── timeline.html
    ├── report.html
    └── trace.tsv
```

---

## Tool choices and rationale

| Step | Tool | Version | Why |
|------|------|---------|-----|
| QC | fastp | 1.0.0 | Fastest adapter trimmer; auto-detects adapters |
| Assembly | Unicycler | 0.5.1 | Best for short-read-only phage assembly; tests multiple k-mers, better repeat resolution than SPAdes `--phage` |
| Quality | CheckV | 1.0.3 | Gold standard for phage completeness; SELECT_CONTIG picks Complete > High-quality > longest |
| Taxonomy | taxmyPHAGE | 0.3.0 | ICTV-compliant assignment; outputs genus/family/order used by downstream scripts |
| Annotation | Pharokka | 1.9.1 | PHROGS-based functional annotation + tRNA (tRNAscan-SE) + CRISPR (MinCED) |
| Re-annotation | Phold | 1.2.5 | Structural homology via ESM + Foldseek; annotates ~49% of CDSs vs ~35% with Pharokka alone. `phold.gbk` is the final annotation used by all downstream steps |
| Lifestyle | VIBRANT | 1.2.1 | Mechanistic evidence (integrase, CI repressor absence); preferred over BACPHLIP which only returns a probability |
| Comparative | FastANI + BLASTn | 2.3 / 2.16 | FastANI for quick ANI; BLASTn-based formula replicates VIRIDIC species/genus demarcation without the web server |
| Synteny | Clinker | 0.0.29 | Standard for phage comparative genomics figures |
| Alignment | MAFFT | 7.525 | L-INS-i for high accuracy on short marker proteins |
| Trimming | trimAl | 1.5.0 | `-automated1` mode; removes poorly aligned columns |
| Phylogeny | IQ-TREE2 | 2.3.6 | ModelFinder (MFP) + ultrafast bootstrap (BB=1000) + SH-aLRT (alrt=1000) |
| Genome map | pycirclize | custom | Publication-quality circular figure; PHROGS colors, GC content, GC skew, coverage |

---

## Execution profiles

Combine one container profile with one executor profile:

```bash
# Local workstation with Apptainer
-profile apptainer,local

# HPC cluster (SLURM) with Apptainer
-profile apptainer,slurm

# Local with Docker
-profile docker,local

# Conda (no containers required)
-profile conda,local
```

---

## Annotation color scheme

Both `genome_map.py` and `clinker_phrogs_colors.py` use the same colorblind-safe palette (Okabe-Ito + Wong):

| PHROGS category | Color |
|-----------------|-------|
| Head and packaging | #0072B2 (blue) |
| Tail | #009E73 (green) |
| DNA, RNA and nucleotide metabolism | #E69F00 (orange) |
| Lysis | #D55E00 (red-orange) |
| Connector | #CC79A7 (pink) |
| Transcription regulation | #56B4E9 (sky blue) |
| Moron / AMG / host takeover | #F0E442 (yellow) |
| Other | #999999 (grey) |
| Unknown function | #DDDDDD (light grey) |

---

## Taxonomy-driven automation

Two Python scripts run automatically based on taxmyPHAGE output, with no manual input:

**`comparative_prep.py`** — reads `Summary_taxonomy.tsv`, determines the best ICTV rank (genus → family → order, requiring ≥5 species in the MSL), downloads all reference genomes, runs FastANI and VIRIDIC-equivalent similarity, and prepares the top-N GBK files for Clinker.

**`phylo_prep.py`** — same taxonomy logic, downloads marker protein sequences for all species in the clade, auto-selects an outgroup from the nearest clade outside the query taxon, and appends the query phage's own proteins from `phold.gbk`.

---

## Re-running individual phases

Nextflow's `-resume` re-uses cached results from completed processes. To rerun only one phase, use skip flags to disable others:

```bash
# Redo only phylogeny (e.g. with different markers)
nextflow run main.nf -params-file config.yaml -profile apptainer,local -resume \
  --skip_comparative true \
  --phylo_markers "terminase large subunit,baseplate"
```
