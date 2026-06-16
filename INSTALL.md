# Installation

## Requirements

| Dependency | Minimum version | Notes |
|------------|-----------------|-------|
| Java | 11 | Required by Nextflow |
| Nextflow | 23.04 | |
| Apptainer / Singularity | 1.1 | Recommended for HPC |
| Docker | 20.10 | Alternative to Apptainer |
| Conda / Mamba | — | Fallback; slower than containers |

At least **one** container runtime (Apptainer or Docker) is strongly recommended. All tools run inside pre-built BioContainers images — no manual installation of bioinformatics software is required.

---

## 1. Install Nextflow

```bash
# Install to ~/bin (add to PATH if needed)
curl -s https://get.nextflow.io | bash
chmod +x nextflow
mv nextflow ~/bin/

# Verify
nextflow -version   # should print >= 23.04
```

---

## 2. Install a container runtime

### Apptainer (recommended for HPC / Linux)

```bash
# Debian / Ubuntu
sudo apt-get install -y apptainer

# RHEL / Rocky / AlmaLinux
sudo dnf install -y apptainer

# Verify
apptainer --version
```

> On shared HPC clusters Apptainer is usually already available as a module:
> `module load apptainer` or `module load singularity`

### Docker (local workstations)

Follow the official guide at https://docs.docker.com/engine/install/ then:

```bash
# Allow running without sudo (log out and back in after)
sudo usermod -aG docker $USER

# Verify
docker run hello-world
```

---

## 3. Install the pipeline

```bash
git clone https://github.com/LymF/phagos.git
cd phagos

# Or copy directly to the server
scp -r phage-nf/ user@server:/path/to/
```

No Python packages need to be installed on the host — the `bin/` scripts run inside their respective containers.

---

## 4. Download databases

All databases can be placed anywhere; update the paths in `config.yaml` accordingly.

### Pharokka database

```bash
mkdir -p /databases/pharokka_db
docker run --rm -v /databases/pharokka_db:/db \
  quay.io/biocontainers/pharokka:1.9.1--pyhdfd78af_0 \
  install_databases.py -o /db

# or with Apptainer
apptainer exec docker://quay.io/biocontainers/pharokka:1.9.1--pyhdfd78af_0 \
  install_databases.py -o /databases/pharokka_db
```

### Phold database

```bash
mkdir -p /databases/phold_db
docker run --rm -v /databases/phold_db:/db \
  quay.io/biocontainers/phold:1.2.5--pyhdfd78af_0 \
  phold install-db -o /db
```

### CheckV database

```bash
mkdir -p /databases/checkv_db
docker run --rm -v /databases:/db \
  quay.io/biocontainers/checkv:1.0.3--pyhdfd78af_0 \
  checkv download_database /db/checkv_db
```

### taxmyPHAGE database

```bash
mkdir -p /databases/taxmyphage
docker run --rm -v /databases/taxmyphage:/db \
  quay.io/biocontainers/taxmyphage:0.3.0--pyhdfd78af_0 \
  taxmyphage setup --db /db
```

### VIBRANT database

```bash
mkdir -p /databases/vibrant-1.0.1
docker run --rm -v /databases/vibrant-1.0.1:/db \
  quay.io/biocontainers/vibrant:1.2.1--hdfd78af_1 \
  download-db.sh /db
```

### ICTV Master Species List (MSL)

Download the latest MSL Excel file from the ICTV website:

```bash
# Example — check https://ictv.global/msl for the current release
wget -O /databases/ICTV_Master_Species_List_2025_MSL41.v1.xlsx \
  "https://ictv.global/filebrowser/api/public/dl/..."
```

Update `msl:` in `config.yaml` with the full path to the downloaded file.

---

## 5. Configure the pipeline

```bash
cp config.yaml my_phage.yaml
```

Edit `my_phage.yaml` and set at minimum:

```yaml
sample:     "MyPhage"
reads_r1:   "/path/to/reads_R1.fq.gz"
reads_r2:   "/path/to/reads_R2.fq.gz"
outdir:     "results_MyPhage"
phage_name: "My Phage Name"

pharokka_db:   "/databases/pharokka_db"
phold_db:      "/databases/phold_db"
checkv_db:     "/databases/checkv_db"
taxmyphage_db: "/databases/taxmyphage"
vibrant_db:    "/databases/vibrant-1.0.1"
msl:           "/databases/ICTV_Master_Species_List_2025_MSL41.v1.xlsx"

ncbi_email:   "your@email.com"
```

---

## 6. Run the pipeline

```bash
# Local workstation with Apptainer
nextflow run main.nf \
  -params-file my_phage.yaml \
  -profile apptainer,local

# HPC (SLURM) with Apptainer
nextflow run main.nf \
  -params-file my_phage.yaml \
  -profile apptainer,slurm

# Docker (local)
nextflow run main.nf \
  -params-file my_phage.yaml \
  -profile docker,local
```

---

## 7. Custom container (genome_map + phylo_prep)

Two processes (`GENOME_MAP` and `PHYLO_PREP`) use a custom image that bundles
`pycirclize`, `matplotlib`, `biopython`, `pandas`, and `openpyxl`.

The image `ghcr.io/lymf/vapor-genome-map:1.0` is pre-built and public.
If you need to rebuild it (e.g. to add packages):

```dockerfile
# Dockerfile
FROM python:3.11-slim
RUN pip install --no-cache-dir \
    pycirclize matplotlib biopython pandas openpyxl requests numpy
```

```bash
docker build -t ghcr.io/<your-org>/vapor-genome-map:1.0 .
docker push ghcr.io/<your-org>/vapor-genome-map:1.0
```

Then update `nextflow.config`:

```groovy
withName: 'GENOME_MAP|PHYLO_PREP' {
    container = 'ghcr.io/<your-org>/vapor-genome-map:1.0'
}
```

---

## 8. Conda fallback (no containers)

If neither Apptainer nor Docker is available, use the conda profile.  
You must create environments for each tool manually — see each tool's documentation.

A minimal working set:

```bash
mamba create -n phage_pipeline \
  fastp unicycler checkv bwa-mem2 samtools \
  taxmyphage pharokka phold vibrant \
  fastani blast mafft trimal iqtree clinker-py \
  biopython pandas openpyxl pycirclize matplotlib requests -c bioconda -c conda-forge

mamba activate phage_pipeline
nextflow run main.nf -params-file my_phage.yaml -profile conda,local
```

> Note: not all tool versions in `nextflow.config` may be available in conda simultaneously. Containers are the recommended approach.

---

## Troubleshooting

**`ERROR: No such container: quay.io/biocontainers/...`**  
Apptainer/Docker needs internet access on first run to pull images. On air-gapped clusters, pre-pull and convert:

```bash
apptainer pull pharokka.sif docker://quay.io/biocontainers/pharokka:1.9.1--pyhdfd78af_0
# then point nextflow to local .sif files via process.container in nextflow.config
```

**`WARN: Task memory limit exceeded`**  
Phold is the most memory-intensive step (32 GB default). On smaller machines:

```bash
nextflow run main.nf -params-file my_phage.yaml -profile apptainer,local \
  --max_memory 64.GB
```

**`ERROR: Could not find embedded JSON data in Clinker HTML`**  
The Clinker HTML format has changed between versions. Make sure `clinker_phrogs_colors.py` is using the version matching `nextflow.config` (`clinker-py:0.0.29`).

**FastANI / BLAST step is slow**  
These steps scale with the number of reference genomes in the MSL taxon. For large families (>200 species), the VIRIDIC-equivalent step can take several hours. Skip it with `--skip_viridic true` if only FastANI results are needed:

```bash
nextflow run main.nf -params-file my_phage.yaml -profile apptainer,local \
  --skip_viridic true
```

> Note: `--skip_viridic` is passed as a CLI override; it is not in `config.yaml` by default.
