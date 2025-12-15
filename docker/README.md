# epi2hit + biallelic_py Docker Container

A containerized analysis pipeline for epigenetic hits identification (epi2hit) combined with biallelic inactivation analysis.

## Overview

This Docker container provides a portable environment for running two-part genomic analysis:

1. **epi2hit**: Identifies epigenetic hits by integrating methylation, ATAC-seq, and gene expression data
2. **biallelic_py**: Analyzes biallelic inactivation patterns using the eCpGs identified in part 1

## Recommended Data Structure

Organize your data like this to simplify Docker mounting:

```
/your/project/
├── data/                          # Single data directory (mounted as /data in container)
│   ├── input/                     # Input files for epi2hit
│   │   ├── chr8_methylation.tsv.gz
│   │   ├── Corces_PRAD_Peaks_hg19_chr8.bed
│   │   ├── chromHMM_chr8.bed
│   │   ├── loops.bedpe.gz
│   │   ├── expression.tsv.gz
│   │   └── ...other input files
│   ├── gencode/                   # Reference files
│   │   ├── gencode_genes.hg19.nohead.nochr.bed
│   │   ├── gencode.v19.protein_coding_gene.minimal_chr.bed.gz
│   │   └── gencode.v19.annotations_chr8.tsv.gz
│   ├── cBio_prad_tcga_2018_mutations.txt.gz
│   ├── cBio_prad_tcga_2018_cn_segs_chr8.txt.gz
│   └── ...other reference files
└── results/                       # Results directory (mounted as /results in container)
    └── (will contain epi2hit outputs and biallelic_py results)
```

## Quick Start

### 1. Build the Docker Image

```bash
docker build -t epi2hit:latest .
```

### 2. Run epi2hit Analysis

```bash
# Set your data and results directories
DATA_DIR="/your/project/data"
RESULTS_DIR="/your/project/results"
mkdir -p "${RESULTS_DIR}"

# Run the epi2hit analysis
docker run -v "${DATA_DIR}:/data" -v "${RESULTS_DIR}:/results" \
  epi2hit:latest run-epi2hit \
    --data-dir /data \
    --output-dir /results \
    --methylation-file chr8_methylation.tsv.gz \
    --atac-file Corces_PRAD_Peaks_hg19_chr8.bed \
    --gene-ref-file gencode.v19.protein_coding_gene.minimal_chr.bed.gz \
    --chromhmm-file chromHMM_chr8.bed \
    --loops-file loops.bedpe.gz \
    --expression-file expression.tsv.gz \
    --correlation-threshold 0.6 \
    --num-cores 4 \
    --methylation-threshold 0.8
```

### 3. Generate Metadata and Run biallelic_py

After epi2hit completes, run both commands (they must be sequential):

```bash
# Step 1: Generate metadata.yaml
docker run -v "${DATA_DIR}:/data" -v "${RESULTS_DIR}:/results" \
  epi2hit:latest generate-metadata \
    --ecpgs /results/results/eCpGs.tsv \
    --genes-ref /data/gencode/gencode_genes.hg19.nohead.nochr.bed \
    --exons-ref /data/gencode/gencode.v19.annotations_chr8.tsv.gz \
    --mutations-ref /data/cBio_prad_tcga_2018_mutations.txt.gz \
    --cn-segs /data/cBio_prad_tcga_2018_cn_segs_chr8.txt.gz \
    --output-yaml /results/metadata.yaml \
    --title "PRAD TCGA epi2hit biallelic analysis"

# Step 2: Run biallelic_py analysis
# CRITICAL: Must mount the SAME data directory so biallelic_py can access reference files
docker run -v "${DATA_DIR}:/data" -v "${RESULTS_DIR}:/results" \
  epi2hit:latest biallelic_inactivation /results/metadata.yaml
```

## Detailed Usage

### Part 1: epi2hit Analysis

The `run-epi2hit` command executes the full epi2hit pipeline with customizable parameters.

#### Required Arguments

- `--data-dir DIR`: Directory containing all input data files
- `--output-dir DIR`: Directory where results will be saved

#### Input Files (in data-dir)

- `chr8_methylation.tsv.gz`: Methylation data
- `Corces_PRAD_Peaks_hg19_chr8.bed`: ATAC-seq peaks
- `gencode.v19.protein_coding_gene.minimal_chr.bed.gz`: Gene reference
- `chromHMM_chr8.bed`: ChromHMM state annotations
- `loops.bedpe.gz`: Chromosome loop annotations
- `expression.tsv.gz`: Gene expression data

You can override any filename using the `--*-file` arguments (see full options below).

#### Optional Arguments

```
--methylation-file FILE         Methylation file (default: chr8_methylation.tsv.gz)
--atac-file FILE               ATAC-seq file (default: Corces_PRAD_Peaks_hg19_chr8.bed)
--gene-ref-file FILE           Gene reference (default: gencode.v19.protein_coding_gene.minimal_chr.bed.gz)
--chromhmm-file FILE           ChromHMM file (default: chromHMM_chr8.bed)
--loops-file FILE              Loops file (default: loops.bedpe.gz)
--expression-file FILE         Expression file (default: expression.tsv.gz)
--correlation-threshold NUM    Correlation threshold for probe merging (default: 0.6)
--num-cores NUM                Number of CPU cores to use (default: 4)
--methylation-threshold NUM    Methylation threshold for final selection (default: 0.8)
```

#### Example: Run with Custom Parameters

```bash
docker run -v /data:/data -v /results:/results \
  epi2hit:latest run-epi2hit \
    --data-dir /data \
    --output-dir /results \
    --correlation-threshold 0.7 \
    --num-cores 8 \
    --methylation-threshold 0.75
```

#### Output Files

The analysis generates the following in `results/results/`:

- `eCpGs.tsv` - **Main output**: Filtered epigenetic hits with methylation patterns
- `probes_with_loops.tsv` - Probes annotated with chromosome loops
- `probes_with_loops_linreg.tsv` - Probes passing linear regression filtering
- `chromHMM_counts.png` - ChromHMM state distribution visualization
- `TCGA_gene_ann_me.pkl` - Pickled methylation-expression dictionary
- `forest_linreg/` - Forest plots from linear regression analysis
- `peaks/` - Methylation peak visualizations
- `methylation_*.pdf` - Methylation category visualizations

### Part 2: Generate Metadata for biallelic_py

The `generate-metadata` command creates a YAML configuration file for biallelic_py analysis.

#### Required Arguments

- `--ecpgs FILE`: Path to eCpGs.tsv (output from epi2hit)
- `--genes-ref FILE`: Path to genes reference BED file
- `--exons-ref FILE`: Path to exons reference TSV/BED file
- `--mutations-ref FILE`: Path to mutations MAF file
- `--cn-segs FILE`: Path to copy number segments file

#### Optional Arguments

```
--output-yaml FILE   Output metadata.yaml path (default: metadata.yaml)
--title TEXT         Analysis title (default: 'Biallelic Inactivation Analysis')
```

#### Example

```bash
docker run -v /data:/data -v /results:/results \
  epi2hit:latest generate-metadata \
    --ecpgs /results/eCpGs.tsv \
    --genes-ref /data/gencode/genes.hg19.nohead.nochr.bed \
    --exons-ref /data/gencode/exons.tsv.gz \
    --mutations-ref /data/cBio_mutations.txt.gz \
    --cn-segs /data/cn_segs_chr8.txt.gz \
    --output-yaml /results/metadata.yaml \
    --title "PRAD TCGA epi2hit"
```

### Part 3: Run biallelic_py

Execute biallelic_py using the generated metadata:

```bash
docker run -v /results:/results \
  epi2hit:latest biallelic_inactivation /results/metadata.yaml
```

This command is passed directly to the biallelic_py CLI.

## Using Test Scripts as Templates

Two example test scripts are provided to help you get started:

### `test_docker.sh` - epi2hit Analysis
```bash
#!/bin/bash
resdir="$(pwd)/results"
datadir="/Users/lgq442/src/ppcg/epi2hit/epi2hit/data"
mkdir -p "${resdir}"

docker run -v "${datadir}:/data" -v "${resdir}:/results" \
  epi2hit:latest run-epi2hit \
    --data-dir /data \
    --output-dir /results \
    --methylation-file chr8_methylation.tsv.gz \
    --atac-file Corces_PRAD_Peaks_hg19_chr8.bed \
    --gene-ref-file gencode.v19.protein_coding_gene.minimal_chr.bed.gz \
    --chromhmm-file chromHMM_chr8.bed \
    --loops-file loops.bedpe.gz \
    --expression-file expression.tsv.gz \
    --correlation-threshold 0.6 \
    --num-cores 4 \
    --methylation-threshold 0.8
```

**To use it:**
1. Update `datadir` to point to your data directory
2. Run: `bash test_docker.sh`

### `test_docker_bi.sh` - Metadata + biallelic_py Analysis
```bash
#!/bin/bash
datadir="/Users/lgq442/src/ppcg/epi2hit/epi2hit"
resdir="$(pwd)/results"

# Step 1: Generate metadata.yaml
echo "Generating metadata.yaml..."
docker run -v "${datadir}:/data" -v "${resdir}:/results" \
  epi2hit:latest generate-metadata \
    --ecpgs /results/results/eCpGs.tsv \
    --genes-ref /data/gencode/gencode_genes.hg19.nohead.nochr.bed \
    --exons-ref /data/gencode/gencode.v19.annotations_chr8.tsv.gz \
    --mutations-ref /data/cBio_prad_tcga_2018_mutations.txt.gz \
    --cn-segs /data/cBio_prad_tcga_2018_cn_segs_chr8.txt.gz \
    --output-yaml /results/metadata.yaml \
    --title "PRAD TCGA epi2hit biallelic analysis"

# Step 2: Run biallelic_py
echo "Running biallelic_py analysis..."
docker run -v "${datadir}:/data" -v "${resdir}:/results" \
  epi2hit:latest biallelic_inactivation /results/metadata.yaml
```

**To use it:**
1. Ensure epi2hit analysis has completed (results in `./results/results/eCpGs.tsv`)
2. Update `datadir` to point to your data directory
3. Run: `bash test_docker_bi.sh`

## Complete Workflow Example

```bash
# 1. Prepare directories
DATA=/home/user/epi2hit_data
RESULTS=/home/user/epi2hit_results
mkdir -p $RESULTS

# 2. Build image
docker build -t epi2hit:latest .

# 3. Run epi2hit analysis
docker run -v $DATA:/data -v $RESULTS:/results \
  epi2hit:latest run-epi2hit \
    --data-dir /data \
    --output-dir /results \
    --methylation-file chr8_methylation.tsv.gz \
    --atac-file Corces_PRAD_Peaks_hg19_chr8.bed \
    --gene-ref-file gencode.v19.protein_coding_gene.minimal_chr.bed.gz \
    --chromhmm-file chromHMM_chr8.bed \
    --loops-file loops.bedpe.gz \
    --expression-file expression.tsv.gz \
    --num-cores 8

# 4. Generate metadata for biallelic_py
docker run -v $DATA:/data -v $RESULTS:/results \
  epi2hit:latest generate-metadata \
    --ecpgs /results/results/eCpGs.tsv \
    --genes-ref /data/gencode/gencode_genes.hg19.nohead.nochr.bed \
    --exons-ref /data/gencode/gencode.v19.annotations_chr8.tsv.gz \
    --mutations-ref /data/cBio_prad_tcga_2018_mutations.txt.gz \
    --cn-segs /data/cBio_prad_tcga_2018_cn_segs_chr8.txt.gz \
    --output-yaml /results/metadata.yaml

# 5. Run biallelic_py analysis (must mount same data directory)
docker run -v $DATA:/data -v $RESULTS:/results \
  epi2hit:latest biallelic_inactivation /results/metadata.yaml
```

## System Requirements

- Docker installed and running
- Sufficient disk space for input data and results
- For fast processing: allocate multiple CPU cores (--cpus option in docker run)
- At least 4GB RAM

### Recommended Docker Run Options

```bash
docker run \
  --cpus=8 \
  --memory=16g \
  -v $DATA:/data \
  -v $RESULTS:/results \
  epi2hit:latest run-epi2hit ...
```

## Important Notes

### Output Directory Structure

When you run `epi2hit run-epi2hit --output-dir /results`, the actual results are saved to `/results/results/`:
- The `--output-dir` specifies a base directory
- epi2hit internally creates a `results/` subdirectory
- So eCpGs.tsv ends up at `/results/results/eCpGs.tsv`

This is why in the metadata generation step, you use:
```bash
--ecpgs /results/results/eCpGs.tsv
```

### Mounting Strategy for Two-Part Analysis

**Critical**: When running the two-part analysis (epi2hit + biallelic_py), you MUST:

1. **Part 1 (epi2hit)**: Mount data and results directories
   ```bash
   docker run -v /data:/data -v /results:/results \
     epi2hit:latest run-epi2hit --data-dir /data --output-dir /results
   ```

2. **Part 2 (generate metadata)**: Mount the SAME data and results directories
   ```bash
   docker run -v /data:/data -v /results:/results \
     epi2hit:latest generate-metadata ...
   ```

3. **Part 3 (biallelic_py)**: Mount the SAME data directory (required!)
   ```bash
   # biallelic_py needs to read reference files at paths specified in metadata.yaml
   # So data must be mounted at the same location
   docker run -v /data:/data -v /results:/results \
     epi2hit:latest biallelic_inactivation /results/metadata.yaml
   ```

The key is that all three commands must mount the data at `/data` so paths in metadata.yaml remain valid.

## Troubleshooting

### Check Available Commands

```bash
docker run epi2hit:latest help
```

### View Logs

Add verbosity or redirect output:

```bash
docker run -v $DATA:/data -v $RESULTS:/results \
  epi2hit:latest run-epi2hit \
    --data-dir /data \
    --output-dir /results 2>&1 | tee analysis.log
```

### File Not Found Errors

**For epi2hit analysis**: Ensure input files are in the mounted data directory with correct relative names.

**For biallelic_py**: The most common error is that biallelic_py cannot find reference files. This happens when:
1. The data directory is not mounted in the biallelic_py command
2. Paths in metadata.yaml don't match the mounted locations
3. The data directory is mounted at a different path than specified in metadata.yaml

Verify with:
```bash
docker run -v /data:/data epi2hit:latest ls -la /data/gencode/
```

### Example: Scattered Data Directories

If your reference files are in multiple locations, you can mount them separately:

```bash
docker run \
  -v /reference/gencode:/gencode \
  -v /reference/cBio:/cBio \
  -v /input/methylation:/meth \
  -v /results:/results \
  epi2hit:latest run-epi2hit \
    --data-dir /meth \
    --output-dir /results \
    --gene-ref-file /gencode/genes.bed
```

But for the full workflow, use a single mounted data directory as shown above.

## Package Versions

- epi2hit: 0.1.2
- biallelic_py: v0.1.0_epi2hit
- bedtools: latest (from apt-get)
- Python: 3.10

## Additional Information

For issues or questions:
- epi2hit documentation: (check package documentation)
- biallelic_py repository: https://bitbucket.org/weischenfeldt/biallelic_py
- Docker documentation: https://docs.docker.com/

## License

See the respective licenses for epi2hit and biallelic_py packages.
