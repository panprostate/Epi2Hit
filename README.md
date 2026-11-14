## Recurrent mechanisms of biallelic epigenetic inactivation reveal new putative tumour suppressor genes in prostate cancer

Kiriy D, … Pan Prostate Cancer Group (PPCG), Weischenfeldt J., ***Nature Communications*** (2026)
https://www.nature.com/articles/s41467-026-72182-5

> If you use this code or the **Epi2Hit** workflow, please cite the paper above.
> 

### Epi2Hit Workflow
This repository contains the code and data necessary to reproduce the results, figures and tables of the publication.  
The workflow uses already processed data (i.e. processed copy-number profiles, DNA-methylation profiles). The pipelines used to process raw data are described in a companion manuscript by the PPCG Technical Working Group.

![Epi2Hit Workflow](https://github.com/panprostate/Epi2Hit/blob/main/img/workflow.png)

**Epi2Hit** is an integrative computational method that leverages DNA sequencing, epigenomic profiling, gene expression, chromatin accessibility, and Hi-C data to identify biallelic inactivation of tumour suppressor genes. This workflow identifies cases where one allele undergoes DNA methylation in promoter and enhancer regions, while the other allele is affected by genomic loss. 

## Repository Structure

- **analysis/**: Contains the primary code for the workflow.
  - `Epi2Hit.ipynb`: Jupyter Notebook illustrating the steps of the workflow, including function calls and intermediate plots.
  - `TCGA_run.ipynb`: Jupyter Notebook that runs the full eCpG discovery pipeline for TCGA samples on chromosome 8.  
    It:
    - loads chr8 methylation, ATAC, chromHMM, loop (BEDPE), and TCGA expression data,
    - runs probe→gene annotation, filtering, chromatin-state assignment, probe merging, loop annotation, and methylation–expression linear regression,
    - produces the eCpG calls and writes all output files (e.g. `probes_with_loops.tsv`, `probes_with_loops_linreg.tsv`, `TCGA_gene_ann_me.pkl`, `methylation_labeled.tsv`, `methylation_peaks.tsv`, `eCpGs.tsv`),
    - and generates the plots described in the manuscript (chromHMM state barplot, merged-region correlation heatmaps, forest plots, and per-probe methylation peak plots).
  - `Epi2Hit.py`: Python script with functions used in the workflow.
  
- **data/**: Contains some publicly available data sets used in the analysis.


- **docs/**: Contains documentation, including [performance.md](docs/performance.md) with runtime benchmarks, hardware specifications, and scalability notes.

## Running Epi2Hit

### Docker

We provide a Docker container that bundles the Epi2Hit workflow (and optional downstream biallelic analysis) with pinned software dependencies to facilitate portable and reproducible execution.

➡️ See **[docker/README.md](docker/README.md)** for full build/run instructions, recommended directory layout, and end-to-end example commands.


## Performance and Scalability

- **Design:** Modular Python-based workflow with multiprocessing support.  
- **Runtime on PPCG (~1,000 samples, ~450,000 probes):**  
  - Jobs were typically submitted with a memory request of up to 40 GB RAM.
  - Preprocessing (input loading, QC): ~30 minutes  
  - Main analysis: ~4 hours (using 16 cores)  
 
- **Scaling:** Runtime increases approximately linearly with sample size. Applicable to other tumor types and multi-omic cohorts.  


## Citation

If you use **Epi2Hit** in your work, please cite:

Kiriy et al. Recurrent mechanisms of biallelic epigenetic inactivation reveal new putative tumour suppressor genes in prostate cancer
