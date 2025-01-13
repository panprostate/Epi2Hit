# Epi2Hit workflow

### This repository contains the code and data necessariy to reproduce the results, figures and tables of the publication "Recurrent mechanisms of biallelic epigenetic inactivation reveal new putative tumour suppressor genes in prostate cancer" by Kiriy D et al. The worklow uses already processed data (i.e. processed copy-number profiles, DNA-methylation profiles). The pipelines used to process raw data are described in a companion manuscript by the PPCG Technical Working Group.

![Epi2Hit Workflow](https://github.com/panprostate/Epi2Hit/blob/main/img/workflow.png)

**Epi2Hit** is an integrative computational method that leverages DNA sequencing, epigenomic profiling, gene expression, chromatin accessibility, and Hi-C data to identify biallelic inactivation of tumor suppressor genes. This workflow identifies cases where one allele undergoes DNA methylation in promoter and enhancer regions, while the other allele is affected by genomic loss. 

## Repository Structure

- **analysis/**: Contains the primary code for the workflow.
  - `Epi2Hit.ipynb`: Jupyter Notebook illustrating the steps of the workflow, including function calls and intermediate plots.
  - `Epi2Hit.py`: Python script with functions used in the workflow.
  
- **data/**: Contains some publicly available data sets used in the analysis.

