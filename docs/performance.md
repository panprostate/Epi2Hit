# Performance and Scalability

This document provides details about runtime, resource requirements, and scalability of the Epi2Hit workflow.

## Modular Python-based Design
Epi2Hit is implemented as a modular Python pipeline. The workflow supports both single-core and multi-core execution, allowing users to adapt runtime to available computational resources.

## Multiprocessing Support
- Parallelization is implemented for probe-level and sample-level operations.  
- Users can specify the number of cores/threads to use when running the workflow.  
- Memory usage scales with the number of samples processed in parallel.

## Hardware Used
The benchmark runs were performed on the PPCG cohort (~1,030 tumor samples, ~480,000 probes).  

- **Cluster environment:** HPC server (Slurm-managed)  
- **CPUs:** 4 × AMD Opteron 6380 @ 2.5 GHz (64 logical cores total, 32 physical cores)  
- **Memory:** 528 GB RAM (44 GB swap)  
- **Operating system:** Debian GNU/Linux 9 (stretch)  

## Approximate Runtime
- **Preprocessing (input loading, QC):** ~30 minutes  
- **Main analysis (Epi2Hit workflow):** ~4 hours (using 16 cores)   

## Scalability
- Runtime scales approximately linearly with the number of samples.  
- Multiprocessing across probes and samples enables efficient use of available cores.  
- The workflow is not limited to prostate cancer: it can be applied to other tumor types and multi-omic cohorts with minimal modification.

