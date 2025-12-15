#!/bin/bash
set -e

# Docker entrypoint script for epi2hit + biallelic_py analysis pipeline

show_help() {
    cat << EOF
epi2hit Analysis Pipeline Docker Container

USAGE:
    entrypoint [COMMAND] [OPTIONS]

COMMANDS:
    run-epi2hit         Run the epi2hit analysis pipeline
    generate-metadata   Generate metadata.yaml for biallelic_py
    help                Show this help message

RUN-EPI2HIT OPTIONS:
    --data-dir DIR              Path to input data directory (required)
    --output-dir DIR            Path where results will be saved (required)
    --methylation-file FILE     Methylation file (default: chr8_methylation.tsv.gz)
    --atac-file FILE            ATAC-seq file (default: Corces_PRAD_Peaks_hg19_chr8.bed)
    --gene-ref-file FILE        Gene reference file (default: gencode.v19.protein_coding_gene.minimal_chr.bed.gz)
    --chromhmm-file FILE        ChromHMM file (default: chromHMM_chr8.bed)
    --loops-file FILE           Loops file (default: loops.bedpe.gz)
    --expression-file FILE      Expression file (default: expression.tsv.gz)
    --correlation-threshold NUM Correlation threshold (default: 0.6)
    --num-cores NUM             Number of CPU cores (default: 4)
    --methylation-threshold NUM Methylation threshold (default: 0.8)

GENERATE-METADATA OPTIONS:
    --ecpgs FILE                eCpGs.tsv file from epi2hit (required)
    --genes-ref FILE            Genes reference file (required)
    --exons-ref FILE            Exons reference file (required)
    --mutations-ref FILE        Mutations/sample file (required)
    --cn-segs FILE              Copy number segments file (required)
    --output-yaml FILE          Output metadata.yaml (default: metadata.yaml)
    --title TEXT                Analysis title (default: 'Biallelic Inactivation Analysis')

EXAMPLES:
    # Run epi2hit analysis
    docker run -v /data:/data epi2hit:latest run-epi2hit \\
        --data-dir /data \\
        --output-dir /results

    # Generate metadata for biallelic_py
    docker run -v /data:/data epi2hit:latest generate-metadata \\
        --ecpgs /results/eCpGs.tsv \\
        --genes-ref /reference/genes.bed \\
        --exons-ref /reference/exons.tsv.gz \\
        --mutations-ref /reference/mutations.maf.gz \\
        --cn-segs /data/cn_segments.txt \\
        --output-yaml /results/metadata.yaml

    # Run biallelic_py (use the standard command)
    docker run -v /data:/data epi2hit:latest \\
        biallelic_inactivation /results/metadata.yaml

EOF
}

if [[ $# -eq 0 ]] || [[ "$1" == "help" ]] || [[ "$1" == "--help" ]] || [[ "$1" == "-h" ]]; then
    show_help
    exit 0
fi

case "$1" in
    run-epi2hit)
        shift
        python /analysis/epi_run_analysis.py "$@"
        ;;
    generate-metadata)
        shift
        python /analysis/generate_metadata.py "$@"
        ;;
    *)
        # Pass through any other command directly (e.g., biallelic_inactivation)
        exec "$@"
        ;;
esac
