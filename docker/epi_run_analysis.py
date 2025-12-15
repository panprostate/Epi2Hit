#!/usr/bin/env python
"""
Entry point script for epi2hit analysis pipeline.
Parameterizes the hardcoded paths from epi_run_test.py to allow for flexible input/output directories.
"""

import argparse
import os
from pathlib import Path

from epi2hit.core import (
    read_dataset,
    find_open_chromatin_probes,
    probes2genes,
    prelim_filters,
    pomerantz_state,
    process_probes,
    resolve_duplicates,
    annotate_probes_with_loops,
    prepare_TCGA_dict,
    lin_reg,
    slice_probes_with_loops_by_linreg,
    methylation,
    filter_probes_by_methylation,
)


def main():
    parser = argparse.ArgumentParser(
        description="Run the epi2hit analysis pipeline for epigenetic hits identification."
    )

    # Input/Output directories
    parser.add_argument(
        "--data-dir",
        type=Path,
        required=True,
        help="Path to directory containing input data files"
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Path to directory where results will be saved"
    )

    # Input files
    parser.add_argument(
        "--methylation-file",
        type=Path,
        default="chr8_methylation.tsv.gz",
        help="Methylation data file (relative to data-dir or absolute path)"
    )
    parser.add_argument(
        "--atac-file",
        type=Path,
        default="Corces_PRAD_Peaks_hg19_chr8.bed",
        help="ATAC-seq peaks BED file (relative to data-dir or absolute path)"
    )
    parser.add_argument(
        "--gene-ref-file",
        type=Path,
        default="gencode.v19.protein_coding_gene.minimal_chr.bed.gz",
        help="Gene reference BED file (relative to data-dir or absolute path)"
    )
    parser.add_argument(
        "--chromhmm-file",
        type=Path,
        default="chromHMM_chr8.bed",
        help="ChromHMM state file (relative to data-dir or absolute path)"
    )
    parser.add_argument(
        "--loops-file",
        type=Path,
        default="loops.bedpe.gz",
        help="Chromosome loops BEDPE file (relative to data-dir or absolute path)"
    )
    parser.add_argument(
        "--expression-file",
        type=Path,
        default="expression.tsv.gz",
        help="Gene expression file (relative to data-dir or absolute path)"
    )

    # Analysis parameters
    parser.add_argument(
        "--correlation-threshold",
        type=float,
        default=0.6,
        help="Correlation threshold for merging probes (default: 0.6)"
    )
    parser.add_argument(
        "--num-cores",
        type=int,
        default=4,
        help="Number of CPU cores to use (default: 4)"
    )
    parser.add_argument(
        "--methylation-threshold",
        type=float,
        default=0.8,
        help="Methylation threshold for final eCpG selection (default: 0.8)"
    )

    args = parser.parse_args()

    # Resolve file paths (allow both relative to data-dir and absolute paths)
    def resolve_path(file_arg):
        file_path = Path(file_arg)
        if file_path.is_absolute():
            return file_path
        else:
            return args.data_dir / file_path

    meth_file = resolve_path(args.methylation_file)
    atac_file = resolve_path(args.atac_file)
    gene_ref_file = resolve_path(args.gene_ref_file)
    chromhmm_file = resolve_path(args.chromhmm_file)
    loops_file = resolve_path(args.loops_file)
    expr_file = resolve_path(args.expression_file)

    # Create output directory
    results_dir = args.output_dir / "results"
    results_dir.mkdir(parents=True, exist_ok=True)

    # Define output files
    probes_with_loops_tsv = results_dir / "probes_with_loops.tsv"
    probes_with_loops_linreg_tsv = results_dir / "probes_with_loops_linreg.tsv"
    tcga_dict_pkl = results_dir / "TCGA_gene_ann_me.pkl"
    chromhmm_counts_png = results_dir / "chromHMM_counts.png"
    ecpg_tsv = results_dir / "eCpGs.tsv"
    forest_dir = results_dir / "forest_linreg"
    peaks_dir = results_dir / "peaks"
    corr_plots_dir = results_dir / "corr_plots"

    print(f"Starting epi2hit analysis pipeline...")
    print(f"Data directory: {args.data_dir}")
    print(f"Output directory: {results_dir}")
    print(f"Using {args.num_cores} CPU cores")

    # Verify input files exist
    input_files = {
        "methylation": meth_file,
        "ATAC": atac_file,
        "gene reference": gene_ref_file,
        "ChromHMM": chromhmm_file,
        "loops": loops_file,
        "expression": expr_file,
    }

    for name, filepath in input_files.items():
        if not filepath.exists():
            raise FileNotFoundError(f"{name} file not found: {filepath}")

    # 1) Load methylation and ATAC data
    print("\n[1/13] Loading methylation and ATAC data...")
    me_full = read_dataset(str(meth_file))
    atac = read_dataset(str(atac_file), sep='\t', index_col=None)

    # 2) Restrict to probes in open chromatin
    print("[2/13] Finding probes in open chromatin...")
    me_open = find_open_chromatin_probes(me_full, atac)

    # 3) Map probes to nearest genes
    print("[3/13] Mapping probes to genes...")
    me_gb = probes2genes(me_open, str(gene_ref_file))

    # 4) Apply filters
    print("[4/13] Applying stromal/immune/epiCMIT filters...")
    me_filtered = prelim_filters(me_gb)

    # 5) Add Pomerantz / chromHMM states
    print("[5/13] Adding Pomerantz/chromHMM states...")
    me_with_pom = pomerantz_state(
        me_filtered,
        reference_bed_path=str(chromhmm_file),
        output_png_path=str(chromhmm_counts_png),
    )

    # Change to results directory for correlation plots
    print("[6/13] Processing probes (merging correlated probes)...")
    os.chdir(results_dir)

    final_combined_df = process_probes(
        me_filtered2=me_with_pom,
        correlation_threshold=args.correlation_threshold,
        num_cores=args.num_cores,corr_plots_dir=str(corr_plots_dir)
    )

    # 7) Fix duplicate indices
    print("[7/13] Resolving duplicate indices...")
    final_combined_df = resolve_duplicates(final_combined_df)

    # 8) Annotate with loops
    print("[8/13] Annotating probes with chromosome loops...")
    final_combined_df2 = annotate_probes_with_loops(
        final_combined_df,
        loop_bedpe_path=str(loops_file),
        output_path=str(probes_with_loops_tsv),
    )

    # 9) Build methylation–expression dictionary
    print("[9/13] Building methylation-expression dictionary...")
    tcga_dict = prepare_TCGA_dict(
        me_with_loops=final_combined_df2,
        exp_path=str(expr_file),
        output_path=str(tcga_dict_pkl),
    )

    # 10) Run linear regression + forest plot
    print("[10/13] Running linear regression and generating forest plots...")
    res = lin_reg(
        tcga_dict,
        plot=True,
        plot_dir=str(forest_dir),
    )

    # 11) Keep only probes with significant hits
    print("[11/13] Filtering probes by linear regression results...")
    sliced_probes = slice_probes_with_loops_by_linreg(
        linreg_results=res,
        probes_path=str(probes_with_loops_tsv),
        output_path=str(probes_with_loops_linreg_tsv),
    )

    # 12) Categorize methylation + detect peaks
    print("[12/13] Categorizing methylation patterns and detecting peaks...")
    labeled_methyl, peak_df = methylation(
        sliced_probes,
        path=False,
        plot=True,
        save=True,
        output_dir=str(results_dir),
        base_name="methylation",
        peaks_dir=str(peaks_dir),
    )

    # 13) Final eCpG selection
    print("[13/13] Final eCpG selection...")
    eCpGs = filter_probes_by_methylation(
        labeled_methyl,
        threshold=args.methylation_threshold,
        save=True,
        output_path=str(ecpg_tsv),
    )

    print(f"\n✓ Pipeline completed successfully!")
    print(f"eCpGs results saved to: {ecpg_tsv}")
    print(f"\nTo proceed with biallelic_py analysis, use the output eCpGs file:")
    print(f"  python /analysis/generate_metadata.py --ecpgs {ecpg_tsv} --output-yaml metadata.yaml")


if __name__ == "__main__":
    main()
