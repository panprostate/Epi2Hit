#!/usr/bin/env python
"""
Helper script to generate metadata.yaml for biallelic_py analysis.
Takes paths to reference files and the eCpGs output from epi2hit analysis.
"""

import argparse
import yaml
from pathlib import Path
from datetime import datetime


def main():
    parser = argparse.ArgumentParser(
        description="Generate metadata.yaml for biallelic_py analysis pipeline.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
NOTES:
  - File paths should be accessible within the container
  - If running in Docker with 'docker run -v /host/data:/data', use paths like /data/reference/genes.bed
  - eCpGs file is typically /results/eCpGs.tsv (from epi2hit output)
  - Reference files should be in /data or accessible via mounted volumes
        """
    )

    # eCpGs file (output from epi2hit)
    parser.add_argument(
        "--ecpgs",
        type=Path,
        required=True,
        help="Path to eCpGs.tsv file (output from epi2hit analysis)"
    )

    # Reference files
    parser.add_argument(
        "--genes-ref",
        type=Path,
        required=True,
        help="Path to genes reference BED file"
    )
    parser.add_argument(
        "--exons-ref",
        type=Path,
        required=True,
        help="Path to exons reference TSV file"
    )
    parser.add_argument(
        "--mutations-ref",
        type=Path,
        required=True,
        help="Path to mutations/sample donors MAF file"
    )

    # Input data files
    parser.add_argument(
        "--cn-segs",
        type=Path,
        required=True,
        help="Path to copy number segments file"
    )

    # Output
    parser.add_argument(
        "--output-yaml",
        type=Path,
        default=Path("metadata.yaml"),
        help="Path where metadata.yaml will be saved (default: metadata.yaml)"
    )

    # Optional metadata
    parser.add_argument(
        "--title",
        type=str,
        default="Biallelic Inactivation Analysis",
        help="Title for the analysis (default: 'Biallelic Inactivation Analysis')"
    )

    args = parser.parse_args()

    # Verify all input files exist
    input_files = {
        "eCpGs": args.ecpgs,
        "genes reference": args.genes_ref,
        "exons reference": args.exons_ref,
        "mutations reference": args.mutations_ref,
        "copy number segments": args.cn_segs,
    }

    for name, filepath in input_files.items():
        if not filepath.exists():
            raise FileNotFoundError(f"{name} file not found: {filepath}")

    # Create metadata dictionary
    # Use absolute paths as provided by the user (they should match container mount paths)
    metadata = {
        "title": args.title,
        "date": datetime.now().strftime("%d/%m/%Y"),
        "ref": {
            "genes": {
                "path": str(args.genes_ref),
                "format_driver": "bed",
            },
            "exons": {
                "path": str(args.exons_ref),
                "format_driver": "bed",
            },
            "sample_donors": {
                "path": str(args.mutations_ref),
                "format_driver": "maf",
            },
        },
        "input": [
            {
                "path": str(args.cn_segs),
                "type": "scna",
                "format_driver": "simple_segments",
                "extra_driver_args": {},
            },
            {
                "path": str(args.ecpgs),
                "type": "meth",
                "format_driver": "methylation",
                "extra_driver_args": {},
            },
        ],
        "analyses": [
            {"name": "write_aberrations"},
            {"name": "write_sample_donor"},
            {"name": "annotate_meth"},
            {"name": "homozygous_inactivations"},
            {"name": "summary_biallelic"},
        ],
    }

    # Write metadata to YAML file
    output_path = args.output_yaml.resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, 'w') as f:
        yaml.dump(metadata, f, default_flow_style=False, sort_keys=False)

    print(f"✓ Metadata file generated successfully: {output_path}")
    print(f"\nTo run biallelic_py analysis:")
    print(f"  biallelic_inactivation {output_path}")


if __name__ == "__main__":
    main()
