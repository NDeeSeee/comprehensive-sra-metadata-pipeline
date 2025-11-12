#!/usr/bin/env python3
"""
Temporary script to process tissue metadata through the full workflow:
1. Filter comprehensive_metadata_paired.tsv by excluding specified BioProjects
2. Copy the classified version from classification directory
3. Convert to XLSX workbook
4. Run srr-to-geo to generate sample lists

Usage:
    python3 process_tissue_metadata.py <tissue_name> [--exclude PRJNA123 PRJNA456]
"""

import argparse
import csv
import subprocess
import sys
from pathlib import Path


def filter_metadata_by_bioproject(input_file, output_file, excluded_projects):
    """Filter TSV file to exclude specified BioProjects."""
    excluded_set = set(excluded_projects)

    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(infile, delimiter='\t')
        fieldnames = reader.fieldnames

        writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()

        total_lines = 0
        kept_lines = 0
        excluded_lines = 0

        for row in reader:
            total_lines += 1
            if row['BioProject'] not in excluded_set:
                writer.writerow(row)
                kept_lines += 1
            else:
                excluded_lines += 1

    print(f"Filtering complete:")
    print(f"  Total input rows: {total_lines}")
    print(f"  Rows kept: {kept_lines}")
    print(f"  Rows excluded: {excluded_lines}")

    return kept_lines


def run_command(cmd, description):
    """Run a shell command and print output."""
    print(f"\n{description}...")
    print(f"Command: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.stdout:
        print(result.stdout)
    if result.stderr:
        print(result.stderr, file=sys.stderr)

    if result.returncode != 0:
        print(f"Error: Command failed with exit code {result.returncode}")
        sys.exit(1)

    return result


def main():
    parser = argparse.ArgumentParser(
        description="Process tissue metadata through filtering, classification, and conversion workflow"
    )
    parser.add_argument(
        "tissue_name",
        help="Tissue name (e.g., 'Bones_and_Joints', 'Melanoma_of_the_skin')"
    )
    parser.add_argument(
        "--exclude",
        nargs="+",
        default=[],
        help="BioProject IDs to exclude (e.g., PRJNA223420 PRJNA1144979)"
    )
    parser.add_argument(
        "--base-dir",
        default="/data/salomonis-archive/FASTQs/NCI-R01/POSEIDON/ValeriiGitRepo/cancer_analysis_output",
        help="Base directory for cancer analysis output"
    )
    parser.add_argument(
        "--sra-dir",
        default="/data/salomonis-archive/FASTQs/NCI-R01/POSEIDON/SRAMetadataFiles",
        help="Directory for SRA metadata files"
    )
    parser.add_argument(
        "--skip-filter",
        action="store_true",
        help="Skip filtering step (use if already done)"
    )

    args = parser.parse_args()

    # Setup paths
    base_dir = Path(args.base_dir)
    tissue_dir = base_dir / args.tissue_name
    metadata_dir = tissue_dir / "metadata"
    classification_dir = tissue_dir / "classification"
    sra_dir = Path(args.sra_dir)

    # Validate directories
    if not metadata_dir.exists():
        print(f"Error: Metadata directory not found: {metadata_dir}")
        sys.exit(1)

    # Define files
    input_metadata = metadata_dir / "comprehensive_metadata_paired.tsv"
    filtered_metadata = metadata_dir / "comprehensive_metadata_paired_selected.tsv"
    classified_metadata = classification_dir / "comprehensive_metadata_paired_selected_classified.tsv"

    # Step 1: Filter metadata (optional)
    if not args.skip_filter:
        if not input_metadata.exists():
            print(f"Error: Input metadata file not found: {input_metadata}")
            sys.exit(1)

        if args.exclude:
            print(f"Step 1: Filtering metadata to exclude BioProjects: {', '.join(args.exclude)}")
            kept_count = filter_metadata_by_bioproject(
                input_metadata,
                filtered_metadata,
                args.exclude
            )

            if kept_count == 0:
                print("Warning: No samples remaining after filtering!")
                sys.exit(1)
        else:
            print("No BioProjects to exclude, skipping filtering step.")
            print("Use --exclude to specify BioProjects to filter out.")
    else:
        print("Skipping filtering step (--skip-filter specified)")

    # Step 2: Check for classified metadata
    if not classified_metadata.exists():
        print(f"\nError: Classified metadata not found: {classified_metadata}")
        print("Please run cancer_classification.py first!")
        sys.exit(1)

    # Step 3: Copy classified metadata to SRA directory
    print(f"\nStep 2: Copying classified metadata to SRA directory")
    # Format tissue name for output file (e.g., Bones_and_Joints -> Bones+Joints)
    output_name = args.tissue_name.replace("_", "+").replace(" ", "+")
    target_tsv = sra_dir / f"{output_name}_target.tsv"

    run_command(
        ["cp", "-f", str(classified_metadata), str(target_tsv)],
        f"Copying to {target_tsv}"
    )

    # Step 4: Convert TSV to XLSX
    print(f"\nStep 3: Converting TSV to XLSX workbook")
    tsv_to_workbook_script = sra_dir / "tsv_to_workbook.py"

    if not tsv_to_workbook_script.exists():
        print(f"Error: tsv_to_workbook.py not found in {sra_dir}")
        sys.exit(1)

    run_command(
        ["python3", str(tsv_to_workbook_script), str(target_tsv)],
        "Running tsv_to_workbook.py"
    )

    # Step 5: Rename XLSX file
    print(f"\nStep 4: Renaming XLSX file")
    target_xlsx = sra_dir / f"{output_name}_target.xlsx"
    final_xlsx = sra_dir / f"{output_name}.xlsx"

    if target_xlsx.exists():
        run_command(
            ["mv", "-f", str(target_xlsx), str(final_xlsx)],
            f"Renaming to {final_xlsx.name}"
        )

    # Step 6: Run srr-to-geo
    print(f"\nStep 5: Running srr-to-geo to generate sample lists")
    run_command(
        ["srr-to-geo", str(final_xlsx)],
        "Running srr-to-geo"
    )

    print(f"\n{'='*60}")
    print(f"SUCCESS! Workflow completed for {args.tissue_name}")
    print(f"{'='*60}")
    print(f"Output files:")
    print(f"  - TSV: {target_tsv}")
    print(f"  - XLSX: {final_xlsx}")
    print(f"  - Sample lists created in POSEIDON directory structure")


if __name__ == "__main__":
    main()
