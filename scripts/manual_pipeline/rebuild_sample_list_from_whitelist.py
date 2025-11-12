#!/usr/bin/env python3
"""
Rebuild sample_list.txt from whitelist by querying NCBI SRA for run accessions.

Usage:
    python rebuild_sample_list_from_whitelist.py <tissue_name> [--output <path>]

Example:
    python rebuild_sample_list_from_whitelist.py Gallbladder
    python rebuild_sample_list_from_whitelist.py "Bones+Joints" --output /path/to/output.txt
"""

import argparse
import sys
import time
import subprocess
from pathlib import Path
from typing import Dict, List, Set

# Hardcoded paths
WHITELIST_PATH = Path("/data/salomonis-archive/FASTQs/NCI-R01/POSEIDON/White_list_samples.txt")
POSEIDON_DIR = Path("/data/salomonis-archive/FASTQs/NCI-R01/POSEIDON")


def load_whitelist(tissue_name: str) -> List[str]:
    """Load sample IDs from whitelist for specified tissue."""
    if not WHITELIST_PATH.exists():
        print(f"ERROR: Whitelist not found: {WHITELIST_PATH}")
        sys.exit(1)

    samples = []
    with open(WHITELIST_PATH, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            parts = line.split('\t')
            if len(parts) >= 2:
                sample_id = parts[0].strip()
                tissue = parts[1].strip()
                if tissue == tissue_name and sample_id:
                    samples.append(sample_id)

    print(f"Found {len(samples)} samples for {tissue_name} in whitelist")
    return samples


def get_srr_runs_for_sample(sample_id: str) -> List[str]:
    """Query NCBI SRA to get SRR/ERR run accessions for a sample ID.

    Uses esearch + efetch runinfo to get runs from NCBI Entrez API.
    """
    try:
        # Query NCBI using esearch | efetch with runinfo CSV format
        cmd = f'esearch -db sra -query {sample_id} | efetch -format runinfo'

        result = subprocess.run(
            cmd,
            shell=True,
            capture_output=True,
            text=True,
            timeout=30
        )

        if result.returncode != 0:
            print(f"  WARNING: NCBI query failed for {sample_id}")
            return []

        # Parse CSV output (first column is Run accession)
        lines = result.stdout.strip().split('\n')
        if len(lines) < 2:  # No data rows (only header or empty)
            return []

        runs = []
        for line in lines[1:]:  # Skip header
            if not line.strip():
                continue

            # Extract first column (Run accession)
            parts = line.split(',')
            if parts:
                run_id = parts[0].strip()
                # Only include SRR/ERR/DRR accessions
                if run_id.startswith(('SRR', 'ERR', 'DRR')):
                    runs.append(run_id)

        return sorted(list(set(runs)))

    except subprocess.TimeoutExpired:
        print(f"  WARNING: NCBI query timeout for {sample_id}")
        return []
    except Exception as e:
        print(f"  WARNING: Failed to query {sample_id}: {e}")
        return []


def format_sample_line(sample_id: str, runs: List[str]) -> str:
    """Format sample line for sample_list.txt.

    Format: SAMPLE_ID\tRUN1_1.fastq.gz,RUN2_1.fastq.gz\tRUN1_2.fastq.gz,RUN2_2.fastq.gz
    """
    if not runs:
        return None

    # Build R1 and R2 lists
    r1_files = [f"{run}_1.fastq.gz" for run in runs]
    r2_files = [f"{run}_2.fastq.gz" for run in runs]

    return f"{sample_id}\t{','.join(r1_files)}\t{','.join(r2_files)}"


def rebuild_sample_list(tissue_name: str, output_path: Path = None):
    """Rebuild sample_list.txt for a tissue from whitelist."""

    # Determine output path
    if output_path is None:
        tissue_dir = POSEIDON_DIR / "Tumors" / tissue_name
        output_path = tissue_dir / "sample_list.txt"

    # Load whitelist samples for this tissue
    samples = load_whitelist(tissue_name)

    if not samples:
        print(f"ERROR: No samples found for {tissue_name} in whitelist")
        sys.exit(1)

    print(f"\nQuerying NCBI for {len(samples)} samples...")
    print("This may take several minutes...")

    # Query NCBI for each sample
    results = []
    failed = []

    for i, sample_id in enumerate(samples, 1):
        print(f"[{i}/{len(samples)}] Querying {sample_id}...", end=' ', flush=True)

        runs = get_srr_runs_for_sample(sample_id)

        if runs:
            print(f"✓ Found {len(runs)} run(s)")
            line = format_sample_line(sample_id, runs)
            if line:
                results.append(line)
        else:
            print("✗ No runs found")
            failed.append(sample_id)

        # Rate limiting - be nice to NCBI
        if i % 3 == 0:
            time.sleep(1)

    print(f"\n{'='*60}")
    print(f"Results: {len(results)} samples with runs, {len(failed)} failed")
    print(f"{'='*60}")

    if failed:
        print(f"\nFailed samples ({len(failed)}):")
        for sample_id in failed[:10]:
            print(f"  - {sample_id}")
        if len(failed) > 10:
            print(f"  ... and {len(failed) - 10} more")

    if not results:
        print("\nERROR: No samples successfully processed")
        sys.exit(1)

    # Write output
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Backup existing file if it exists
    if output_path.exists():
        backup_path = output_path.with_suffix('.txt.backup')
        print(f"\nBacking up existing file to: {backup_path}")
        import shutil
        shutil.copy2(output_path, backup_path)

    print(f"\nWriting {len(results)} samples to: {output_path}")
    with open(output_path, 'w') as f:
        for line in results:
            f.write(line + '\n')

    print(f"✓ Done! sample_list.txt written with {len(results)} samples")


def main():
    parser = argparse.ArgumentParser(
        description="Rebuild sample_list.txt from whitelist by querying NCBI SRA",
        epilog="Example: python rebuild_sample_list_from_whitelist.py Gallbladder"
    )
    parser.add_argument('tissue_name', help='Tissue name (must match whitelist exactly, e.g., "Gallbladder")')
    parser.add_argument('--output', '-o', type=Path, help='Output path (default: Tumors/<tissue>/sample_list.txt)')

    args = parser.parse_args()

    # Check if NCBI E-utilities are available
    try:
        subprocess.run(['esearch', '-version'], capture_output=True, timeout=5)
    except FileNotFoundError:
        print("ERROR: NCBI E-utilities (esearch/efetch) not found in PATH")
        print("Install with: conda install -c bioconda entrez-direct")
        sys.exit(1)

    rebuild_sample_list(args.tissue_name, args.output)


if __name__ == "__main__":
    main()