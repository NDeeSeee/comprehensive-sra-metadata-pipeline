#!/usr/bin/env python3
"""
POSEIDON FASTQ Workflow
Downloads SRA files and converts to FASTQ using sample_list.txt
Maintains the same functionality as the original bash script with improved error handling
"""

import os
import sys
import subprocess
import logging
import argparse
import time
import re
from pathlib import Path
from collections import defaultdict
import gzip
import shutil
import contextlib

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class SRAWorkflow:
    """Handles SRA download and FASTQ conversion workflow"""
    
    def __init__(self, cancer_dir, no_wait: bool = False, poll_interval_sec: int = 60):
        self.cancer_dir = Path(cancer_dir)
        self.sample_list_path = self.cancer_dir / "sample_list.txt"
        self.logs_dir = self.cancer_dir / "logs"
        self.samples = {}  # Will store sample_id -> {srr_ids: [...], status: ...}
        self.no_wait = no_wait
        self.poll_interval_sec = poll_interval_sec
        self.submitted_jobs = {}  # srr_id -> job_id
        
        # Create logs directory if needed
        self.logs_dir.mkdir(exist_ok=True)
        
    def parse_sample_list(self):
        """Parse sample_list.txt and extract SRR/ERR IDs for each sample"""
        if not self.sample_list_path.exists():
            logger.error(f"sample_list.txt not found in {self.cancer_dir}")
            sys.exit(1)
            
        with open(self.sample_list_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                    
                parts = line.split('\t')
                if len(parts) >= 3:
                    sample_id = parts[0]
                    # Extract SRR/ERR IDs from columns 2 and 3
                    srr_ids = self._extract_srr_ids(parts[1], parts[2])
                    self.samples[sample_id] = {
                        'srr_ids': srr_ids,
                        'col2': parts[1],
                        'col3': parts[2],
                        'status': 'PENDING'
                    }
        
        logger.info(f"Found {len(self.samples)} samples")
        
    def _extract_srr_ids(self, col2, col3):
        """Extract unique SRR/ERR IDs from the two columns"""
        ids = set()
        # Combine both columns and split by comma
        combined = f"{col2},{col3}"
        for item in combined.split(','):
            # Remove _1.fastq.gz or _2.fastq.gz suffixes
            cleaned = item.replace('_1.fastq.gz', '').replace('_2.fastq.gz', '')
            # Check if it matches SRR/ERR pattern
            if cleaned.startswith(('SRR', 'ERR')) and cleaned[3:].isdigit():
                ids.add(cleaned)
        return sorted(list(ids))
    
    def load_modules(self):
        """Load required modules (sratoolkit and aspera)"""
        # Best-effort; ignore failures since 'module' may be a shell function
        try:
            subprocess.run('module load sratoolkit/2.10.4', shell=True, check=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            subprocess.run('module load aspera/3.9.1', shell=True, check=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except Exception:
            pass
    
    def download_sra_files(self):
        """Step 1: Download SRA files for all samples"""
        logger.info("=" * 50)
        logger.info("STEP 1: Download SRA files (if needed)")
        logger.info("=" * 50)
        
        for sample_id, sample_data in self.samples.items():
            # Check if sample-level FASTQs already exist
            sample_r1 = self.cancer_dir / f"{sample_id}_1.fastq.gz"
            sample_r2 = self.cancer_dir / f"{sample_id}_2.fastq.gz"
            
            if sample_r1.exists() and sample_r2.exists():
                logger.info(f"  ✓ Sample {sample_id} FASTQs exist, skipping downloads")
                sample_data['status'] = 'DONE'
                continue
            
            # Download each SRR file for this sample
            for srr_id in sample_data['srr_ids']:
                self._download_single_sra(srr_id, sample_id)
    
    def _download_single_sra(self, srr_id, sample_id):
        """Download a single SRA file using prefetch"""
        sra_file = self.cancer_dir / f"{srr_id}.sra"
        srr_r1 = self.cancer_dir / f"{srr_id}_1.fastq.gz"
        srr_r2 = self.cancer_dir / f"{srr_id}_2.fastq.gz"
        
        # Skip if SRR-level FASTQs already exist
        if srr_r1.exists() and srr_r2.exists():
            logger.info(f"    ✓ {srr_id} FASTQs exist, skipping prefetch")
            return
        
        # Skip if SRA file already exists
        if sra_file.exists():
            logger.info(f"    ✓ {srr_id}.sra already exists")
            return
        
        # Download using prefetch
        logger.info(f"    → Downloading {srr_id}...")
        log_out = self.logs_dir / f"prefetch_{srr_id}.out.txt"
        log_err = self.logs_dir / f"prefetch_{srr_id}.err.txt"
        
        try:
            with open(log_out, 'w') as out, open(log_err, 'w') as err:
                result = subprocess.run(
                    ['prefetch', srr_id, '-X', '35000000'],
                    stdout=out, stderr=err, cwd=self.cancer_dir
                )
            
            # Move .sra file from subdirectory if created
            srr_dir = self.cancer_dir / srr_id
            if srr_dir.exists():
                for sra in srr_dir.glob('*.sra'):
                    sra.rename(self.cancer_dir / sra.name)
                try:
                    # Remove empty directory if possible
                    if not any(srr_dir.iterdir()):
                        srr_dir.rmdir()
                except Exception as e:
                    logger.debug(f"Could not remove {srr_dir}: {e}")
            
            # Check if download succeeded
            if sra_file.exists():
                logger.info(f"    ✓ {srr_id}.sra downloaded successfully")
            else:
                # Check for dbGaP access issues
                with open(log_err, 'r') as f:
                    error_content = f.read().lower()
                    if any(x in error_content for x in ['dbgap', 'unauthorized', 'permission denied', '403']):
                        logger.warning(f"    ✗ {srr_id}.sra download failed: dbGaP/permission required")
                        status_file = self.logs_dir / f"prefetch_{srr_id}.status"
                        status_file.write_text("DBGaP_REQUIRED")
                    else:
                        logger.error(f"    ✗ {srr_id}.sra download failed")
                        
        except Exception as e:
            logger.error(f"Error downloading {srr_id}: {e}")
    
    def convert_sra_to_fastq(self):
        """Step 2: Convert SRA files to FASTQ format"""
        logger.info("=" * 50)
        logger.info("STEP 2: Convert SRA to FASTQ")
        logger.info("=" * 50)
        
        # Path to the conversion script (fdump.sh)
        fdump_script = Path("/data/salomonis-archive/FASTQs/NCI-R01/POSEIDON/ValeriiGitRepo/scripts/manual_pipeline/fdump.sh")
        
        jobs_to_wait = {}
        for sample_id, sample_data in self.samples.items():
            # Skip if sample already complete
            if sample_data['status'] == 'DONE':
                continue
            
            for srr_id in sample_data['srr_ids']:
                sra_file = self.cancer_dir / f"{srr_id}.sra"
                srr_r1 = self.cancer_dir / f"{srr_id}_1.fastq.gz"
                srr_r2 = self.cancer_dir / f"{srr_id}_2.fastq.gz"
                
                # Skip if already converted or if marked as dbGaP
                status_file = self.logs_dir / f"prefetch_{srr_id}.status"
                if status_file.exists() and "DBGaP_REQUIRED" in status_file.read_text():
                    continue
                
                if srr_r1.exists() and srr_r2.exists():
                    continue
                
                if sra_file.exists():
                    logger.info(f"→ Submitting conversion job for {srr_id}.sra")
                    try:
                        # Capture bsub output from fdump.sh to parse Job ID
                        proc = subprocess.run(
                            ['bash', str(fdump_script), str(sra_file)],
                            cwd=self.cancer_dir,
                            check=True,
                            capture_output=True,
                            text=True
                        )
                        job_id = self._parse_bsub_job_id(proc.stdout + (proc.stderr or ''))
                        if job_id:
                            self.submitted_jobs[srr_id] = job_id
                            jobs_to_wait[srr_id] = job_id
                            logger.info(f"  ✓ Submitted as Job <{job_id}>")
                        else:
                            logger.warning("  ! Could not parse Job ID; will proceed without waiting for this SRR")
                    except subprocess.CalledProcessError as e:
                        logger.error(f"Failed to submit conversion for {srr_id}: {e}")
        
        # Optionally wait for jobs to finish and then clean up SRAs
        if not self.no_wait and jobs_to_wait:
            self._wait_for_jobs_and_cleanup(jobs_to_wait)

    def _parse_bsub_job_id(self, text: str):
        match = re.search(r"Job\s*<(\d+)>", text or "")
        return match.group(1) if match else None

    def _bjobs_status(self, job_id: str):
        try:
            proc = subprocess.run(
                ['bjobs', '-noheader', '-o', 'stat', job_id],
                capture_output=True, text=True, check=False
            )
            out = (proc.stdout or '').strip()
            err = (proc.stderr or '').strip()
            if proc.returncode != 0 and ('not found' in err.lower() or 'job' in err.lower() and 'not' in err.lower()):
                return 'UNKNOWN'  # Might be finished/cleaned from queue
            return out if out else 'UNKNOWN'
        except Exception:
            return 'UNKNOWN'

    def _gzip_test(self, *paths: Path) -> bool:
        try:
            cmd = ['gzip', '-t'] + [str(p) for p in paths]
            r = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            return r.returncode == 0
        except Exception:
            # Fallback: try opening with Python gzip
            try:
                for p in paths:
                    with gzip.open(p, 'rb') as f:
                        while f.read(1024 * 1024):
                            pass
                return True
            except Exception:
                return False

    def _wait_for_jobs_and_cleanup(self, jobs_to_wait: dict):
        logger.info("Waiting for conversion jobs to finish...")
        remaining = dict(jobs_to_wait)
        while remaining:
            finished = []
            for srr_id, job_id in remaining.items():
                status = self._bjobs_status(job_id)
                if status in ('DONE', 'EXIT', 'ZOMBIE', 'UNKNOWN'):
                    # Check outputs regardless of status
                    srr_r1 = self.cancer_dir / f"{srr_id}_1.fastq.gz"
                    srr_r2 = self.cancer_dir / f"{srr_id}_2.fastq.gz"
                    sra_file = self.cancer_dir / f"{srr_id}.sra"
                    if srr_r1.exists() and srr_r2.exists() and self._gzip_test(srr_r1, srr_r2):
                        if sra_file.exists():
                            try:
                                sra_file.unlink()
                                logger.info(f"  ✓ Cleaned {srr_id}.sra after successful conversion")
                            except Exception as e:
                                logger.warning(f"  ! Could not remove {srr_id}.sra: {e}")
                        finished.append(srr_id)
                    else:
                        if status == 'DONE':
                            logger.warning(f"  ! {srr_id} job DONE but FASTQs missing or invalid; keeping .sra")
                        finished.append(srr_id)
                # else still running (e.g., RUN, PEND)
            # Remove finished from remaining
            for s in finished:
                remaining.pop(s, None)
            if remaining:
                time.sleep(self.poll_interval_sec)
        logger.info("All submitted conversion jobs processed.")
    
    def merge_fastq_files(self):
        """Step 3: Merge SRR-level FASTQs into sample-level FASTQs"""
        logger.info("=" * 50)
        logger.info("STEP 3: Merge SRR FASTQs into sample-level FASTQs")
        logger.info("=" * 50)
        
        for sample_id, sample_data in self.samples.items():
            sample_r1 = self.cancer_dir / f"{sample_id}_1.fastq.gz"
            sample_r2 = self.cancer_dir / f"{sample_id}_2.fastq.gz"
            
            # Skip if already exists and valid
            if sample_r1.exists() and sample_r2.exists():
                if self._gzip_test(sample_r1, sample_r2):
                    continue
                else:
                    logger.warning(f"    ! Found corrupt/incomplete sample FASTQs for {sample_id}; removing to re-merge")
                    with contextlib.suppress(Exception):
                        sample_r1.unlink()
                    with contextlib.suppress(Exception):
                        sample_r2.unlink()
            
            # Check if all SRR FASTQs are ready
            srr_files_r1 = []
            srr_files_r2 = []
            all_ready = True
            
            for srr_id in sample_data['srr_ids']:
                srr_r1 = self.cancer_dir / f"{srr_id}_1.fastq.gz"
                srr_r2 = self.cancer_dir / f"{srr_id}_2.fastq.gz"
                
                if srr_r1.exists() and srr_r2.exists() and self._gzip_test(srr_r1, srr_r2):
                    srr_files_r1.append(srr_r1)
                    srr_files_r2.append(srr_r2)
                else:
                    all_ready = False
                    break
            
            if all_ready and srr_files_r1:
                if len(srr_files_r1) == 1:
                    # Single SRR: create hard link or copy
                    logger.info(f"→ Linking {sample_id} to single SRR file")
                    try:
                        os.link(srr_files_r1[0], sample_r1)
                        os.link(srr_files_r2[0], sample_r2)
                    except OSError:
                        # If hard link fails, copy the file
                        shutil.copy2(srr_files_r1[0], sample_r1)
                        shutil.copy2(srr_files_r2[0], sample_r2)
                else:
                    # Multiple SRRs: merge them atomically and log inputs
                    logger.info(f"→ Merging {len(srr_files_r1)} SRR files for {sample_id}")
                    logger.info("    Inputs: %s", ",".join([p.stem.replace('_1','').replace('_2','') for p in srr_files_r1]))
                    tmp1 = self.cancer_dir / f".{sample_id}_1.fastq.gz.tmp"
                    tmp2 = self.cancer_dir / f".{sample_id}_2.fastq.gz.tmp"
                    with contextlib.suppress(Exception):
                        tmp1.unlink()
                    with contextlib.suppress(Exception):
                        tmp2.unlink()
                    ok1 = self._merge_gzip_files(srr_files_r1, tmp1)
                    ok2 = self._merge_gzip_files(srr_files_r2, tmp2)
                    if ok1 and ok2 and self._gzip_test(tmp1, tmp2):
                        tmp1.replace(sample_r1)
                        tmp2.replace(sample_r2)
                    else:
                        logger.error(f"    ✗ Merge failed for {sample_id}; cleaning partial outputs")
                        with contextlib.suppress(Exception):
                            tmp1.unlink()
                        with contextlib.suppress(Exception):
                            tmp2.unlink()
    
    def _merge_gzip_files(self, input_files, output_file) -> bool:
        """Merge multiple gzipped files into one; return True on success."""
        try:
            with gzip.open(output_file, 'wb') as outf:
                for input_file in input_files:
                    with gzip.open(input_file, 'rb') as inf:
                        shutil.copyfileobj(inf, outf)
            return True
        except Exception as e:
            logger.error(f"    ✗ Error during merge into {output_file}: {e}")
            return False

    def cleanup_empty_srr_dirs(self):
        """Remove empty accession directories (SRR*/ERR*) left by prefetch."""
        for entry in self.cancer_dir.iterdir():
            if entry.is_dir() and (entry.name.startswith('SRR') or entry.name.startswith('ERR')):
                try:
                    if not any(entry.iterdir()):
                        entry.rmdir()
                        logger.info(f"  ✓ Removed empty directory {entry}")
                except Exception:
                    pass
    
    def generate_status_report(self):
        """Step 4: Generate sample status report"""
        logger.info("=" * 50)
        logger.info("STEP 4: Generate sample status snapshot")
        logger.info("=" * 50)
        
        status_file = self.cancer_dir / "sample_list.with_status.txt"
        
        # Update status for each sample
        for sample_id, sample_data in self.samples.items():
            sample_r1 = self.cancer_dir / f"{sample_id}_1.fastq.gz"
            sample_r2 = self.cancer_dir / f"{sample_id}_2.fastq.gz"
            
            if sample_r1.exists() and sample_r2.exists() and self._gzip_test(sample_r1, sample_r2):
                sample_data['status'] = 'DONE'
            else:
                # Check individual SRR statuses
                has_dbgap = False
                has_pending = False
                has_missing = False
                
                for srr_id in sample_data['srr_ids']:
                    srr_r1 = self.cancer_dir / f"{srr_id}_1.fastq.gz"
                    srr_r2 = self.cancer_dir / f"{srr_id}_2.fastq.gz"
                    sra_file = self.cancer_dir / f"{srr_id}.sra"
                    status_path = self.logs_dir / f"prefetch_{srr_id}.status"
                    
                    if srr_r1.exists() and srr_r2.exists():
                        continue
                    elif sra_file.exists():
                        has_pending = True
                    elif status_path.exists() and "DBGaP_REQUIRED" in status_path.read_text():
                        has_dbgap = True
                    else:
                        has_missing = True
                
                if has_dbgap:
                    sample_data['status'] = 'DBGaP_REQUIRED'
                elif has_missing and not has_pending:
                    sample_data['status'] = 'MISSING'
                else:
                    sample_data['status'] = 'PENDING'
        
        # Write status file
        with open(status_file, 'w') as f:
            for sample_id, sample_data in self.samples.items():
                f.write(f"{sample_id}\t{sample_data['col2']}\t{sample_data['col3']}\t{sample_data['status']}\n")
        
        # Print summary
        status_counts = defaultdict(int)
        for sample_data in self.samples.values():
            status_counts[sample_data['status']] += 1
        
        logger.info("\nSummary (by sample status):")
        for status in ['DONE', 'PENDING', 'DBGaP_REQUIRED', 'MISSING']:
            count = status_counts.get(status, 0)
            logger.info(f"  {status}: {count}")
        
        # Final counts
        done_count = sum(1 for s in self.samples.values() if s['status'] == 'DONE')
        logger.info(f"\nExpected sample-level FASTQ files: {len(self.samples) * 2}")
        logger.info(f"Current sample-level FASTQ files: {done_count * 2}")
    
    def run(self):
        """Execute the complete workflow"""
        logger.info("=" * 50)
        logger.info("POSEIDON FASTQ WORKFLOW")
        logger.info("=" * 50)
        logger.info(f"Cancer directory: {self.cancer_dir}")
        logger.info(f"Sample list: {self.sample_list_path}")
        logger.info("=" * 50)
        
        # Change to cancer directory
        os.chdir(self.cancer_dir)
        
        # Execute workflow steps
        self.parse_sample_list()
        self.load_modules()
        self.download_sra_files()
        self.convert_sra_to_fastq()
        self.merge_fastq_files()
        self.generate_status_report()
        
        logger.info("=" * 50)
        logger.info("WORKFLOW COMPLETE")
        logger.info("=" * 50)
        logger.info("Check job status with: bjobs")
        logger.info(f"Monitor logs in: {self.logs_dir}/")
        # Final directory cleanup pass
        self.cleanup_empty_srr_dirs()


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description="POSEIDON FASTQ Workflow - Downloads and processes SRA files",
        usage="%(prog)s <cancer_directory>",
        epilog="Example: %(prog)s /data/salomonis-archive/FASTQs/NCI-R01/POSEIDON/Tumors/Tongue"
    )
    parser.add_argument('cancer_directory', help='Path to cancer directory containing sample_list.txt')
    parser.add_argument('--no-wait', action='store_true', help='Do not wait for conversion jobs; submit and exit')
    
    args = parser.parse_args()
    
    # Validate directory exists
    if not os.path.isdir(args.cancer_directory):
        logger.error(f"Directory not found: {args.cancer_directory}")
        sys.exit(1)
    
    # Run workflow
    workflow = SRAWorkflow(args.cancer_directory, no_wait=args.no_wait)
    workflow.run()


if __name__ == "__main__":
    main()