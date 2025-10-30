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

    def _sample_has_bam(self, sample_id: str) -> bool:
        """Return True if a BAM for this sample exists (sample_id*.bam)."""
        try:
            for p in self.cancer_dir.glob(f"{sample_id}*.bam"):
                if p.is_file() and p.stat().st_size > 0:
                    return True
        except Exception:
            pass
        return False
        
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
        
        # Path to the conversion script (submit_fastq_dump_jobs.sh)
        fdump_script = Path("/data/salomonis-archive/FASTQs/NCI-R01/POSEIDON/ValeriiGitRepo/scripts/manual_pipeline/submit_fastq_dump_jobs.sh")
        
        jobs_to_wait = {}
        for sample_id, sample_data in self.samples.items():
            
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
                        # Capture bsub output from submit_fastq_I_jobs.sh to parse Job ID
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
        # Global SRA cleanup pass: remove any .sra with verified FASTQs
        self.cleanup_converted_sras_global()
        # Refresh status after conversion stage
        try:
            self.generate_status_report()
        except Exception as e:
            logger.warning(f"Failed to update status report: {e}")

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

    def _is_job_active(self, job_id: str) -> bool:
        status = self._bjobs_status(job_id)
        return status in ('PEND', 'RUN', 'PSUSP', 'USUSP', 'SSUSP')

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
            # Periodic status refresh each poll cycle
            try:
                self.generate_status_report()
            except Exception:
                pass
            if remaining:
                time.sleep(self.poll_interval_sec)
        logger.info("All submitted conversion jobs processed.")
    
    # merge_fastq_files disabled: SRR-level FASTQs are final outputs

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

    def cleanup_converted_sras_global(self):
        """Delete any lingering .sra when both paired FASTQs exist and validate."""
        # Build unique SRR list from samples
        srr_ids = set()
        for s in self.samples.values():
            for sid in s['srr_ids']:
                srr_ids.add(sid)
        removed = 0
        for srr_id in sorted(srr_ids):
            sra_file = self.cancer_dir / f"{srr_id}.sra"
            if not sra_file.exists():
                continue
            r1 = self.cancer_dir / f"{srr_id}_1.fastq.gz"
            r2 = self.cancer_dir / f"{srr_id}_2.fastq.gz"
            if r1.exists() and r2.exists() and self._gzip_test(r1, r2):
                with contextlib.suppress(Exception):
                    sra_file.unlink()
                    removed += 1
        if removed:
            logger.info(f"  ✓ Global cleanup removed {removed} converted .sra files")
    
    def generate_status_report(self):
        """Generate sample status report (4 columns, exact format)."""
        logger.info("=" * 50)
        logger.info("STEP 3: Generate sample status snapshot")
        logger.info("=" * 50)
        
        status_file = self.cancer_dir / "sample_list.with_status.txt"
        lines = []
        for sample_id, sample_data in self.samples.items():
            srr_ids = sample_data['srr_ids']
            r1_list = ",".join([f"{sid}_1.fastq.gz" for sid in srr_ids])
            r2_list = ",".join([f"{sid}_2.fastq.gz" for sid in srr_ids])

            # Determine status by strict priority
            status = None
            # 1) DBGaP_REQUIRED
            for sid in srr_ids:
                status_path = self.logs_dir / f"prefetch_{sid}.status"
                if status_path.exists() and "DBGaP_REQUIRED" in status_path.read_text():
                    status = 'DBGaP_REQUIRED'
                    break
            
            # Check FASTQ completeness once, reuse it
            all_fastqs_ok = True
            for sid in srr_ids:
                r1 = self.cancer_dir / f"{sid}_1.fastq.gz"
                r2 = self.cancer_dir / f"{sid}_2.fastq.gz"
                if not (r1.exists() and r2.exists() and self._gzip_test(r1, r2)):
                    all_fastqs_ok = False
                    break

            # 2) BAM_DONE (only if all fastqs are good AND bam exists)
            if status is None and all_fastqs_ok:
                if any((p.is_file() and p.stat().st_size > 0) for p in self.cancer_dir.glob(f"{sample_id}*.bam")):
                    status = 'BAM_DONE'

            # 3) FASTQ_DONE
            if status is None and all_fastqs_ok:
                status = 'FASTQ_DONE'
  
            # 4) PENDING (FASTQ -> BAM, job ID <...>) — not tracked here
            # 5) PENDING (SRA -> FASTQ, job ID <...>) if any active conversion job
            if status is None:
                active_job_id = None
                for sid in srr_ids:
                    job_id = self.submitted_jobs.get(sid)
                    if job_id and self._is_job_active(job_id):
                        active_job_id = job_id
                        break
                if active_job_id:
                    status = f"PENDING (SRA -> FASTQ, job ID <{active_job_id}>)"
            # 6) PENDING (SRA downloading, job ID <...>) — not applicable
            # 7) PENDING fallback
            if status is None:
                status = 'PENDING'

            lines.append(f"{sample_id}\t{r1_list}\t{r2_list}\t{status}")

        with open(status_file, 'w') as f:
            f.write("\n".join(lines) + "\n")
    
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
        # Refresh status after download stage
        try:
            self.generate_status_report()
        except Exception:
            pass
        self.convert_sra_to_fastq()
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