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
        self.job_poll_cycles = {}  # srr_id -> number of poll cycles observed
        # Minimal color support
        # Enable color when terminal supports it or FORCE_COLOR is set (NO_COLOR disables)
        try:
            force_color = bool(os.environ.get("FORCE_COLOR"))
            no_color = bool(os.environ.get("NO_COLOR"))
            tty_out = sys.stdout.isatty() if hasattr(sys.stdout, "isatty") else False
            tty_err = sys.stderr.isatty() if hasattr(sys.stderr, "isatty") else False
            self._use_color = (force_color or tty_out or tty_err) and not no_color
        except Exception:
            self._use_color = False
        self._C_GREEN = "\033[92m"
        self._C_YELLOW = "\033[93m"
        self._C_RED = "\033[91m"
        self._C_CYAN = "\033[96m"
        self._C_MAGENTA = "\033[95m"
        self._C_RESET = "\033[0m"

        # Create logs directory if needed
        self.logs_dir.mkdir(exist_ok=True)

    def _c(self, text: str, color: str) -> str:
        return f"{color}{text}{self._C_RESET}" if self._use_color else text
        
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
        """Return True if a BAM for this sample exists (sample_id*.bam).

        Checks both the cancer directory root and an optional 'bams/' subdirectory.
        """
        try:
            # Check root
            for p in self.cancer_dir.glob(f"{sample_id}*.bam"):
                if p.is_file() and p.stat().st_size > 0:
                    return True
            # Check bams/ subdir
            bams_dir = self.cancer_dir / "bams"
            if bams_dir.is_dir():
                for p in bams_dir.glob(f"{sample_id}*.bam"):
                    if p.is_file() and p.stat().st_size > 0:
                        return True
        except Exception:
            pass
        return False

    def _has_star_progress_dirs(self, sample_id: str) -> bool:
        """Return True if any STAR output/progress directories exist for this sample.

        This detects work started by external STAR jobs, e.g. directories like:
        <sample_id>__STARpass1, <sample_id>__STARgenome, <sample_id>__STARtmp.
        """
        try:
            for pattern in (
                f"{sample_id}__STARpass1",
                f"{sample_id}__STARgenome",
                f"{sample_id}__STARtmp",
            ):
                for p in self.cancer_dir.glob(pattern):
                    if p.is_dir():
                        # Consider as progress if directory exists (non-empty check optional)
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
        logger.info(self._c("STEP 1: Download SRA files (if needed)", self._C_CYAN))
        logger.info("=" * 50)
        
        for sample_id, sample_data in self.samples.items():
            # Skip entire sample if BAM already exists
            if self._sample_has_bam(sample_id):
                logger.info(self._c(f"✓ {sample_id}: BAM exists; skipping SRA downloads for all SRRs", self._C_GREEN))
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
            logger.info(self._c(f"    ✓ {srr_id} FASTQs exist, skipping prefetch", self._C_GREEN))
            return
        
        # Skip if SRA file already exists
        if sra_file.exists():
            logger.info(self._c(f"    ✓ {srr_id}.sra already exists", self._C_GREEN))
            return
        
        # Download using prefetch
        logger.info(self._c(f"    → Downloading {srr_id}...", self._C_CYAN))
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
                logger.info(self._c(f"    ✓ {srr_id}.sra downloaded successfully", self._C_GREEN))
            else:
                # Check for dbGaP access issues
                with open(log_err, 'r') as f:
                    error_content = f.read().lower()
                    if any(x in error_content for x in ['dbgap', 'unauthorized', 'permission denied', '403']):
                        logger.warning(self._c(f"    ✗ {srr_id}.sra download failed: dbGaP/permission required", self._C_YELLOW))
                        status_file = self.logs_dir / f"prefetch_{srr_id}.status"
                        status_file.write_text("DBGaP_REQUIRED")
                    else:
                        logger.error(self._c(f"    ✗ {srr_id}.sra download failed", self._C_RED))
                        
        except Exception as e:
            logger.error(f"Error downloading {srr_id}: {e}")
    
    def convert_sra_to_fastq(self):
        """Step 2: Convert SRA files to FASTQ format"""
        logger.info("=" * 50)
        logger.info(self._c("STEP 2: Convert SRA to FASTQ", self._C_CYAN))
        logger.info("=" * 50)
        
        # Path to the conversion script (submit_fastq_dump_jobs.sh)
        fdump_script = Path("/data/salomonis-archive/FASTQs/NCI-R01/POSEIDON/ValeriiGitRepo/scripts/manual_pipeline/submit_fastq_dump_jobs.sh")
        
        jobs_to_wait = {}
        total_srrs = 0
        skipped_fastq = 0
        skipped_dbgap = 0
        submitted = 0
        for sample_id, sample_data in self.samples.items():
            # Skip conversion entirely if BAM already exists for this sample
            if self._sample_has_bam(sample_id):
                logger.info(self._c(f"✓ {sample_id}: BAM exists; skipping SRA->FASTQ conversion", self._C_GREEN))
                continue
            
            for srr_id in sample_data['srr_ids']:
                total_srrs += 1
                sra_file = self.cancer_dir / f"{srr_id}.sra"
                srr_r1 = self.cancer_dir / f"{srr_id}_1.fastq.gz"
                srr_r2 = self.cancer_dir / f"{srr_id}_2.fastq.gz"
                
                # Skip if already converted or if marked as dbGaP
                status_file = self.logs_dir / f"prefetch_{srr_id}.status"
                if status_file.exists() and "DBGaP_REQUIRED" in status_file.read_text():
                    skipped_dbgap += 1
                    continue
                
                if srr_r1.exists() and srr_r2.exists():
                    skipped_fastq += 1
                    continue
                
                if sra_file.exists():
                    logger.info(self._c(f"→ Submitting conversion job for {srr_id}.sra", self._C_CYAN))
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
                            logger.info(self._c(f"  ✓ Submitted as Job <{job_id}>", self._C_GREEN))
                            submitted += 1
                        else:
                            logger.warning(self._c("  ! Could not parse Job ID; will proceed without waiting for this SRR", self._C_YELLOW))
                    except subprocess.CalledProcessError as e:
                        logger.error(self._c(f"Failed to submit conversion for {srr_id}: {e}", self._C_RED))
                # If neither FASTQs nor .sra exist, nothing to do for this SRR
        
        logger.info(self._c(f"Summary: SRRs={total_srrs}, skipped_fastq={skipped_fastq}, skipped_dbgap={skipped_dbgap}, submitted={submitted}", self._C_MAGENTA))
        # Optionally wait for jobs to finish
        if not self.no_wait and jobs_to_wait:
            self._wait_for_jobs_and_cleanup(jobs_to_wait)
            # After actual conversions, thorough cleanup using gzip test
            self.cleanup_converted_sras_global()
        else:
            logger.info(self._c("No conversions needed; skipping wait.", self._C_CYAN))
            # Fast cleanup when nothing was submitted (size-only check)
            self.cleanup_converted_sras_lightweight()
        # Refresh status after conversion stage (quiet)
        try:
            self.generate_status_report(log_header=False)
        except Exception as e:
            logger.warning(f"Failed to update status report: {e}")

        # If BAMs exist, perform final cleanup of SRAs and FASTQs for those samples
        try:
            self.cleanup_artifacts_for_completed_samples()
        except Exception as e:
            logger.warning(self._c(f"Final cleanup skipped due to error: {e}", self._C_YELLOW))

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
        """Validate gzip files quickly without long blocking.

        - Prefer `gzip -t` with a timeout.
        - Fallback: read a small chunk via Python gzip to check header.
        """
        try:
            cmd = ['gzip', '-t'] + [str(p) for p in paths]
            r = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=20)
            return r.returncode == 0
        except subprocess.TimeoutExpired:
            return False
        except Exception:
            try:
                for p in paths:
                    with gzip.open(p, 'rb') as f:
                        f.read(4096)
                return True
            except Exception:
                return False

    def _wait_for_jobs_and_cleanup(self, jobs_to_wait: dict):
        logger.info(self._c("Waiting for conversion jobs to finish...", self._C_CYAN))
        remaining = dict(jobs_to_wait)
        while remaining:
            finished = []
            for srr_id, job_id in remaining.items():
                # Fast path: if FASTQs are already present and valid, finish regardless of scheduler status
                srr_r1 = self.cancer_dir / f"{srr_id}_1.fastq.gz"
                srr_r2 = self.cancer_dir / f"{srr_id}_2.fastq.gz"
                sra_file = self.cancer_dir / f"{srr_id}.sra"
                if srr_r1.exists() and srr_r2.exists() and self._gzip_test(srr_r1, srr_r2):
                    if sra_file.exists():
                        with contextlib.suppress(Exception):
                            sra_file.unlink()
                            logger.info(self._c(f"  ✓ Cleaned {srr_id}.sra after successful conversion", self._C_GREEN))
                    finished.append(srr_id)
                    continue

                status = self._bjobs_status(job_id)
                # Treat common terminal/unknown states as finished
                if status in ('DONE', 'EXIT', 'ZOMBIE', 'ZOMBI', 'UNKNOWN', 'UNKWN'):
                    # Check outputs regardless of status
                    if srr_r1.exists() and srr_r2.exists() and self._gzip_test(srr_r1, srr_r2):
                        if sra_file.exists():
                            try:
                                sra_file.unlink()
                                logger.info(f"  ✓ Cleaned {srr_id}.sra after successful conversion")
                            except Exception as e:
                                logger.warning(f"  ! Could not remove {srr_id}.sra: {e}")
                        finished.append(srr_id)
                    else:
                        if status in ('DONE', 'EXIT'):
                            logger.warning(self._c(f"  ! {srr_id} job {status} but FASTQs missing/invalid; keeping .sra", self._C_YELLOW))
                        elif status in ('ZOMBIE', 'ZOMBI', 'UNKNOWN', 'UNKWN'):
                            logger.warning(self._c(f"  ! {srr_id} job {status}; FASTQs missing; treating as finished", self._C_YELLOW))
                        finished.append(srr_id)
                # else still running (e.g., RUN, PEND)
                else:
                    # Detect stuck jobs: if status stays non-terminal for many cycles, break out with warning
                    cycles = self.job_poll_cycles.get(srr_id, 0) + 1
                    self.job_poll_cycles[srr_id] = cycles
                    if cycles >= 30:  # ~30 minutes at default poll interval
                        logger.warning(self._c(f"  ! {srr_id} appears stuck in status {status} for {cycles} cycles; marking finished without FASTQs", self._C_YELLOW))
                        finished.append(srr_id)
            # Remove finished from remaining
            for s in finished:
                remaining.pop(s, None)
            # Periodic status refresh each poll cycle
            try:
                logger.info(self._c(f"Remaining conversion jobs: {len(remaining)}", self._C_CYAN))
                self.generate_status_report(log_header=False)
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
                        logger.info(self._c(f"  ✓ Removed empty directory {entry}", self._C_GREEN))
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
            logger.info(self._c(f"  ✓ Global cleanup removed {removed} converted .sra files", self._C_GREEN))

    def cleanup_artifacts_for_completed_samples(self):
        """If a sample has a BAM, remove its SRR .sra and FASTQs to free space.

        This is more aggressive than conversion cleanup and runs after STEP 2.
        """
        total_sra_removed = 0
        total_fastq_removed = 0
        for sample_id, sample_data in self.samples.items():
            if not self._sample_has_bam(sample_id):
                continue
            for srr_id in sample_data['srr_ids']:
                sra_file = self.cancer_dir / f"{srr_id}.sra"
                r1 = self.cancer_dir / f"{srr_id}_1.fastq.gz"
                r2 = self.cancer_dir / f"{srr_id}_2.fastq.gz"
                with contextlib.suppress(Exception):
                    if sra_file.exists():
                        sra_file.unlink()
                        total_sra_removed += 1
                with contextlib.suppress(Exception):
                    if r1.exists():
                        r1.unlink()
                        total_fastq_removed += 1
                with contextlib.suppress(Exception):
                    if r2.exists():
                        r2.unlink()
                        total_fastq_removed += 1
        if total_sra_removed or total_fastq_removed:
            logger.info(self._c(
                f"  ✓ Final cleanup (BAM present): removed {total_sra_removed} .sra and {total_fastq_removed} FASTQ files",
                self._C_GREEN
            ))

    def cleanup_converted_sras_lightweight(self):
        """Remove .sra files when both FASTQs exist and are non-empty (no gzip validation)."""
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
            try:
                if r1.exists() and r2.exists() and r1.stat().st_size > 0 and r2.stat().st_size > 0:
                    with contextlib.suppress(Exception):
                        sra_file.unlink()
                        removed += 1
            except Exception:
                continue
        if removed:
            logger.info(self._c(f"  ✓ Lightweight cleanup removed {removed} converted .sra files", self._C_GREEN))
        else:
            logger.info(self._c("  ✓ Lightweight cleanup: no .sra files to remove", self._C_GREEN))
    
    def generate_status_report(self, log_header: bool = True):
        """Generate sample status report (4 columns, exact format)."""
        if log_header:
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
            
            # Fast checks (existence and size only) for status reporting
            all_fastqs_ok = True
            any_fastq_missing = False
            any_sra_present_missing_fastq = False
            any_missing_both = False
            for sid in srr_ids:
                r1 = self.cancer_dir / f"{sid}_1.fastq.gz"
                r2 = self.cancer_dir / f"{sid}_2.fastq.gz"
                sra = self.cancer_dir / f"{sid}.sra"
                try:
                    r1_ok = r1.exists() and r1.stat().st_size > 0
                    r2_ok = r2.exists() and r2.stat().st_size > 0
                except Exception:
                    r1_ok = r2_ok = False
                if not (r1_ok and r2_ok):
                    all_fastqs_ok = False
                    any_fastq_missing = True
                    if sra.exists():
                        any_sra_present_missing_fastq = True
                    else:
                        any_missing_both = True

            # 2) BAM_DONE (only if all fastqs are good AND bam exists)
            if status is None and all_fastqs_ok:
                if self._sample_has_bam(sample_id):
                    status = 'BAM_DONE'

            # 3) ALIGN_IN_PROGRESS if STAR work directories exist
            if status is None and all_fastqs_ok:
                if self._has_star_progress_dirs(sample_id):
                    status = 'ALIGN_IN_PROGRESS'
                else:
                    # 4) FASTQ_DONE
                    status = 'FASTQ_DONE'
  
            # 5) CONVERTING (SRA -> FASTQ, job ID <...>) if any active conversion job
            if status is None:
                active_job_id = None
                for sid in srr_ids:
                    job_id = self.submitted_jobs.get(sid)
                    if job_id and self._is_job_active(job_id):
                        active_job_id = job_id
                        break
                if active_job_id:
                    status = f"CONVERTING (SRA -> FASTQ, job ID <{active_job_id}>)"

            # 6) NEEDS_CONVERSION if any .sra exists and some FASTQs missing
            if status is None and any_sra_present_missing_fastq:
                status = 'NEEDS_CONVERSION'

            # 7) NEEDS_PREFETCH if no .sra and some FASTQs missing
            if status is None and any_missing_both and any_fastq_missing:
                status = 'NEEDS_PREFETCH'

            # 8) UNKNOWN fallback (should be rare)
            if status is None:
                status = 'UNKNOWN'

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
        # Refresh status after download stage (quiet)
        try:
            self.generate_status_report(log_header=False)
        except Exception:
            pass
        self.convert_sra_to_fastq()
        self.generate_status_report()
        
        logger.info("=" * 50)
        # Minimal colorized completion message for better skimmability
        try:
            class _Ansi:
                GREEN = "\033[92m"
                RESET = "\033[0m"
            def _color(txt: str, code: str) -> str:
                try:
                    return f"{code}{txt}{_Ansi.RESET}" if sys.stdout.isatty() else txt
                except Exception:
                    return txt
            logger.info(_color("WORKFLOW COMPLETE", _Ansi.GREEN))
        except Exception:
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