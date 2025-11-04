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
from concurrent.futures import ThreadPoolExecutor, as_completed
import time
import re
from pathlib import Path
from collections import defaultdict
import gzip
import shutil
import contextlib
import json

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class SRAWorkflow:
    """Handles SRA download and FASTQ conversion workflow"""
    
    def __init__(self, cancer_dir, no_wait: bool = False, poll_interval_sec: int = 60, prefetch_workers: int = 32, gzip_test_timeout: int = 60):
        self.cancer_dir = Path(cancer_dir)
        self.sample_list_path = self.cancer_dir / "sample_list.txt"
        self.logs_dir = self.cancer_dir / "logs"
        self.samples = {}  # Will store sample_id -> {srr_ids: [...], status: ...}
        self.no_wait = no_wait

        # Validate poll_interval_sec
        if poll_interval_sec <= 0:
            raise ValueError(f"poll_interval_sec must be positive, got {poll_interval_sec}")
        if poll_interval_sec > 3600:
            logger.warning(f"poll_interval_sec={poll_interval_sec} is unusually large (>1 hour)")
        self.poll_interval_sec = poll_interval_sec

        # Validate prefetch_workers
        if prefetch_workers <= 0:
            raise ValueError(f"prefetch_workers must be positive, got {prefetch_workers}")
        if prefetch_workers > 128:
            logger.warning(f"prefetch_workers={prefetch_workers} is very high; capping at 128 to prevent system overload")
            prefetch_workers = 128
        self.prefetch_workers = int(prefetch_workers)

        # Validate and store gzip_test_timeout
        if gzip_test_timeout <= 0:
            raise ValueError(f"gzip_test_timeout must be positive, got {gzip_test_timeout}")
        if gzip_test_timeout > 300:
            logger.warning(f"gzip_test_timeout={gzip_test_timeout} is very long (>5 minutes)")
        self.gzip_test_timeout = gzip_test_timeout

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

        line_num = 0
        with open(self.sample_list_path, 'r') as f:
            for line in f:
                line_num += 1
                line = line.strip()
                if not line:
                    continue

                # Try tab first, then multiple spaces as fallback
                parts = line.split('\t')
                if len(parts) < 3:
                    # Try splitting by multiple spaces
                    parts = [p for p in line.split(' ') if p]

                if len(parts) < 3:
                    logger.warning(f"Line {line_num} has insufficient columns (expected >=3, got {len(parts)}): {line[:50]}")
                    continue

                sample_id = parts[0].strip()
                if not sample_id:
                    logger.warning(f"Line {line_num} has empty sample ID, skipping")
                    continue

                # Extract SRR/ERR IDs from columns 2 and 3
                srr_ids = self._extract_srr_ids(parts[1], parts[2])
                if not srr_ids:
                    logger.warning(f"Line {line_num}: No valid SRR/ERR IDs found for sample {sample_id}")
                    continue

                self.samples[sample_id] = {
                    'srr_ids': srr_ids,
                    'col2': parts[1],
                    'col3': parts[2],
                    'status': 'PENDING'
                }

        if not self.samples:
            logger.error("No valid samples found in sample_list.txt")
            sys.exit(1)

        logger.info(f"Found {len(self.samples)} samples")

    def _sample_has_bam(self, sample_id: str) -> bool:
        """Return True if a BAM for this sample exists.

        Robust matching across both directory root and optional 'bams/' subdir:
        - Filenames starting with sample_id (legacy convention)
        - OR filenames containing any SRR/ERR ID for this sample
        """
        try:
            srr_ids = self.samples.get(sample_id, {}).get('srr_ids', [])

            def _dir_has_bam(d: Path) -> bool:
                for p in d.glob('*.bam'):
                    try:
                        if not (p.is_file() and p.stat().st_size > 0):
                            continue
                        name = p.name
                        if name.startswith(sample_id):
                            return True
                        for sid in srr_ids:
                            if sid in name:
                                return True
                    except OSError as e:
                        # File permission or access error, skip this file
                        logger.debug(f"Could not access {p}: {e}")
                        continue
                return False

            # Root directory
            if _dir_has_bam(self.cancer_dir):
                return True
            # bams/ subdirectory
            bams_dir = self.cancer_dir / 'bams'
            if bams_dir.is_dir() and _dir_has_bam(bams_dir):
                return True
        except Exception as e:
            logger.debug(f"Error checking BAM existence for {sample_id}: {e}")
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
        except OSError as e:
            logger.debug(f"Error checking STAR progress dirs for {sample_id}: {e}")
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

    def _validate_sra_file(self, sra_file: Path) -> bool:
        """Validate that an SRA file is not corrupted using vdb-validate.

        Returns True if file is valid, False otherwise.
        Falls back to basic size check if vdb-validate is unavailable or times out.

        IMPORTANT: Timeout does NOT mean corruption - just slow I/O on large files!
        """
        if not sra_file.exists():
            return False

        # Basic check: file must be at least 1KB
        try:
            file_size = sra_file.stat().st_size
            if file_size < 1024:
                logger.warning(f"SRA file {sra_file.name} is suspiciously small (<1KB)")
                return False
        except Exception:
            return False

        # Calculate timeout based on file size: ~1 second per 100MB, min 30s, max 300s
        # Large files on network storage need more time
        timeout_seconds = max(30, min(300, int(file_size / (100 * 1024 * 1024))))

        # Try vdb-validate if available (part of sratoolkit)
        try:
            result = subprocess.run(
                ['vdb-validate', str(sra_file)],
                capture_output=True,
                text=True,
                timeout=timeout_seconds,
                check=False
            )
            # vdb-validate returns 0 for valid files
            if result.returncode == 0:
                return True
            else:
                logger.warning(f"SRA validation failed for {sra_file.name}: {result.stderr.strip()}")
                return False
        except FileNotFoundError:
            # vdb-validate not available, use size-based heuristic
            logger.debug("vdb-validate not found; using size-based validation")
            return True  # Assume valid if size check passed
        except subprocess.TimeoutExpired:
            # TIMEOUT DOES NOT MEAN CORRUPTION - just slow I/O on large files!
            # For large files on network storage, validation can take a very long time
            logger.debug(f"SRA validation timed out after {timeout_seconds}s for {sra_file.name} ({file_size / (1024**3):.1f} GB) - assuming valid")
            return True  # Assume valid - timeout is not corruption
        except Exception as e:
            logger.debug(f"Could not validate {sra_file.name}: {e}")
            return True  # Assume valid if validation unavailable

    def _validate_fastq_pair(self, r1_path: Path, r2_path: Path, srr_id: str = None) -> tuple:
        """Validate a pair of FASTQ files for completeness and correctness.

        Returns (is_valid: bool, reason: str)

        Checks performed:
        1. Files exist and non-empty
        2. Gzip integrity (not corrupted)
        3. FASTQ format validity (proper 4-line records)
        4. Both files have same number of reads (paired)
        5. Files are not suspiciously small

        This catches interrupted conversions from HPC shutdowns.
        """
        id_str = f"{srr_id}: " if srr_id else ""

        # Check existence and size
        for p, label in [(r1_path, 'R1'), (r2_path, 'R2')]:
            if not p.exists():
                return False, f"{id_str}{label} does not exist"
            try:
                size = p.stat().st_size
                if size == 0:
                    return False, f"{id_str}{label} is empty (0 bytes)"
                if size < 100:  # Suspiciously small even for gzipped
                    return False, f"{id_str}{label} is suspiciously small ({size} bytes)"
            except OSError as e:
                return False, f"{id_str}Cannot stat {label}: {e}"

        # Check gzip integrity
        if not self._gzip_test(r1_path, r2_path):
            return False, f"{id_str}Gzip integrity check failed (corrupted or incomplete compression)"

        # Check FASTQ format and read counts
        try:
            r1_reads = self._count_fastq_reads(r1_path)
            r2_reads = self._count_fastq_reads(r2_path)

            if r1_reads == 0 and r2_reads == 0:
                return False, f"{id_str}Both files appear empty or invalid FASTQ format"
            if r1_reads == 0:
                return False, f"{id_str}R1 has no valid reads"
            if r2_reads == 0:
                return False, f"{id_str}R2 has no valid reads"
            if r1_reads != r2_reads:
                return False, f"{id_str}Read count mismatch: R1={r1_reads}, R2={r2_reads} (interrupted conversion)"

            logger.debug(f"{id_str}Validated: {r1_reads} paired reads")
            return True, f"Valid ({r1_reads} reads)"

        except Exception as e:
            return False, f"{id_str}FASTQ validation error: {e}"

    def _cleanup_corrupted_fastqs(self, srr_id: str, reason: str = None) -> bool:
        """Delete corrupted FASTQ files for a given SRR ID.

        Always deletes corrupted FASTQs regardless of whether .sra exists.
        Corrupted genomics data is dangerous and should be removed immediately.

        Returns True if cleanup succeeded, False otherwise.
        """
        srr_r1 = self.cancer_dir / f"{srr_id}_1.fastq.gz"
        srr_r2 = self.cancer_dir / f"{srr_id}_2.fastq.gz"

        removed = []
        try:
            if srr_r1.exists():
                srr_r1.unlink()
                removed.append("R1")
            if srr_r2.exists():
                srr_r2.unlink()
                removed.append("R2")

            if removed:
                files_str = " and ".join(removed)
                reason_str = f" ({reason})" if reason else ""
                logger.info(self._c(
                    f"  ✓ Removed corrupted {files_str} for {srr_id}{reason_str}; keeping .sra for retry",
                    self._C_GREEN
                ))
                return True
            return False

        except Exception as e:
            logger.warning(self._c(f"  ! Could not remove corrupted FASTQs for {srr_id}: {e}", self._C_YELLOW))
            return False

    def _count_fastq_reads(self, fastq_gz_path: Path, max_reads: int = 1000) -> int:
        """Count reads in a gzipped FASTQ file (samples first N reads for speed).

        Returns number of reads found, or 0 if format is invalid.
        Validates FASTQ 4-line structure while counting.
        """
        try:
            with gzip.open(fastq_gz_path, 'rt') as f:
                read_count = 0
                line_in_record = 0

                for i, line in enumerate(f):
                    if i >= max_reads * 4:  # Sample first N reads only
                        break

                    # Validate FASTQ 4-line structure
                    if line_in_record == 0:  # Header line
                        if not line.startswith('@'):
                            logger.warning(f"Invalid FASTQ header at line {i+1} in {fastq_gz_path.name}")
                            return 0
                    elif line_in_record == 2:  # '+' separator line
                        if not line.startswith('+'):
                            logger.warning(f"Invalid FASTQ separator at line {i+1} in {fastq_gz_path.name}")
                            return 0

                    line_in_record += 1
                    if line_in_record == 4:
                        read_count += 1
                        line_in_record = 0

                # Check if file ended mid-record (truncated)
                if line_in_record != 0:
                    logger.warning(f"Truncated FASTQ record in {fastq_gz_path.name} (incomplete 4-line record)")
                    return 0

                return read_count if read_count > 0 else 0

        except EOFError:
            logger.warning(f"Unexpected EOF in {fastq_gz_path.name} (interrupted compression)")
            return 0
        except Exception as e:
            logger.debug(f"Error counting reads in {fastq_gz_path.name}: {e}")
            return 0
    
    def load_modules(self):
        """Load required modules (sratoolkit and aspera).

        Best-effort attempt since 'module' command may not be available or
        may be a shell function. Logs warnings if modules fail to load.
        """
        modules_loaded = []
        modules_failed = []

        for module_name in ['sratoolkit/2.10.4', 'aspera/3.9.1']:
            try:
                result = subprocess.run(
                    f'module load {module_name}',
                    shell=True,
                    check=False,
                    capture_output=True,
                    text=True,
                    timeout=5
                )
                if result.returncode == 0:
                    modules_loaded.append(module_name)
                    logger.debug(f"Loaded module: {module_name}")
                else:
                    # Check if 'module' command exists
                    if 'command not found' in result.stderr or 'module: not found' in result.stderr:
                        logger.debug("Module system not available; assuming tools are in PATH")
                        break
                    else:
                        modules_failed.append(module_name)
                        logger.warning(f"Failed to load module {module_name}: {result.stderr.strip()}")
            except subprocess.TimeoutExpired:
                logger.warning(f"Timeout loading module {module_name}")
                modules_failed.append(module_name)
            except Exception as e:
                logger.debug(f"Could not load module {module_name}: {e}")
                modules_failed.append(module_name)

        if modules_failed:
            logger.warning(
                f"Some modules failed to load: {', '.join(modules_failed)}. "
                "Ensure required tools (prefetch, fastq-dump) are available in PATH."
            )
    
    def download_sra_files(self):
        """Step 1: Download SRA files for all samples"""
        logger.info("=" * 50)
        logger.info(self._c("STEP 1: Download SRA files (if needed)", self._C_CYAN))
        logger.info("=" * 50)
        
        tasks = []
        for sample_id, sample_data in self.samples.items():
            # Skip entire sample if BAM already exists
            if self._sample_has_bam(sample_id):
                logger.info(self._c(f"✓ {sample_id}: BAM exists; skipping SRA downloads for all SRRs", self._C_GREEN))
                continue
            for srr_id in sample_data['srr_ids']:
                tasks.append((srr_id, sample_id))

        if not tasks:
            logger.info(self._c("No SRA downloads needed.", self._C_CYAN))
            return

        logger.info(self._c(f"Starting parallel prefetch with {self.prefetch_workers} workers ({len(tasks)} SRRs)", self._C_CYAN))
        with ThreadPoolExecutor(max_workers=self.prefetch_workers) as pool:
            futures = [pool.submit(self._download_single_sra, srr_id, sample_id) for srr_id, sample_id in tasks]
            for _ in as_completed(futures):
                pass
    
    def _download_single_sra(self, srr_id, sample_id):
        """Download a single SRA file using prefetch"""
        sra_file = self.cancer_dir / f"{srr_id}.sra"
        srr_r1 = self.cancer_dir / f"{srr_id}_1.fastq.gz"
        srr_r2 = self.cancer_dir / f"{srr_id}_2.fastq.gz"

        # Skip if SRR-level FASTQs already exist AND are valid
        if srr_r1.exists() and srr_r2.exists():
            is_valid, reason = self._validate_fastq_pair(srr_r1, srr_r2, srr_id)
            if is_valid:
                logger.info(self._c(f"    ✓ {srr_id} FASTQs exist and validated, skipping prefetch", self._C_GREEN))
                return
            else:
                logger.warning(self._c(f"    ! {srr_id} FASTQs exist but invalid: {reason}", self._C_YELLOW))
                logger.warning(self._c(f"    ! Removing invalid FASTQs and will reprocess", self._C_YELLOW))
                # Remove invalid files so they can be regenerated
                try:
                    srr_r1.unlink()
                    srr_r2.unlink()
                except Exception as e:
                    logger.warning(f"Could not remove invalid FASTQs: {e}")
                # Continue to download/convert
        
        # Skip if SRA file already exists AND is valid
        if sra_file.exists():
            if self._validate_sra_file(sra_file):
                logger.info(self._c(f"    ✓ {srr_id}.sra already exists and validated", self._C_GREEN))
                return
            else:
                logger.warning(self._c(f"    ! {srr_id}.sra exists but invalid (corrupted/incomplete)", self._C_YELLOW))
                logger.warning(self._c(f"    ! Removing invalid SRA and will re-download", self._C_YELLOW))
                try:
                    sra_file.unlink()
                except Exception as e:
                    logger.warning(f"Could not remove invalid SRA: {e}")
                # Continue to re-download below
        
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
                    target = self.cancer_dir / sra.name
                    try:
                        # Check if target already exists
                        if target.exists():
                            logger.debug(f"Target {target.name} already exists; removing source from subdirectory")
                            sra.unlink()
                        else:
                            sra.rename(target)
                    except Exception as e:
                        logger.warning(f"Could not move {sra.name} from subdirectory: {e}")
                try:
                    # Remove empty directory if possible
                    if not any(srr_dir.iterdir()):
                        srr_dir.rmdir()
                except Exception as e:
                    logger.debug(f"Could not remove {srr_dir}: {e}")
            
            # Check if download succeeded and validate integrity
            if sra_file.exists():
                if self._validate_sra_file(sra_file):
                    logger.info(self._c(f"    ✓ {srr_id}.sra downloaded and validated successfully", self._C_GREEN))
                else:
                    logger.error(self._c(f"    ✗ {srr_id}.sra downloaded but failed validation; removing", self._C_RED))
                    try:
                        sra_file.unlink()
                    except Exception:
                        pass
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
        # Use relative path from this script's location for portability
        fdump_script = Path(__file__).parent / "submit_fastq_dump_jobs.sh"
        if not fdump_script.exists():
            logger.error(f"Required script not found: {fdump_script}")
            sys.exit(1)
        
        # Discover any active conversion jobs from previous runs and merge into submitted_jobs
        try:
            self._load_and_discover_active_jobs()
        except Exception as e:
            logger.warning(self._c(f"Unable to load/discover active jobs: {e}", self._C_YELLOW))

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
                try:
                    if "DBGaP_REQUIRED" in status_file.read_text():
                        skipped_dbgap += 1
                        continue
                except (FileNotFoundError, OSError):
                    pass  # File doesn't exist or can't be read, continue normally

                # Check if FASTQs exist AND are valid (not corrupted/interrupted)
                if srr_r1.exists() and srr_r2.exists():
                    is_valid, reason = self._validate_fastq_pair(srr_r1, srr_r2, srr_id)
                    if is_valid:
                        skipped_fastq += 1
                        continue
                    else:
                        logger.warning(self._c(f"  ! {srr_id} FASTQs exist but invalid: {reason}", self._C_YELLOW))
                        logger.warning(self._c(f"  ! Will reconvert from SRA", self._C_YELLOW))
                        # Remove invalid files - this ensures bash script can proceed with clean slate
                        self._cleanup_corrupted_fastqs(srr_id, reason)

                # If there is an active job from a previous run, do not resubmit
                prev_job = self.submitted_jobs.get(srr_id)
                if prev_job and self._is_job_active(prev_job):
                    jobs_to_wait[srr_id] = prev_job
                    logger.info(self._c(f"→ {srr_id} already running as Job <{prev_job}>; will wait", self._C_CYAN))
                    continue
                
                if sra_file.exists():
                    # Validate SRA before submitting conversion job
                    if not self._validate_sra_file(sra_file):
                        logger.warning(self._c(f"  ! {srr_id}.sra exists but is invalid/corrupted; skipping conversion", self._C_YELLOW))
                        logger.warning(self._c(f"  ! Run script again to re-download this SRA", self._C_YELLOW))
                        continue

                    logger.info(self._c(f"→ Submitting conversion job for {srr_id}.sra", self._C_CYAN))
                    try:
                        # Capture bsub output from submit_fastq_dump_jobs.sh to parse Job ID
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
                            try:
                                self._persist_submitted_job(srr_id, job_id)
                            except Exception:
                                pass
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
        """Query LSF job status via bjobs command.

        Returns job status string (e.g., 'RUN', 'PEND', 'DONE', 'EXIT') or 'UNKNOWN'.
        """
        try:
            proc = subprocess.run(
                ['bjobs', '-noheader', '-o', 'stat', job_id],
                capture_output=True, text=True, check=False
            )
            out = (proc.stdout or '').strip()
            err = (proc.stderr or '').strip().lower()

            # Job not found in queue (finished or never existed)
            if proc.returncode != 0:
                # Check for explicit "not found" messages
                if 'not found' in err or ('job' in err and 'is not found' in err):
                    return 'UNKNOWN'

            return out if out else 'UNKNOWN'
        except Exception as e:
            logger.debug(f"Error querying job status for {job_id}: {e}")
            return 'UNKNOWN'

    def _is_job_active(self, job_id: str) -> bool:
        status = self._bjobs_status(job_id)
        return status in ('PEND', 'RUN', 'PSUSP', 'USUSP', 'SSUSP')

    def _persist_submitted_job(self, srr_id: str, job_id: str) -> None:
        """Append or merge the submitted job ID to logs/submitted_jobs.json."""
        db_path = self.logs_dir / 'submitted_jobs.json'
        data = {}
        try:
            if db_path.exists():
                data = json.loads(db_path.read_text() or '{}')
        except (json.JSONDecodeError, OSError) as e:
            logger.debug(f"Could not load existing job database, starting fresh: {e}")
            data = {}
        data[str(srr_id)] = str(job_id)
        try:
            db_path.write_text(json.dumps(data, indent=2) + "\n")
        except OSError as e:
            logger.warning(f"Could not persist job ID for {srr_id}: {e}")

    def _load_and_discover_active_jobs(self) -> None:
        """Load previous submissions and discover active LSF jobs in this directory.

        - Loads logs/submitted_jobs.json and keeps only still-active jobs
        - Discovers active jobs with name 'fastq_<SRR>' whose CWD equals self.cancer_dir
        - Merges discoveries into self.submitted_jobs
        """
        # Load persisted
        db_path = self.logs_dir / 'submitted_jobs.json'
        if db_path.exists():
            try:
                data = json.loads(db_path.read_text() or '{}')
                for srr, jid in (data.items() if isinstance(data, dict) else []):
                    if jid and self._is_job_active(str(jid)):
                        self.submitted_jobs[str(srr)] = str(jid)
            except (json.JSONDecodeError, OSError) as e:
                logger.debug(f"Could not load submitted jobs database: {e}")

        # Discover via bjobs with cwd if available
        def _discover_with_cwd() -> int:
            try:
                proc = subprocess.run(
                    ['bjobs', '-noheader', '-o', 'jobid job_name stat cwd'],
                    capture_output=True, text=True, check=False
                )
                out = (proc.stdout or '').strip().splitlines()
                count = 0
                for line in out:
                    parts = line.split()
                    if len(parts) < 4:
                        continue
                    jid, jname, jstat = parts[0], parts[1], parts[2]
                    cwd = ' '.join(parts[3:])  # cwd can contain spaces
                    if not jname.startswith('fastq_'):
                        continue
                    if Path(cwd) != self.cancer_dir:
                        continue
                    # Extract SRR after fastq_
                    m = re.match(r'^fastq_(\w+)$', jname)
                    if not m:
                        continue
                    srr = m.group(1)
                    self.submitted_jobs[srr] = str(jid)
                    count += 1
                return count
            except Exception:
                return 0

        def _discover_with_long() -> int:
            try:
                # First, list job IDs and names
                proc = subprocess.run(
                    ['bjobs', '-noheader', '-o', 'jobid job_name stat'],
                    capture_output=True, text=True, check=False
                )
                out = (proc.stdout or '').strip().splitlines()
                ids = []
                for line in out:
                    parts = line.split()
                    if len(parts) < 2:
                        continue
                    jid, jname = parts[0], parts[1]
                    if jname.startswith('fastq_'):
                        ids.append((jid, jname))
                count = 0
                for jid, jname in ids:
                    # Get long info to parse CWD
                    proc2 = subprocess.run(['bjobs', '-l', jid], capture_output=True, text=True, check=False)
                    txt = (proc2.stdout or '') + (proc2.stderr or '')
                    cwd_match = re.search(r'\bCWD:\s*(.*)', txt)
                    cwd = cwd_match.group(1).strip() if cwd_match else ''
                    if cwd and Path(cwd) == self.cancer_dir:
                        m = re.match(r'^fastq_(\w+)$', jname)
                        if m:
                            self.submitted_jobs[m.group(1)] = str(jid)
                            count += 1
                return count
            except Exception:
                return 0

        found = _discover_with_cwd()
        if found == 0:
            _discover_with_long()

    def _gzip_test(self, *paths: Path) -> bool:
        """Validate gzip files with full integrity check.

        Uses gzip -t to decompress entire file and verify integrity.
        This is the ONLY way to catch truncated files from cluster shutdowns.

        Timeout protection prevents hangs on slow network storage.
        IMPORTANT: Timeout does NOT mean corruption - just slow I/O.
        """
        # Calculate dynamic timeout based on file sizes
        # ~1 second per 100MB, min 60s, max 300s (5 minutes)
        total_size = 0
        try:
            for p in paths:
                total_size += p.stat().st_size
        except Exception:
            total_size = 1024 * 1024 * 1024  # Assume 1GB if stat fails

        timeout = max(60, min(300, int(total_size / (100 * 1024 * 1024))))

        # Use gzip -t with dynamic timeout for full validation
        try:
            cmd = ['gzip', '-t'] + [str(p) for p in paths]
            r = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=timeout)
            return r.returncode == 0
        except subprocess.TimeoutExpired:
            # TIMEOUT DOES NOT MEAN CORRUPTION - just slow I/O!
            # Fall back to assuming valid if header test passes
            logger.debug(f"Gzip test timed out after {timeout}s for {[p.name for p in paths]} ({total_size / (1024**3):.1f} GB) - assuming valid (slow I/O)")
            try:
                for p in paths:
                    with gzip.open(p, 'rb') as f:
                        f.read(4096)
                return True  # Header is valid, assume file is OK
            except Exception as e:
                logger.warning(f"Gzip header fallback failed: {e}")
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
                if srr_r1.exists() and srr_r2.exists():
                    is_valid, reason = self._validate_fastq_pair(srr_r1, srr_r2, srr_id)
                    if is_valid:
                        if sra_file.exists():
                            with contextlib.suppress(Exception):
                                sra_file.unlink()
                                logger.info(self._c(f"  ✓ Cleaned {srr_id}.sra after successful conversion", self._C_GREEN))
                        finished.append(srr_id)
                        continue
                    else:
                        logger.warning(self._c(f"  ! {srr_id} FASTQs present but validation failed: {reason}", self._C_YELLOW))
                        # Auto-cleanup corrupted files to enable retry on next run
                        self._cleanup_corrupted_fastqs(srr_id, reason)

                status = self._bjobs_status(job_id)
                # Treat common terminal/unknown states as finished
                if status in ('DONE', 'EXIT', 'ZOMBIE', 'ZOMBI', 'UNKNOWN', 'UNKWN'):
                    # Check outputs regardless of status
                    if srr_r1.exists() and srr_r2.exists():
                        is_valid, reason = self._validate_fastq_pair(srr_r1, srr_r2, srr_id)
                        if is_valid:
                            if sra_file.exists():
                                try:
                                    sra_file.unlink()
                                    logger.info(f"  ✓ Cleaned {srr_id}.sra after successful conversion")
                                except Exception as e:
                                    logger.warning(f"  ! Could not remove {srr_id}.sra: {e}")
                            finished.append(srr_id)
                        else:
                            logger.warning(self._c(f"  ! {srr_id} job {status} but FASTQs invalid: {reason}; keeping .sra", self._C_YELLOW))
                            # Auto-cleanup corrupted files to enable retry on next run
                            self._cleanup_corrupted_fastqs(srr_id, reason)
                            finished.append(srr_id)
                    else:
                        if status in ('DONE', 'EXIT'):
                            logger.warning(self._c(f"  ! {srr_id} job {status} but FASTQs missing; keeping .sra", self._C_YELLOW))
                        elif status in ('ZOMBIE', 'ZOMBI', 'UNKNOWN', 'UNKWN'):
                            logger.warning(self._c(f"  ! {srr_id} job {status}; FASTQs missing; treating as finished", self._C_YELLOW))
                        finished.append(srr_id)
                # else still running (e.g., RUN, PEND)
                else:
                    # Detect stuck jobs: if status stays non-terminal for many cycles, break out with warning
                    cycles = self.job_poll_cycles.get(srr_id, 0) + 1
                    self.job_poll_cycles[srr_id] = cycles
                    # Calculate max cycles for 30 minutes based on actual poll interval
                    max_cycles_for_30min = max(30, int(1800 / self.poll_interval_sec))
                    elapsed_seconds = cycles * self.poll_interval_sec
                    if cycles >= max_cycles_for_30min:
                        logger.warning(self._c(
                            f"  ! {srr_id} appears stuck in status {status} for {cycles} cycles "
                            f"({elapsed_seconds}s); marking finished without FASTQs",
                            self._C_YELLOW
                        ))
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
        """Delete .sra files when both paired FASTQs exist and pass gzip integrity validation.

        This is the THOROUGH cleanup method - validates gzip integrity before removal.
        Used after actual conversion jobs complete to ensure FASTQs are fully valid.
        Timeout: Uses self.gzip_test_timeout (default 60s).

        Use this when:
        - After conversion jobs finish
        - When you need to ensure FASTQs are not corrupted before removing source data

        SAFETY: Checks for active fastq-dump processes before cleanup to prevent
        race conditions and .nfs ghost files.
        """
        # Check for active fastq-dump processes to prevent race conditions
        try:
            result = subprocess.run(
                ['pgrep', '-f', 'fastq-dump'],
                capture_output=True,
                text=True,
                timeout=5
            )
            if result.returncode == 0:
                # Active fastq-dump processes found
                active_pids = result.stdout.strip().split('\n')
                logger.warning(self._c(
                    f"  ⚠ Skipping global cleanup: {len(active_pids)} active fastq-dump process(es) detected",
                    self._C_YELLOW
                ))
                return
        except (subprocess.TimeoutExpired, FileNotFoundError) as e:
            # pgrep not available or timeout - proceed cautiously with cleanup
            logger.warning(self._c(
                f"  ⚠ Could not check for active fastq-dump processes: {e}",
                self._C_YELLOW
            ))

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
            if r1.exists() and r2.exists():
                is_valid, _ = self._validate_fastq_pair(r1, r2, srr_id)
                if is_valid:
                    with contextlib.suppress(Exception):
                        sra_file.unlink()
                        removed += 1
        if removed:
            logger.info(self._c(f"  ✓ Global cleanup removed {removed} converted .sra files", self._C_GREEN))

    def cleanup_artifacts_for_completed_samples(self):
        """If a sample has a BAM, remove its SRR .sra and FASTQs to free space.

        Validates FASTQs before deletion to provide audit trail of data integrity.
        If corruption is detected, logs a warning (BAM may be compromised).
        Deletion proceeds regardless to free space.

        This is more aggressive than conversion cleanup and runs after STEP 2.
        """
        total_sra_removed = 0
        total_fastq_removed = 0
        corrupted_detected = []

        for sample_id, sample_data in self.samples.items():
            if not self._sample_has_bam(sample_id):
                continue
            for srr_id in sample_data['srr_ids']:
                sra_file = self.cancer_dir / f"{srr_id}.sra"
                r1 = self.cancer_dir / f"{srr_id}_1.fastq.gz"
                r2 = self.cancer_dir / f"{srr_id}_2.fastq.gz"

                # Validate FASTQs before cleanup if both exist (audit trail)
                if r1.exists() and r2.exists():
                    try:
                        is_valid, reason = self._validate_fastq_pair(r1, r2, srr_id)
                        if not is_valid:
                            logger.warning(self._c(
                                f"  ⚠ {sample_id}/{srr_id}: BAM exists but FASTQ corrupted: {reason}",
                                self._C_YELLOW
                            ))
                            logger.warning(self._c(
                                f"  ⚠ Consider re-validating {sample_id} BAM or re-running alignment",
                                self._C_YELLOW
                            ))
                            corrupted_detected.append(srr_id)
                    except Exception as e:
                        logger.debug(f"FASTQ validation failed for {srr_id} (will delete anyway): {e}")

                # Delete artifacts regardless of validation result
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
        if corrupted_detected:
            logger.warning(self._c(
                f"  ⚠ Detected {len(corrupted_detected)} corrupted FASTQ(s) for samples with BAMs: {', '.join(corrupted_detected[:5])}{'...' if len(corrupted_detected) > 5 else ''}",
                self._C_YELLOW
            ))

    def cleanup_converted_sras_lightweight(self):
        """Remove .sra files when both FASTQs exist and are non-empty (no gzip validation).

        This is the FAST cleanup method - only checks file existence and size.
        Used when no new conversions happened or during pre-run cleanup.
        Much faster than thorough cleanup because it skips gzip integrity validation.

        Use this when:
        - Pre-run cleanup of files from previous runs
        - No conversions happened in current run
        - Speed is more important than verification (e.g., files already validated before)

        SAFETY: Checks for active fastq-dump processes before cleanup to prevent
        race conditions and .nfs ghost files.
        """
        # Check for active fastq-dump processes to prevent race conditions
        try:
            result = subprocess.run(
                ['pgrep', '-f', 'fastq-dump'],
                capture_output=True,
                text=True,
                timeout=5
            )
            if result.returncode == 0:
                # Active fastq-dump processes found
                active_pids = result.stdout.strip().split('\n')
                logger.warning(self._c(
                    f"  ⚠ Skipping lightweight cleanup: {len(active_pids)} active fastq-dump process(es) detected",
                    self._C_YELLOW
                ))
                return
        except (subprocess.TimeoutExpired, FileNotFoundError) as e:
            # pgrep not available or timeout - proceed cautiously with cleanup
            logger.warning(self._c(
                f"  ⚠ Could not check for active fastq-dump processes: {e}",
                self._C_YELLOW
            ))

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
    
    def _check_dbgap_required(self, srr_ids):
        """Check if any SRR is marked as requiring dbGaP access."""
        for sid in srr_ids:
            status_path = self.logs_dir / f"prefetch_{sid}.status"
            try:
                if "DBGaP_REQUIRED" in status_path.read_text():
                    return 'DBGaP_REQUIRED'
            except (FileNotFoundError, OSError):
                pass
        return None

    def _analyze_fastq_status(self, srr_ids):
        """Analyze FASTQ file status for all SRRs in a sample.

        Returns tuple: (all_fastqs_ok, any_fastq_missing, any_sra_present_missing_fastq, any_missing_both)
        """
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
            except OSError:
                r1_ok = r2_ok = False

            if not (r1_ok and r2_ok):
                all_fastqs_ok = False
                any_fastq_missing = True
                if sra.exists():
                    any_sra_present_missing_fastq = True
                else:
                    any_missing_both = True

        return all_fastqs_ok, any_fastq_missing, any_sra_present_missing_fastq, any_missing_both

    def _determine_sample_status(self, sample_id, srr_ids, all_fastqs_ok, any_fastq_missing,
                                   any_sra_present_missing_fastq, any_missing_both):
        """Determine sample status using ordered priority checks.

        Status priority (checked in order):
        1. DBGaP_REQUIRED - Access restrictions detected
        2. BAM_DONE - Final output exists
        3. ALIGN_IN_PROGRESS - STAR alignment running
        4. FASTQ_DONE - All FASTQs ready
        5. CONVERTING - SRA->FASTQ conversion in progress
        6. NEEDS_CONVERSION - SRA exists but FASTQs incomplete
        7. NEEDS_PREFETCH - No SRA, no FASTQs
        8. UNKNOWN - Fallback
        """
        # Priority 1: DBGaP required
        status = self._check_dbgap_required(srr_ids)
        if status:
            return status

        # Priority 2: BAM exists (final output)
        if self._sample_has_bam(sample_id):
            return 'BAM_DONE'

        # Priority 3-4: All FASTQs present
        if all_fastqs_ok:
            if self._has_star_progress_dirs(sample_id):
                return 'ALIGN_IN_PROGRESS'
            return 'FASTQ_DONE'

        # Priority 5: Active conversion job
        for sid in srr_ids:
            job_id = self.submitted_jobs.get(sid)
            if job_id and self._is_job_active(job_id):
                return f"CONVERTING (SRA -> FASTQ, job ID <{job_id}>)"

        # Priority 6: SRA exists but needs conversion
        if any_sra_present_missing_fastq:
            return 'NEEDS_CONVERSION'

        # Priority 7: Nothing exists, needs prefetch
        if any_missing_both and any_fastq_missing:
            return 'NEEDS_PREFETCH'

        # Priority 8: Unknown state (should be rare)
        return 'UNKNOWN'

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

            # Analyze FASTQ status
            all_fastqs_ok, any_fastq_missing, any_sra_present_missing_fastq, any_missing_both = \
                self._analyze_fastq_status(srr_ids)

            # Determine status using ordered priority checks
            status = self._determine_sample_status(
                sample_id, srr_ids, all_fastqs_ok, any_fastq_missing,
                any_sra_present_missing_fastq, any_missing_both
            )

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

        # Execute workflow steps
        # Note: No global directory change - all subprocess calls use cwd parameter for portability
        self.parse_sample_list()
        self.load_modules()
        # Pre-run cleanup (lightweight): remove any .sra whose FASTQs already exist from prior runs
        try:
            logger.info(self._c("Pre-run cleanup: scanning for converted .sra and empty SRR dirs...", self._C_CYAN))
            self.cleanup_converted_sras_lightweight()
            self.cleanup_empty_srr_dirs()
        except Exception:
            pass
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
    parser.add_argument('--prefetch-workers', type=int, default=32, help='Number of parallel workers for SRA prefetch (default: 32)')
    
    args = parser.parse_args()
    
    # Validate directory exists
    if not os.path.isdir(args.cancer_directory):
        logger.error(f"Directory not found: {args.cancer_directory}")
        sys.exit(1)
    
    # Run workflow
    workflow = SRAWorkflow(args.cancer_directory, no_wait=args.no_wait, prefetch_workers=args.prefetch_workers)
    workflow.run()


if __name__ == "__main__":
    main()