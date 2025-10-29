#!/usr/bin/env python3
"""
FASTQ Download Script by Category
Downloads FASTQ files for SRR IDs organized by classification categories
"""

import os
import sys
import subprocess
import argparse
import pandas as pd
import time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('fastq_download.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

class FASTQDownloader:
    def __init__(self, base_dir="fastq_downloads", max_workers=4):
        self.base_dir = Path(base_dir)
        self.max_workers = max_workers
        self.categories = ["Tumor", "Pre-malignant", "Normal", "Cell_line", "Unknown"]
        
        # Create base directory structure
        self.create_directory_structure()
        
    def create_directory_structure(self):
        """Create organized directory structure for FASTQ files"""
        logger.info(f"Creating directory structure in {self.base_dir}")
        
        for category in self.categories:
            category_dir = self.base_dir / category
            category_dir.mkdir(parents=True, exist_ok=True)
            
            # Create subdirectories for each category
            (category_dir / "raw").mkdir(exist_ok=True)
            (category_dir / "processed").mkdir(exist_ok=True)
            (category_dir / "logs").mkdir(exist_ok=True)
            
        logger.info("Directory structure created successfully")
    
    def download_srr(self, srr_id, category, retry_count=3):
        """Download FASTQ files for a single SRR ID"""
        category_dir = self.base_dir / category / "raw"
        log_file = self.base_dir / category / "logs" / f"{srr_id}.log"
        
        for attempt in range(retry_count):
            try:
                logger.info(f"Downloading {srr_id} to {category} (attempt {attempt + 1})")
                
                # Use prefetch to download SRA file first
                prefetch_cmd = [
                    "prefetch",
                    "--output-directory", str(category_dir),
                    "--progress",
                    srr_id
                ]
                
                with open(log_file, 'a') as log:
                    log.write(f"Attempt {attempt + 1}: Starting prefetch for {srr_id}\n")
                    result = subprocess.run(prefetch_cmd, capture_output=True, text=True, timeout=1800)
                    
                    if result.returncode != 0:
                        log.write(f"Prefetch failed: {result.stderr}\n")
                        # Check if it's a dbGaP authorization issue
                        if "dbGaP" in result.stderr or "Access denied" in result.stderr:
                            log.write(f"dbGaP authorization required for {srr_id}\n")
                            logger.warning(f"dbGaP authorization required for {srr_id} - skipping")
                            return False
                        raise subprocess.CalledProcessError(result.returncode, prefetch_cmd, result.stderr)
                    
                    log.write(f"Prefetch successful for {srr_id}\n")
                
                # Convert SRA to FASTQ
                sra_file = category_dir / srr_id / f"{srr_id}.sra"
                if sra_file.exists():
                    fastq_dump_cmd = [
                        "fastq-dump",
                        "--outdir", str(category_dir),
                        "--gzip",
                        "--split-files",
                        "--skip-technical",
                        "--readids",
                        "--read-filter", "pass",
                        str(sra_file)
                    ]
                    
                    with open(log_file, 'a') as log:
                        log.write(f"Converting SRA to FASTQ for {srr_id}\n")
                        result = subprocess.run(fastq_dump_cmd, capture_output=True, text=True, timeout=1800)
                        
                        if result.returncode != 0:
                            log.write(f"Fastq-dump failed: {result.stderr}\n")
                            raise subprocess.CalledProcessError(result.returncode, fastq_dump_cmd, result.stderr)
                        
                        log.write(f"FASTQ conversion successful for {srr_id}\n")
                    
                    # Clean up SRA file to save space
                    subprocess.run(["rm", "-rf", str(category_dir / srr_id)], check=False)
                    
                    logger.info(f"Successfully downloaded {srr_id} to {category}")
                    return True
                else:
                    logger.error(f"SRA file not found for {srr_id}")
                    return False
                    
            except subprocess.TimeoutExpired:
                logger.warning(f"Timeout downloading {srr_id}, attempt {attempt + 1}")
                if attempt < retry_count - 1:
                    time.sleep(30)  # Wait before retry
                continue
            except Exception as e:
                logger.error(f"Error downloading {srr_id}: {e}")
                if attempt < retry_count - 1:
                    time.sleep(30)  # Wait before retry
                continue
        
        logger.error(f"Failed to download {srr_id} after {retry_count} attempts")
        return False
    
    def download_category(self, category, srr_list, max_downloads=None):
        """Download FASTQ files for all SRR IDs in a category"""
        if not srr_list:
            logger.info(f"No SRR IDs found for category {category}")
            return
        
        if max_downloads:
            srr_list = srr_list[:max_downloads]
            logger.info(f"Limiting downloads to {max_downloads} samples for {category}")
        
        logger.info(f"Starting downloads for {category}: {len(srr_list)} samples")
        
        # Save SRR list for this category
        srr_list_file = self.base_dir / category / f"{category.lower()}_srr_list.txt"
        with open(srr_list_file, 'w') as f:
            for srr_id in srr_list:
                f.write(f"{srr_id}\n")
        logger.info(f"Saved SRR list to {srr_list_file}")
        
        # Download with parallel processing
        successful_downloads = 0
        failed_downloads = 0
        
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            future_to_srr = {
                executor.submit(self.download_srr, srr_id, category): srr_id 
                for srr_id in srr_list
            }
            
            for future in as_completed(future_to_srr):
                srr_id = future_to_srr[future]
                try:
                    success = future.result()
                    if success:
                        successful_downloads += 1
                    else:
                        failed_downloads += 1
                except Exception as e:
                    logger.error(f"Exception downloading {srr_id}: {e}")
                    failed_downloads += 1
        
        logger.info(f"Category {category} complete: {successful_downloads} successful, {failed_downloads} failed")
        return successful_downloads, failed_downloads

def load_classified_metadata(metadata_file):
    """Load classified metadata and organize SRR IDs by category"""
    logger.info(f"Loading classified metadata from {metadata_file}")
    
    df = pd.read_csv(metadata_file, sep='\t')
    logger.info(f"Loaded {len(df)} samples")
    
    # Group by classification category
    category_srr = {}
    
    for category in ["Tumor", "Pre-malignant", "Normal", "Cell_line", "Unknown"]:
        if category == "Cell_line":
            # Handle both "Cell line" and "Cell_line" variations
            mask = (df['top_label'] == "Cell line") | (df['top_label'] == "Cell_line")
        else:
            mask = df['top_label'] == category
        
        srr_ids = df[mask]['run_accession'].dropna().tolist()
        category_srr[category] = srr_ids
        logger.info(f"{category}: {len(srr_ids)} samples")
    
    return category_srr

def main():
    parser = argparse.ArgumentParser(description="Download FASTQ files organized by classification category")
    parser.add_argument("-i", "--input", required=True, help="Input classified metadata TSV file")
    parser.add_argument("-o", "--output-dir", default="fastq_downloads", help="Output directory for FASTQ files")
    parser.add_argument("-w", "--workers", type=int, default=4, help="Number of parallel download workers")
    parser.add_argument("--max-per-category", type=int, help="Maximum downloads per category (for testing)")
    parser.add_argument("--categories", nargs="+", default=["Tumor", "Pre-malignant", "Normal", "Cell_line"], 
                       help="Categories to download (default: all except Unknown)")
    
    args = parser.parse_args()
    
    # Load and organize SRR IDs by category
    category_srr = load_classified_metadata(args.input)
    
    # Initialize downloader
    downloader = FASTQDownloader(args.output_dir, args.workers)
    
    # Download each category
    total_successful = 0
    total_failed = 0
    
    for category in args.categories:
        if category in category_srr:
            successful, failed = downloader.download_category(
                category, 
                category_srr[category], 
                args.max_per_category
            )
            total_successful += successful
            total_failed += failed
        else:
            logger.warning(f"Category {category} not found in metadata")
    
    logger.info(f"Download complete: {total_successful} successful, {total_failed} failed")
    
    # Create summary report
    summary_file = Path(args.output_dir) / "download_summary.txt"
    with open(summary_file, 'w') as f:
        f.write("FASTQ Download Summary\n")
        f.write("=====================\n\n")
        f.write(f"Total successful downloads: {total_successful}\n")
        f.write(f"Total failed downloads: {total_failed}\n\n")
        f.write("Category breakdown:\n")
        for category in args.categories:
            if category in category_srr:
                f.write(f"  {category}: {len(category_srr[category])} samples\n")
    
    logger.info(f"Summary saved to {summary_file}")

if __name__ == "__main__":
    main()