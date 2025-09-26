#!/usr/bin/env python3
"""
Complete Cancer Analysis Pipeline
Integrates cancer type search, metadata collection, classification, and FASTQ download
"""

import os
import sys
import subprocess
import argparse
import pandas as pd
import time
from pathlib import Path
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('cancer_analysis_pipeline.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

class CancerAnalysisPipeline:
    def __init__(self, cancer_type, output_dir="cancer_analysis_output"):
        self.cancer_type = cancer_type
        self.output_dir = Path(output_dir)
        self.cancer_type_dir = self.output_dir / cancer_type.replace(" ", "_").replace("/", "_")
        
        # Create directory structure
        self.create_directory_structure()
        
    def create_directory_structure(self):
        """Create organized directory structure for the analysis"""
        logger.info(f"Creating directory structure in {self.cancer_type_dir}")
        
        # Main directories
        (self.cancer_type_dir / "metadata").mkdir(parents=True, exist_ok=True)
        (self.cancer_type_dir / "classification").mkdir(exist_ok=True)
        (self.cancer_type_dir / "fastq_downloads").mkdir(exist_ok=True)
        (self.cancer_type_dir / "logs").mkdir(exist_ok=True)
        
        # FASTQ subdirectories by category
        for category in ["Tumor", "Pre-malignant", "Normal", "Cell_line", "Unknown"]:
            fastq_dir = self.cancer_type_dir / "fastq_downloads" / category
            fastq_dir.mkdir(exist_ok=True)
            (fastq_dir / "raw").mkdir(exist_ok=True)
            (fastq_dir / "processed").mkdir(exist_ok=True)
            (fastq_dir / "logs").mkdir(exist_ok=True)
        
        logger.info("Directory structure created successfully")
    
    def step1_search_cancer_type(self, max_results=1000):
        """Step 1: Search for cancer type and extract SRR IDs"""
        logger.info(f"Step 1: Searching for '{self.cancer_type}' cancer type")
        
        srr_list_file = self.cancer_type_dir / "metadata" / f"{self.cancer_type.replace(' ', '_')}_srr_list.txt"
        
        # Run cancer type search
        search_cmd = [
            "python", "scripts/cancer_type_search.py",
            "-c", self.cancer_type,
            "-o", str(srr_list_file),
            "-m", str(max_results)
        ]
        
        logger.info(f"Running: {' '.join(search_cmd)}")
        result = subprocess.run(search_cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            logger.error(f"Cancer type search failed: {result.stderr}")
            return False
        
        logger.info(f"Cancer type search completed. SRR list saved to {srr_list_file}")
        return True
    
    def step2_collect_metadata(self, srr_list_file):
        """Step 2: Collect comprehensive metadata for SRR IDs"""
        logger.info("Step 2: Collecting comprehensive metadata")
        
        metadata_file = self.cancer_type_dir / "metadata" / "comprehensive_metadata.tsv"
        
        # Run comprehensive metadata collection
        metadata_cmd = [
            "bash", "scripts/comprehensive_metadata_pipeline.sh",
            str(srr_list_file),
            str(metadata_file)
        ]
        
        logger.info(f"Running: {' '.join(metadata_cmd)}")
        result = subprocess.run(metadata_cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            logger.error(f"Metadata collection failed: {result.stderr}")
            return False
        
        logger.info(f"Metadata collection completed. Results saved to {metadata_file}")
        return True
    
    def step3_classify_samples(self, metadata_file):
        """Step 3: Classify samples into categories"""
        logger.info("Step 3: Classifying samples into categories")
        
        classified_file = self.cancer_type_dir / "classification" / "classified_metadata.tsv"
        
        # Run cancer classification
        classify_cmd = [
            "python", "scripts/cancer_classification.py",
            "-i", str(metadata_file),
            "-o", str(classified_file)
        ]
        
        logger.info(f"Running: {' '.join(classify_cmd)}")
        result = subprocess.run(classify_cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            logger.error(f"Classification failed: {result.stderr}")
            return False
        
        logger.info(f"Classification completed. Results saved to {classified_file}")
        return True
    
    def step4_download_fastq(self, classified_file, max_per_category=None):
        """Step 4: Download FASTQ files organized by category"""
        logger.info("Step 4: Downloading FASTQ files by category")
        
        fastq_dir = self.cancer_type_dir / "fastq_downloads"
        
        # Run FASTQ download
        download_cmd = [
            "python", "scripts/download_fastq_by_category.py",
            "-i", str(classified_file),
            "-o", str(fastq_dir),
            "-w", "4"  # 4 parallel workers
        ]
        
        if max_per_category:
            download_cmd.extend(["--max-per-category", str(max_per_category)])
        
        logger.info(f"Running: {' '.join(download_cmd)}")
        result = subprocess.run(download_cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            logger.error(f"FASTQ download failed: {result.stderr}")
            return False
        
        logger.info(f"FASTQ download completed. Files saved to {fastq_dir}")
        return True
    
    def run_complete_pipeline(self, max_results=1000, max_per_category=None, skip_download=False):
        """Run the complete cancer analysis pipeline"""
        logger.info(f"Starting complete pipeline for '{self.cancer_type}'")
        start_time = time.time()
        
        # Step 1: Search for cancer type
        if not self.step1_search_cancer_type(max_results):
            logger.error("Pipeline failed at Step 1: Cancer type search")
            return False
        
        srr_list_file = self.cancer_type_dir / "metadata" / f"{self.cancer_type.replace(' ', '_')}_srr_list.txt"
        
        # Step 2: Collect metadata
        if not self.step2_collect_metadata(srr_list_file):
            logger.error("Pipeline failed at Step 2: Metadata collection")
            return False
        
        metadata_file = self.cancer_type_dir / "metadata" / "comprehensive_metadata.tsv"
        
        # Step 3: Classify samples
        if not self.step3_classify_samples(metadata_file):
            logger.error("Pipeline failed at Step 3: Classification")
            return False
        
        classified_file = self.cancer_type_dir / "classification" / "classified_metadata.tsv"
        
        # Step 4: Download FASTQ files (optional)
        if not skip_download:
            if not self.step4_download_fastq(classified_file, max_per_category):
                logger.error("Pipeline failed at Step 4: FASTQ download")
                return False
        
        # Create pipeline summary
        self.create_pipeline_summary(classified_file, start_time)
        
        logger.info(f"Complete pipeline finished successfully for '{self.cancer_type}'")
        return True
    
    def create_pipeline_summary(self, classified_file, start_time):
        """Create a summary report of the pipeline results"""
        logger.info("Creating pipeline summary")
        
        # Load classification results
        df = pd.read_csv(classified_file, sep='\t')
        
        # Count samples by category
        category_counts = df['top_label'].value_counts()
        
        # Create summary report
        summary_file = self.cancer_type_dir / "pipeline_summary.txt"
        with open(summary_file, 'w') as f:
            f.write(f"Cancer Analysis Pipeline Summary\n")
            f.write("=" * 50 + "\n\n")
            f.write(f"Cancer Type: {self.cancer_type}\n")
            f.write(f"Analysis Date: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Total Runtime: {time.time() - start_time:.2f} seconds\n\n")
            
            f.write("Sample Classification Results:\n")
            f.write("-" * 30 + "\n")
            for category, count in category_counts.items():
                percentage = (count / len(df)) * 100
                f.write(f"{category:15}: {count:4} samples ({percentage:5.1f}%)\n")
            
            f.write(f"\nTotal Samples: {len(df)}\n")
            
            f.write("\nDirectory Structure:\n")
            f.write("-" * 20 + "\n")
            f.write(f"Metadata: {self.cancer_type_dir}/metadata/\n")
            f.write(f"Classification: {self.cancer_type_dir}/classification/\n")
            f.write(f"FASTQ Downloads: {self.cancer_type_dir}/fastq_downloads/\n")
            f.write(f"Logs: {self.cancer_type_dir}/logs/\n")
        
        logger.info(f"Pipeline summary saved to {summary_file}")

def main():
    parser = argparse.ArgumentParser(description="Complete cancer analysis pipeline")
    parser.add_argument("-c", "--cancer-type", required=True, help="Cancer type to analyze")
    parser.add_argument("-o", "--output-dir", default="cancer_analysis_output", help="Output directory")
    parser.add_argument("-m", "--max-results", type=int, default=1000, help="Maximum SRR results to search")
    parser.add_argument("--max-per-category", type=int, help="Maximum FASTQ downloads per category")
    parser.add_argument("--skip-download", action="store_true", help="Skip FASTQ download step")
    parser.add_argument("--test", action="store_true", help="Run in test mode with limited results")
    
    args = parser.parse_args()
    
    # Adjust parameters for test mode
    if args.test:
        args.max_results = 50
        args.max_per_category = 5
        logger.info("Running in test mode with limited results")
    
    # Initialize and run pipeline
    pipeline = CancerAnalysisPipeline(args.cancer_type, args.output_dir)
    
    success = pipeline.run_complete_pipeline(
        max_results=args.max_results,
        max_per_category=args.max_per_category,
        skip_download=args.skip_download
    )
    
    if success:
        logger.info("Pipeline completed successfully!")
        sys.exit(0)
    else:
        logger.error("Pipeline failed!")
        sys.exit(1)

if __name__ == "__main__":
    main()