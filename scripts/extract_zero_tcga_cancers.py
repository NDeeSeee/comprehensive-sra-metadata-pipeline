#!/usr/bin/env python3
"""
Extract cancer types with TCGA=0 from cancers_poseidon_rna_seq.csv
"""

import csv
import re
from pathlib import Path

def clean_cancer_type_name(cancer_type):
    """Clean and normalize cancer type names for use in file paths and searches"""
    # Remove extra spaces and clean up
    cancer_type = cancer_type.strip()
    
    # Replace problematic characters for file paths
    clean_name = re.sub(r'[^\w\s-]', '', cancer_type)  # Remove special chars except spaces and hyphens
    clean_name = re.sub(r'\s+', '_', clean_name)  # Replace spaces with underscores
    clean_name = re.sub(r'_+', '_', clean_name)  # Replace multiple underscores with single
    clean_name = clean_name.strip('_')  # Remove leading/trailing underscores
    
    return clean_name.lower()

def extract_zero_tcga_cancers(csv_file):
    """Extract cancer types with TCGA=0"""
    zero_tcga_cancers = []
    
    with open(csv_file, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        
        for i, row in enumerate(reader):
            # Skip header rows
            if i < 2:
                continue
                
            # Skip empty rows
            if len(row) < 4:
                continue
                
            # Check if TCGA column (index 3) is 0 and cancer type exists
            tcga_value = row[3].strip()
            cancer_type = row[2].strip() if len(row) > 2 else ""
            
            if tcga_value == "0" and cancer_type:
                other_value = row[4].strip() if len(row) > 4 else "0"
                source = row[5].strip() if len(row) > 5 else ""
                
                zero_tcga_cancers.append({
                    'original_name': cancer_type,
                    'clean_name': clean_cancer_type_name(cancer_type),
                    'tcga': tcga_value,
                    'other': other_value,
                    'source': source,
                    'row_number': i + 1
                })
    
    return zero_tcga_cancers

def main():
    csv_file = "data/cancers_poseidon_rna_seq.csv"
    
    if not Path(csv_file).exists():
        print(f"Error: {csv_file} not found")
        return
    
    zero_tcga_cancers = extract_zero_tcga_cancers(csv_file)
    
    print("Cancer types with TCGA=0:")
    print("=" * 50)
    
    for i, cancer in enumerate(zero_tcga_cancers, 1):
        print(f"{i:2d}. {cancer['original_name']}")
        print(f"    Clean name: {cancer['clean_name']}")
        print(f"    TCGA: {cancer['tcga']}, Other: {cancer['other']}, Source: {cancer['source']}")
        print()
    
    print(f"Total cancer types with TCGA=0: {len(zero_tcga_cancers)}")
    
    # Save to file for batch processing
    output_file = "data/zero_tcga_cancers.txt"
    with open(output_file, 'w') as f:
        for cancer in zero_tcga_cancers:
            f.write(f"{cancer['original_name']}\n")
    
    print(f"Cancer types saved to: {output_file}")

if __name__ == "__main__":
    main()