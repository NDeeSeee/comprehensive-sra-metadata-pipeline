#!/usr/bin/env python3
"""
Cancer Classification Script
Classifies samples into: Cell line, Tumor, Pre-malignant, Normal, Unknown
Based on metadata fields for downstream analysis filtering
"""

import pandas as pd
import re
import argparse
import pathlib

def classify_cancer_type(row):
    """
    Classify sample based on metadata using decision tree logic
    """
    # Initialize classification
    top_label = "Unknown"
    is_cell_line = "no"
    is_bulk_sorted = "no"
    is_control = "no"
    adjacent_normal = "no"
    barrett_grade = "unknown"
    tissue_origin = ""
    
    # Get relevant fields (case-insensitive)
    cell_line = str(row.get('cell_line', '')).lower()
    disease = str(row.get('disease', '')).lower()
    tissue_type = str(row.get('tissue_type', '')).lower()
    sample_title = str(row.get('sample_title', '')).lower()
    study_title = str(row.get('study_title', '')).lower()
    experiment_title = str(row.get('experiment_title', '')).lower()
    source_name = str(row.get('source_name', '')).lower()
    cell_type = str(row.get('cell_type', '')).lower()
    disease_state = str(row.get('disease_state', '')).lower()
    histological_type = str(row.get('histological_type', '')).lower()
    treatment = str(row.get('treatment', '')).lower()
    genotype = str(row.get('genotype', '')).lower()
    phenotype = str(row.get('phenotype', '')).lower()
    
    # Combine all text fields for analysis
    all_text = f"{cell_line} {disease} {tissue_type} {sample_title} {study_title} {experiment_title} {source_name} {cell_type} {disease_state} {histological_type} {treatment} {genotype} {phenotype}"
    
    # Decision Tree Logic
    
    # 1. Check for Cell Line
    cell_line_keywords = ['cell line', 'culture', 'passaged', 'cellline', 'cl', 'line']
    if any(keyword in all_text for keyword in cell_line_keywords):
        top_label = "Cell line"
        is_cell_line = "yes"
    
    # 2. If not cell line, check for cancer/tumor
    elif any(keyword in all_text for keyword in [
        'carcinoma', 'cancer', 'tumor', 'tumour', 'malignant', 'adenocarcinoma',
        'squamous cell carcinoma', 'scc', 'invasive', 'metastatic', 'neoplasm',
        'esophageal adenocarcinoma', 'eac', 'esophageal cancer'
    ]):
        top_label = "Tumor"
    
    # 3. Check for Pre-malignant conditions
    elif any(keyword in all_text for keyword in [
        'barrett', 'metaplasia', 'dysplasia', 'lgd', 'hgd', 'indefinite',
        'pre-malignant', 'premalignant', 'precursor', 'atypia'
    ]):
        top_label = "Pre-malignant"
        
        # Determine Barrett's grade
        if 'lgd' in all_text or 'low grade dysplasia' in all_text:
            barrett_grade = "LGD"
        elif 'hgd' in all_text or 'high grade dysplasia' in all_text:
            barrett_grade = "HGD"
        elif 'indefinite' in all_text:
            barrett_grade = "indefinite"
        elif 'no dysplasia' in all_text:
            barrett_grade = "no dysplasia"
    
    # 4. Check for Normal/Healthy
    elif any(keyword in all_text for keyword in [
        'normal', 'healthy', 'control', 'non-diseased', 'squamous epithelium',
        'healthy donor', 'baseline', 'wild type', 'wt'
    ]):
        top_label = "Normal"
        is_control = "yes"
    
    # Additional classifications
    
    # Check for bulk sorted/purified populations
    if any(keyword in all_text for keyword in [
        'sorted', 'purified', 'isolated', 'enriched', 'facs', 'magnetic',
        'cd4', 'cd8', 't cell', 'b cell', 'monocyte', 'macrophage', 'fibroblast', 'CAFS'
    ]):
        is_bulk_sorted = "yes"
    
    # Check for adjacent normal
    if any(keyword in all_text for keyword in [
        'adjacent', 'marginal normal', 'paired normal', 'surrounding normal',
        'normal adjacent', 'adjacent normal'
    ]):
        adjacent_normal = "yes"
        if top_label == "Unknown":
            top_label = "Normal"
    
    # Determine tissue origin
    if 'esophagus' in all_text or 'esophageal' in all_text:
        tissue_origin = "esophagus"
    elif 'stomach' in all_text or 'gastric' in all_text:
        tissue_origin = "stomach"
    elif 'intestine' in all_text or 'intestinal' in all_text:
        tissue_origin = "intestine"
    elif 'lung' in all_text or 'pulmonary' in all_text:
        tissue_origin = "lung"
    else:
        tissue_origin = "unknown"
    
    return {
        'top_label': top_label,
        'is_cell_line': is_cell_line,
        'is_bulk_sorted': is_bulk_sorted,
        'is_control': is_control,
        'adjacent_normal': adjacent_normal,
        'barrett_grade': barrett_grade,
        'tissue_origin': tissue_origin
    }

def main():
    parser = argparse.ArgumentParser(description='Classify samples for cancer analysis')
    parser.add_argument('-i', '--input', required=True, help='Input TSV file with metadata')
    parser.add_argument('-o', '--output', required=True, help='Output TSV file with classifications')
    args = parser.parse_args()
    
    # Load metadata
    print(f"Loading metadata from: {args.input}")
    df = pd.read_csv(args.input, sep='\t', low_memory=False)
    print(f"Loaded {len(df)} samples")
    
    # Apply classification
    print("Applying cancer classification...")
    classifications = []
    
    for idx, row in df.iterrows():
        classification = classify_cancer_type(row)
        classifications.append(classification)
        
        if (idx + 1) % 50 == 0:
            print(f"Processed {idx + 1}/{len(df)} samples")
    
    # Add classification columns to dataframe
    classification_df = pd.DataFrame(classifications)
    result_df = pd.concat([df, classification_df], axis=1)
    
    # Reorder columns to put classification fields at the front
    classification_cols = ['top_label', 'is_cell_line', 'is_bulk_sorted', 'is_control', 
                          'adjacent_normal', 'barrett_grade', 'tissue_origin']
    
    # Get other columns
    other_cols = [col for col in result_df.columns if col not in classification_cols]
    
    # Reorder
    result_df = result_df[classification_cols + other_cols]
    
    # Save results
    result_df.to_csv(args.output, sep='\t', index=False)
    
    # Print summary
    print(f"\nClassification Summary:")
    print(f"Total samples: {len(result_df)}")
    print(f"\nTop Label Distribution:")
    print(result_df['top_label'].value_counts())
    print(f"\nCell Line Samples: {len(result_df[result_df['is_cell_line'] == 'yes'])}")
    print(f"Bulk Sorted Samples: {len(result_df[result_df['is_bulk_sorted'] == 'yes'])}")
    print(f"Control Samples: {len(result_df[result_df['is_control'] == 'yes'])}")
    print(f"Adjacent Normal Samples: {len(result_df[result_df['adjacent_normal'] == 'yes'])}")
    
    print(f"\nBarrett's Grade Distribution:")
    print(result_df['barrett_grade'].value_counts())
    
    print(f"\nTissue Origin Distribution:")
    print(result_df['tissue_origin'].value_counts())
    
    print(f"\nResults saved to: {args.output}")

if __name__ == "__main__":
    main()