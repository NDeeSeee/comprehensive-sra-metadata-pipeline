#!/usr/bin/env python3
"""
Cancer Classification Script
Classifies samples into tumor/normal/cell line categories using comprehensive metadata
Integrates cell line reference database for accurate cell line detection
"""

import pandas as pd
import re
import argparse
import pathlib

def load_cell_line_reference(cell_line_path="data/cell_line_model_from_audrey.csv"):
    """Load cell line names from reference database"""
    cell_line_names = set()
    
    try:
        df = pd.read_csv(cell_line_path)
        print(f"Loaded {len(df)} cell lines from reference database")
        
        # Extract cell line names (both original and stripped)
        for _, row in df.iterrows():
            # Original cell line name
            if pd.notna(row.get('CellLineName')):
                name = str(row['CellLineName']).strip()
                if name and len(name) >= 3:  # Minimum length to avoid false matches
                    cell_line_names.add(name.lower())
            
            # Stripped cell line name
            if pd.notna(row.get('StrippedCellLineName')):
                stripped = str(row['StrippedCellLineName']).strip()
                if stripped and len(stripped) >= 3:
                    cell_line_names.add(stripped.lower())
        
        print(f"Extracted {len(cell_line_names)} unique cell line names")
        
    except Exception as e:
        print(f"Warning: Could not load cell line reference: {e}")
        print("Falling back to generic cell line detection")
        return None  # Return None to indicate failure
    
    return cell_line_names

def classify_cancer_type(row, cell_line_names=None):
    """
    Enhanced classification using comprehensive metadata fields and cell line reference
    """
    # Initialize classification
    top_label = "Unknown"
    is_cell_line = "no"
    is_bulk_sorted = "no"
    is_control = "no"
    adjacent_normal = "no"
    barrett_grade = "unknown"
    tissue_origin = ""
    
    # Get all relevant fields (case-insensitive)
    fields_to_check = [
        'cell_line', 'disease', 'tissue_type', 'sample_title', 'study_title', 
        'experiment_title', 'source_name', 'cell_type', 'disease_state', 
        'histological_type', 'treatment', 'genotype', 'phenotype',
        'SampleName', 'Histological_Type', 'Body_Site', 'Tumor', 'Analyte_Type', 
        'Disease', 'experiment_title', 'study_title', 'sample_title',
        'age', 'altitude', 'environment_biome', 'host', 'host_body_site',
        'isolate', 'location', 'sex', 'strain', 'temperature', 'dev_stage',
        'environment_feature', 'environment_material', 'environmental_medium',
        'environmental_sample', 'host_genotype', 'host_phenotype',
        'biomaterial_provider', 'organism_part', 'sampling_site', 'analyte_type',
        'body_site', 'disease_staging', 'is_tumor', 'subject_is_affected',
        'individual', 'replicate', 'experimental_factor'
    ]
    
    # Combine all text fields for analysis
    all_text_parts = []
    for field in fields_to_check:
        value = str(row.get(field, '')).lower()
        if value and value != 'nan':
            all_text_parts.append(value)
    
    all_text = ' '.join(all_text_parts)
    
    # Decision Tree Logic - Check in priority order
    
    # 1. Check explicit tumor flags first
    tumor_flags = ['yes', 'true', '1']
    if str(row.get('Tumor', '')).lower() in tumor_flags or str(row.get('is_tumor', '')).lower() in tumor_flags:
        top_label = "Tumor"
    
    # 2. Check for Cell Line using reference database (if available)
    elif cell_line_names and len(cell_line_names) > 0:
        # Check for known cell line names with word boundary matching
        for cell_line_name in cell_line_names:
            if len(cell_line_name) >= 3:  # Minimum length
                # Use word boundary matching to avoid false positives
                if re.search(rf'\b{re.escape(cell_line_name)}\b', all_text):
                    top_label = "Cell line"
                    is_cell_line = "yes"
                    break
        
        # If no specific cell line found, check for generic cell line keywords
        if top_label != "Cell line":
            if any(keyword in all_text for keyword in ['cell line', 'culture', 'passaged', 'cellline', 'cl', 'line']):
                top_label = "Cell line"
                is_cell_line = "yes"
    
    # 3. If no cell line reference available, check for generic cell line keywords
    elif cell_line_names is None:
        if any(keyword in all_text for keyword in ['cell line', 'culture', 'passaged', 'cellline', 'cl', 'line']):
            top_label = "Cell line"
            is_cell_line = "yes"
    
    # 4. Check for Pre-malignant conditions
    if top_label == "Unknown":
        if any(keyword in all_text for keyword in [
            'barrett', 'metaplasia', 'dysplasia', 'lgd', 'hgd', 'indefinite',
            'pre-malignant', 'premalignant', 'precursor', 'atypia', 'hyperplasia'
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
    
    # 5. Check for Normal/Healthy
    if top_label == "Unknown":
        if any(keyword in all_text for keyword in [
            'normal', 'healthy', 'control', 'non-diseased', 'squamous epithelium',
            'healthy donor', 'baseline', 'wild type', 'wt', 'benign'
        ]):
            top_label = "Normal"
            is_control = "yes"
    
    # 6. Check for cancer/tumor keywords (generalizable) - only if still unknown
    if top_label == "Unknown":
        if any(keyword in all_text for keyword in [
            'carcinoma', 'cancer', 'tumor', 'tumour', 'malignant', 'adenocarcinoma',
            'squamous cell carcinoma', 'scc', 'invasive', 'metastatic', 'neoplasm',
            'esophageal adenocarcinoma', 'eac', 'esophageal cancer', 'eac1416',
            'sarcoma', 'lymphoma', 'leukemia', 'melanoma', 'glioblastoma',
            'blastoma', 'carcinoid', 'adenoma', 'papilloma', 'fibroma'
        ]):
            top_label = "Tumor"
    
    # Additional classifications
    
    # Check for bulk sorted/purified populations
    if any(keyword in all_text for keyword in [
        'sorted', 'purified', 'isolated', 'enriched', 'facs', 'magnetic',
        'cd4', 'cd8', 't cell', 'b cell', 'monocyte', 'macrophage'
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
    
    # Determine tissue origin (generalizable)
    tissue_text = f"{all_text} {str(row.get('Body_Site', '')).lower()}"
    if 'esophagus' in tissue_text or 'esophageal' in tissue_text:
        tissue_origin = "esophagus"
    elif 'stomach' in tissue_text or 'gastric' in tissue_text:
        tissue_origin = "stomach"
    elif 'intestine' in tissue_text or 'intestinal' in tissue_text or 'colon' in tissue_text:
        tissue_origin = "intestine"
    elif 'lung' in tissue_text or 'pulmonary' in tissue_text:
        tissue_origin = "lung"
    elif 'liver' in tissue_text or 'hepatic' in tissue_text:
        tissue_origin = "liver"
    elif 'breast' in tissue_text or 'mammary' in tissue_text:
        tissue_origin = "breast"
    elif 'brain' in tissue_text or 'cerebral' in tissue_text:
        tissue_origin = "brain"
    elif 'ovary' in tissue_text or 'ovarian' in tissue_text:
        tissue_origin = "ovary"
    elif 'prostate' in tissue_text or 'prostatic' in tissue_text:
        tissue_origin = "prostate"
    elif 'pancreas' in tissue_text or 'pancreatic' in tissue_text:
        tissue_origin = "pancreas"
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
    parser = argparse.ArgumentParser(description='Cancer classification with cell line reference integration')
    parser.add_argument('-i', '--input', required=True, help='Input TSV file with comprehensive metadata')
    parser.add_argument('-o', '--output', required=True, help='Output TSV file with classifications')
    parser.add_argument('--cell-line-ref', default='data/cell_line_model_from_audrey.csv',
                       help='Path to cell line reference database')
    args = parser.parse_args()
    
    # Load cell line reference
    print("Loading cell line reference database...")
    cell_line_names = load_cell_line_reference(args.cell_line_ref)
    
    # Load metadata
    print(f"Loading comprehensive metadata from: {args.input}")
    df = pd.read_csv(args.input, sep='\t', low_memory=False)
    print(f"Loaded {len(df)} samples with {len(df.columns)} columns")
    
    # Apply classification
    print("Applying cancer classification with cell line reference...")
    classifications = []
    
    for idx, row in df.iterrows():
        classification = classify_cancer_type(row, cell_line_names)
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
    print(f"\nCancer Classification Summary:")
    print(f"Total samples: {len(result_df)}")
    print(f"Total columns: {len(result_df.columns)}")
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