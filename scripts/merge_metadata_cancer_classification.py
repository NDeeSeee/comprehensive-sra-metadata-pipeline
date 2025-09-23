#!/usr/bin/env python3
"""
Enhanced merge script specifically for cancer classification
Includes additional fields needed for tumor/normal/cell line classification
"""

import argparse, json, pathlib, re
import pandas as pd
import xml.etree.ElementTree as ET

def load_runinfo(path):
    """Load SRA RunInfo CSV with error handling"""
    try:
        df = pd.read_csv(path, low_memory=False)
        return df
    except Exception as e:
        print(f"Warning: Could not load RunInfo: {e}")
        return pd.DataFrame()

def load_ena(path):
    """Load ENA filereport TSV with duplicate header handling"""
    try:
        lines = open(path, 'r').read().strip().splitlines()
        if not lines:
            return pd.DataFrame()
        
        # Keep first header; drop repeated headers
        header = lines[0]
        rows = [header] + [ln for ln in lines[1:] if ln != header]
        
        from io import StringIO
        return pd.read_csv(StringIO("\n".join(rows)), sep="\t")
    except Exception as e:
        print(f"Warning: Could not load ENA data: {e}")
        return pd.DataFrame()

def parse_biosample_jsonl(path):
    """Enhanced BioSample JSON parser with clinical data extraction"""
    records = []
    try:
        with open(path, 'r') as fh:
            for line in fh:
                line = line.strip()
                if not line: continue
                try:
                    obj = json.loads(line)
                except json.JSONDecodeError:
                    continue
                
                # Handle different JSON structures
                biosamples = []
                if isinstance(obj, dict):
                    if 'BioSampleSet' in obj:
                        bs = obj['BioSampleSet'].get('BioSample', [])
                        if isinstance(bs, dict): biosamples = [bs]
                        else: biosamples = bs
                    elif 'BioSample' in obj:
                        biosamples = [obj['BioSample']]
                
                for bs in biosamples:
                    acc = bs.get('accession')
                    attrs = {}
                    
                    # Extract attributes
                    at_root = bs.get('Attributes', {}).get('Attribute', [])
                    if isinstance(at_root, dict): at_root = [at_root]
                    
                    for a in at_root:
                        name = a.get('@attribute_name') or a.get('attribute_name') or a.get('harmonized_name')
                        val = a.get('#text') or a.get('value') or a.get('text') or ''
                        if name:
                            key = name.strip().lower().replace(' ', '_').replace('-', '_')
                            attrs[key] = val
                    
                    # Capture title/description
                    title = bs.get('Title')
                    descr = bs.get('Description')
                    if title: attrs['biosample_title'] = title
                    if descr: attrs['biosample_description'] = descr
                    
                    rec = {'BioSample': acc}
                    rec.update(attrs)
                    records.append(rec)
    except Exception as e:
        print(f"Warning: Could not parse BioSample data: {e}")
    
    if records:
        return pd.DataFrame(records).drop_duplicates(subset=['BioSample'])
    else:
        return pd.DataFrame()

def parse_bioproject_jsonl(path):
    """Parse BioProject JSON data"""
    records = []
    try:
        with open(path, 'r') as fh:
            for line in fh:
                line = line.strip()
                if not line: continue
                try:
                    obj = json.loads(line)
                except json.JSONDecodeError:
                    continue
                
                # Extract BioProject information
                if isinstance(obj, dict):
                    project_data = {}
                    if 'Project' in obj:
                        project = obj['Project']
                        project_data['BioProject'] = project.get('ProjectID', {}).get('ArchiveID', {}).get('#text', '')
                        project_data['project_title'] = project.get('ProjectDescr', {}).get('Title', '')
                        project_data['project_description'] = project.get('ProjectDescr', {}).get('Description', '')
                        project_data['project_type'] = project.get('ProjectType', {}).get('ProjectTypeSubmission', '')
                        
                        # Extract study information
                        if 'Study' in project:
                            study = project['Study']
                            project_data['study_title'] = study.get('Descriptor', {}).get('StudyTitle', '')
                            project_data['study_description'] = study.get('Descriptor', {}).get('StudyAbstract', '')
                        
                        records.append(project_data)
    except Exception as e:
        print(f"Warning: Could not parse BioProject data: {e}")
    
    if records:
        return pd.DataFrame(records).drop_duplicates(subset=['BioProject'])
    else:
        return pd.DataFrame()

def parse_sra_xml_jsonl(path):
    """Parse SRA XML data converted to JSON"""
    records = []
    try:
        with open(path, 'r') as fh:
            for line in fh:
                line = line.strip()
                if not line: continue
                try:
                    obj = json.loads(line)
                    if obj:
                        records.append(obj)
                except json.JSONDecodeError:
                    continue
    except Exception as e:
        print(f"Warning: Could not parse SRA XML data: {e}")
    
    if records:
        return pd.DataFrame(records)
    else:
        return pd.DataFrame()

def parse_geo_tsv(path):
    """Parse GEO metadata"""
    records = []
    try:
        with open(path, 'r') as fh:
            for line in fh:
                line = line.strip()
                if not line: continue
                try:
                    obj = json.loads(line)
                    if obj:
                        records.append(obj)
                except json.JSONDecodeError:
                    continue
    except Exception as e:
        print(f"Warning: Could not parse GEO data: {e}")
    
    if records:
        return pd.DataFrame(records)
    else:
        return pd.DataFrame()

def parse_ffq_jsonl(path):
    """Parse ffq JSON data"""
    records = []
    try:
        with open(path, 'r') as fh:
            for line in fh:
                line = line.strip()
                if not line: continue
                try:
                    obj = json.loads(line)
                    if obj:
                        records.append(obj)
                except json.JSONDecodeError:
                    continue
    except Exception as e:
        print(f"Warning: Could not parse ffq data: {e}")
    
    if records:
        return pd.DataFrame(records)
    else:
        return pd.DataFrame()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--indir", required=True, help="Directory with raw/ files")
    ap.add_argument("-o", "--outfile", required=True, help="Output TSV")
    args = ap.parse_args()
    indir = pathlib.Path(args.indir)

    print("Loading data sources for cancer classification...")
    
    # Load all data sources with proper file existence checks
    runinfo = load_runinfo(indir / "runinfo.csv") if (indir / "runinfo.csv").exists() else pd.DataFrame()
    ena = load_ena(indir / "ena_read_run.tsv") if (indir / "ena_read_run.tsv").exists() else pd.DataFrame()
    
    # Check for both possible filenames (original and extracted)
    biosample_file = indir / "biosample_extracted.jsonl" if (indir / "biosample_extracted.jsonl").exists() else indir / "biosample.jsonl"
    bios = parse_biosample_jsonl(biosample_file) if biosample_file.exists() else pd.DataFrame()
    
    bioproject_file = indir / "bioproject_extracted.jsonl" if (indir / "bioproject_extracted.jsonl").exists() else indir / "bioproject.jsonl"
    bioproj = parse_bioproject_jsonl(bioproject_file) if bioproject_file.exists() else pd.DataFrame()
    
    sra_xml = parse_sra_xml_jsonl(indir / "sra_xml.jsonl") if (indir / "sra_xml.jsonl").exists() else pd.DataFrame()
    geo = parse_geo_tsv(indir / "geo_metadata.tsv") if (indir / "geo_metadata.tsv").exists() else pd.DataFrame()
    ffq = parse_ffq_jsonl(indir / "ffq.jsonl") if (indir / "ffq.jsonl").exists() else pd.DataFrame()

    print(f"Loaded: RunInfo({len(runinfo)}), ENA({len(ena)}), BioSample({len(bios)}), BioProject({len(bioproj)}), SRA-XML({len(sra_xml)}), GEO({len(geo)}), ffq({len(ffq)})")

    # Validate that we have at least one data source
    if ena.empty and runinfo.empty:
        print("ERROR: No primary data sources found (ENA and RunInfo both empty)")
        print("Available files in directory:")
        for f in indir.iterdir():
            print(f"  - {f.name}")
        return

    # Start with ENA data as primary source (most comprehensive)
    if not ena.empty:
        merged = ena.copy()
        print(f"Using ENA as primary source: {len(merged)} samples")
    elif not runinfo.empty:
        # If ENA is empty, start with RunInfo
        merged = runinfo.copy()
        merged = merged.rename(columns={'Run': 'run_accession'})
        print(f"Using RunInfo as primary source: {len(merged)} samples")
    else:
        merged = pd.DataFrame()
    
    # Merge with RunInfo (if not already used as primary)
    if not runinfo.empty and not ena.empty:
        runinfo = runinfo.rename(columns={'Run': 'run_accession'})
        merged = merged.merge(runinfo, on='run_accession', how='left', suffixes=('', '_sra'))
    
    # Merge with BioSample
    if not bios.empty and 'BioSample' in merged.columns:
        merged = merged.merge(bios, on='BioSample', how='left', suffixes=('', '_bs'))
    
    # Merge with BioProject
    if not bioproj.empty and 'BioProject' in merged.columns:
        merged = merged.merge(bioproj, on='BioProject', how='left', suffixes=('', '_bp'))
    
    # Merge with SRA XML
    if not sra_xml.empty and 'run_accession' in sra_xml.columns:
        merged = merged.merge(sra_xml, on='run_accession', how='left', suffixes=('', '_xml'))
    
    # Merge with GEO (more efficient approach)
    if not geo.empty and not merged.empty:
        # Create a GEO lookup dictionary for efficient matching
        geo_lookup = {}
        for _, geo_row in geo.iterrows():
            for col, val in geo_row.items():
                if isinstance(val, str) and val.startswith('GSE'):
                    geo_lookup[val] = geo_row.to_dict()
        
        # Add GEO columns to merged data
        geo_cols = []
        for col in geo.columns:
            geo_cols.append(f'geo_{col}')
            merged[f'geo_{col}'] = ''
        
        # Match GEO data efficiently
        if 'study_title' in merged.columns:
            for idx, row in merged.iterrows():
                study_title = str(row.get('study_title', ''))
                for gse_id, geo_data in geo_lookup.items():
                    if gse_id in study_title:
                        for col, val in geo_data.items():
                            merged.at[idx, f'geo_{col}'] = val
                        break
    
    # Merge with ffq
    if not ffq.empty and 'run_accession' in ffq.columns:
        merged = merged.merge(ffq, on='run_accession', how='left', suffixes=('', '_ffq'))

    # Reorder important columns for cancer classification
    cancer_classification_cols = [
        'run_accession', 'Experiment', 'SampleName', 'BioSample', 'SRAStudy', 'BioProject',
        'study_accession', 'secondary_study_accession', 'experiment_accession', 'sample_accession',
        'secondary_sample_accession',
        'LibraryLayout', 'LibraryStrategy', 'LibrarySource', 'Platform', 'Model',
        'instrument_model', 'library_layout', 'library_strategy', 'library_source',
        'read_count', 'base_count', 'first_public', 'last_updated',
        # Cancer classification specific fields
        'age', 'altitude', 'cell_line', 'cell_type', 'disease', 'environment_biome',
        'host', 'host_body_site', 'isolate', 'location', 'sex', 'strain', 'temperature', 'tissue_type',
        'treatment', 'genotype', 'phenotype', 'source_name', 'biomaterial_provider', 'organism_part',
        'sampling_site', 'analyte_type', 'body_site', 'histological_type', 'disease_staging',
        'is_tumor', 'subject_is_affected', 'individual', 'replicate', 'experimental_factor',
        'broker_name', 'center_name', 'experiment_title', 'library_name', 'library_selection',
        'scientific_name', 'collection_date', 'study_title', 'sample_title',
        'submitted_ftp', 'fastq_ftp', 'dev_stage', 'environment_feature', 'environment_material',
        'environmental_medium', 'environmental_sample', 'host_genotype', 'host_phenotype',
        'ReleaseDate', 'LoadDate', 'spots', 'bases', 'spots_with_mates', 'avgLength', 'size_MB',
        'AssemblyName', 'download_path', 'LibraryName', 'LibrarySelection', 'InsertSize', 'InsertDev',
        'Study_Pubmed_id', 'ProjectID', 'Sample', 'SampleType', 'TaxID', 'ScientificName',
        'g1k_pop_code', 'source', 'g1k_analysis_group', 'Subject_ID', 'Sex', 'Disease', 'Tumor',
        'Affection_Status', 'Analyte_Type', 'Histological_Type', 'Body_Site', 'CenterName',
        'Submission', 'dbgap_study_accession', 'Consent', 'RunHash', 'ReadHash'
    ]
    
    prefer = [c for c in cancer_classification_cols if c in merged.columns]
    others = [c for c in merged.columns if c not in prefer]
    merged = merged[prefer + others]

    # Remove duplicate columns (keep first occurrence)
    merged = merged.loc[:, ~merged.columns.duplicated()]

    merged.to_csv(args.outfile, sep="\t", index=False)
    print(f"Wrote: {args.outfile}  (rows={len(merged)}, columns={len(merged.columns)})")
    print(f"Cancer classification fields available: {len([c for c in merged.columns if c in cancer_classification_cols])}")

if __name__ == "__main__":
    main()