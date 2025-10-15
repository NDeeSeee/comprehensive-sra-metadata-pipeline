#!/usr/bin/env python3
"""
MAXIMUM metadata merge script
Merges ALL data sources with ALL fields for comprehensive metadata table
Target: 200+ columns
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

def load_ena_maximum(path):
    """Load ENA filereport with ALL fields and handle duplicates"""
    try:
        lines = open(path, 'r').read().strip().splitlines()
        if not lines:
            return pd.DataFrame()
        
        # Keep first header; drop repeated headers
        header = lines[0]
        rows = [header] + [ln for ln in lines[1:] if ln != header]
        
        from io import StringIO
        df = pd.read_csv(StringIO("\n".join(rows)), sep="\t")
        print(f"Loaded ENA data with {len(df.columns)} columns")
        return df
    except Exception as e:
        print(f"Warning: Could not load ENA data: {e}")
        return pd.DataFrame()

def parse_biosample_maximum(path):
    """Enhanced BioSample JSON parser with comprehensive attribute extraction"""
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
                    
                    # Extract ALL attributes comprehensively
                    at_root = bs.get('Attributes', {}).get('Attribute', [])
                    if isinstance(at_root, dict): at_root = [at_root]
                    
                    for a in at_root:
                        name = a.get('@attribute_name') or a.get('attribute_name') or a.get('harmonized_name')
                        val = a.get('#text') or a.get('value') or a.get('text') or ''
                        if name:
                            # Clean and normalize attribute names
                            key = name.strip().lower().replace(' ', '_').replace('-', '_').replace('(', '').replace(')', '')
                            attrs[f'biosample_{key}'] = val
                    
                    # Capture ALL metadata fields
                    for field in ['Title', 'Description', 'Organism', 'Taxonomy_ID', 'BioSampleModel']:
                        if field in bs:
                            attrs[f'biosample_{field.lower()}'] = bs[field]
                    
                    # Extract Links and References
                    if 'Links' in bs:
                        links = bs['Links']
                        if 'Link' in links:
                            link_list = links['Link'] if isinstance(links['Link'], list) else [links['Link']]
                            for link in link_list:
                                if 'target' in link and 'label' in link:
                                    attrs[f'biosample_link_{link["label"].lower().replace(" ", "_")}'] = link['target']
                    
                    rec = {'BioSample': acc}
                    rec.update(attrs)
                    records.append(rec)
    except Exception as e:
        print(f"Warning: Could not parse BioSample data: {e}")
    
    if records:
        df = pd.DataFrame(records).drop_duplicates(subset=['BioSample'])
        print(f"Loaded BioSample data with {len(df.columns)} columns")
        return df
    else:
        return pd.DataFrame()

def parse_bioproject_maximum(path):
    """Enhanced BioProject JSON parser with comprehensive project data"""
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
                
                # Extract ALL BioProject information
                if isinstance(obj, dict):
                    project_data = {}
                    if 'Project' in obj:
                        project = obj['Project']
                        
                        # Project ID and basic info
                        project_data['BioProject'] = project.get('ProjectID', {}).get('ArchiveID', {}).get('#text', '')
                        project_data['bioproject_id'] = project.get('ProjectID', {}).get('ID', {}).get('#text', '')
                        
                        # Project description
                        descr = project.get('ProjectDescr', {})
                        project_data['bioproject_title'] = descr.get('Title', '')
                        project_data['bioproject_description'] = descr.get('Description', '')
                        project_data['bioproject_type'] = project.get('ProjectType', {}).get('ProjectTypeSubmission', '')
                        
                        # Extract study information comprehensively
                        if 'Study' in project:
                            study = project['Study']
                            study_descr = study.get('Descriptor', {})
                            project_data['bioproject_study_title'] = study_descr.get('StudyTitle', '')
                            project_data['bioproject_study_description'] = study_descr.get('StudyAbstract', '')
                            project_data['bioproject_study_type'] = study_descr.get('StudyType', '')
                            
                            # Study attributes
                            if 'StudyAttributes' in study:
                                attrs = study['StudyAttributes']
                                if 'StudyAttribute' in attrs:
                                    attr_list = attrs['StudyAttribute'] if isinstance(attrs['StudyAttribute'], list) else [attrs['StudyAttribute']]
                                    for attr in attr_list:
                                        name = attr.get('@attribute_name', '')
                                        val = attr.get('#text', '')
                                        if name:
                                            key = f'bioproject_study_{name.lower().replace(" ", "_")}'
                                            project_data[key] = val
                        
                        # Extract organism information
                        if 'ProjectType' in project:
                            proj_type = project['ProjectType']
                            if 'Target' in proj_type:
                                target = proj_type['Target']
                                if 'Organism' in target:
                                    org = target['Organism']
                                    project_data['bioproject_target_organism'] = org.get('#text', '')
                                    project_data['bioproject_target_taxonomy_id'] = org.get('@taxonomy_id', '')
                        
                        records.append(project_data)
    except Exception as e:
        print(f"Warning: Could not parse BioProject data: {e}")
    
    if records:
        df = pd.DataFrame(records).drop_duplicates(subset=['BioProject'])
        print(f"Loaded BioProject data with {len(df.columns)} columns")
        return df
    else:
        return pd.DataFrame()

def parse_sra_xml_maximum(path):
    """Enhanced SRA XML parser with comprehensive metadata extraction"""
    records = []
    try:
        with open(path, 'r') as fh:
            for line in fh:
                line = line.strip()
                if not line: continue
                try:
                    obj = json.loads(line)
                    if obj:
                        # Add prefix to distinguish from other sources
                        prefixed_obj = {}
                        for key, value in obj.items():
                            prefixed_obj[f'sra_xml_{key.lower()}'] = value
                        records.append(prefixed_obj)
                except json.JSONDecodeError:
                    continue
    except Exception as e:
        print(f"Warning: Could not parse SRA XML data: {e}")
    
    if records:
        df = pd.DataFrame(records)
        print(f"Loaded SRA XML data with {len(df.columns)} columns")
        return df
    else:
        return pd.DataFrame()

def parse_geo_maximum(path):
    """Enhanced GEO metadata parser"""
    records = []
    try:
        with open(path, 'r') as fh:
            for line in fh:
                line = line.strip()
                if not line: continue
                try:
                    obj = json.loads(line)
                    if obj:
                        # Add prefix and clean up
                        prefixed_obj = {}
                        for key, value in obj.items():
                            clean_key = key.strip().lower().replace(' ', '_').replace('-', '_')
                            prefixed_obj[f'geo_{clean_key}'] = value
                        records.append(prefixed_obj)
                except json.JSONDecodeError:
                    continue
    except Exception as e:
        print(f"Warning: Could not parse GEO data: {e}")
    
    if records:
        df = pd.DataFrame(records)
        print(f"Loaded GEO data with {len(df.columns)} columns")
        return df
    else:
        return pd.DataFrame()

def parse_ffq_maximum(path):
    """Enhanced ffq JSON parser with comprehensive file information"""
    records = []
    try:
        with open(path, 'r') as fh:
            for line in fh:
                line = line.strip()
                if not line: continue
                try:
                    obj = json.loads(line)
                    if obj:
                        # Extract comprehensive file information
                        ffq_data = {}
                        for run_id, run_data in obj.items():
                            ffq_data['run_accession'] = run_id
                            
                            # Extract all available fields
                            for key, value in run_data.items():
                                if isinstance(value, dict):
                                    for sub_key, sub_value in value.items():
                                        ffq_data[f'ffq_{key}_{sub_key}'] = sub_value
                                else:
                                    ffq_data[f'ffq_{key}'] = value
                            
                            records.append(ffq_data)
                except json.JSONDecodeError:
                    continue
    except Exception as e:
        print(f"Warning: Could not parse ffq data: {e}")
    
    if records:
        df = pd.DataFrame(records)
        print(f"Loaded ffq data with {len(df.columns)} columns")
        return df
    else:
        return pd.DataFrame()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--indir", required=True, help="Directory with raw/ files")
    ap.add_argument("-o", "--outfile", required=True, help="Output TSV")
    args = ap.parse_args()
    indir = pathlib.Path(args.indir)

    print("=== MAXIMUM METADATA MERGE ===")
    print("Loading ALL data sources with ALL fields...")
    
    # Load all data sources
    runinfo = load_runinfo(indir / "runinfo.csv") if (indir / "runinfo.csv").exists() else pd.DataFrame()
    ena = load_ena_maximum(indir / "ena_read_run.tsv") if (indir / "ena_read_run.tsv").exists() else pd.DataFrame()
    bios = parse_biosample_maximum(indir / "biosample.jsonl") if (indir / "biosample.jsonl").exists() else pd.DataFrame()
    bioproj = parse_bioproject_maximum(indir / "bioproject.jsonl") if (indir / "bioproject.jsonl").exists() else pd.DataFrame()
    sra_xml = parse_sra_xml_maximum(indir / "sra_xml.jsonl") if (indir / "sra_xml.jsonl").exists() else pd.DataFrame()
    geo = parse_geo_maximum(indir / "geo_metadata.tsv") if (indir / "geo_metadata.tsv").exists() else pd.DataFrame()
    ffq = parse_ffq_maximum(indir / "ffq.jsonl") if (indir / "ffq.jsonl").exists() else pd.DataFrame()

    print(f"Loaded: RunInfo({len(runinfo)}), ENA({len(ena)}), BioSample({len(bios)}), BioProject({len(bioproj)}), SRA-XML({len(sra_xml)}), GEO({len(geo)}), ffq({len(ffq)})")

    # Validate that we have at least one data source
    if ena.empty and runinfo.empty:
        print("ERROR: No primary data sources found (ENA and RunInfo both empty)")
        return

    # Start with ENA data as primary source (most comprehensive)
    if not ena.empty:
        merged = ena.copy()
        print(f"Using ENA as primary source: {len(merged)} samples with {len(merged.columns)} columns")
    elif not runinfo.empty:
        merged = runinfo.copy()
        merged = merged.rename(columns={'Run': 'run_accession'})
        print(f"Using RunInfo as primary source: {len(merged)} samples")
    else:
        merged = pd.DataFrame()
    
    # Merge with RunInfo (if not already used as primary)
    if not runinfo.empty and not ena.empty:
        runinfo = runinfo.rename(columns={'Run': 'run_accession'})
        merged = merged.merge(runinfo, on='run_accession', how='left', suffixes=('', '_sra'))
        print(f"After RunInfo merge: {len(merged.columns)} columns")
    
    # Merge with BioSample
    if not bios.empty and 'BioSample' in merged.columns:
        merged = merged.merge(bios, on='BioSample', how='left', suffixes=('', '_bs'))
        print(f"After BioSample merge: {len(merged.columns)} columns")
    
    # Merge with BioProject
    if not bioproj.empty and 'BioProject' in merged.columns:
        merged = merged.merge(bioproj, on='BioProject', how='left', suffixes=('', '_bp'))
        print(f"After BioProject merge: {len(merged.columns)} columns")
    
    # Merge with SRA XML
    if not sra_xml.empty and 'run_accession' in sra_xml.columns:
        merged = merged.merge(sra_xml, on='run_accession', how='left', suffixes=('', '_xml'))
        print(f"After SRA XML merge: {len(merged.columns)} columns")
    
    # Merge with GEO (efficient approach)
    if not geo.empty and not merged.empty:
        # Create a GEO lookup dictionary
        geo_lookup = {}
        for _, geo_row in geo.iterrows():
            for col, val in geo_row.items():
                if isinstance(val, str) and val.startswith('GSE'):
                    geo_lookup[val] = geo_row.to_dict()
        
        # Add GEO columns to merged data
        for col in geo.columns:
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
        print(f"After GEO merge: {len(merged.columns)} columns")
    
    # Merge with ffq
    if not ffq.empty and 'run_accession' in ffq.columns:
        merged = merged.merge(ffq, on='run_accession', how='left', suffixes=('', '_ffq'))
        print(f"After ffq merge: {len(merged.columns)} columns")

    # Remove duplicate columns (keep first occurrence)
    merged = merged.loc[:, ~merged.columns.duplicated()]
    
    # Reorder important columns to the front
    important_cols = [
        'run_accession', 'experiment_accession', 'sample_accession', 'study_accession',
        'secondary_sample_accession', 'secondary_study_accession', 'BioSample', 'BioProject',
        'experiment_title', 'sample_title', 'study_title', 'sample_description',
        'library_name', 'library_strategy', 'library_source', 'library_layout',
        'instrument_model', 'read_count', 'base_count', 'scientific_name',
        'collection_date', 'first_public', 'last_updated'
    ]
    
    prefer = [c for c in important_cols if c in merged.columns]
    others = [c for c in merged.columns if c not in prefer]
    merged = merged[prefer + others]

    merged.to_csv(args.outfile, sep="\t", index=False)
    
    print(f"\n=== MAXIMUM METADATA MERGE COMPLETE ===")
    print(f"Final result: {len(merged)} samples with {len(merged.columns)} columns")
    print(f"Saved to: {args.outfile}")
    
    # Show column categories
    print(f"\nColumn breakdown:")
    print(f"  - ENA fields: {len([c for c in merged.columns if not c.startswith(('biosample_', 'bioproject_', 'geo_', 'sra_xml_', 'ffq_'))])}")
    print(f"  - BioSample fields: {len([c for c in merged.columns if c.startswith('biosample_')])}")
    print(f"  - BioProject fields: {len([c for c in merged.columns if c.startswith('bioproject_')])}")
    print(f"  - GEO fields: {len([c for c in merged.columns if c.startswith('geo_')])}")
    print(f"  - SRA XML fields: {len([c for c in merged.columns if c.startswith('sra_xml_')])}")
    print(f"  - ffq fields: {len([c for c in merged.columns if c.startswith('ffq_')])}")

if __name__ == "__main__":
    main()