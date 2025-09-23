#!/usr/bin/env python3
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
    """Enhanced BioSample JSON parser with better error handling"""
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

    print("Loading data sources...")
    
    # Load all data sources
    runinfo = load_runinfo(indir / "runinfo.csv")
    ena = load_ena(indir / "ena_read_run.tsv")
    bios = parse_biosample_jsonl(indir / "biosample.jsonl")
    bioproj = parse_bioproject_jsonl(indir / "bioproject.jsonl")
    sra_xml = parse_sra_xml_jsonl(indir / "sra_xml.jsonl")
    geo = parse_geo_tsv(indir / "geo_metadata.tsv")
    ffq = parse_ffq_jsonl(indir / "ffq.jsonl")

    print(f"Loaded: RunInfo({len(runinfo)}), ENA({len(ena)}), BioSample({len(bios)}), BioProject({len(bioproj)}), SRA-XML({len(sra_xml)}), GEO({len(geo)}), ffq({len(ffq)})")

    # Start with ENA data as primary source
    if not ena.empty:
        merged = ena.copy()
    else:
        merged = pd.DataFrame()
    
    # Merge with RunInfo
    if not runinfo.empty:
        runinfo = runinfo.rename(columns={'Run': 'run_accession'})
        merged = merged.merge(runinfo, on='run_accession', how='outer', suffixes=('', '_sra'))
    
    # Merge with BioSample
    if not bios.empty and 'BioSample' in merged.columns:
        merged = merged.merge(bios, on='BioSample', how='left', suffixes=('', '_bs'))
    
    # Merge with BioProject
    if not bioproj.empty and 'BioProject' in merged.columns:
        merged = merged.merge(bioproj, on='BioProject', how='left', suffixes=('', '_bp'))
    
    # Merge with SRA XML
    if not sra_xml.empty and 'run_accession' in sra_xml.columns:
        merged = merged.merge(sra_xml, on='run_accession', how='left', suffixes=('', '_xml'))
    
    # Merge with GEO
    if not geo.empty:
        # Try to match GEO data with study information
        if 'study_title' in merged.columns:
            for idx, row in merged.iterrows():
                study_title = str(row.get('study_title', ''))
                geo_match = geo[geo.apply(lambda x: any(gse in study_title for gse in x.values() if isinstance(gse, str) and gse.startswith('GSE')), axis=1)]
                if not geo_match.empty:
                    for col in geo_match.columns:
                        merged.at[idx, f'geo_{col}'] = geo_match.iloc[0][col]
    
    # Merge with ffq
    if not ffq.empty and 'run_accession' in ffq.columns:
        merged = merged.merge(ffq, on='run_accession', how='left', suffixes=('', '_ffq'))

    # Reorder important columns to the front
    prefer = [c for c in [
        'run_accession','Experiment','SampleName','BioSample','SRAStudy','BioProject',
        'study_accession','secondary_study_accession','experiment_accession','sample_accession',
        'secondary_sample_accession',
        'LibraryLayout','LibraryStrategy','LibrarySource','Platform','Model',
        'instrument_model','library_layout','library_strategy','library_source',
        'read_count','base_count','first_public','last_updated',
        'age','altitude','cell_line','cell_type','disease','environment_biome',
        'host','host_body_site','isolate','location','organism','phenotype',
        'sex','strain','temperature','tissue_type','treatment',
        'biosample_title','biosample_description','project_title','project_description'
    ] if c in merged.columns]
    
    others = [c for c in merged.columns if c not in prefer]
    merged = merged[prefer + others]

    # Remove duplicate columns (keep first occurrence)
    merged = merged.loc[:, ~merged.columns.duplicated()]

    merged.to_csv(args.outfile, sep="\t", index=False)
    print(f"Wrote: {args.outfile}  (rows={len(merged)}, columns={len(merged.columns)})")

if __name__ == "__main__":
    main()