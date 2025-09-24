#!/usr/bin/env python3
"""
Cancer Type Search Script
Searches SRA/ENA databases for cancer types and extracts associated SRR IDs
"""

import argparse
import requests
import json
import time
import re
from pathlib import Path
import pandas as pd
from typing import List, Dict, Set

class CancerTypeSearcher:
    def __init__(self):
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'Cancer-Type-Metadata-Collector/1.0'
        })
        
    def normalize_cancer_type(self, cancer_type: str) -> List[str]:
        """
        Normalize cancer type name to generate multiple search variations
        """
        cancer_type = cancer_type.strip().lower()
        
        # Common cancer type mappings and variations
        variations = [cancer_type]
        
        # Add common synonyms and variations
        synonyms = {
            'esophageal adenocarcinoma': ['esophageal adenocarcinoma', 'eac', 'esophagus adenocarcinoma', 'oesophageal adenocarcinoma'],
            'lung squamous cell carcinoma': ['lung squamous cell carcinoma', 'lung scc', 'pulmonary squamous cell carcinoma', 'lung squamous carcinoma'],
            'breast cancer': ['breast cancer', 'breast carcinoma', 'mammary cancer'],
            'colorectal cancer': ['colorectal cancer', 'colon cancer', 'rectal cancer', 'colorectal carcinoma'],
            'prostate cancer': ['prostate cancer', 'prostate carcinoma', 'prostatic cancer'],
            'ovarian cancer': ['ovarian cancer', 'ovarian carcinoma', 'ovary cancer'],
            'pancreatic cancer': ['pancreatic cancer', 'pancreatic adenocarcinoma', 'pancreas cancer'],
            'liver cancer': ['liver cancer', 'hepatocellular carcinoma', 'hcc', 'hepatic cancer'],
            'kidney cancer': ['kidney cancer', 'renal cell carcinoma', 'rcc', 'renal cancer'],
            'bladder cancer': ['bladder cancer', 'bladder carcinoma', 'urothelial carcinoma'],
            'brain cancer': ['brain cancer', 'glioma', 'glioblastoma', 'brain tumor'],
            'melanoma': ['melanoma', 'skin cancer', 'cutaneous melanoma'],
            'leukemia': ['leukemia', 'leukaemia', 'blood cancer'],
            'lymphoma': ['lymphoma', 'hodgkin lymphoma', 'non-hodgkin lymphoma']
        }
        
        # Check if input matches any known cancer types
        for key, values in synonyms.items():
            if any(v in cancer_type for v in values) or cancer_type in values:
                variations.extend(values)
                break
        
        # Add generic cancer terms
        if 'cancer' not in cancer_type and 'carcinoma' not in cancer_type and 'tumor' not in cancer_type:
            variations.extend([f"{cancer_type} cancer", f"{cancer_type} carcinoma", f"{cancer_type} tumor"])
        
        # Remove duplicates and return
        return list(set(variations))
    
    def search_sra_studies(self, cancer_type: str, max_results: int = 1000) -> List[Dict]:
        """
        Search SRA database for studies matching cancer type
        """
        print(f"Searching SRA for: {cancer_type}")
        
        # Construct search query
        search_terms = self.normalize_cancer_type(cancer_type)
        query_parts = []
        
        for term in search_terms[:5]:  # Limit to avoid overly complex queries
            query_parts.append(f'"{term}"')
        
        query = " OR ".join(query_parts)
        
        # SRA E-utilities API search
        base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        
        # Search for studies
        search_url = f"{base_url}esearch.fcgi"
        search_params = {
            'db': 'sra',
            'term': query,
            'retmax': max_results,
            'retmode': 'json',
            'sort': 'relevance'
        }
        
        try:
            response = self.session.get(search_url, params=search_params, timeout=30)
            response.raise_for_status()
            data = response.json()
            
            if 'esearchresult' in data and 'idlist' in data['esearchresult']:
                study_ids = data['esearchresult']['idlist']
                print(f"Found {len(study_ids)} studies for '{cancer_type}'")
                return study_ids
            else:
                print(f"No studies found for '{cancer_type}'")
                return []
                
        except Exception as e:
            print(f"Error searching SRA: {e}")
            return []
    
    def get_study_runs(self, study_ids: List[str]) -> Set[str]:
        """
        Get all run accessions (SRR IDs) from study IDs using SRA API
        """
        all_runs = set()
        
        for study_id in study_ids:
            print(f"Fetching runs for study {study_id}...")
            
            try:
                # Use SRA API to get run info
                url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
                params = {
                    'db': 'sra',
                    'id': study_id,
                    'rettype': 'runinfo',
                    'retmode': 'text'
                }
                
                response = self.session.get(url, params=params, timeout=30)
                response.raise_for_status()
                
                lines = response.text.strip().split('\n')
                if len(lines) > 1:  # Has header + data
                    for line in lines[1:]:
                        if line.strip():
                            parts = line.split(',')
                            if len(parts) > 0 and parts[0].startswith('SRR'):
                                all_runs.add(parts[0])
                
                time.sleep(0.5)  # Rate limiting
                
            except Exception as e:
                print(f"Error fetching runs for study {study_id}: {e}")
                continue
        
        return all_runs
    
    def search_ena_studies(self, cancer_type: str, max_results: int = 1000) -> Set[str]:
        """
        Search ENA database for studies matching cancer type
        """
        print(f"Searching ENA for: {cancer_type}")
        
        search_terms = self.normalize_cancer_type(cancer_type)
        all_runs = set()
        
        for term in search_terms[:3]:  # Limit to avoid rate limiting
            try:
                # ENA API search
                url = "https://www.ebi.ac.uk/ena/portal/api/search"
                params = {
                    'dataPortal': 'ena',
                    'result': 'read_run',
                    'query': f'study_title="{term}" OR sample_title="{term}" OR experiment_title="{term}"',
                    'fields': 'run_accession,study_accession,sample_accession,experiment_accession,study_title,sample_title,experiment_title',
                    'limit': max_results,
                    'format': 'json'
                }
                
                response = self.session.get(url, params=params, timeout=30)
                response.raise_for_status()
                data = response.json()
                
                # Handle different response formats
                entries = []
                if isinstance(data, dict):
                    if 'entries' in data:
                        entries = data['entries']
                    elif 'data' in data:
                        entries = data['data']
                elif isinstance(data, list):
                    entries = data
                
                for entry in entries:
                    if isinstance(entry, dict) and 'run_accession' in entry:
                        all_runs.add(entry['run_accession'])
                
                print(f"Found {len(entries)} ENA entries for '{term}'")
                time.sleep(1)  # Rate limiting
                
            except Exception as e:
                print(f"Error searching ENA for '{term}': {e}")
                continue
        
        return all_runs
    
    def search_cancer_type(self, cancer_type: str, max_results: int = 1000) -> Set[str]:
        """
        Main search function that combines SRA and ENA searches
        """
        print(f"Searching for cancer type: {cancer_type}")
        print(f"Maximum results per source: {max_results}")
        
        all_runs = set()
        
        # Search SRA
        study_ids = self.search_sra_studies(cancer_type, max_results)
        if study_ids:
            sra_runs = self.get_study_runs(study_ids)
            all_runs.update(sra_runs)
            print(f"SRA search yielded {len(sra_runs)} runs")
        
        # Search ENA
        ena_runs = self.search_ena_studies(cancer_type, max_results)
        all_runs.update(ena_runs)
        print(f"ENA search yielded {len(ena_runs)} runs")
        
        print(f"Total unique runs found: {len(all_runs)}")
        return all_runs

def main():
    parser = argparse.ArgumentParser(description='Search for cancer types and extract SRR IDs')
    parser.add_argument('-c', '--cancer-type', required=True, 
                       help='Cancer type to search for (e.g., "esophageal adenocarcinoma")')
    parser.add_argument('-o', '--output', required=True,
                       help='Output file for SRR IDs')
    parser.add_argument('-m', '--max-results', type=int, default=1000,
                       help='Maximum results per source (default: 1000)')
    parser.add_argument('--test', action='store_true',
                       help='Run in test mode with limited results')
    
    args = parser.parse_args()
    
    if args.test:
        args.max_results = 50
        print("Running in test mode (limited to 50 results per source)")
    
    searcher = CancerTypeSearcher()
    
    # Search for cancer type
    runs = searcher.search_cancer_type(args.cancer_type, args.max_results)
    
    if not runs:
        print(f"No runs found for cancer type: {args.cancer_type}")
        return
    
    # Write results
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_path, 'w') as f:
        for run in sorted(runs):
            f.write(f"{run}\n")
    
    print(f"Found {len(runs)} runs for '{args.cancer_type}'")
    print(f"SRR IDs saved to: {output_path}")
    
    # Show some sample runs
    sample_runs = sorted(list(runs))[:10]
    print(f"Sample runs: {', '.join(sample_runs)}")

if __name__ == "__main__":
    main()