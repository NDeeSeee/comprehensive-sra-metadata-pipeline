#!/usr/bin/env python3
"""
Simple Cancer Type Search Script
Uses ENA API to search for cancer types and extract SRR IDs
"""

import argparse
import requests
import time
from pathlib import Path
from typing import Set

def search_ena_for_cancer_type(cancer_type: str, max_results: int = 1000) -> Set[str]:
    """
    Search ENA database for cancer type using multiple search strategies
    """
    session = requests.Session()
    session.headers.update({
        'User-Agent': 'Cancer-Type-Metadata-Collector/1.0'
    })
    
    all_runs = set()
    
    # Create search variations
    search_terms = [
        cancer_type.lower(),
        cancer_type.lower().replace('cancer', 'carcinoma'),
        cancer_type.lower().replace('carcinoma', 'cancer'),
        cancer_type.lower().replace(' ', '_'),
        cancer_type.lower().replace(' ', '-')
    ]
    
    # Remove duplicates while preserving order
    search_terms = list(dict.fromkeys(search_terms))
    
    print(f"Searching ENA for cancer type: {cancer_type}")
    print(f"Search variations: {search_terms}")
    
    for term in search_terms:
        print(f"Searching for: '{term}'")
        
        try:
            # ENA API search for read runs
            url = "https://www.ebi.ac.uk/ena/portal/api/search"
            params = {
                'dataPortal': 'ena',
                'result': 'read_run',
                'query': f'study_title="{term}" OR sample_title="{term}" OR experiment_title="{term}" OR study_description="{term}" OR sample_description="{term}"',
                'fields': 'run_accession,study_accession,sample_accession,experiment_accession,study_title,sample_title,experiment_title,study_description,sample_description',
                'limit': max_results,
                'format': 'json'
            }
            
            response = session.get(url, params=params, timeout=30)
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
            
            runs_found = 0
            for entry in entries:
                if isinstance(entry, dict) and 'run_accession' in entry:
                    all_runs.add(entry['run_accession'])
                    runs_found += 1
            
            print(f"  Found {runs_found} runs for '{term}'")
            time.sleep(1)  # Rate limiting
            
        except Exception as e:
            print(f"  Error searching for '{term}': {e}")
            continue
    
    return all_runs

def main():
    parser = argparse.ArgumentParser(description='Simple cancer type search using ENA API')
    parser.add_argument('-c', '--cancer-type', required=True, 
                       help='Cancer type to search for (e.g., "esophageal adenocarcinoma")')
    parser.add_argument('-o', '--output', required=True,
                       help='Output file for SRR IDs')
    parser.add_argument('-m', '--max-results', type=int, default=1000,
                       help='Maximum results per search term (default: 1000)')
    parser.add_argument('--test', action='store_true',
                       help='Run in test mode with limited results')
    
    args = parser.parse_args()
    
    if args.test:
        args.max_results = 50
        print("Running in test mode (limited to 50 results per search term)")
    
    # Search for cancer type
    runs = search_ena_for_cancer_type(args.cancer_type, args.max_results)
    
    if not runs:
        print(f"No runs found for cancer type: {args.cancer_type}")
        return
    
    # Write results
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_path, 'w') as f:
        for run in sorted(runs):
            f.write(f"{run}\n")
    
    print(f"\nFound {len(runs)} total runs for '{args.cancer_type}'")
    print(f"SRR IDs saved to: {output_path}")
    
    # Show some sample runs
    sample_runs = sorted(list(runs))[:10]
    print(f"Sample runs: {', '.join(sample_runs)}")

if __name__ == "__main__":
    main()