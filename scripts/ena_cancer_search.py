#!/usr/bin/env python3
"""
ENA Cancer Type Search Script
Uses ENA filereport API to search for cancer types
"""

import argparse
import requests
import time
from pathlib import Path
from typing import Set, List
import re

def search_ena_filereport(cancer_type: str, max_results: int = 1000) -> Set[str]:
    """
    Search ENA using filereport API with cancer type terms
    """
    session = requests.Session()
    session.headers.update({
        'User-Agent': 'Cancer-Type-Metadata-Collector/1.0'
    })
    
    all_runs = set()
    
    # Create search variations - simpler terms work better
    search_terms = [
        cancer_type.lower(),
        cancer_type.lower().replace('cancer', 'carcinoma'),
        cancer_type.lower().replace('carcinoma', 'cancer'),
        cancer_type.lower().split()[0] if ' ' in cancer_type else cancer_type.lower(),  # First word only
    ]
    
    # Remove duplicates while preserving order
    search_terms = list(dict.fromkeys(search_terms))
    
    print(f"Searching ENA filereport for cancer type: {cancer_type}")
    print(f"Search variations: {search_terms}")
    
    for term in search_terms:
        print(f"Searching for: '{term}'")
        
        try:
            # Use ENA filereport API with broader search
            url = "https://www.ebi.ac.uk/ena/portal/api/filereport"
            params = {
                'result': 'read_run',
                'fields': 'run_accession,study_accession,sample_accession,experiment_accession,study_title,sample_title,experiment_title,study_description,sample_description',
                'format': 'tsv',
                'download': 'true',
                'limit': max_results
            }
            
            # Try different search strategies
            search_strategies = [
                f'study_title="{term}"',
                f'sample_title="{term}"',
                f'experiment_title="{term}"',
                f'study_description="{term}"',
                f'sample_description="{term}"',
                term,  # Simple term search
            ]
            
            for strategy in search_strategies:
                try:
                    params['query'] = strategy
                    response = session.get(url, params=params, timeout=30)
                    response.raise_for_status()
                    
                    lines = response.text.strip().split('\n')
                    if len(lines) > 1:  # Has header + data
                        runs_found = 0
                        for line in lines[1:]:
                            if line.strip():
                                parts = line.split('\t')
                                if len(parts) > 0 and parts[0].startswith('SRR'):
                                    all_runs.add(parts[0])
                                    runs_found += 1
                        
                        if runs_found > 0:
                            print(f"  Strategy '{strategy}': Found {runs_found} runs")
                            break  # Found results, move to next term
                    
                    time.sleep(0.5)  # Rate limiting
                    
                except Exception as e:
                    print(f"  Strategy '{strategy}' failed: {e}")
                    continue
            
            time.sleep(1)  # Rate limiting between terms
            
        except Exception as e:
            print(f"  Error searching for '{term}': {e}")
            continue
    
    return all_runs

def search_ena_broad(cancer_type: str, max_results: int = 1000) -> Set[str]:
    """
    Broader search using ENA API with multiple approaches
    """
    session = requests.Session()
    session.headers.update({
        'User-Agent': 'Cancer-Type-Metadata-Collector/1.0'
    })
    
    all_runs = set()
    
    # Extract key terms from cancer type
    terms = cancer_type.lower().split()
    key_terms = [term for term in terms if len(term) > 3]  # Skip short words
    
    print(f"Broad ENA search for terms: {key_terms}")
    
    for term in key_terms:
        print(f"Searching for term: '{term}'")
        
        try:
            # Try different ENA endpoints
            endpoints = [
                {
                    'url': 'https://www.ebi.ac.uk/ena/portal/api/filereport',
                    'params': {
                        'result': 'read_run',
                        'fields': 'run_accession,study_title,sample_title',
                        'format': 'tsv',
                        'download': 'true',
                        'limit': max_results
                    }
                }
            ]
            
            for endpoint in endpoints:
                try:
                    # Try exact match
                    endpoint['params']['query'] = f'study_title="{term}"'
                    response = session.get(endpoint['url'], params=endpoint['params'], timeout=30)
                    
                    if response.status_code == 200:
                        lines = response.text.strip().split('\n')
                        if len(lines) > 1:
                            runs_found = 0
                            for line in lines[1:]:
                                if line.strip():
                                    parts = line.split('\t')
                                    if len(parts) > 0 and parts[0].startswith('SRR'):
                                        all_runs.add(parts[0])
                                        runs_found += 1
                            
                            if runs_found > 0:
                                print(f"  Found {runs_found} runs for '{term}'")
                                break
                    
                    time.sleep(0.5)
                    
                except Exception as e:
                    print(f"  Endpoint failed for '{term}': {e}")
                    continue
            
            time.sleep(1)
            
        except Exception as e:
            print(f"  Error with term '{term}': {e}")
            continue
    
    return all_runs

def main():
    parser = argparse.ArgumentParser(description='ENA cancer type search')
    parser.add_argument('-c', '--cancer-type', required=True, 
                       help='Cancer type to search for (e.g., "esophageal adenocarcinoma")')
    parser.add_argument('-o', '--output', required=True,
                       help='Output file for SRR IDs')
    parser.add_argument('-m', '--max-results', type=int, default=1000,
                       help='Maximum results per search (default: 1000)')
    parser.add_argument('--test', action='store_true',
                       help='Run in test mode with limited results')
    parser.add_argument('--broad', action='store_true',
                       help='Use broader search strategy')
    
    args = parser.parse_args()
    
    if args.test:
        args.max_results = 50
        print("Running in test mode (limited to 50 results per search)")
    
    # Search for cancer type
    if args.broad:
        runs = search_ena_broad(args.cancer_type, args.max_results)
    else:
        runs = search_ena_filereport(args.cancer_type, args.max_results)
    
    if not runs:
        print(f"No runs found for cancer type: {args.cancer_type}")
        print("Try using --broad flag for broader search")
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