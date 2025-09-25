#!/usr/bin/env python3
"""
Working Cancer Type Search Script
Uses SRA API to search for cancer types and extract SRR IDs
"""

import argparse
import requests
import time
import re
from pathlib import Path
from typing import Set, List

def search_sra_for_cancer_type(cancer_type: str, max_results: int = 1000) -> Set[str]:
    """
    Search SRA database for cancer type using E-utilities API
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
        cancer_type.lower().replace(' ', '+'),
        cancer_type.lower().replace(' ', '_'),
    ]
    
    # Remove duplicates while preserving order
    search_terms = list(dict.fromkeys(search_terms))
    
    print(f"Searching SRA for cancer type: {cancer_type}")
    print(f"Search variations: {search_terms}")
    
    for term in search_terms:
        print(f"Searching for: '{term}'")
        
        try:
            # Step 1: Search for studies
            search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            search_params = {
                'db': 'sra',
                'term': term,
                'retmax': max_results,
                'retmode': 'json',
                'sort': 'relevance'
            }
            
            response = session.get(search_url, params=search_params, timeout=30)
            response.raise_for_status()
            data = response.json()
            
            if 'esearchresult' in data and 'idlist' in data['esearchresult']:
                study_ids = data['esearchresult']['idlist']
                print(f"  Found {len(study_ids)} studies for '{term}'")
                
                # Step 2: Get run info for each study
                for study_id in study_ids[:10]:  # Limit to first 10 studies
                    try:
                        # Get run info for this study
                        fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
                        fetch_params = {
                            'db': 'sra',
                            'id': study_id,
                            'rettype': 'runinfo',
                            'retmode': 'text'
                        }
                        
                        fetch_response = session.get(fetch_url, params=fetch_params, timeout=30)
                        fetch_response.raise_for_status()
                        
                        lines = fetch_response.text.strip().split('\n')
                        if len(lines) > 1:  # Has header + data
                            for line in lines[1:]:
                                if line.strip():
                                    parts = line.split(',')
                                    if len(parts) > 0 and parts[0].startswith('SRR'):
                                        all_runs.add(parts[0])
                        
                        time.sleep(0.5)  # Rate limiting
                        
                    except Exception as e:
                        print(f"    Error fetching runs for study {study_id}: {e}")
                        continue
            
            time.sleep(1)  # Rate limiting between terms
            
        except Exception as e:
            print(f"  Error searching for '{term}': {e}")
            continue
    
    return all_runs

def search_with_existing_data(cancer_type: str, max_results: int = 1000) -> Set[str]:
    """
    Use existing SRR data for demonstration purposes
    """
    print(f"DEMO: Using existing data for cancer type: {cancer_type}")
    print("This demonstrates how the system would work with real database searches.")
    
    # Use existing SRR list for demo
    existing_srr_file = Path("data/srr_list.txt")
    
    if existing_srr_file.exists():
        with open(existing_srr_file, 'r') as f:
            all_srrs = [line.strip() for line in f if line.strip()]
        
        # Simulate different results based on cancer type
        if "esophageal" in cancer_type.lower() or "esophagus" in cancer_type.lower():
            demo_srrs = all_srrs[:max_results]
            print(f"Found {len(demo_srrs)} SRR IDs for esophageal cancer (demo)")
        elif "lung" in cancer_type.lower():
            # Simulate fewer results for lung cancer
            demo_srrs = all_srrs[:min(30, max_results)]
            print(f"Found {len(demo_srrs)} SRR IDs for lung cancer (demo)")
        elif "breast" in cancer_type.lower():
            # Simulate different results for breast cancer
            demo_srrs = all_srrs[:min(25, max_results)]
            print(f"Found {len(demo_srrs)} SRR IDs for breast cancer (demo)")
        else:
            # Default for other cancer types
            demo_srrs = all_srrs[:min(20, max_results)]
            print(f"Found {len(demo_srrs)} SRR IDs for {cancer_type} (demo)")
        
        return set(demo_srrs)
    else:
        print("No existing SRR data found for demo")
        return set()

def main():
    parser = argparse.ArgumentParser(description='Working cancer type search')
    parser.add_argument('-c', '--cancer-type', required=True, 
                       help='Cancer type to search for (e.g., "esophageal adenocarcinoma")')
    parser.add_argument('-o', '--output', required=True,
                       help='Output file for SRR IDs')
    parser.add_argument('-m', '--max-results', type=int, default=1000,
                       help='Maximum results (default: 1000)')
    parser.add_argument('--test', action='store_true',
                       help='Run in test mode with limited results')
    parser.add_argument('--demo', action='store_true',
                       help='Use demo mode with existing data')
    
    args = parser.parse_args()
    
    if args.test:
        args.max_results = 50
        print("Running in test mode (limited to 50 results)")
    
    # Search for cancer type
    if args.demo:
        runs = search_with_existing_data(args.cancer_type, args.max_results)
    else:
        runs = search_sra_for_cancer_type(args.cancer_type, args.max_results)
    
    if not runs:
        print(f"No runs found for cancer type: {args.cancer_type}")
        print("Try using --demo flag for demonstration with existing data")
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