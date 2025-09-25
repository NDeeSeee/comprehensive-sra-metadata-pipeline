#!/usr/bin/env python3
"""
Real Cancer Type Search Script
Actually searches for different cancer types instead of using demo data
"""

import argparse
import requests
import time
import re
from pathlib import Path
from typing import Set, List

def search_sra_real(cancer_type: str, max_results: int = 1000) -> Set[str]:
    """
    Real SRA search for cancer types using E-utilities API
    """
    session = requests.Session()
    session.headers.update({
        'User-Agent': 'Cancer-Type-Metadata-Collector/1.0'
    })
    
    all_runs = set()
    
    # Create search variations for the specific cancer type
    search_terms = create_cancer_search_terms(cancer_type)
    
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
                for study_id in study_ids[:5]:  # Limit to first 5 studies
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

def create_cancer_search_terms(cancer_type: str) -> List[str]:
    """
    Create specific search terms for each cancer type
    """
    cancer_type = cancer_type.lower().strip()
    
    # Define specific search terms for each cancer type
    cancer_terms = {
        'gallbladder': [
            'gallbladder cancer',
            'gallbladder carcinoma',
            'gallbladder tumor',
            'gallbladder neoplasm',
            'gallbladder adenocarcinoma'
        ],
        'large cell carcinoma of the lung and bronchus': [
            'large cell lung cancer',
            'large cell carcinoma lung',
            'LCLC',
            'large cell neuroendocrine carcinoma',
            'lung large cell carcinoma'
        ],
        'renal pelvis': [
            'renal pelvis cancer',
            'renal pelvis carcinoma',
            'renal pelvis tumor',
            'upper urinary tract cancer',
            'renal pelvis adenocarcinoma'
        ],
        'small cell carcinoma of the lung and bronchus': [
            'small cell lung cancer',
            'SCLC',
            'small cell carcinoma lung',
            'oat cell carcinoma',
            'small cell neuroendocrine carcinoma'
        ],
        'salivary gland': [
            'salivary gland cancer',
            'salivary gland carcinoma',
            'salivary gland tumor',
            'parotid cancer',
            'submandibular cancer'
        ],
        'anus, anal canal & anorectum': [
            'anal cancer',
            'anal carcinoma',
            'anal canal cancer',
            'anorectal cancer',
            'anal squamous cell carcinoma'
        ],
        'small intestine': [
            'small intestine cancer',
            'small bowel cancer',
            'small intestine carcinoma',
            'jejunal cancer',
            'ileal cancer'
        ],
        'bones and joints': [
            'bone cancer',
            'osteosarcoma',
            'chondrosarcoma',
            'ewing sarcoma',
            'bone sarcoma'
        ],
        'meningioma of the brain/other nervous system': [
            'meningioma',
            'brain meningioma',
            'meningeal tumor',
            'meningioma brain',
            'meningeal neoplasm'
        ]
    }
    
    # Check if we have specific terms for this cancer type
    for key, terms in cancer_terms.items():
        if key in cancer_type:
            return terms
    
    # Default fallback - create generic terms
    base_terms = [
        cancer_type,
        cancer_type.replace('cancer', 'carcinoma'),
        cancer_type.replace('carcinoma', 'cancer'),
        cancer_type.replace(' ', '+'),
        cancer_type.replace(' ', '_')
    ]
    
    return list(set(base_terms))

def main():
    parser = argparse.ArgumentParser(description='Real cancer type search')
    parser.add_argument('-c', '--cancer-type', required=True, 
                       help='Cancer type to search for (e.g., "gallbladder cancer")')
    parser.add_argument('-o', '--output', required=True,
                       help='Output file for SRR IDs')
    parser.add_argument('-m', '--max-results', type=int, default=1000,
                       help='Maximum results (default: 1000)')
    parser.add_argument('--test', action='store_true',
                       help='Run in test mode with limited results')
    
    args = parser.parse_args()
    
    if args.test:
        args.max_results = 50
        print("Running in test mode (limited to 50 results)")
    
    # Search for cancer type
    runs = search_sra_real(args.cancer_type, args.max_results)
    
    if not runs:
        print(f"No runs found for cancer type: {args.cancer_type}")
        print("This could mean:")
        print("1. No studies exist for this cancer type")
        print("2. Search terms need refinement")
        print("3. API issues")
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