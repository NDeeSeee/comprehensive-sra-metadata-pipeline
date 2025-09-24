#!/usr/bin/env python3
"""
Demo Cancer Type Search Script
Demonstrates the cancer type search functionality using existing data
"""

import argparse
from pathlib import Path
from typing import Set

def demo_cancer_search(cancer_type: str, max_results: int = 1000) -> Set[str]:
    """
    Demo function that simulates cancer type search using existing SRR data
    """
    print(f"DEMO: Searching for cancer type: {cancer_type}")
    print("This is a demonstration of how the cancer type search would work.")
    print("In a real implementation, this would search SRA/ENA databases.")
    
    # For demo purposes, use existing SRR list if it contains esophageal data
    existing_srr_file = Path("data/srr_list.txt")
    
    if existing_srr_file.exists():
        print(f"Using existing SRR data from: {existing_srr_file}")
        
        with open(existing_srr_file, 'r') as f:
            all_srrs = [line.strip() for line in f if line.strip()]
        
        # Simulate filtering based on cancer type
        if "esophageal" in cancer_type.lower() or "esophagus" in cancer_type.lower():
            # Use all existing SRRs for esophageal demo
            demo_srrs = all_srrs[:max_results]
            print(f"Found {len(demo_srrs)} SRR IDs for esophageal cancer (demo)")
        else:
            # For other cancer types, simulate finding fewer results
            demo_srrs = all_srrs[:min(50, max_results)]
            print(f"Found {len(demo_srrs)} SRR IDs for {cancer_type} (demo)")
        
        return set(demo_srrs)
    else:
        print("No existing SRR data found for demo")
        return set()

def main():
    parser = argparse.ArgumentParser(description='Demo cancer type search')
    parser.add_argument('-c', '--cancer-type', required=True, 
                       help='Cancer type to search for (e.g., "esophageal adenocarcinoma")')
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
    runs = demo_cancer_search(args.cancer_type, args.max_results)
    
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