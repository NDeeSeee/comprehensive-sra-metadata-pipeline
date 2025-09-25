#!/usr/bin/env bash
# Batch processing script for cancer types with TCGA=0
# Runs cancer type search for each cancer type and organizes outputs

set -euo pipefail

# Configuration
BASE_OUTPUT_DIR="./output/zero_tcga_cancers"
CANCER_TYPES_FILE="data/zero_tcga_cancers.txt"
TEST_MODE=1
MAX_RESULTS=1000

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to clean cancer type name for directory naming
clean_name() {
    echo "$1" | sed 's/[^a-zA-Z0-9_-]/_/g' | sed 's/__*/_/g' | sed 's/^_\|_$//g' | tr '[:upper:]' '[:lower:]'
}

# Function to run cancer type search
run_cancer_search() {
    local cancer_type="$1"
    local clean_name="$2"
    local output_dir="$3"
    
    echo -e "${BLUE}Processing: $cancer_type${NC}"
    echo -e "${YELLOW}Output directory: $output_dir${NC}"
    
    # Create output directory
    mkdir -p "$output_dir"
    
    # Run the cancer type search
    if [[ ${TEST_MODE} -eq 1 ]]; then
        echo -e "${YELLOW}Running in TEST mode${NC}"
        scripts/collect_cancer_metadata.sh -c "$cancer_type" -o "$output_dir" --test --max-results "$MAX_RESULTS"
    else
        echo -e "${YELLOW}Running in FULL mode${NC}"
        scripts/collect_cancer_metadata.sh -c "$cancer_type" -o "$output_dir" --max-results "$MAX_RESULTS"
    fi
    
    local exit_code=$?
    
    if [[ $exit_code -eq 0 ]]; then
        echo -e "${GREEN}✅ SUCCESS: $cancer_type completed successfully${NC}"
        
        # Check if results were generated
        if [[ -f "$output_dir/cancer_type_srr_list.txt" ]]; then
            local srr_count=$(wc -l < "$output_dir/cancer_type_srr_list.txt")
            echo -e "${GREEN}   Found $srr_count SRR IDs${NC}"
        fi
        
        if [[ -f "$output_dir/classified_metadata.tsv" ]]; then
            local metadata_count=$(wc -l < "$output_dir/classified_metadata.tsv")
            echo -e "${GREEN}   Generated metadata with $metadata_count rows${NC}"
        fi
        
    else
        echo -e "${RED}❌ FAILED: $cancer_type failed with exit code $exit_code${NC}"
    fi
    
    echo -e "${BLUE}----------------------------------------${NC}"
    echo
}

# Function to generate summary report
generate_summary_report() {
    local summary_file="$BASE_OUTPUT_DIR/SUMMARY_REPORT.md"
    
    echo "# Zero TCGA Cancer Types - Processing Summary" > "$summary_file"
    echo "" >> "$summary_file"
    echo "**Date**: $(date)" >> "$summary_file"
    echo "**Test Mode**: $TEST_MODE" >> "$summary_file"
    echo "**Max Results**: $MAX_RESULTS" >> "$summary_file"
    echo "" >> "$summary_file"
    
    echo "## Processing Results" >> "$summary_file"
    echo "" >> "$summary_file"
    echo "| Cancer Type | Status | SRR IDs | Metadata Rows | Output Directory |" >> "$summary_file"
    echo "|-------------|--------|---------|---------------|-----------------|" >> "$summary_file"
    
    while IFS= read -r cancer_type; do
        [[ -z "$cancer_type" ]] && continue
        
        clean_name=$(clean_name "$cancer_type")
        output_dir="$BASE_OUTPUT_DIR/$clean_name"
        
        status="❌ Failed"
        srr_count="0"
        metadata_count="0"
        
        if [[ -f "$output_dir/cancer_type_srr_list.txt" ]]; then
            srr_count=$(wc -l < "$output_dir/cancer_type_srr_list.txt")
            status="✅ Success"
        fi
        
        if [[ -f "$output_dir/classified_metadata.tsv" ]]; then
            metadata_count=$(wc -l < "$output_dir/classified_metadata.tsv")
        fi
        
        echo "| $cancer_type | $status | $srr_count | $metadata_count | $clean_name |" >> "$summary_file"
        
    done < "$CANCER_TYPES_FILE"
    
    echo "" >> "$summary_file"
    echo "## Directory Structure" >> "$summary_file"
    echo "" >> "$summary_file"
    echo "\`\`\`" >> "$summary_file"
    tree "$BASE_OUTPUT_DIR" 2>/dev/null || find "$BASE_OUTPUT_DIR" -type d | sort >> "$summary_file"
    echo "\`\`\`" >> "$summary_file"
    
    echo -e "${GREEN}Summary report generated: $summary_file${NC}"
}

# Main execution
main() {
    echo -e "${BLUE}========================================${NC}"
    echo -e "${BLUE}Zero TCGA Cancer Types - Batch Processing${NC}"
    echo -e "${BLUE}========================================${NC}"
    echo
    
    # Check if cancer types file exists
    if [[ ! -f "$CANCER_TYPES_FILE" ]]; then
        echo -e "${RED}Error: $CANCER_TYPES_FILE not found${NC}"
        echo -e "${YELLOW}Please run: /usr/bin/python3 scripts/extract_zero_tcga_cancers.py${NC}"
        exit 1
    fi
    
    # Create base output directory
    mkdir -p "$BASE_OUTPUT_DIR"
    
    # Count total cancer types
    total_count=$(grep -c . "$CANCER_TYPES_FILE")
    echo -e "${YELLOW}Total cancer types to process: $total_count${NC}"
    echo -e "${YELLOW}Test mode: $TEST_MODE${NC}"
    echo -e "${YELLOW}Max results per cancer type: $MAX_RESULTS${NC}"
    echo
    
    # Process each cancer type
    current=0
    while IFS= read -r cancer_type; do
        [[ -z "$cancer_type" ]] && continue
        
        current=$((current + 1))
        clean_name=$(clean_name "$cancer_type")
        output_dir="$BASE_OUTPUT_DIR/$clean_name"
        
        echo -e "${BLUE}[$current/$total_count] Processing: $cancer_type${NC}"
        
        # Run the cancer search
        run_cancer_search "$cancer_type" "$clean_name" "$output_dir"
        
        # Add delay between requests to be respectful to APIs
        if [[ $current -lt $total_count ]]; then
            echo -e "${YELLOW}Waiting 5 seconds before next cancer type...${NC}"
            sleep 5
        fi
        
    done < "$CANCER_TYPES_FILE"
    
    echo -e "${GREEN}========================================${NC}"
    echo -e "${GREEN}Batch processing completed!${NC}"
    echo -e "${GREEN}========================================${NC}"
    
    # Generate summary report
    generate_summary_report
    
    echo -e "${YELLOW}All outputs saved under: $BASE_OUTPUT_DIR${NC}"
    echo -e "${YELLOW}Check the SUMMARY_REPORT.md for detailed results${NC}"
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --full)
            TEST_MODE=0
            shift
            ;;
        --max-results)
            MAX_RESULTS="$2"
            shift 2
            ;;
        --help)
            echo "Usage: $0 [--full] [--max-results N]"
            echo "  --full        Run in full mode (not test mode)"
            echo "  --max-results Set maximum results per cancer type (default: 1000)"
            echo "  --help        Show this help message"
            exit 0
            ;;
        *)
            echo "Unknown argument: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Run main function
main