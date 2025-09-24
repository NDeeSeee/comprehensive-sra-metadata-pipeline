#!/usr/bin/env bash
# Comprehensive cancer metadata collection pipeline
# Collects metadata for cancer types from multiple sources

set -euo pipefail

# Default parameters
CANCER_TYPE=""
OUTDIR="./output/cancer_type_collection"
WITH_FFQ=0
WITH_GEO=0
WITH_XML=0
MAX_RESULTS=1000
TEST_MODE=0
SKIP_SEARCH=0
SKIP_MERGE=0
SKIP_CLASSIFY=0

# Help function
show_help() {
    cat << EOF
Cancer Metadata Collection Pipeline

USAGE:
    $0 -c "cancer type" [OPTIONS]

REQUIRED:
    -c, --cancer-type    Cancer type to search for (e.g., "esophageal adenocarcinoma")

OPTIONS:
    -o, --output-dir      Output directory (default: ./output/cancer_type_collection)
    --with-geo            Include GEO metadata collection
    --with-xml            Include SRA XML metadata collection  
    --with-ffq            Include ffq metadata collection
    --max-results         Maximum results per source (default: 1000)
    --test                Run in test mode (limited results)
    --skip-search         Skip cancer type search (use existing SRR list)
    --skip-merge          Skip metadata merging step
    --skip-classify       Skip cancer classification step
    -h, --help            Show this help message

EXAMPLES:
    # Basic collection for esophageal adenocarcinoma
    $0 -c "esophageal adenocarcinoma"
    
    # Comprehensive collection with all metadata sources
    $0 -c "lung squamous cell carcinoma" --with-geo --with-xml --with-ffq
    
    # Test mode with limited results
    $0 -c "breast cancer" --test --max-results 100
    
    # Skip search and use existing SRR list
    $0 -c "pancreatic cancer" --skip-search

WORKFLOW:
    1. Search for cancer type in SRA/ENA databases
    2. Extract SRR IDs from matching studies
    3. Collect comprehensive metadata from multiple sources
    4. Merge all metadata into unified dataset
    5. Apply cancer classification (tumor/normal/cell line)

OUTPUT FILES:
    - cancer_type_srr_list.txt     List of SRR IDs found
    - raw/                         Raw metadata from each source
    - ultimate_metadata.tsv        Merged metadata
    - classified_metadata.tsv      Metadata with cancer classifications

EOF
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    -c|--cancer-type)
      CANCER_TYPE="$2"
      shift 2
      ;;
    -o|--output-dir)
      OUTDIR="$2"
      shift 2
      ;;
    --with-geo)
      WITH_GEO=1
      shift
      ;;
    --with-xml)
      WITH_XML=1
      shift
      ;;
    --with-ffq)
      WITH_FFQ=1
      shift
      ;;
    --max-results)
      MAX_RESULTS="$2"
      shift 2
      ;;
    --test)
      TEST_MODE=1
      shift
      ;;
    --skip-search)
      SKIP_SEARCH=1
      shift
      ;;
    --skip-merge)
      SKIP_MERGE=1
      shift
      ;;
    --skip-classify)
      SKIP_CLASSIFY=1
      shift
      ;;
    -h|--help)
      show_help
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      echo "Use --help for usage information" >&2
      exit 1
      ;;
  esac
done

# Validate required parameters
if [[ -z "${CANCER_TYPE}" ]]; then
  echo "ERROR: Cancer type is required" >&2
  echo "Use --help for usage information" >&2
  exit 1
fi

# Create output directory
mkdir -p "${OUTDIR}"

echo "=========================================="
echo "Cancer Metadata Collection Pipeline"
echo "=========================================="
echo "Cancer type: ${CANCER_TYPE}"
echo "Output directory: ${OUTDIR}"
echo "Test mode: ${TEST_MODE}"
echo "Max results: ${MAX_RESULTS}"
echo "=========================================="

# Step 1: Search for cancer type (unless skipped)
if [[ ${SKIP_SEARCH} -eq 0 ]]; then
  echo ""
  echo "STEP 1: Searching for cancer type"
  echo "----------------------------------"
  
  SRR_LIST="${OUTDIR}/cancer_type_srr_list.txt"
  
  if [[ ${TEST_MODE} -eq 1 ]]; then
    /usr/bin/python3 scripts/demo_cancer_search.py -c "${CANCER_TYPE}" -o "${SRR_LIST}" --test
  else
    /usr/bin/python3 scripts/demo_cancer_search.py -c "${CANCER_TYPE}" -o "${SRR_LIST}" --max-results "${MAX_RESULTS}"
  fi
  
  if [[ ! -f "${SRR_LIST}" || ! -s "${SRR_LIST}" ]]; then
    echo "ERROR: No SRR IDs found for cancer type: ${CANCER_TYPE}"
    exit 1
  fi
  
  SRR_COUNT=$(wc -l < "${SRR_LIST}")
  echo "Found ${SRR_COUNT} SRR IDs for cancer type: ${CANCER_TYPE}"
else
  echo ""
  echo "STEP 1: Skipped (using existing SRR list)"
  echo "----------------------------------------"
  SRR_LIST="${OUTDIR}/cancer_type_srr_list.txt"
  if [[ ! -f "${SRR_LIST}" ]]; then
    echo "ERROR: SRR list file not found: ${SRR_LIST}"
    exit 1
  fi
  SRR_COUNT=$(wc -l < "${SRR_LIST}")
  echo "Using existing SRR list with ${SRR_COUNT} entries"
fi

# Step 2: Collect metadata
echo ""
echo "STEP 2: Collecting metadata from multiple sources"
echo "--------------------------------------------------"

# Build command for metadata collection
METADATA_CMD="scripts/build_metadata_by_cancer_type.sh -c \"${CANCER_TYPE}\" -o \"${OUTDIR}\""

if [[ ${WITH_GEO} -eq 1 ]]; then
  METADATA_CMD="${METADATA_CMD} --with-geo"
fi

if [[ ${WITH_XML} -eq 1 ]]; then
  METADATA_CMD="${METADATA_CMD} --with-xml"
fi

if [[ ${WITH_FFQ} -eq 1 ]]; then
  METADATA_CMD="${METADATA_CMD} --with-ffq"
fi

if [[ ${TEST_MODE} -eq 1 ]]; then
  METADATA_CMD="${METADATA_CMD} --test"
fi

METADATA_CMD="${METADATA_CMD} --max-results ${MAX_RESULTS}"

# Execute metadata collection
eval "${METADATA_CMD}"

# Step 3: Merge metadata (unless skipped)
if [[ ${SKIP_MERGE} -eq 0 ]]; then
  echo ""
  echo "STEP 3: Merging metadata from all sources"
  echo "------------------------------------------"
  
  ULTIMATE_METADATA="${OUTDIR}/ultimate_metadata.tsv"
  /usr/bin/python3 scripts/merge_metadata_cancer_classification.py -i "${OUTDIR}/raw" -o "${ULTIMATE_METADATA}"
  
  if [[ -f "${ULTIMATE_METADATA}" ]]; then
    METADATA_ROWS=$(wc -l < "${ULTIMATE_METADATA}")
    echo "Merged metadata saved: ${ULTIMATE_METADATA} (${METADATA_ROWS} rows)"
  else
    echo "ERROR: Failed to create merged metadata file"
    exit 1
  fi
else
  echo ""
  echo "STEP 3: Skipped (metadata merging)"
  echo "----------------------------------"
fi

# Step 4: Apply cancer classification (unless skipped)
if [[ ${SKIP_CLASSIFY} -eq 0 ]]; then
  echo ""
  echo "STEP 4: Applying cancer classification"
  echo "---------------------------------------"
  
  ULTIMATE_METADATA="${OUTDIR}/ultimate_metadata.tsv"
  CLASSIFIED_METADATA="${OUTDIR}/classified_metadata.tsv"
  
  if [[ -f "${ULTIMATE_METADATA}" ]]; then
    /usr/bin/python3 scripts/cancer_classification.py -i "${ULTIMATE_METADATA}" -o "${CLASSIFIED_METADATA}"
    
    if [[ -f "${CLASSIFIED_METADATA}" ]]; then
      echo "Cancer classification completed: ${CLASSIFIED_METADATA}"
    else
      echo "ERROR: Failed to create classified metadata file"
      exit 1
    fi
  else
    echo "ERROR: Ultimate metadata file not found: ${ULTIMATE_METADATA}"
    exit 1
  fi
else
  echo ""
  echo "STEP 4: Skipped (cancer classification)"
  echo "---------------------------------------"
fi

# Summary
echo ""
echo "=========================================="
echo "PIPELINE COMPLETED SUCCESSFULLY"
echo "=========================================="
echo "Cancer type: ${CANCER_TYPE}"
echo "SRR IDs found: ${SRR_COUNT}"
echo "Output directory: ${OUTDIR}"
echo ""
echo "Output files:"
echo "  - cancer_type_srr_list.txt     (${SRR_COUNT} SRR IDs)"
if [[ ${SKIP_MERGE} -eq 0 ]]; then
  echo "  - ultimate_metadata.tsv        (merged metadata)"
fi
if [[ ${SKIP_CLASSIFY} -eq 0 ]]; then
  echo "  - classified_metadata.tsv      (with cancer classifications)"
fi
echo "  - raw/                          (raw metadata from each source)"
echo ""
echo "Next steps:"
echo "  1. Review the classified_metadata.tsv file"
echo "  2. Filter samples based on cancer classifications"
echo "  3. Use the metadata for downstream analysis"
echo "=========================================="