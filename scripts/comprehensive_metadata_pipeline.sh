#!/usr/bin/env bash
# Comprehensive SRA Metadata Pipeline
# Downloads maximum metadata from all available sources and classifies cancer samples
# 
# Usage: ./comprehensive_metadata_pipeline.sh -i srr_list.txt [-o output_dir] [--classify]
#
# Features:
# - Downloads ALL 192 ENA fields (maximum metadata)
# - Collects BioSample, BioProject, GEO, SRA XML, ffq data
# - Merges all sources into comprehensive metadata table
# - Optional cancer classification (tumor/normal/cell line)
#
# Requirements:
# - Entrez Direct (esearch, efetch)
# - curl, jq, python3, pandas
# - ffq (optional, for file information)

set -euo pipefail

# Default values
OUTDIR="./output/comprehensive_metadata"
SRR_LIST=""
CLASSIFY=false

# Parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i) SRR_LIST="$2"; shift 2;;
    -o) OUTDIR="$2"; shift 2;;
    --classify) CLASSIFY=true; shift;;
    *) echo "Usage: $0 -i srr_list.txt [-o output_dir] [--classify]" >&2; exit 1;;
  esac
done

if [[ -z "${SRR_LIST}" || ! -f "${SRR_LIST}" ]]; then
  echo "ERROR: SRR list file required" >&2
  echo "Usage: $0 -i srr_list.txt [-o output_dir] [--classify]" >&2
  exit 1
fi

echo "=== COMPREHENSIVE SRA METADATA PIPELINE ==="
echo "Input: ${SRR_LIST}"
echo "Output: ${OUTDIR}"
echo "Cancer Classification: ${CLASSIFY}"
echo ""

# Create output directory
mkdir -p "${OUTDIR}"/{raw,logs}

# File paths
RUNINFO_CSV="${OUTDIR}/raw/runinfo.csv"
ENA_TSV="${OUTDIR}/raw/ena_read_run.tsv"
BIOSAMPLE_JSONL="${OUTDIR}/raw/biosample.jsonl"
BIOPROJECT_JSONL="${OUTDIR}/raw/bioproject.jsonl"
GEO_TSV="${OUTDIR}/raw/geo_metadata.tsv"
SRA_XML="${OUTDIR}/raw/sra_xml.jsonl"
FFQ_JSONL="${OUTDIR}/raw/ffq.jsonl"

# Clear files
: > "${RUNINFO_CSV}"
: > "${ENA_TSV}"
: > "${BIOSAMPLE_JSONL}"
: > "${BIOPROJECT_JSONL}"
: > "${GEO_TSV}"
: > "${SRA_XML}"
: > "${FFQ_JSONL}"

# Check dependencies
need() { command -v "$1" >/dev/null 2>&1 || { echo "Missing: $1" >&2; exit 1; }; }
need esearch; need efetch; need curl; need jq; need python3

# Clear proxy settings
unset https_proxy http_proxy HTTP_PROXY HTTPS_PROXY
export https_proxy="" http_proxy="" HTTP_PROXY="" HTTPS_PROXY=""

echo "[1/7] Fetching SRA RunInfo (via EDirect) …"
while read -r SRR; do
  [[ -z "${SRR}" ]] && continue
  esearch -db sra -query "${SRR}" | efetch -format runinfo >> "${RUNINFO_CSV}" || \
    echo "[WARN] RunInfo failed for ${SRR}" >&2
done < "${SRR_LIST}"

echo "[2/7] Fetching MAXIMUM ENA filereport (ALL 192 fields) …"
while read -r SRR; do
  [[ -z "${SRR}" ]] && continue
  curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${SRR}&result=read_run&fields=all&format=tsv&download=true" \
    >> "${ENA_TSV}" || echo "[WARN] ENA filereport failed for ${SRR}" >&2
done < "${SRR_LIST}"

echo "[3/7] Fetching BioSample JSON (clinical metadata) …"
BIOSAMPLES=$(awk -F',' 'NR==1{for(i=1;i<=NF;i++){if($i ~ /BioSample/i) col=i}} NR>1 && col && $col != "" && $col ~ /^SAM/{print $col}' "${RUNINFO_CSV}" | sort -u)
echo "Found BioSample accessions: ${BIOSAMPLES}"
for SAMN in ${BIOSAMPLES}; do
  [[ -z "${SAMN}" ]] && continue
  if [[ "${SAMN}" =~ ^SAM[NED] ]]; then
    echo "Fetching BioSample: ${SAMN}"
    esearch -db biosample -query "${SAMN}" | efetch -format json >> "${BIOSAMPLE_JSONL}" || \
      echo "[WARN] BioSample fetch failed for ${SAMN}" >&2
  fi
done

echo "[4/7] Fetching BioProject metadata …"
BIOPROJECTS=$(awk -F',' 'NR==1{for(i=1;i<=NF;i++){if($i ~ /BioProject/i) col=i}} NR>1 && col && $col != "" && $col ~ /^PRJ/{print $col}' "${RUNINFO_CSV}" | sort -u)
echo "Found BioProject accessions: ${BIOPROJECTS}"
for PRJ in ${BIOPROJECTS}; do
  [[ -z "${PRJ}" ]] && continue
  if [[ "${PRJ}" =~ ^PRJ ]]; then
    echo "Fetching BioProject: ${PRJ}"
    esearch -db bioproject -query "${PRJ}" | efetch -format json >> "${BIOPROJECT_JSONL}" || \
      echo "[WARN] BioProject fetch failed for ${PRJ}" >&2
  fi
done

echo "[5/7] Fetching SRA XML format (detailed metadata) …"
while read -r SRR; do
  [[ -z "${SRR}" ]] && continue
  echo "Fetching SRA XML for: ${SRR}"
  esearch -db sra -query "${SRR}" | efetch -format xml | python3 -c "
import sys, json, xml.etree.ElementTree as ET
try:
    root = ET.fromstring(sys.stdin.read())
    for run in root.findall('.//RUN'):
        run_data = {}
        for elem in run.iter():
            if elem.text and elem.text.strip():
                run_data[elem.tag] = elem.text.strip()
        print(json.dumps(run_data))
except Exception as e:
    print(f'Error parsing XML: {e}', file=sys.stderr)
    pass
" >> "${SRA_XML}" || echo "[WARN] SRA XML failed for ${SRR}" >&2
done < "${SRR_LIST}"

echo "[6/7] Fetching GEO metadata …"
GEO_ACCESSIONS=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++){if($i ~ /study_title/i) col=i}} NR>1 && col && $col ~ /GSE[0-9]+/{gsub(/.*GSE([0-9]+).*/, "GSE\\1", $col); print $col}' "${ENA_TSV}" | sort -u)
GEO_ACCESSIONS="$GEO_ACCESSIONS $(awk -F',' 'NR==1{for(i=1;i<=NF;i++){if($i ~ /Study_Pubmed_id/i) col=i}} NR>1 && col && $col ~ /GSE[0-9]+/{gsub(/.*GSE([0-9]+).*/, "GSE\\1", $col); print $col}' "${RUNINFO_CSV}" | sort -u)"
echo "Found GEO accessions: ${GEO_ACCESSIONS}"
for GSE in ${GEO_ACCESSIONS}; do
  [[ -z "${GSE}" ]] && continue
  if [[ "${GSE}" =~ ^GSE ]]; then
    echo "Fetching GEO: ${GSE}"
    curl -s "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${GSE}&targ=all&view=data&form=text" | \
      python3 -c "
import sys, json, re
geo_data = {}
for line in sys.stdin:
    if line.startswith('!') and '=' in line:
        key, value = line[1:].strip().split('=', 1)
        geo_data[key.strip()] = value.strip()
if geo_data:
    print(json.dumps(geo_data))
" >> "${GEO_TSV}" || echo "[WARN] GEO fetch failed for ${GSE}" >&2
  fi
done

echo "[7/7] Fetching ffq JSON (comprehensive file info) …"
if command -v ffq >/dev/null 2>&1; then
  while read -r SRR; do
    [[ -z "${SRR}" ]] && continue
    echo "Fetching ffq for: ${SRR}"
    ffq "${SRR}" >> "${FFQ_JSONL}" || echo "[WARN] ffq failed for ${SRR}" >&2
  done < "${SRR_LIST}"
else
  echo "[SKIP] ffq not found; install with: pip install ffq" >&2
fi

echo ""
echo "=== DATA COLLECTION COMPLETE ==="
echo "Data source summary:"
echo "  - RunInfo: $(wc -l < "${RUNINFO_CSV}") lines"
echo "  - ENA (ALL fields): $(wc -l < "${ENA_TSV}") lines"
echo "  - BioSample: $(wc -l < "${BIOSAMPLE_JSONL}") lines"
echo "  - BioProject: $(wc -l < "${BIOPROJECT_JSONL}") lines"
echo "  - GEO: $(wc -l < "${GEO_TSV}") lines"
echo "  - SRA XML: $(wc -l < "${SRA_XML}") lines"
echo "  - ffq: $(wc -l < "${FFQ_JSONL}") lines"

echo ""
echo "[8/8] Merging comprehensive metadata …"
python3 scripts/merge_metadata_maximum.py -i "${OUTDIR}/raw" -o "${OUTDIR}/comprehensive_metadata.tsv"

echo ""
echo "=== METADATA MERGE COMPLETE ==="
echo "Final result: $(wc -l < "${OUTDIR}/comprehensive_metadata.tsv") samples with $(head -1 "${OUTDIR}/comprehensive_metadata.tsv" | tr '\t' '\n' | wc -l) columns"
echo "Saved to: ${OUTDIR}/comprehensive_metadata.tsv"

# Optional cancer classification
if [[ "${CLASSIFY}" == "true" ]]; then
  echo ""
  echo "[9/9] Applying cancer classification …"
  python3 scripts/cancer_classification.py -i "${OUTDIR}/comprehensive_metadata.tsv" -o "${OUTDIR}/classified_metadata.tsv"
  echo "Cancer classification complete: ${OUTDIR}/classified_metadata.tsv"
fi

echo ""
echo "=== PIPELINE COMPLETE ==="
echo "Results:"
echo "  - Comprehensive metadata: ${OUTDIR}/comprehensive_metadata.tsv"
if [[ "${CLASSIFY}" == "true" ]]; then
  echo "  - Classified metadata: ${OUTDIR}/classified_metadata.tsv"
fi
echo ""
echo "Usage examples:"
echo "  # Basic metadata collection:"
echo "  $0 -i data/srr_list.txt"
echo ""
echo "  # With cancer classification:"
echo "  $0 -i data/srr_list.txt --classify"
echo ""
echo "  # Custom output directory:"
echo "  $0 -i data/srr_list.txt -o my_results --classify"