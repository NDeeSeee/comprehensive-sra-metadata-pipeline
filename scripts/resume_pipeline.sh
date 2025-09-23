#!/usr/bin/env bash
# Resume pipeline from where it left off
# Usage: ./resume_pipeline.sh -i srr_list.txt -o output_dir [--with-geo] [--with-xml] [--with-ffq]

set -euo pipefail

OUTDIR="./meta_out"
SRR_LIST=""
WITH_FFQ=0
WITH_GEO=0
WITH_XML=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i) SRR_LIST="$2"; shift 2;;
    -o) OUTDIR="$2"; shift 2;;
    --with-ffq) WITH_FFQ=1; shift;;
    --with-geo) WITH_GEO=1; shift;;
    --with-xml) WITH_XML=1; shift;;
    *) echo "Unknown arg: $1"; exit 1;;
  esac
done

if [[ -z "${SRR_LIST}" || ! -f "${SRR_LIST}" ]]; then
  echo "Usage: $0 -i srr_list.txt [-o outdir] [--with-geo] [--with-xml] [--with-ffq]" >&2
  exit 1
fi

# Check if output directory exists
if [[ ! -d "${OUTDIR}" ]]; then
  echo "Error: Output directory ${OUTDIR} does not exist. Run the main pipeline first." >&2
  exit 1
fi

# Clear proxy settings
unset https_proxy http_proxy HTTP_PROXY HTTPS_PROXY
export https_proxy="" http_proxy="" HTTP_PROXY="" HTTPS_PROXY=""

# File paths
RUNINFO_CSV="${OUTDIR}/raw/runinfo.csv"
ENA_TSV="${OUTDIR}/raw/ena_read_run.tsv"
BIOSAMPLE_JSONL="${OUTDIR}/raw/biosample.jsonl"
FFQ_JSONL="${OUTDIR}/raw/ffq.jsonl"
GEO_TSV="${OUTDIR}/raw/geo_metadata.tsv"
SRA_XML="${OUTDIR}/raw/sra_xml.jsonl"
BIOPROJECT_JSONL="${OUTDIR}/raw/bioproject.jsonl"

# Function to find missing SRRs
find_missing_srrs() {
    local source_file="$1"
    local srr_list="$2"
    
    if [[ ! -f "${source_file}" ]]; then
        echo "Source file ${source_file} not found, processing all SRRs"
        cat "${srr_list}"
        return
    fi
    
    # Get processed SRRs from the source file
    local processed_srrs=""
    if [[ "${source_file}" == "${RUNINFO_CSV}" ]]; then
        processed_srrs=$(awk -F',' 'NR>1 && $1 != "" {print $1}' "${source_file}" 2>/dev/null || true)
    elif [[ "${source_file}" == "${ENA_TSV}" ]]; then
        processed_srrs=$(awk -F'\t' 'NR>1 && $1 != "" {print $1}' "${source_file}" 2>/dev/null || true)
    elif [[ "${source_file}" == "${FFQ_JSONL}" ]]; then
        processed_srrs=$(grep -o '"run_accession":"[^"]*"' "${source_file}" 2>/dev/null | sed 's/"run_accession":"\([^"]*\)"/\1/' || true)
    fi
    
    # Find missing SRRs
    if [[ -n "${processed_srrs}" ]]; then
        comm -23 <(sort "${srr_list}") <(echo "${processed_srrs}" | sort)
    else
        cat "${srr_list}"
    fi
}

echo "🔍 Analyzing incomplete samples..."

# Find missing SRRs for each source
MISSING_RUNINFO=$(find_missing_srrs "${RUNINFO_CSV}" "${SRR_LIST}")
MISSING_ENA=$(find_missing_srrs "${ENA_TSV}" "${SRR_LIST}")
MISSING_FFQ=$(find_missing_srrs "${FFQ_JSONL}" "${SRR_LIST}")

echo "📊 Missing samples analysis:"
echo "  RunInfo missing: $(echo "${MISSING_RUNINFO}" | wc -l) samples"
echo "  ENA missing: $(echo "${MISSING_ENA}" | wc -l) samples"
echo "  ffq missing: $(echo "${MISSING_FFQ}" | wc -l) samples"

# Resume RunInfo collection
if [[ -n "${MISSING_RUNINFO}" ]]; then
    echo "[RESUME 1/7] Fetching missing SRA RunInfo..."
    echo "${MISSING_RUNINFO}" | while read -r SRR; do
        [[ -z "${SRR}" ]] && continue
        echo "Processing: ${SRR}"
        esearch -db sra -query "${SRR}" | efetch -format runinfo >> "${RUNINFO_CSV}" || \
            echo "[WARN] RunInfo failed for ${SRR}" >&2
    done
fi

# Resume ENA collection
if [[ -n "${MISSING_ENA}" ]]; then
    echo "[RESUME 2/7] Fetching missing ENA filereport..."
    ENA_FIELDS="run_accession,experiment_accession,sample_accession,study_accession,secondary_sample_accession,secondary_study_accession,broker_name,center_name,experiment_title,library_name,library_layout,library_selection,library_strategy,library_source,instrument_model,read_count,base_count,first_public,last_updated,scientific_name,collection_date,study_title,sample_title,submitted_ftp,fastq_ftp,age,altitude,cell_line,cell_type,dev_stage,disease,environment_biome,environment_feature,environment_material,environmental_medium,environmental_sample,host,host_body_site,host_genotype,host_phenotype,isolate,location,sex,strain,temperature,tissue_type"
    
    echo "${MISSING_ENA}" | while read -r SRR; do
        [[ -z "${SRR}" ]] && continue
        echo "Processing: ${SRR}"
        curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${SRR}&result=read_run&fields=${ENA_FIELDS}&format=tsv&download=true" \
            >> "${ENA_TSV}" || echo "[WARN] ENA filereport failed for ${SRR}" >&2
    done
fi

# Resume BioSample collection
if [[ -n "${MISSING_RUNINFO}" ]]; then
    echo "[RESUME 3/7] Fetching missing BioSample data..."
    BIOSAMPLES=$(awk -F',' 'NR==1{for(i=1;i<=NF;i++){if($i ~ /BioSample/i) col=i}} NR>1 && col{print $col}' "${RUNINFO_CSV}" | sort -u)
    for SAMN in ${BIOSAMPLES}; do
        [[ -z "${SAMN}" ]] && continue
        if [[ "${SAMN}" =~ ^SAM[NED] ]]; then
            # Check if already processed
            if ! grep -q "${SAMN}" "${BIOSAMPLE_JSONL}" 2>/dev/null; then
                echo "Processing BioSample: ${SAMN}"
                esearch -db biosample -query "${SAMN}" | efetch -format json >> "${BIOSAMPLE_JSONL}" || \
                    echo "[WARN] BioSample fetch failed for ${SAMN}" >&2
            fi
        fi
    done
fi

# Resume BioProject collection
if [[ -n "${MISSING_RUNINFO}" ]]; then
    echo "[RESUME 4/7] Fetching missing BioProject data..."
    BIOPROJECTS=$(awk -F',' 'NR==1{for(i=1;i<=NF;i++){if($i ~ /BioProject/i) col=i}} NR>1 && col && $col != "" && $col ~ /^PRJ/{print $col}' "${RUNINFO_CSV}" | sort -u)
    for PRJ in ${BIOPROJECTS}; do
        [[ -z "${PRJ}" ]] && continue
        if [[ "${PRJ}" =~ ^PRJ ]]; then
            # Check if already processed
            if ! grep -q "${PRJ}" "${BIOPROJECT_JSONL}" 2>/dev/null; then
                echo "Processing BioProject: ${PRJ}"
                esearch -db bioproject -query "${PRJ}" | efetch -format json >> "${BIOPROJECT_JSONL}" || \
                    echo "[WARN] BioProject fetch failed for ${PRJ}" >&2
            fi
        fi
    done
fi

# Resume SRA XML collection
if [[ ${WITH_XML} -eq 1 ]] && [[ -n "${MISSING_RUNINFO}" ]]; then
    echo "[RESUME 5/7] Fetching missing SRA XML data..."
    echo "${MISSING_RUNINFO}" | while read -r SRR; do
        [[ -z "${SRR}" ]] && continue
        echo "Processing SRA XML: ${SRR}"
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
except:
    pass
" >> "${SRA_XML}" || echo "[WARN] SRA XML failed for ${SRR}" >&2
    done
fi

# Resume GEO collection
if [[ ${WITH_GEO} -eq 1 ]] && [[ -n "${MISSING_ENA}" ]]; then
    echo "[RESUME 6/7] Fetching missing GEO data..."
    GEO_ACCESSIONS=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++){if($i ~ /study_title/i) col=i}} NR>1 && col && $col ~ /GSE[0-9]+/{gsub(/.*GSE([0-9]+).*/, "GSE\\1", $col); print $col}' "${ENA_TSV}" | sort -u)
    GEO_ACCESSIONS="$GEO_ACCESSIONS $(awk -F',' 'NR==1{for(i=1;i<=NF;i++){if($i ~ /Study_Pubmed_id/i) col=i}} NR>1 && col && $col ~ /GSE[0-9]+/{gsub(/.*GSE([0-9]+).*/, "GSE\\1", $col); print $col}' "${RUNINFO_CSV}" | sort -u)"
    for GSE in ${GEO_ACCESSIONS}; do
        [[ -z "${GSE}" ]] && continue
        if [[ "${GSE}" =~ ^GSE ]]; then
            # Check if already processed
            if ! grep -q "${GSE}" "${GEO_TSV}" 2>/dev/null; then
                echo "Processing GEO: ${GSE}"
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
        fi
    done
fi

# Resume ffq collection
if [[ ${WITH_FFQ} -eq 1 ]] && [[ -n "${MISSING_FFQ}" ]]; then
    if command -v ffq >/dev/null 2>&1; then
        echo "[RESUME 7/7] Fetching missing ffq data..."
        echo "${MISSING_FFQ}" | while read -r SRR; do
            [[ -z "${SRR}" ]] && continue
            echo "Processing ffq: ${SRR}"
            ffq "${SRR}" >> "${FFQ_JSONL}" || echo "[WARN] ffq failed for ${SRR}" >&2
        done
    else
        echo "[SKIP] ffq not found; install with: pip install ffq" >&2
    fi
fi

echo "✅ Resume completed!"
echo "Next: python3 merge_metadata_enhanced.py -i ${OUTDIR}/raw -o ${OUTDIR}/ultimate_metadata.tsv"