#!/usr/bin/env bash
# Enhanced metadata collection specifically for cancer classification
# Collects additional fields needed for tumor/normal/cell line classification

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

mkdir -p "${OUTDIR}"/{raw,logs}

RUNINFO_CSV="${OUTDIR}/raw/runinfo.csv"
ENA_TSV="${OUTDIR}/raw/ena_read_run.tsv"
BIOSAMPLE_JSONL="${OUTDIR}/raw/biosample.jsonl"
FFQ_JSONL="${OUTDIR}/raw/ffq.jsonl"
GEO_TSV="${OUTDIR}/raw/geo_metadata.tsv"
SRA_XML="${OUTDIR}/raw/sra_xml.jsonl"
BIOPROJECT_JSONL="${OUTDIR}/raw/bioproject.jsonl"

: > "${RUNINFO_CSV}"
: > "${ENA_TSV}"
: > "${BIOSAMPLE_JSONL}"
: > "${FFQ_JSONL}"
: > "${GEO_TSV}"
: > "${SRA_XML}"
: > "${BIOPROJECT_JSONL}"

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

echo "[2/7] Fetching enhanced ENA filereport with cancer classification fields …"
# Enhanced ENA fields specifically for cancer classification
ENA_FIELDS="run_accession,experiment_accession,sample_accession,study_accession,secondary_sample_accession,secondary_study_accession,broker_name,center_name,experiment_title,library_name,library_layout,library_selection,library_strategy,library_source,instrument_model,read_count,base_count,first_public,last_updated,scientific_name,collection_date,study_title,sample_title,submitted_ftp,fastq_ftp,age,altitude,cell_line,cell_type,dev_stage,disease,environment_biome,environment_feature,environment_material,environmental_medium,environmental_sample,host,host_body_site,host_genotype,host_phenotype,isolate,location,sex,strain,temperature,tissue_type,treatment,genotype,phenotype,source_name,biomaterial_provider,organism_part,sampling_site,analyte_type,body_site,histological_type,disease_staging,is_tumor,subject_is_affected,individual,replicate,experimental_factor"

while read -r SRR; do
  [[ -z "${SRR}" ]] && continue
  curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${SRR}&result=read_run&fields=${ENA_FIELDS}&format=tsv&download=true" \
    >> "${ENA_TSV}" || echo "[WARN] ENA filereport failed for ${SRR}" >&2
done < "${SRR_LIST}"

echo "[3/7] Fetching enhanced BioSample JSON with clinical data …"
BIOSAMPLES=$(awk -F',' 'NR==1{for(i=1;i<=NF;i++){if($i ~ /BioSample/i) col=i}} NR>1 && col{print $col}' "${RUNINFO_CSV}" | sort -u)
for SAMN in ${BIOSAMPLES}; do
  [[ -z "${SAMN}" ]] && continue
  if [[ "${SAMN}" =~ ^SAM[NED] ]]; then
    esearch -db biosample -query "${SAMN}" | efetch -format json >> "${BIOSAMPLE_JSONL}" || \
      echo "[WARN] BioSample fetch failed for ${SAMN}" >&2
  fi
done

echo "[4/7] Fetching BioProject metadata …"
BIOPROJECTS=$(awk -F',' 'NR==1{for(i=1;i<=NF;i++){if($i ~ /BioProject/i) col=i}} NR>1 && col && $col != "" && $col ~ /^PRJ/{print $col}' "${RUNINFO_CSV}" | sort -u)
for PRJ in ${BIOPROJECTS}; do
  [[ -z "${PRJ}" ]] && continue
  if [[ "${PRJ}" =~ ^PRJ ]]; then
    esearch -db bioproject -query "${PRJ}" | efetch -format json >> "${BIOPROJECT_JSONL}" || \
      echo "[WARN] BioProject fetch failed for ${PRJ}" >&2
  fi
done

if [[ ${WITH_XML} -eq 1 ]]; then
  echo "[5/7] Fetching SRA XML format (detailed metadata) …"
  while read -r SRR; do
    [[ -z "${SRR}" ]] && continue
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
  done < "${SRR_LIST}"
fi

if [[ ${WITH_GEO} -eq 1 ]]; then
  echo "[6/7] Fetching GEO metadata …"
  GEO_ACCESSIONS=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++){if($i ~ /study_title/i) col=i}} NR>1 && col && $col ~ /GSE[0-9]+/{gsub(/.*GSE([0-9]+).*/, "GSE\\1", $col); print $col}' "${ENA_TSV}" | sort -u)
  GEO_ACCESSIONS="$GEO_ACCESSIONS $(awk -F',' 'NR==1{for(i=1;i<=NF;i++){if($i ~ /Study_Pubmed_id/i) col=i}} NR>1 && col && $col ~ /GSE[0-9]+/{gsub(/.*GSE([0-9]+).*/, "GSE\\1", $col); print $col}' "${RUNINFO_CSV}" | sort -u)"
  for GSE in ${GEO_ACCESSIONS}; do
    [[ -z "${GSE}" ]] && continue
    if [[ "${GSE}" =~ ^GSE ]]; then
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
fi

if [[ ${WITH_FFQ} -eq 1 ]]; then
  if command -v ffq >/dev/null 2>&1; then
    echo "[7/7] Fetching ffq JSON (fixed execution) …"
    while read -r SRR; do
      [[ -z "${SRR}" ]] && continue
      ffq "${SRR}" >> "${FFQ_JSONL}" || echo "[WARN] ffq failed for ${SRR}" >&2
    done < "${SRR_LIST}"
  else
    echo "[SKIP] ffq not found; install with: pip install ffq" >&2
  fi
fi

echo "Enhanced cancer classification metadata saved under: ${OUTDIR}/raw"
echo "Next: python3 merge_metadata_cancer_classification.py -i ${OUTDIR}/raw -o ${OUTDIR}/ultimate_metadata.tsv"