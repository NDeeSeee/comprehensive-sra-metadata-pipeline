#!/usr/bin/env bash
# Wrapper to invoke the Python cohort builder (discovery → QC → metadata → merge)
# Replaces legacy per-accession curl/esearch loops.
#
# Usage:
#   ./comprehensive_metadata_pipeline.sh -c "Anus" -o cancer_analysis_output/Anus/metadata [-y data/cancer_terms.yml] [-m 1000]
# Notes:
#   - The -i flag from the legacy pipeline is accepted but ignored (deprecated).

set -euo pipefail

OUTDIR=""
CANCER_LABEL=""
YAML="data/cancer_terms.yml"
MAX_RESULTS="1000"

while [[ $# -gt 0 ]]; do
  case "$1" in
    -o) OUTDIR="$2"; shift 2;;
    -c) CANCER_LABEL="$2"; shift 2;;
    -y) YAML="$2"; shift 2;;
    -m) MAX_RESULTS="$2"; shift 2;;
    -i) # deprecated - was SRR list; discovery now handled in Python
        shift 2;;
    *) echo "Usage: $0 -c <cancer_label> -o <output_dir> [-y cancer_terms.yml] [-m max_results]" >&2; exit 1;;
  esac
done

if [[ -z "${OUTDIR}" || -z "${CANCER_LABEL}" ]]; then
  echo "ERROR: Both -c <cancer_label> and -o <output_dir> are required" >&2
  exit 1
fi

mkdir -p "${OUTDIR}"

echo "=== COHORT BUILDER WRAPPER ==="
echo "Cancer label: ${CANCER_LABEL}"
echo "Output dir:  ${OUTDIR}"
echo "YAML terms:  ${YAML}"
echo "Max results: ${MAX_RESULTS}"

echo "Running cohort builder ..."
python scripts/automated_pipeline/cohort_builder.py \
  -c "${CANCER_LABEL}" \
  -o "${OUTDIR}" \
  -y "${YAML}" \
  -m "${MAX_RESULTS}"
RC=$?

if [[ ${RC} -ne 0 ]]; then
  echo "Cohort builder exited with code ${RC}" >&2
  exit ${RC}
fi

echo "Cohort builder completed. Outputs saved to ${OUTDIR}"
