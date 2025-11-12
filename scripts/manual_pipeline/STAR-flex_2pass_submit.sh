#!/bin/bash
# ================================================================
# STAR 2-pass alignment LSF submission script (flexible layout)
# - Supports sample list with either 2 columns (SAMPLE FASTQ1)
#   or 3 columns (SAMPLE FASTQ1 FASTQ2)
# - Chooses single- vs paired-end automatically based on FASTQ2 presence
# ================================================================
set -euo pipefail

# Create output directories
mkdir -p bams logs

# Paths
ROOT=$PWD
SAMPLE_LIST=${1:-"${ROOT}/sample_list.txt"}

if [[ ! -f "${SAMPLE_LIST}" ]]; then
  echo "Sample list not found: ${SAMPLE_LIST}" >&2
  exit 1
fi

# Submit one job per line (skip empty lines and lines starting with #)
while IFS=$'\t' read -r SAMPLE FQ1 FQ2 || [[ -n "${SAMPLE}" ]]; do
  # Handle whitespace-delimited as well
  if [[ -z "${FQ1:-}" ]]; then
    # try splitting by whitespace if tabs not used
    read -r SAMPLE FQ1 FQ2 <<<"${SAMPLE}"
    [[ -z "${FQ1:-}" ]] && continue
  fi
  [[ -z "${SAMPLE}" ]] && continue
  [[ "${SAMPLE}" =~ ^# ]] && continue

  # Normalize to files relative to current dir if not absolute
  F1_PATH="${FQ1}"
  F2_PATH="${FQ2:-}"

  # Decide paired vs single based on existence of FASTQ2
  if [[ -n "${F2_PATH}" && -f "${F2_PATH}" ]]; then
    bsub -W 12:00 -n 2 -M 128000 \
         -R "rusage[mem=16000] span[hosts=1]" \
         -J "align_${SAMPLE}" \
         -o "logs/STAR2pass_${SAMPLE}.out" \
         -e "logs/STAR2pass_${SAMPLE}.err" \
         "${ROOT}/run_star-flex.sh" "${SAMPLE}" "${F1_PATH}" "${F2_PATH}"
  else
    bsub -W 12:00 -n 2 -M 128000 \
         -R "rusage[mem=16000] span[hosts=1]" \
         -J "align_${SAMPLE}" \
         -o "logs/STAR2pass_${SAMPLE}.out" \
         -e "logs/STAR2pass_${SAMPLE}.err" \
         "${ROOT}/run_star-flex.sh" "${SAMPLE}" "${F1_PATH}"
  fi

done < "${SAMPLE_LIST}"

# End of script
