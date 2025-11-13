#!/bin/bash

# Flexible STAR 2-pass alignment runner (supports single- and paired-end)
# Usage: ./run_star-flex.sh <SAMPLE> <FASTQ1> [FASTQ2]

set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <SAMPLE> <FASTQ1> [FASTQ2]" >&2
  exit 1
fi

SAMPLE=$1
FASTQ1=$2
FASTQ2=${3:-}

# Load STAR and common tools (adjust for your environment)
module load STAR/2.7.10b || true
module load samtools || true
module load bedtools || true

# Paths (copied from existing run_star-new.sh)
GENOME_DIR=/data/salomonis2/Genomes/STAR-2.7.10b-Index-GRCH38/Grch38-STAR-index
GENOME=/data/salomonis2/Genomes/STAR-2.7.10b-Index-GRCH38/Grch38_r85.all.fa
GTF=/data/salomonis2/Genomes/STAR-2.7.10b-Index-GRCH38/Homo_sapiens.GRCh38.85.gtf
ROOT_DIR=$PWD

mkdir -p "${ROOT_DIR}/bams"

if [[ ! -f "${FASTQ1}" ]]; then
  echo "FASTQ1 not found: ${FASTQ1}" >&2
  exit 2
fi

# Determine layout: paired if FASTQ2 provided AND exists; otherwise single-end
IS_PAIRED=0
if [[ -n "${FASTQ2}" && -f "${FASTQ2}" ]]; then
  IS_PAIRED=1
fi

echo "Running STAR 2-pass alignment for sample: ${SAMPLE} ("$( (( IS_PAIRED )) && echo paired-end || echo single-end )")"

# Build readFilesIn argument
if (( IS_PAIRED )); then
  READS_ARGS=("${FASTQ1}" "${FASTQ2}")
else
  READS_ARGS=("${FASTQ1}")
fi

STAR \
  --runThreadN 8 \
  --genomeDir "${GENOME_DIR}" \
  --readFilesIn "${READS_ARGS[@]}" \
  --readFilesCommand "gunzip -c" \
  --outFileNamePrefix "${ROOT_DIR}/${SAMPLE}_" \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMunmapped Within \
  --outSAMattributes NH HI NM MD AS XS \
  --outSAMstrandField intronMotif \
  --twopassMode Basic \
  --limitBAMsortRAM 200000000000 \
  --outFilterMultimapScoreRange 1 \
  --outFilterMultimapNmax 20 \
  --outFilterMismatchNmax 10 \
  --outFilterMatchNminOverLread 0.33 \
  --outFilterScoreMinOverLread 0.33 \
  --alignIntronMax 500000 \
  --alignMatesGapMax 1000000 \
  --alignSJDBoverhangMin 1 \
  --sjdbGTFfile "${GTF}" \
  --sjdbOverhang 100

# Move BAM to final location
if [[ -f "${ROOT_DIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam" ]]; then
  mv "${ROOT_DIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam" "${ROOT_DIR}/bams/${SAMPLE}.bam"
fi

# Optional cleanup of large intermediates
rm -f "${ROOT_DIR}/${SAMPLE}_SJ.out.tab" \
      "${ROOT_DIR}/${SAMPLE}_Log.final.out" \
      "${ROOT_DIR}/${SAMPLE}_Log.out" \
      "${ROOT_DIR}/${SAMPLE}_Log.progress.out" || true

echo "${SAMPLE} STAR 2-pass complete"
