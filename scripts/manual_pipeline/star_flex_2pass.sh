#!/bin/bash
# Unified STAR 2-pass (flex) submitter/runner
# - If invoked with a single existing file: treats it as sample list and submits jobs (LSF)
# - If invoked with >=2 args: treats as a direct run for one sample: <SAMPLE> <FASTQ1> [FASTQ2]

set -euo pipefail

# Default resources (aligned with LSF submission)
THREADS=${STAR_THREADS:-2}

# Reference data (matches existing run_star-flex.sh)
GENOME_DIR=/data/salomonis2/Genomes/STAR-2.7.10b-Index-GRCH38/Grch38-STAR-index
GENOME=/data/salomonis2/Genomes/STAR-2.7.10b-Index-GRCH38/Grch38_r85.all.fa
GTF=/data/salomonis2/Genomes/STAR-2.7.10b-Index-GRCH38/Homo_sapiens.GRCh38.85.gtf

run_one_sample() {
  local sample="$1" fq1="$2" fq2="${3:-}"

  # Load modules if available (non-fatal if module not present)
  module load STAR/2.7.10b >/dev/null 2>&1 || true
  module load samtools >/dev/null 2>&1 || true
  module load bedtools >/dev/null 2>&1 || true

  local root_dir
  root_dir="$PWD"
  mkdir -p "${root_dir}/bams"

  # Skip recomputation if final BAM already exists and is non-empty
  local final_bam="${root_dir}/bams/${sample}.bam"
  if [[ -e "$final_bam" ]]; then
    if [[ -s "$final_bam" ]]; then
      echo "Output BAM already exists and is non-empty: $final_bam — skipping."
      return 0
    else
      echo "Zero-size BAM detected at $final_bam — removing and recalculating."
      rm -f "$final_bam" || true
    fi
  fi

  if [[ ! -f "$fq1" ]]; then
    echo "FASTQ1 not found: $fq1" >&2
    exit 2
  fi
  local is_paired=0
  if [[ -n "$fq2" ]]; then
    if [[ ! -f "$fq2" ]]; then
      echo "FASTQ2 provided but not found: $fq2" >&2
      exit 3
    fi
    is_paired=1
  fi

  echo "Running STAR 2-pass for sample: ${sample} ($([[ $is_paired -eq 1 ]] && echo paired-end || echo single-end))"

  # Build read arguments
  local -a reads_args
  if [[ $is_paired -eq 1 ]]; then
    reads_args=("$fq1" "$fq2")
  else
    reads_args=("$fq1")
  fi

  # Use decompression only for gzip-compressed inputs
  local -a star_cmd
  star_cmd=(
    STAR
    --runThreadN "${THREADS}"
    --genomeDir "${GENOME_DIR}"
    --readFilesIn "${reads_args[@]}"
    --outFileNamePrefix "${root_dir}/${sample}_"
    --outSAMtype BAM SortedByCoordinate
    --outSAMunmapped Within
    --outSAMattributes NH HI NM MD AS XS
    --outSAMstrandField intronMotif
    --twopassMode Basic
    --limitBAMsortRAM 200000000000
    --outFilterMultimapScoreRange 1
    --outFilterMultimapNmax 20
    --outFilterMismatchNmax 10
    --outFilterMatchNminOverLread 0.33
    --outFilterScoreMinOverLread 0.33
    --alignIntronMax 500000
    --alignMatesGapMax 1000000
    --alignSJDBoverhangMin 1
    --sjdbGTFfile "${GTF}"
    --sjdbOverhang 100
  )

  if [[ "$fq1" == *.gz ]] || ([[ $is_paired -eq 1 ]] && [[ "$fq2" == *.gz ]]); then
    star_cmd+=( --readFilesCommand "gunzip -c" )
  fi

  "${star_cmd[@]}"

  # Move BAM to final location
  if [[ -f "${root_dir}/${sample}_Aligned.sortedByCoord.out.bam" ]]; then
    mv "${root_dir}/${sample}_Aligned.sortedByCoord.out.bam" "${root_dir}/bams/${sample}.bam"
  fi

  # Keep Log.final.out; clean large intermediates
  rm -f "${root_dir}/${sample}_SJ.out.tab" \
        "${root_dir}/${sample}_Log.out" \
        "${root_dir}/${sample}_Log.progress.out" || true

  echo "${sample} STAR 2-pass complete"
}

submit_from_list() {
  local sample_list="$1"
  if [[ ! -f "$sample_list" ]]; then
    echo "Sample list not found: $sample_list" >&2
    exit 1
  fi
  mkdir -p bams logs

  # Submit one job per line; support tab- or whitespace-delimited lines; skip comments/empties
  while IFS=$'\t' read -r SAMPLE FQ1 FQ2 || [[ -n "${SAMPLE}" ]]; do
    if [[ -z "${FQ1:-}" ]]; then
      # Fallback: split by whitespace if tabs not used
      read -r SAMPLE FQ1 FQ2 <<<"${SAMPLE}"
      [[ -z "${FQ1:-}" ]] && continue
    fi
    [[ -z "${SAMPLE}" ]] && continue
    [[ "${SAMPLE}" =~ ^# ]] && continue

    # Submit; the run mode will validate file existence and pairedness
    bsub -W 12:00 -n "${THREADS}" -M 128000 \
         -R "rusage[mem=16000] span[hosts=1]" \
         -J "align_${SAMPLE}" \
         -o "logs/STAR2pass_${SAMPLE}.out" \
         -e "logs/STAR2pass_${SAMPLE}.err" \
         "$0" "${SAMPLE}" "${FQ1}" ${FQ2:+"${FQ2}"}
  done < "$sample_list"
}

main() {
  if [[ $# -eq 1 && -f "$1" ]]; then
    submit_from_list "$1"
    exit 0
  fi
  if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <sample_list.tsv>  OR  $0 <SAMPLE> <FASTQ1> [FASTQ2]" >&2
    exit 1
  fi
  run_one_sample "$@"
}

main "$@"
