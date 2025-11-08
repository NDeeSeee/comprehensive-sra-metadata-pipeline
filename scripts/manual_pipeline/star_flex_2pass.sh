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

  # Support comma-separated lists of FASTQs (from sample_list.with_status.txt)
  IFS=',' read -r -a fq1_list <<<"$fq1"
  # Trim whitespace tokens
  local cleaned_fq1_list=()
  for t in "${fq1_list[@]}"; do
    t="${t##[[:space:]]}"
    t="${t%%[[:space:]]}"
    [[ -n "$t" ]] && cleaned_fq1_list+=("$t")
  done

  local is_paired=0
  local cleaned_fq2_list=()
  if [[ -n "$fq2" ]]; then
    IFS=',' read -r -a fq2_list <<<"$fq2"
    for t in "${fq2_list[@]}"; do
      t="${t##[[:space:]]}"
      t="${t%%[[:space:]]}"
      [[ -n "$t" ]] && cleaned_fq2_list+=("$t")
    done
    if (( ${#cleaned_fq2_list[@]} > 0 )); then
      is_paired=1
    fi
  fi

  # Validate file lists and existence
  if (( ${#cleaned_fq1_list[@]} == 0 )); then
    echo "FASTQ1 list is empty for sample ${sample}" >&2
    exit 2
  fi
  if (( is_paired )) && (( ${#cleaned_fq2_list[@]} != ${#cleaned_fq1_list[@]} )); then
    echo "Mismatched R1/R2 counts for ${sample}: R1=${#cleaned_fq1_list[@]} R2=${#cleaned_fq2_list[@]}" >&2
    exit 3
  fi
  local i
  for (( i=0; i<${#cleaned_fq1_list[@]}; i++ )); do
    local p1="${cleaned_fq1_list[$i]}"
    if [[ ! -f "$p1" ]]; then
      echo "FASTQ1 not found: $p1" >&2
      exit 2
    fi
    if (( is_paired )); then
      local p2="${cleaned_fq2_list[$i]}"
      if [[ ! -f "$p2" ]]; then
        echo "FASTQ2 not found: $p2" >&2
        exit 3
      fi
    fi
  done

  echo "Running STAR 2-pass for sample: ${sample} ($([[ $is_paired -eq 1 ]] && echo paired-end || echo single-end))"

  # Build read arguments; STAR expects multiple files as separate args
  local -a reads_args
  if (( is_paired )); then
    # Interleave mates by STAR convention: all R1 then all R2
    reads_args=("${cleaned_fq1_list[@]}" "${cleaned_fq2_list[@]}")
  else
    reads_args=("${cleaned_fq1_list[@]}")
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

  # Determine compression: enable gz decompression only if ALL files end with .gz
  local all_gz=1
  for f in "${reads_args[@]}"; do
    [[ "$f" == *.gz ]] || { all_gz=0; break; }
  done
  if (( all_gz )); then
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
  local list_dir
  list_dir="$(cd "$(dirname "$sample_list")" && pwd)"

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
         -cwd "${list_dir}" \
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
