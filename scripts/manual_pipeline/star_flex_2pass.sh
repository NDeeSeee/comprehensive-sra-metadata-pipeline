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

  # Optional status tracking file passed via environment (set by submitter)
  local status_file="${STAR_STATUS_FILE:-}"

  update_status() {
    local new_status="$1"
    if [[ -n "$status_file" && -f "$status_file" ]]; then
      local tmp="${status_file}.tmp.$$"
      awk -v s="$sample" -v st="$new_status" 'BEGIN{FS=OFS="\t"} { if ($1==s) { $4=st } print $0 }' "$status_file" > "$tmp" && mv "$tmp" "$status_file" || true
    fi
  }

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
      update_status "BAM_DONE"
      return 0
    else
      echo "Zero-size BAM detected at $final_bam — removing and recalculating."
      rm -f "$final_bam" || true
    fi
  fi

  # Mark status as in-progress at start
  update_status "BAM_IN_PROGRESS"

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

  # Run STAR; on failure, mark error and exit
  if ! "${star_cmd[@]}"; then
    echo "STAR failed for sample ${sample}" >&2
    update_status "BAM_ERROR"
    exit 4
  fi

  # Move BAM to final location
  if [[ -f "${root_dir}/${sample}_Aligned.sortedByCoord.out.bam" ]]; then
    mv "${root_dir}/${sample}_Aligned.sortedByCoord.out.bam" "${root_dir}/bams/${sample}.bam"
  fi

  # Validate final BAM and set status accordingly
  if [[ -s "${root_dir}/bams/${sample}.bam" ]]; then
    update_status "BAM_DONE"
  else
    echo "Alignment finished but BAM missing or zero-size for ${sample}" >&2
    update_status "BAM_ERROR"
    exit 5
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
  local wait_interval
  wait_interval=${STAR_WAIT_INTERVAL_SEC:-60}

  # Helper: fetch current columns for a sample from the status file
  _read_sample_fields() {
    local sfile="$1" sid="$2"
    awk -v s="$2" 'BEGIN{FS=OFS="\t"} $1==s {print $2"\t"$3"\t"$4; exit}' "$1"
  }

  # Helper: wait until FASTQ_DONE, or skip on DBGaP_REQUIRED; returns 0 when done, 1 to skip
  _wait_for_fastq_done() {
    local sfile="$1" sid="$2"
    while true; do
      local line
      line="$(_read_sample_fields "$sfile" "$sid")"
      # If not found or no status col, proceed immediately
      if [[ -z "$line" ]]; then
        echo "${sid}: no status info; proceeding"
        return 0
      fi
      local r1 r2 st
      IFS=$'\t' read -r r1 r2 st <<<"$line"
      if [[ "$st" == "FASTQ_DONE" ]]; then
        return 0
      fi
      if [[ "$st" == "DBGaP_REQUIRED" ]]; then
        echo "${sid}: DBGaP_REQUIRED; skipping"
        return 1
      fi
      echo "${sid}: waiting for FASTQ_DONE (current=${st:-unknown})..."
      sleep "$wait_interval"
    done
  }

  # Submit one job per line; support tab- or whitespace-delimited lines; skip comments/empties
  # Read up to 4 columns (SAMPLE, R1, R2, STATUS). STATUS may be used to gate submission.
  while IFS=$'\t' read -r SAMPLE FQ1 FQ2 _STATUS || [[ -n "${SAMPLE}" ]]; do
    if [[ -z "${FQ1:-}" ]]; then
      # Fallback: split by whitespace if tabs not used
      read -r SAMPLE FQ1 FQ2 _STATUS <<<"${SAMPLE}"
      [[ -z "${FQ1:-}" ]] && continue
    fi
    [[ -z "${SAMPLE}" ]] && continue
    [[ "${SAMPLE}" =~ ^# ]] && continue

    # If status file has statuses, wait for FASTQ_DONE; skip DBGaP_REQUIRED
    if ! _wait_for_fastq_done "$sample_list" "$SAMPLE"; then
      continue
    fi

    # Re-read the latest R1/R2 fields right before submission (they may have changed)
    {
      latest_line="$(_read_sample_fields "$sample_list" "$SAMPLE")"
      if [[ -n "$latest_line" ]]; then
        IFS=$'\t' read -r FQ1 FQ2 _STATUS <<<"$latest_line"
      fi
    }

    echo "Submit ${SAMPLE}: status=FASTQ_DONE"

    # Submit; the run mode will validate file existence and pairedness. Pass status file via env
    STAR_STATUS_FILE="${sample_list}" \
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
