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
  # Returns: R1, R2, STATUS (for paired-end) or R1, "", STATUS (for single-end)
  _read_sample_fields() {
    local sfile="$1" sid="$2"
    awk -v s="$sid" 'BEGIN{FS=OFS="\t"}
      $1==s {
        if (NF == 4 && $3 == "") {
          # Single-end: col2=R1, col3=empty, col4=status; use placeholder for R2
          print $2, "-", $4
        } else if (NF == 4) {
          # Paired-end: col2=R1, col3=R2, col4=status
          print $2, $3, $4
        } else if (NF == 3) {
          # Legacy single-end format: col2=R1, col3=status
          print $2, "-", $3
        } else {
          # Unknown format, output as-is
          print $2, $3, $4
        }
        exit
      }' "$sfile"
  }

  # Helper: check if sample is ready for submission; returns 0 if ready, 1 to skip, 2 to wait
  _check_sample_status() {
    local sfile="$1" sid="$2"
    local line
    line="$(_read_sample_fields "$sfile" "$sid")"
    # If not found or no status col, proceed immediately
    if [[ -z "$line" ]]; then
      echo "${sid}: no status info; proceeding"
      return 0
    fi
    local r1 r2 st
    IFS=$'\t' read -r r1 r2 st <<<"$line"

    # Skip if already done or if DBGaP required
    if [[ "$st" == "BAM_DONE" ]]; then
      echo "${sid}: BAM_DONE; skipping"
      return 1
    fi
    if [[ "$st" == "DBGaP_REQUIRED" ]]; then
      echo "${sid}: DBGaP_REQUIRED; skipping"
      return 1
    fi

    # Ready if FASTQ_DONE
    if [[ "$st" == "FASTQ_DONE" ]]; then
      return 0
    fi

    # Otherwise wait
    return 2
  }

  # First pass: submit all samples that are ready (FASTQ_DONE and not BAM_DONE)
  echo "=== First pass: submitting all ready samples ==="
  local submitted_count=0
  local skipped_count=0
  local pending_count=0

  # CRITICAL FIX: Read entire file into array BEFORE modifying it in the loop
  local -a sample_lines
  mapfile -t sample_lines < "$sample_list"

  for line in "${sample_lines[@]}"; do
    # Parse line into variables
    IFS=$'\t' read -r SAMPLE FQ1 FQ2 _STATUS <<<"$line" || continue
    if [[ -z "${FQ1:-}" ]]; then
      # Fallback: split by whitespace if tabs not used
      read -r SAMPLE FQ1 FQ2 _STATUS <<<"${SAMPLE}"
      [[ -z "${FQ1:-}" ]] && continue
    fi
    [[ -z "${SAMPLE}" ]] && continue
    [[ "${SAMPLE}" =~ ^# ]] && continue

    # Check status: 0=ready, 1=skip, 2=wait
    local status_code=0
    _check_sample_status "$sample_list" "$SAMPLE" || status_code=$?

    if [[ $status_code -eq 1 ]]; then
      # Skip (BAM_DONE or DBGaP_REQUIRED)
      skipped_count=$((skipped_count + 1))
      continue
    elif [[ $status_code -eq 2 ]]; then
      # Not ready yet, will check in wait loop
      pending_count=$((pending_count + 1))
      continue
    fi

    # Re-read the latest R1/R2 fields right before submission (they may have changed)
    local latest_line
    latest_line="$(_read_sample_fields "$sample_list" "$SAMPLE")"
    if [[ -n "$latest_line" ]]; then
      IFS=$'\t' read -r FQ1 FQ2 _STATUS <<<"$latest_line"
      # Convert placeholder "-" back to empty string for single-end
      [[ "$FQ2" == "-" ]] && FQ2=""
    fi

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

    # Mark as in progress immediately after submission to prevent re-submission
    local tmp="${sample_list}.tmp.$$"
    awk -v s="$SAMPLE" 'BEGIN{FS=OFS="\t"} { if ($1==s) { $4="BAM_IN_PROGRESS" } print $0 }' "$sample_list" > "$tmp" && mv "$tmp" "$sample_list" || true

    submitted_count=$((submitted_count + 1))
  done

  echo "First pass complete: submitted=${submitted_count}, skipped=${skipped_count}, pending=${pending_count}"

  # Wait loop: continuously check for samples that become ready
  if [[ $pending_count -eq 0 ]]; then
    echo "No pending samples to monitor. Exiting."
    return 0
  fi

  echo "=== Monitoring for $pending_count pending samples that become ready ==="
  while true; do
    local any_pending=0
    local any_submitted=0

    # Read file into array to avoid modification-during-read issues
    local -a sample_lines_inner
    mapfile -t sample_lines_inner < "$sample_list"

    for line in "${sample_lines_inner[@]}"; do
      IFS=$'\t' read -r SAMPLE FQ1 FQ2 _STATUS <<<"$line" || continue
      if [[ -z "${FQ1:-}" ]]; then
        read -r SAMPLE FQ1 FQ2 _STATUS <<<"${SAMPLE}"
        [[ -z "${FQ1:-}" ]] && continue
      fi
      [[ -z "${SAMPLE}" ]] && continue
      [[ "${SAMPLE}" =~ ^# ]] && continue

      local status_code=0
      _check_sample_status "$sample_list" "$SAMPLE" >/dev/null || status_code=$?

      if [[ $status_code -eq 1 ]]; then
        # Skip (already done or excluded)
        continue
      elif [[ $status_code -eq 2 ]]; then
        # Still waiting
        any_pending=1
        continue
      fi

      # Newly ready! Submit it
      local latest_line
      latest_line="$(_read_sample_fields "$sample_list" "$SAMPLE")"
      if [[ -n "$latest_line" ]]; then
        IFS=$'\t' read -r FQ1 FQ2 _STATUS <<<"$latest_line"
        # Convert placeholder "-" back to empty string for single-end
        [[ "$FQ2" == "-" ]] && FQ2=""
      fi

      echo "Submit ${SAMPLE}: status=FASTQ_DONE (newly ready)"
      STAR_STATUS_FILE="${sample_list}" \
      bsub -W 12:00 -n "${THREADS}" -M 128000 \
           -R "rusage[mem=16000] span[hosts=1]" \
           -J "align_${SAMPLE}" \
           -o "logs/STAR2pass_${SAMPLE}.out" \
           -e "logs/STAR2pass_${SAMPLE}.err" \
           -cwd "${list_dir}" \
           "$0" "${SAMPLE}" "${FQ1}" ${FQ2:+"${FQ2}"}

      # Mark as in progress immediately after submission to prevent re-submission
      local tmp="${sample_list}.tmp.$$"
      awk -v s="$SAMPLE" 'BEGIN{FS=OFS="\t"} { if ($1==s) { $4="BAM_IN_PROGRESS" } print $0 }' "$sample_list" > "$tmp" && mv "$tmp" "$sample_list" || true

      any_submitted=$((any_submitted + 1))
    done

    # If nothing pending and nothing submitted in this round, we're done
    if [[ $any_pending -eq 0 ]]; then
      echo "All samples processed or skipped. Exiting."
      break
    fi

    if [[ $any_submitted -eq 0 ]]; then
      echo "Waiting for samples to become ready... (checking every ${wait_interval}s)"
      sleep "$wait_interval"
    fi
  done
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
