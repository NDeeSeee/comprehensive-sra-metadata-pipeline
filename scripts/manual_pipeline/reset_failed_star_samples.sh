#!/bin/bash
# Reset samples that had STAR failures so they can be reprocessed
# Enhanced with job status verification for BAM_IN_PROGRESS samples

if [ $# -ne 1 ]; then
  echo "Usage: $0 <tumor_directory>"
  echo "Example: $0 /data/salomonis-archive/FASTQs/NCI-R01/POSEIDON/Tumors/LungLargeCell"
  exit 1
fi

TUMOR_DIR="$1"
cd "$TUMOR_DIR" || exit 1

echo "Working directory: $PWD"
echo ""

# Function to check LSF job status for a sample
check_job_status() {
  local sample="$1"

  # Check if bjobs is available
  if ! command -v bjobs &>/dev/null; then
    echo "BJOBS_NOT_AVAILABLE"
    return
  fi

  # Get detailed job info with full job names (no truncation)
  # Job names are like: align_SAMN34721927
  local job_output
  job_output=$(bjobs -a -l 2>/dev/null | grep -A 5 "Job Name <align_$sample>")

  if [[ -z "$job_output" ]]; then
    echo "NOT_FOUND"
    return
  fi

  # Extract status from the line containing "Status <RUN>" or "Status <DONE>", etc.
  local status
  status=$(echo "$job_output" | grep -oP 'Status <\K[^>]+' | head -n1)

  echo "${status:-UNKNOWN}"
}

# Check if status file exists
if [ ! -f sample_list.with_status.txt ]; then
  echo "ERROR: sample_list.with_status.txt not found"
  exit 1
fi

echo "=== Analyzing samples for cleanup ==="
echo ""

# Arrays to track samples
declare -a cleanup_samples=()
declare -a skip_samples=()
declare -a done_anomalies=()

# Parse status file and categorize samples
while IFS=$'\t' read -r sample run_id srr status; do
  # Skip if no status column (3-column format) or empty
  [[ -z "$status" ]] && continue

  case "$status" in
    BAM_ERROR)
      echo "  $sample: BAM_ERROR → will clean up"
      cleanup_samples+=("$sample")
      ;;

    BAM_IN_PROGRESS)
      echo -n "  $sample: BAM_IN_PROGRESS → checking job status..."
      job_stat=$(check_job_status "$sample")
      echo " $job_stat"

      case "$job_stat" in
        RUN|PEND)
          echo "    → Job is active, skipping"
          skip_samples+=("$sample")
          ;;
        DONE)
          echo "    → Job completed but status not updated (manual review needed), skipping cleanup"
          done_anomalies+=("$sample")
          skip_samples+=("$sample")
          ;;
        EXIT|NOT_FOUND)
          echo "    → Job failed or not found, will clean up"
          cleanup_samples+=("$sample")
          ;;
        BJOBS_NOT_AVAILABLE)
          echo "    → WARNING: bjobs command not available, will clean up"
          cleanup_samples+=("$sample")
          ;;
        *)
          echo "    → Unknown status '$job_stat', will clean up to be safe"
          cleanup_samples+=("$sample")
          ;;
      esac
      ;;

    BAM_DONE)
      # Don't touch completed samples
      skip_samples+=("$sample")
      ;;

    *)
      # Unknown status - skip to be safe
      skip_samples+=("$sample")
      ;;
  esac
done < sample_list.with_status.txt

echo ""
echo "=== Summary ==="
echo "Samples to clean up: ${#cleanup_samples[@]}"
echo "Samples to skip: ${#skip_samples[@]}"
[[ ${#done_anomalies[@]} -gt 0 ]] && echo "⚠ Anomalies (DONE jobs with BAM_IN_PROGRESS): ${#done_anomalies[@]}"
echo ""

# Exit if nothing to clean
if [[ ${#cleanup_samples[@]} -eq 0 ]]; then
  echo "No samples need cleanup. Exiting."
  exit 0
fi

# Perform targeted cleanup
echo "=== Removing STAR output files for failed samples ==="
for sample in "${cleanup_samples[@]}"; do
  echo "  Cleaning $sample..."
  # Remove STAR outputs
  rm -f "${sample}_Aligned.sortedByCoord.out.bam" 2>/dev/null
  rm -f "${sample}_Log.out" "${sample}_Log.progress.out" "${sample}_Log.final.out" 2>/dev/null
  rm -f "${sample}_SJ.out.tab" 2>/dev/null
  rm -rf "${sample}__STAR"* 2>/dev/null

  # Remove prefetch artifacts
  rm -f "${sample}.sra.vdbcache" 2>/dev/null
done
echo "✓ Removed STAR output files for ${#cleanup_samples[@]} samples"

# Remove SRR subdirectories (can't easily tie to specific samples, so remove all incomplete ones)
echo ""
echo "=== Removing SRR subdirectories ==="
rm -rf SRR*/ ERR*/ DRR*/ 2>/dev/null
echo "✓ Removed SRR/ERR/DRR subdirectories"

# Update status file - reset only cleaned samples
echo ""
echo "=== Updating status file ==="

# Backup original
cp sample_list.with_status.txt sample_list.with_status.txt.backup

# Build a temp file with samples to clean
printf '%s\n' "${cleanup_samples[@]}" > /tmp/cleanup_samples.$$.txt

# Reset status for cleaned samples only (two-pass awk)
awk 'BEGIN{FS=OFS="\t"}
NR==FNR { cleanup[$1]=1; next }
{
  if ($1 in cleanup && ($4 == "BAM_ERROR" || $4 == "BAM_IN_PROGRESS")) {
    print $1, $2, $3
  } else {
    print $0
  }
}' /tmp/cleanup_samples.$$.txt sample_list.with_status.txt > sample_list.with_status.txt.tmp

rm -f /tmp/cleanup_samples.$$.txt
mv sample_list.with_status.txt.tmp sample_list.with_status.txt

echo "✓ Reset ${#cleanup_samples[@]} sample statuses"
echo "  Backup saved: sample_list.with_status.txt.backup"

# Report anomalies if any
if [[ ${#done_anomalies[@]} -gt 0 ]]; then
  echo ""
  echo "⚠ WARNING: The following samples have jobs marked DONE but status is BAM_IN_PROGRESS:"
  for sample in "${done_anomalies[@]}"; do
    echo "    $sample"
  done
  echo "  These may need manual verification to update status to BAM_DONE"
fi

echo ""
echo "=== Cleanup complete! ==="
echo "Cleaned ${#cleanup_samples[@]} samples, skipped ${#skip_samples[@]} samples"
echo ""
echo "Next steps:"
echo "  1. Re-run FASTQ workflow to regenerate corrupted files (if needed):"
echo "     fastq-workflow $TUMOR_DIR"
echo "  2. Then run STAR alignment:"
echo "     bash .../star_flex_2pass.sh sample_list.with_status.txt"
