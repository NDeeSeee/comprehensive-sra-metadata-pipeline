#!/bin/bash
# Reset samples that had STAR failures so they can be reprocessed

if [ $# -ne 1 ]; then
  echo "Usage: $0 <tumor_directory>"
  echo "Example: $0 /data/salomonis-archive/FASTQs/NCI-R01/POSEIDON/Tumors/LungLargeCell"
  exit 1
fi

TUMOR_DIR="$1"
cd "$TUMOR_DIR" || exit 1

echo "Working directory: $PWD"
echo ""

# 1. Remove STAR debris
echo "=== Removing STAR output files ==="
rm -f SAMN*_Aligned.sortedByCoord.out.bam
rm -f SAMN*_Log.out SAMN*_Log.progress.out SAMN*_Log.final.out
rm -f SAMN*_SJ.out.tab
rm -rf SAMN*__STAR*
echo "✓ Removed STAR output files"

# 2. Remove prefetch cache files
echo ""
echo "=== Removing .vdbcache files ==="
rm -f *.sra.vdbcache
echo "✓ Removed .vdbcache files"

# 3. Remove SRR subdirectories (from incomplete prefetch)
echo ""
echo "=== Removing SRR subdirectories ==="
rm -rf SRR*/
rm -rf ERR*/
rm -rf DRR*/
echo "✓ Removed SRR/ERR/DRR subdirectories"

# 4. Update status file: BAM_IN_PROGRESS/BAM_ERROR -> (empty status)
if [ -f sample_list.with_status.txt ]; then
  echo ""
  echo "=== Resetting status for failed samples ==="

  # Backup original
  cp sample_list.with_status.txt sample_list.with_status.txt.backup

  # Count samples to reset
  bam_error=$(grep -c "BAM_ERROR" sample_list.with_status.txt 2>/dev/null || echo 0)
  bam_in_progress=$(grep -c "BAM_IN_PROGRESS" sample_list.with_status.txt 2>/dev/null || echo 0)

  # Reset to 3-column format (no status) so workflow will recheck
  awk 'BEGIN{FS=OFS="\t"} {
    if ($4 == "BAM_ERROR" || $4 == "BAM_IN_PROGRESS") {
      # Check if corresponding FASTQs exist
      print $1, $2, $3
    } else {
      print $0
    }
  }' sample_list.with_status.txt > sample_list.with_status.txt.tmp

  mv sample_list.with_status.txt.tmp sample_list.with_status.txt

  echo "✓ Reset $((bam_error + bam_in_progress)) sample statuses"
  echo "  - BAM_ERROR: $bam_error"
  echo "  - BAM_IN_PROGRESS: $bam_in_progress"
  echo "  - Backup saved: sample_list.with_status.txt.backup"
fi

echo ""
echo "=== Cleanup complete! ==="
echo "Next steps:"
echo "  1. Re-run FASTQ workflow to regenerate corrupted files:"
echo "     fastq-workflow $TUMOR_DIR"
echo "  2. Then run STAR alignment:"
echo "     bash .../star_flex_2pass.sh sample_list.with_status.txt"
