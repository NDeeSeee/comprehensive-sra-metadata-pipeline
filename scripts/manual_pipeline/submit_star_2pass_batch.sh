#!/bin/bash
# ================================================================
# STAR 2-pass alignment LSF submission script
# Submits one job per sample from a sample list file
# Each job runs run_star-new.sh with appropriate FASTQ files.
# ================================================================

# Create output directories
mkdir -p bams logs

# Paths
ROOT=$PWD
SAMPLE_LIST=${ROOT}/sample_list.txt

# ----------------------------------------------------------------
# Recommended LSF settings:
#   -W 12:00          walltime (12 hours)
#   -n 8              use 8 threads per STAR run
#   -R "rusage[mem=16000]"  request 16 GB per core (128 GB total)
#   -M 128000         total memory limit (matches 128 GB)
#   span[hosts=1]     keep all threads on one host
# ----------------------------------------------------------------
# Adjust -W or memory values if your cluster has tighter limits.
# ----------------------------------------------------------------

while read SAMPLE FASTQ1 FASTQ2; do
    bsub -W 12:00 -n 2 -M 128000 \
         -R "rusage[mem=16000] span[hosts=1]" \
         -J "align_${SAMPLE}" \
         -o logs/STAR2pass_${SAMPLE}.out \
         -e logs/STAR2pass_${SAMPLE}.err \
         "$ROOT/run_star-new.sh $SAMPLE $FASTQ1 $FASTQ2"
done < ${SAMPLE_LIST}

# End of script
