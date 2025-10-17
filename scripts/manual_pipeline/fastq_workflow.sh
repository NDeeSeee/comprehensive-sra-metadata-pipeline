#!/bin/bash

# POSEIDON FASTQ Workflow Script
# Downloads SRA files and converts to FASTQ using sample_list.txt
# Skips existing files automatically

if [ $# -ne 1 ]; then
    echo "Usage: $0 <cancer_directory>"
    echo "Example: $0 /data/salomonis-archive/FASTQs/NCI-R01/POSEIDON/Tumors/Tongue"
    exit 1
fi

CANCER_DIR=$1
SAMPLE_LIST="$CANCER_DIR/sample_list.txt"

if [ ! -f "$SAMPLE_LIST" ]; then
    echo "ERROR: sample_list.txt not found in $CANCER_DIR"
    exit 1
fi

echo "=========================================="
echo "POSEIDON FASTQ WORKFLOW"
echo "=========================================="
echo "Cancer directory: $CANCER_DIR"
echo "Sample list: $SAMPLE_LIST"
echo "=========================================="
echo ""

cd "$CANCER_DIR"

# Extract SRR IDs from sample_list.txt
SRR_IDS=$(tail -n +1 "$SAMPLE_LIST" | awk '{print $2}' | sed 's/_1.fastq.gz//g' | sed 's/_2.fastq.gz//g' | \
          grep -E '^[SE]RR[0-9]+' | tr ',' '\n' | sort | uniq)

NUM_SAMPLES=$(echo "$SRR_IDS" | wc -l)
echo "Found $NUM_SAMPLES SRR IDs to process"

# Create logs directory
mkdir -p logs

# Load modules
module load sratoolkit/2.10.4
module load aspera/3.9.1

echo ""
echo "=========================================="
echo "STEP 1: Download SRA files (if needed)"
echo "=========================================="

for SRR_ID in $SRR_IDS; do
    if [ -f "${SRR_ID}.sra" ]; then
        echo "  ✓ ${SRR_ID}.sra already exists"
    else
        echo "  → Downloading ${SRR_ID}..."
        prefetch ${SRR_ID} -X 35000000
        
        # Move .sra file to current directory
        if [ -d "${SRR_ID}" ]; then
            find "${SRR_ID}" -type f -name '*.sra' -exec mv {} . \;
            rmdir "${SRR_ID}" 2>/dev/null || true
        fi
        
        if [ -f "${SRR_ID}.sra" ]; then
            echo "  ✓ ${SRR_ID}.sra downloaded successfully"
        else
            echo "  ✗ ${SRR_ID}.sra download failed"
        fi
    fi
done

echo ""
echo "=========================================="
echo "STEP 2: Convert SRA to FASTQ"
echo "=========================================="

# Use the enhanced fdump.sh script
echo "→ Submitting FASTQ conversion jobs..."
bash /data/salomonis-archive/FASTQs/NCI-R01/POSEIDON/ValeriiGitRepo/scripts/manual_pipeline/fdump.sh sample_list.txt | bsub

echo ""
echo "=========================================="
echo "WORKFLOW COMPLETE"
echo "=========================================="
echo "Check job status with: bjobs"
echo "Monitor logs in: $CANCER_DIR/logs/"
echo ""
echo "Expected FASTQ files: $((NUM_SAMPLES * 2))"
echo "Current FASTQ files: $(ls -1 *.fastq.gz 2>/dev/null | wc -l)"