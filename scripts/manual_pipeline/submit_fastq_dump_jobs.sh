#!/bin/bash

# Enhanced FASTQ Conversion Script (previously fdump.sh)
# Can work with individual SRA files OR sample_list.txt
# Skips existing FASTQ files automatically

if [ $# -eq 0 ]; then
    echo "Usage: $0 <sra_file_or_sample_list>"
    echo "Examples:"
    echo "  $0 SRR123456.sra                    # Convert single SRA file"
    echo "  $0 sample_list.txt                   # Convert all SRAs from sample list"
    exit 1
fi

INPUT=$1
DIR=$(pwd)

# Check if input is sample_list.txt
if [ "$INPUT" = "sample_list.txt" ] || [ "$(basename "$INPUT")" = "sample_list.txt" ]; then
    echo "Processing sample_list.txt..."
    
    if [ ! -f "$INPUT" ]; then
        echo "ERROR: sample_list.txt not found in current directory"
        exit 1
    fi
    
    # Extract SRR IDs from sample_list.txt
    SRR_IDS=$(tail -n +1 "$INPUT" | awk '{print $2}' | sed 's/_1.fastq.gz//g' | sed 's/_2.fastq.gz//g' | \
              grep -E '^[SE]RR[0-9]+' | tr ',' '\n' | sort | uniq)
    
    echo "Found $(echo "$SRR_IDS" | wc -l) SRR IDs to process"
    
    # Create logs directory
    mkdir -p logs
    
    # Process each SRR ID
    for SRR_ID in $SRR_IDS; do
        # Check if FASTQ files already exist AND pass quick validation
        # Supports both single-end (R1 only) and paired-end (R1 + R2)
        if [ -f "${SRR_ID}_1.fastq.gz" ]; then
            # Detect layout: single-end (R2 absent) or paired-end (R2 present)
            if [ -f "${SRR_ID}_2.fastq.gz" ]; then
                # Paired-end: validate both files
                if [ -s "${SRR_ID}_1.fastq.gz" ] && [ -s "${SRR_ID}_2.fastq.gz" ]; then
                    R1_RECENT=$(find "${SRR_ID}_1.fastq.gz" -mmin -2 2>/dev/null)
                    R2_RECENT=$(find "${SRR_ID}_2.fastq.gz" -mmin -2 2>/dev/null)
                    if [ -z "$R1_RECENT" ] && [ -z "$R2_RECENT" ]; then
                        echo "✓ ${SRR_ID} paired-end FASTQ files already exist (stable), skipping"
                        continue
                    else
                        echo "⚠ ${SRR_ID} paired-end FASTQ files recently modified; will re-convert"
                        rm -f "${SRR_ID}_1.fastq.gz" "${SRR_ID}_2.fastq.gz"
                    fi
                else
                    echo "⚠ ${SRR_ID} paired-end FASTQ files exist but one is empty; will re-convert"
                    rm -f "${SRR_ID}_1.fastq.gz" "${SRR_ID}_2.fastq.gz"
                fi
            else
                # Single-end: validate R1 only
                if [ -s "${SRR_ID}_1.fastq.gz" ]; then
                    R1_RECENT=$(find "${SRR_ID}_1.fastq.gz" -mmin -2 2>/dev/null)
                    if [ -z "$R1_RECENT" ]; then
                        echo "✓ ${SRR_ID} single-end FASTQ file already exists (stable), skipping"
                        continue
                    else
                        echo "⚠ ${SRR_ID} single-end FASTQ file recently modified; will re-convert"
                        rm -f "${SRR_ID}_1.fastq.gz"
                    fi
                else
                    echo "⚠ ${SRR_ID} single-end FASTQ file exists but empty; will re-convert"
                    rm -f "${SRR_ID}_1.fastq.gz"
                fi
            fi
        fi
        
        # Check if SRA file exists
        if [ ! -f "${SRR_ID}.sra" ]; then
            echo "✗ ${SRR_ID}.sra not found, skipping"
            continue
        fi
        
        echo "→ Submitting conversion job for ${SRR_ID}.sra"

        # Submit LSF job for this SRR
        bsub <<EOF
#BSUB -L /bin/bash
#BSUB -W 10:00
#BSUB -n 1
#BSUB -M 32000
#BSUB -e "$DIR/logs/${SRR_ID}_fastqdump.err.txt"
#BSUB -o "$DIR/logs/${SRR_ID}_fastqdump.out.txt"
#BSUB -J fastq_${SRR_ID}

mkdir -p logs

module load sratoolkit/2.10.4

cd "$DIR"

fastq-dump --split-files ${SRR_ID}.sra --origfmt --gzip -O .
FASTQ_DUMP_EXIT=\$?

# Only cleanup if fastq-dump succeeded AND FASTQs exist with content
# Handles both single-end (R1 only) and paired-end (R1+R2)
if [ \$FASTQ_DUMP_EXIT -eq 0 ]; then
    if [ -s "${SRR_ID}_1.fastq.gz" ]; then
        # R1 exists and non-empty. For paired-end, R2 must also exist and be non-empty.
        # For single-end, R2 won't exist.
        if [ ! -f "${SRR_ID}_2.fastq.gz" ] || [ -s "${SRR_ID}_2.fastq.gz" ]; then
            rm -f "${SRR_ID}.sra"
            echo "✓ FASTQ conversion succeeded; removed ${SRR_ID}.sra"
        else
            echo "FASTQ conversion succeeded but R2 exists but empty for ${SRR_ID}; retaining .sra" 1>&2
        fi
    else
        echo "FASTQ conversion succeeded but R1 missing or empty for ${SRR_ID}; retaining .sra" 1>&2
    fi
else
    echo "FASTQ conversion failed (exit \$FASTQ_DUMP_EXIT) for ${SRR_ID}; retaining .sra" 1>&2
fi

EOF
        if [ $? -eq 0 ]; then
            echo "  ✓ Job submitted for ${SRR_ID}"
        else
            echo "  ✗ Job submission failed for ${SRR_ID}"
        fi
    done
    
else
    # Single SRA file processing (original functionality)
    INPUTFILE="$1"
    SAMPLE=$(basename "$INPUTFILE" .sra)

    # Check if FASTQ files already exist AND pass quick validation
    # Supports both single-end (R1 only) and paired-end (R1 + R2)
    if [ -f "${SAMPLE}_1.fastq.gz" ]; then
        # Detect layout: single-end (R2 absent) or paired-end (R2 present)
        if [ -f "${SAMPLE}_2.fastq.gz" ]; then
            # Paired-end: validate both files
            if [ -s "${SAMPLE}_1.fastq.gz" ] && [ -s "${SAMPLE}_2.fastq.gz" ]; then
                R1_RECENT=$(find "${SAMPLE}_1.fastq.gz" -mmin -2 2>/dev/null)
                R2_RECENT=$(find "${SAMPLE}_2.fastq.gz" -mmin -2 2>/dev/null)
                if [ -z "$R1_RECENT" ] && [ -z "$R2_RECENT" ]; then
                    echo "✓ ${SAMPLE} paired-end FASTQ files already exist (stable), skipping"
                    exit 0
                else
                    echo "⚠ ${SAMPLE} paired-end FASTQ files recently modified; will re-convert"
                    rm -f "${SAMPLE}_1.fastq.gz" "${SAMPLE}_2.fastq.gz"
                fi
            else
                echo "⚠ ${SAMPLE} paired-end FASTQ files exist but one is empty; will re-convert"
                rm -f "${SAMPLE}_1.fastq.gz" "${SAMPLE}_2.fastq.gz"
            fi
        else
            # Single-end: validate R1 only
            if [ -s "${SAMPLE}_1.fastq.gz" ]; then
                R1_RECENT=$(find "${SAMPLE}_1.fastq.gz" -mmin -2 2>/dev/null)
                if [ -z "$R1_RECENT" ]; then
                    echo "✓ ${SAMPLE} single-end FASTQ file already exists (stable), skipping"
                    exit 0
                else
                    echo "⚠ ${SAMPLE} single-end FASTQ file recently modified; will re-convert"
                    rm -f "${SAMPLE}_1.fastq.gz"
                fi
            else
                echo "⚠ ${SAMPLE} single-end FASTQ file exists but empty; will re-convert"
                rm -f "${SAMPLE}_1.fastq.gz"
            fi
        fi
    fi
    
    bsub <<EOF
#BSUB -L /bin/bash
#BSUB -W 10:00
#BSUB -n 1
#BSUB -M 32000
#BSUB -e "$DIR/logs/${SAMPLE}_fastqdump.err.txt"
#BSUB -o "$DIR/logs/${SAMPLE}_fastqdump.out.txt"
#BSUB -J fastq_${SAMPLE}

mkdir -p logs

module load sratoolkit/2.10.4

cd "$DIR"

# Use basename only since we cd'd into the directory
# This avoids fastq-dump URL-decoding issues with special chars like '+'
fastq-dump --split-files "${SAMPLE}.sra" --origfmt --gzip -O .
FASTQ_DUMP_EXIT=\$?

# Only cleanup if fastq-dump succeeded AND FASTQs exist with content
# Handles both single-end (R1 only) and paired-end (R1+R2)
if [ \$FASTQ_DUMP_EXIT -eq 0 ]; then
    if [ -s "${SAMPLE}_1.fastq.gz" ]; then
        # R1 exists and non-empty. For paired-end, R2 must also exist and be non-empty.
        # For single-end, R2 won't exist.
        if [ ! -f "${SAMPLE}_2.fastq.gz" ] || [ -s "${SAMPLE}_2.fastq.gz" ]; then
            rm -f "${SAMPLE}.sra"
            echo "✓ FASTQ conversion succeeded; removed ${SAMPLE}.sra"
        else
            echo "FASTQ conversion succeeded but R2 exists but empty for ${SAMPLE}; retaining .sra" 1>&2
        fi
    else
        echo "FASTQ conversion succeeded but R1 missing or empty for ${SAMPLE}; retaining .sra" 1>&2
    fi
else
    echo "FASTQ conversion failed (exit \$FASTQ_DUMP_EXIT) for ${SAMPLE}; retaining .sra" 1>&2
fi

EOF
fi

# Usage examples:
# Single file: ./submit_fastq_dump_jobs.sh SRR123456.sra | bsub
# Sample list: ./submit_fastq_dump_jobs.sh sample_list.txt | bsub