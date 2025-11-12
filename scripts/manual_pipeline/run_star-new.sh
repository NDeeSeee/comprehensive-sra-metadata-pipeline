#!/bin/bash

# Script to run STAR 2-pass alignment for a single sample
# Usage: ./run_star.sh <SAMPLE> <FASTQ1> <FASTQ2>

SAMPLE=$1
FASTQ1=$2
FASTQ2=$3

#!/bin/bash

# Script to run STAR 2-pass alignment for a single sample using STAR ≥2.7
# Usage: ./run_star_modern.sh <SAMPLE> <FASTQ1> <FASTQ2>

SAMPLE=$1
FASTQ1=$2
FASTQ2=$3

# Load updated STAR version (adjust this line based on your module system)
module load STAR/2.7.10b
module load samtools
module load bedtools

# Paths
GENOME_DIR=/data/salomonis2/Genomes/STAR-2.7.10b-Index-GRCH38/Grch38-STAR-index
GENOME=/data/salomonis2/Genomes/STAR-2.7.10b-Index-GRCH38/Grch38_r85.all.fa
GTF=/data/salomonis2/Genomes/STAR-2.7.10b-Index-GRCH38/Homo_sapiens.GRCh38.85.gtf
ROOT_DIR=$PWD

# STAR 2-pass alignment
echo "Running STAR 2-pass alignment for sample: $SAMPLE"

STAR --runThreadN 8 \
     --genomeDir ${GENOME_DIR} \
     --readFilesIn $FASTQ1 $FASTQ2 \
     --readFilesCommand gunzip -c \
     --outFileNamePrefix ${ROOT_DIR}/${SAMPLE}_ \
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
     --sjdbGTFfile ${GTF} \
     --sjdbOverhang 100

# Move BAM to final location
mkdir -p ${ROOT_DIR}/bams
mv ${ROOT_DIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam ${ROOT_DIR}/bams/${SAMPLE}.bam

# Cleanup intermediate STAR files (optional)
echo "$SAMPLE cleanup"
rm -f ${ROOT_DIR}/${SAMPLE}_SJ.out.tab
rm -f ${ROOT_DIR}/${SAMPLE}_Log.final.out
rm -f ${ROOT_DIR}/${SAMPLE}_Log.out
rm -f ${ROOT_DIR}/${SAMPLE}_Log.progress.out

echo "$SAMPLE STAR 2-pass (modern) complete"

#bsub < STAR_2pass_submit.sh
