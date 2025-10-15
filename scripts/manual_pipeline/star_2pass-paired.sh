#!/bin/bash

# Get the input FASTQ file from the command-line argument
FASTQ1=$1
FASTQ2=${FASTQ1/_1./_2.}
SAMPLE=$(basename $FASTQ1 _1.fastq.gz)
DIR=$(pwd)

OUT_DIR="${DIR}/star_output"
LOG_DIR="${DIR}/logs"
BAM_DIR="${DIR}/bams"
GENOME_DIR=/data/salomonis2/Genomes/Star2pass-GRCH38/GenomeRef
GENOME=/data/salomonis2/Genomes/Star2pass-GRCH38/GRCh38.d1.vd1.fa

mkdir -p $OUT_DIR $LOG_DIR $BAM_DIR

cat <<EOF
#BSUB -L /bin/bash
#BSUB -W 160:00
#BSUB -n 4
#BSUB -M 128000
#BSUB -e $LOG_DIR/${SAMPLE}_star_redo_%J.err
#BSUB -o $LOG_DIR/${SAMPLE}_star_redo_%J.out
#BSUB -J $SAMPLE

cd $DIR
module load STAR/2.4.0h

# 1st Pass: Generate the splice junctions database
STAR --genomeDir ${GENOME_DIR} \
     --readFilesIn $FASTQ1 $FASTQ2 \
     --runThreadN 8 \
     --outFilterMultimapScoreRange 1 \
     --outFilterMultimapNmax 20 \
     --outFilterMismatchNmax 10 \
     --alignIntronMax 500000 \
     --alignMatesGapMax 1000000 \
     --sjdbScore 2 \
     --alignSJDBoverhangMin 1 \
     --genomeLoad NoSharedMemory \
     --limitBAMsortRAM 100000000000 \
     --readFilesCommand gunzip -c \
     --outFileNamePrefix ${OUT_DIR}/${SAMPLE}_pass1_ \
     --outFilterMatchNminOverLread 0.33 \
     --outFilterScoreMinOverLread 0.33 \
     --sjdbOverhang 100 \
     --outSAMstrandField intronMotif \
     --outSAMtype None \
     --outSAMmode None

# 2nd Pass: Use the splice junctions from the first pass
mkdir -p ${OUT_DIR}/GenomeRef_${SAMPLE}
STAR --runMode genomeGenerate \
     --genomeDir ${OUT_DIR}/GenomeRef_${SAMPLE} \
     --genomeFastaFiles $GENOME \
     --sjdbOverhang 100 \
     --runThreadN 8 \
     --sjdbFileChrStartEnd ${OUT_DIR}/${SAMPLE}_pass1_SJ.out.tab \
     --outFileNamePrefix ${OUT_DIR}/${SAMPLE}_pass2_

STAR --genomeDir ${OUT_DIR}/GenomeRef_${SAMPLE} \
     --readFilesIn $FASTQ1 $FASTQ2 \
     --runThreadN 8 \
     --outFilterMultimapScoreRange 1 \
     --outFilterMultimapNmax 20 \
     --outFilterMismatchNmax 10 \
     --alignIntronMax 500000 \
     --alignMatesGapMax 1000000 \
     --sjdbScore 2 \
     --alignSJDBoverhangMin 1 \
     --genomeLoad NoSharedMemory \
     --limitBAMsortRAM 100000000000 \
     --readFilesCommand gunzip -c \
     --outFileNamePrefix ${OUT_DIR}/${SAMPLE}_second_ \
     --outFilterMatchNminOverLread 0.33 \
     --outFilterScoreMinOverLread 0.33 \
     --sjdbOverhang 100 \
     --outSAMstrandField intronMotif \
     --outSAMattributes NH HI NM MD AS XS \
     --outSAMunmapped Within \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMheaderHD @HD VN:1.4 

# Cleanup and move final BAM
mv ${OUT_DIR}/${SAMPLE}_second_Aligned.sortedByCoord.out.bam ${BAM_DIR}/${SAMPLE}.bam
rm -r ${OUT_DIR}/GenomeRef_${SAMPLE}
rm ${OUT_DIR}/${SAMPLE}_pass1_SJ.out.tab
EOF


#for i in *1.fastq.gz; do bash star_2pass-paired.sh $i | bsub; done
