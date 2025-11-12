#!/bin/bash

BAM=$1
SAMPLE=$(basename $BAM .txt)
DIR=$(pwd)

cat <<EOF
#BSUB -L /bin/bash
#BSUB -W 1:00
#BSUB -n 1
#BSUB -R "span[ptile=4]"
#BSUB -M 16000
#BSUB -e $DIR/logs/%J.err
#BSUB -o $DIR/logs/%J.out
#BSUB -J $SAMPLE

cd $DIR

python3 /data/salomonis2/software/AltAnalyze/import_scripts/hla.py --i $DIR/$BAM --o $DIR/fastqs

EOF
#for i in *.bam; do ./hla.sh $i | bsub; done
