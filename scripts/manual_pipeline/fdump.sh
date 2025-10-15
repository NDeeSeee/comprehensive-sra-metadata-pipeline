INPUTFILE=$1
SAMPLE=$(basename $INPUTFILE .sra) #removes the .txt from the file name only (parent folder name)
DIR=$(pwd)


cat <<EOF
#BSUB -L /bin/bash
#BSUB -W 10:00
#BSUB -n 1
#BSUB -R "span[ptile=4]"
#BSUB -M 32000
#BSUB -e $DIR/logs/fastqdump.err.txt
#BSUB -o $DIR/logs/fastqdump.out.txt
#BSUB -J SAMPLE

mkdir -p logs

module load sratoolkit/2.10.4
module load aspera/3.9.1 

cd $DIR

fastq-dump --split-files $INPUTFILE --origfmt --gzip  -O .

EOF
#for i in *.sra; do ./fdump.sh $i | bsub; done

