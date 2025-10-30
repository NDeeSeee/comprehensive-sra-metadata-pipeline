DIR=$(pwd)

cat <<EOF
#BSUB -L /bin/bash
#BSUB -W 5:00
#BSUB -n 1
#BSUB -R "span[ptile=4]"
#BSUB -M 8000
#BSUB -e $DIR/logs/%J_sra.err.txt
#BSUB -o $DIR/logs/%J_sra.out.txt
#BSUB -J sra_download

cd $DIR

mkdir -p logs

module load sratoolkit/2.10.4
module load aspera/3.9.1 
module load python3/3.8.6

prefetch --option-file *samples_GEO.txt

# Move .sra files to the main directory
find . -type f -name '*.sra' -exec mv {} . \;


# Remove empty SRR* directories. This will work but also throw an error. I haven't figured out how to fix that yet. 
find . -type d -name 'SRR*' -exec rmdir /{} \;


EOF

# ./sratoolkit_GEO.sh | bsub
