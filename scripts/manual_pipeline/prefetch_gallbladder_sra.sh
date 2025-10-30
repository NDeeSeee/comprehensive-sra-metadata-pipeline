DIR=$(pwd)

cat <<EOF
#BSUB -L /bin/bash
#BSUB -W 50:00
#BSUB -n 1
#BSUB -R "span[ptile=4]"
#BSUB -M 132000
#BSUB -e $DIR/logs/%J_sra.err.txt
#BSUB -o $DIR/logs/%J_sra.out.txt
#BSUB -J sra_download

cd $DIR

mkdir -p logs

module load sratoolkit/2.10.4
module load aspera/3.9.1 
module load python3/3.8.6

prefetch --option-file Gallbladder-SRA_2.txt

# Move .sra files to the main directory
find . -type f -name '*.sra' -exec mv {} . \;

#Python script creates a text file of the GSM sample IDs and the corresponding Fastq file.
#This will be used for the STAR alignment step to make sure samples with multiple SRR files are processed together.
python /data/salomonis2/software/LabShellScripts/DownloadGeoData/GEO_sampleSetup.py ./Gallbladder-SraRunTable.txt

# Remove empty SRR* directories. This will work but also throw an error. I haven't figured out how to fix that yet. 
find . -type d -name 'SRR*' -exec rmdir /{} \;


EOF

# ./sratoolkit.sh | bsub
