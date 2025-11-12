#!/bin/bash

file1=$1  # First read (R1)
sample="${file1%_1.fastq.gz}"  # Sample name without "_1.fastq.gz"
file2="${sample}_2.fastq.gz"  # Corresponding second read (R2)

op_name=$(basename $file1 _1.fastq.gz)  # Extract sample name for output

sif_path="/data/salomonis-archive/BAMs/NCI-R01/TCGA/TCGA-OV/optitype_container.sif"  # Path to Singularity container

DIR=$(pwd)  # Current working directory

# Check if the result file already exists in any timestamped subdirectory
if find processed/${op_name} -name "*_result.tsv" | grep -q .; then
    echo "Results for ${op_name} already exist. Exiting."
    exit 0
fi

cat << EOF
#BSUB -W 8:00
#BSUB -M 32000
#BSUB -n 4
#BSUB -J OptiType_${op_name}_%J
#BSUB -R "span[hosts=1]"
#BSUB -o $DIR/logs/OptiType_${op_name}_%J.out
#BSUB -e $DIR/logs/OptiType_${op_name}_%J.err

## Create logs directory if it doesn't exist
if [ ! -d logs ]; then
    mkdir logs
fi

## Create processed directory for the sample if it doesn't exist
if [ ! -d processed/${op_name} ]; then
    mkdir -p processed/${op_name}
fi

module load singularity/3.7.0

# Run OptiType with both reads provided together
singularity exec -W /mnt -B $(pwd):/mnt ${sif_path} /bin/bash -c \
"cd /mnt && /usr/local/bin/OptiType/OptiTypePipeline.py -i ${file1} ${file2} --rna -v -o processed/${op_name}"

EOF

##for f in fastqs/*_1.fastq.gz; do ./run_optiplex_code3.sh $f | bsub; done