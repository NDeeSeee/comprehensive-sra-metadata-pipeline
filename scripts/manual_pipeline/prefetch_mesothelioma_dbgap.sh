INPUTFILE=$1
SAMPLE=$(basename $INPUTFILE .sra) #removes the .txt from the file name only (parent folder name)
DIR=$(pwd)


cat <<EOF
#BSUB -L /bin/bash
#BSUB -W 50:00
#BSUB -n 4
#BSUB -R "span[ptile=4]"
#BSUB -M 10000
#BSUB -e $DIR/logs/%J_sra.err.txt
#BSUB -o $DIR/logs/%J_sra.out.txt
#BSUB -J sra_download

cd $DIR

mkdir -p logs

module load sratoolkit/3.1.1
module load aspera/3.9.1 

#Bioproject PRJNA671233. Mesothelioma Samples
#generating a cart file will all 100 samples was impossible for some reason. Thus, I've resorted to this stupid, clunky way of downloading these samples.
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18543964
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18543969
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18543972
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18543981
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18543984
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18543987
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18543990
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18543993
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18543996
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18543999
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544005
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544008
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544013
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544018
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544021
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544024
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544027
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544032
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544035
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544039
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544042
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544045
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544048
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544051
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544054
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544057
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544060
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544063
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544066
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544070
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544073
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544076
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544079
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544082
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544085
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544088
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544091
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544095
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544098
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544101
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544104
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544107
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544111
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544116
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544123
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544126
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544129
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544132
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544135
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544138
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544141
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544144
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544149
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544152
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544155
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544158
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544161
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544164
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544167
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544171
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544176
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544179
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544183
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544186
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544189
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544194
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544197
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544200
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544203
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544206
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544209
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544212
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544216
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544219
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544222
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544225
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544231
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544235
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544238
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544241
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544245
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544248
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544251
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544254
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544257
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544262
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544265
prefetch --ngc prj_40747_D41857.ngc --max-size 35000000 SRR18544268



EOF

# ./sratoolkit.sh | bsub
