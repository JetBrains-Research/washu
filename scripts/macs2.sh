#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

WORK_DIR=`pwd`
GENOME=$1
Q=$2
BAM_FILE=$3
NAME=${BAM_FILE%%.bam} # file name without extension
ID=${GENOME}_${Q}_$NAME

# Convert Genome build to macs2 species
[[ $GENOME =~ ^hg[0-9]+$ ]] && SPECIES="hs"
[[ $GENOME =~ ^mm[0-9]+$ ]] && SPECIES="mm"
[[ -z "$SPECIES" ]] && echo "Unknown species for macs: $GENOME" && exit 1

if [ ! -f "${ID}_peaks.bed" ]; then
    echo $(qsub -d $WORK_DIR << ENDINPUT
#!/bin/sh
#PBS -N macs2_$ID
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o $WORK_DIR/qsub/macs2_${GENOME}_${Q}_$NAME.log

# Loading modules. TODO: install macs2 on washu cluster
# module load macs2

/home/oshpynov/miniconda2/bin/macs2 callpeak -t $BAM_FILE -f BAM -g $SPECIES -n $ID -B -q $Q
mv ${ID}_peaks.narrowPeak ${ID}_peaks.bed
ENDINPUT
)
fi