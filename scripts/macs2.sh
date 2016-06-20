#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

WORK_DIR=`pwd`
GENOME=$1
Q=$2
BAM_FILE=$3
NAME=${BAM_FILE%%.bam} # file name without extension
ID=${NAME}_${GENOME}_${Q}

# Convert Genome build to macs2 species
[[ ${GENOME} =~ ^hg[0-9]+$ ]] && SPECIES="hs"
[[ ${GENOME} =~ ^mm[0-9]+$ ]] && SPECIES="mm"
[[ -z "$SPECIES" ]] && echo "Unknown species for macs: $GENOME" && exit 1

if [ ! -f "${ID}_peaks.bed" ]; then
    echo $(qsub << ENDINPUT
#!/bin/sh
#PBS -N macs2_${ID}
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${WORK_DIR}/qsub/${NAME}_macs2_${GENOME}.log

# Loading modules. TODO: install macs2 on washu cluster
# module load macs2

# This is necessary because qsub default working dir is user home
cd ${WORK_DIR}
/home/oshpynov/miniconda2/bin/macs2 callpeak -t ${BAM_FILE} -f BAM -g ${SPECIES} -n ${ID} -B -q ${Q}

# Cleanup
mv ${ID}_peaks.narrowPeak do_not_remove_${ID}_peaks.bed
rm ${ID}*
mv do_not_remove_${ID}_peaks.bed ${ID}_peaks.bed
ENDINPUT
)
fi