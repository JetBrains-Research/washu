#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

WORK_DIR=`pwd`
BAM_FILE=$1
GENOME=$2
NAME=${BAM_FILE%%.bam} # file name without extension

if [ ! -f ${NAME}.tdf ]; then
    echo $(qsub << ENDINPUT
#!/bin/sh
#PBS -N tdf_${ID}
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${NAME}_tdf.log

# Loading modules. TODO: install macs2 on washu cluster
# module load macs2

# This is necessary because qsub default working dir is user home
cd ${WORK_DIR}
/home/oshpynov/IGVTools/igvtools count -z 5 -w 50 -e 0 ${BAM_FILE} ${NAME}.tdf ${GENOME}
ENDINPUT
)
fi