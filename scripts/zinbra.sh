#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

WORK_DIR=`pwd`
GENOME=$1
INDEX_DIR=$2
Q=$3
BAM_FILE=$4
NAME=${BAM_FILE%%.bam} # file name without extension
ID=${NAME}_${GENOME}_${Q}

if [ ! -f "${INDEX_DIR}/${GENOME}.2bit" ]; then
    echo "Downloading sequence 2bit required for zinbra"
    cd ${INDEX_DIR}
    wget http://hgdownload.cse.ucsc.edu/goldenPath/${GENOME}/bigZips/${GENOME}.2bit
    chmod a+r *
    cd ${WORK_DIR}
fi

if [ ! -f "${ID}_peaks.bed" ]; then
    echo $(qsub << ENDINPUT
#!/bin/sh
#PBS -N zinbra_${ID}
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${NAME}_zinbra_${GENOME}.log

# This is necessary because qsub default working dir is user home
cd ${WORK_DIR}

module load java
export _JAVA_OPTIONS="-Xmx20G"
java -jar /home/oshpynov/zinbra/zinbra-0.2.4.jar analyze --input ${BAM_FILE} --reference ${INDEX_DIR}/${GENOME}.2bit --fdr ${Q} --bed ${ID}_peaks.bed

# Cleanup
mv ${ID}_peaks.narrowPeak do_not_remove_${ID}_peaks.bed
rm ${ID}*
mv do_not_remove_${ID}_peaks.bed ${ID}_peaks.bed
ENDINPUT
)
fi