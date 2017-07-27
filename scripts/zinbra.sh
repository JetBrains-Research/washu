#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

# Load technical stuff, not available in qsub emulation
if [ -f "$(dirname $0)/util.sh" ]; then
    source "$(dirname $0)/util.sh"
fi

if [ $# -lt 5 ]; then
    echo "Need 5 parameters! <ZINBRA_JAR_PATH> <WORK_DIR> <GENOME> <CHROM_SIZES> <Q>"
    exit 1
fi

ZINBRA_JAR_PATH=$1
WORK_DIR=$2
GENOME=$3
CHROM_SIZES=$4
Q=$5

echo "Batch zinbra: ${ZINBRA_JAR_PATH} ${WORK_DIR} ${GENOME} ${CHROM_SIZES} ${Q}"
cd ${WORK_DIR}

TASKS=""
for FILE in $(find . -name '*.bam' | sed 's#./##g' | grep -v 'input')
do :
    INPUT=$(python $(dirname $0)/util.py find_input ${WORK_DIR}/${FILE})
    echo "${FILE}: control file: ${INPUT}"

    NAME=${FILE%%.bam} # file name without extension
    ID=${NAME}_${GENOME}_${Q}

    # Submit task
    QSUB_ID=$(qsub << ENDINPUT
#!/bin/sh
#PBS -N zinbra_${ID}
#PBS -l nodes=1:ppn=4,walltime=24:00:00,vmem=64gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${NAME}_zinbra_${GENOME}.log

# This is necessary because qsub default working dir is user home
cd ${WORK_DIR}

module load java
# PROBLEM: vmem is much bigger, however we face with the problem with bigger values:
# There is insufficient memory for the Java Runtime Environment to continue.
export _JAVA_OPTIONS="-Xmx30g"

if [ -f "${INPUT}" ]; then
    echo "${FILE}: control file found: ${INPUT}"
    java -cp ${ZINBRA_JAR_PATH} org.jetbrains.bio.zinbra.ZinbraCLA analyze -t ${FILE} -c ${FILE} \
        --chrom.sizes ${CHROM_SIZES} --fdr ${Q} --bed ${ID}_peaks.bed

else
    echo "${FILE}: no control file"
    java -cp ${ZINBRA_JAR_PATH} org.jetbrains.bio.zinbra.ZinbraCLA analyze -t ${FILE} \
        --chrom.sizes ${CHROM_SIZES} --fdr ${Q} --bed ${ID}_peaks.bed
fi

# Compute Reads in Peaks
bash $(dirname $0)/../reports/rip.sh ${FILE} ${ID}_peaks.bed
ENDINPUT
)
    echo "FILE: ${FILE}; TASK: ${QSUB_ID}"
    TASKS="$TASKS $QSUB_ID"
done
wait_complete ${TASKS}
check_logs

echo "Done. Batch zinbra: ${ZINBRA_JAR_PATH} ${WORK_DIR} ${GENOME} ${CHROM_SIZES} ${Q}"