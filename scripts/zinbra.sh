#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

# Load technical stuff, not available in qsub emulation
if [ -f "$(dirname $0)/util.sh" ]; then
    source "$(dirname $0)/util.sh"
fi

if [ $# -lt 4 ]; then
    echo "Need 4 parameters! <WORK_DIR> <GENOME> <FOLDER> <Q>"
    exit 1
fi
WORK_DIR=$1
GENOME=$2
FOLDER=$3
Q=$4

echo "Batch zinbra: ${WORK_DIR} ${GENOME} ${FOLDER} ${Q}"
cd ${WORK_DIR}

TASKS=""
for FILE in $(find . -name '*.bam' | sed 's#./##g')
do :
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
java -jar /home/oshpynov/zinbra/zinbra-0.2.4.jar analyze --input ${FILE} --reference ${FOLDER}/${GENOME}.2bit --fdr ${Q} --bed ${ID}_peaks.bed
ENDINPUT
)
    echo "FILE: ${FILE}; TASK: ${QSUB_ID}"
    TASKS="$TASKS $QSUB_ID"
done
wait_complete ${TASKS}
check_logs

echo "Done. Batch zinbra: ${WORK_DIR} ${GENOME} ${FOLDER} ${Q}"