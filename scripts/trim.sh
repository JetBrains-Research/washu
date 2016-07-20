#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

# Load technical stuff
source ~/work/washu/scripts/util.sh

WORK_DIR=$1
TRIM=$2

echo "Batch Trim ${TRIM}: ${WORK_DIR}"
cd ${WORK_DIR}

TASKS=""
for FILE in $(find . -type f -name '*.f*q' -printf '%P\n')
do :
    NAME=${FILE%%.f*q} # file name without extension

    # Submit task
    QSUB_ID=$(qsub << ENDINPUT
#!/bin/sh
#PBS -N trim_${NAME}
#PBS -l nodes=1:ppn=8,walltime=4:00:00,vmem=4gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${NAME}_trim_${TRIM}.log

# This is necessary because qsub default working dir is user home
cd ${WORK_DIR}

seqtk trimfq -b 5 ${FILE} > ${NAME}_${TRIM}.fq
ENDINPUT
)
    echo "FILE: ${FILE}; JOB: ${QSUB_ID}"
    TASKS="$TASKS $QSUB_ID"
done

wait_complete ${TASKS}
check_logs
echo "Done. Batch Trim ${TRIM}: ${WORK_DIR}"