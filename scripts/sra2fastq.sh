#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

WORK_DIR=$1

# Load technical stuff
source ~/work/washu/scripts/util.sh

echo "Batch Sra2fastq: ${WORK_DIR}"
cd ${WORK_DIR}

TASKS=""
for FILE in $(find . -type f -name '*.sra' -printf '%P\n')
do :
    NAME=${FILE%%.sra} # file name without extension

    # Submit task
    QSUB_ID=$(qsub << ENDINPUT
#!/bin/sh
#PBS -N sra2fastq_${NAME}
#PBS -l nodes=1:ppn=8,walltime=2:00:00,vmem=8gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${NAME}_sra2fastq.log

# Loading sratoolkit module
module load sratoolkit

# This is necessary because qsub default working dir is user home
cd ${WORK_DIR}
fastq-dump --split-3 --outdir ${WORK_DIR} ${FILE}
ENDINPUT
)
    echo "FILE: ${FILE}; JOB: ${QSUB_ID}"
    TASKS="$TASKS $QSUB_ID"
done

wait_complete ${TASKS}
check_logs
echo "Done. Batch Sra2fastq: ${WORK_DIR}"
