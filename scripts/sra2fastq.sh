#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

if [ $# -lt 1 ]; then
    echo "Need 1 parameter! <WORK_DIR>"
    exit 1
fi
WORK_DIR=$1

# Load technical stuff, not available in qsub emulation
if [ -f "$(dirname $0)/util.sh" ]; then
    source "$(dirname $0)/util.sh"
fi

echo "Batch Sra2fastq: ${WORK_DIR}"
cd ${WORK_DIR}

TASKS=""
for FILE in $(find . -name '*.sra' | sed 's#./##g')
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
