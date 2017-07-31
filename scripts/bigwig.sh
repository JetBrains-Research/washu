#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

# Load technical stuff, not available in qsub emulation
if [ -f "$(dirname $0)/util.sh" ]; then
    source "$(dirname $0)/util.sh"
fi

if [ $# -lt 2 ]; then
    echo "Need 2 parameters! <WORK_DIR> <CHROM_SIZES>"
    exit 1
fi
WORK_DIR=$1
CHROM_SIZES=$2

echo "Batch BigWig: ${WORK_DIR} ${CHROM_SIZES}"
cd ${WORK_DIR}

TASKS=""
for FILE in $(find . -name '*.bam' | sed 's#\./##g' | grep -vE ".tr")
do :
    NAME=${FILE%%.bam} # file name without extension

    # Submit task
    QSUB_ID=$(qsub << ENDINPUT
#!/bin/sh
#PBS -N bw_${NAME}
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${NAME}_bw.log

# This is necessary because qsub default working dir is user home
cd ${WORK_DIR}

module load bedtools2
bash $(dirname $0)/../bam2bw.sh ${FILE} ${CHROM_SIZES}
ENDINPUT
)
    echo "FILE: ${FILE}; TASK: ${QSUB_ID}"
    TASKS="$TASKS $QSUB_ID"
done
wait_complete ${TASKS}
check_logs

echo "Done. Batch BigWig: ${WORK_DIR} ${CHROM_SIZES}"
