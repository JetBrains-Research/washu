#!/usr/bin/env bash
# Script to create _pileup.bed files for BAM alignment to compute reads coverage.
# author oleg.shpynov@jetbrains.com

# Load technical stuff, not available in qsub emulation
if [ -f "$(dirname $0)/util.sh" ]; then
    source "$(dirname $0)/util.sh"
fi

if [ $# -lt 1 ]; then
    echo "Need 1 parameters! <WORK_DIR>"
    exit 1
fi

WORK_DIR=$1

echo "Batch pileup: ${WORK_DIR}"
cd ${WORK_DIR}

PROCESSED=""
TASKS=""

for FILE in $(find . -name '*.bam' | sed 's#./##g' | sort)
do :
    NAME=${FILE%%.bam}
    FILE_BED=${NAME}.bed

    # Submit task
    QSUB_ID=$(qsub << ENDINPUT
#!/bin/sh
#PBS -N pileup_${NAME}
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${NAME}_pileup.log

cd ${WORK_DIR}
module load bedtools2

# To pileup sorted bed
export LC_ALL=C
bedtools bamtobed -i ${FILE} | sort -k1,1 -k3,3n -k2,2n -k6,6 > ${FILE_BED}
ENDINPUT
)
    TASKS="$TASKS $QSUB_ID"
done
wait_complete ${TASKS}
check_logs

echo "Done. Batch pileup: ${WORK_DIR}"