#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

# Load technical stuff, not available in qsub emulation
if [ -f "$(dirname $0)/util.sh" ]; then
    source "$(dirname $0)/util.sh"
fi

if [ $# -lt 1 ]; then
    echo "Need 1 parameter! <WORK_DIR>"
    exit 1
fi
WORK_DIR=$1

echo "Batch Fastqc: ${WORK_DIR}"
cd ${WORK_DIR}

TASKS=""
for FILE in $(find . -regextype posix-extended -regex '.*\.f.*q(\.gz)?' | sed 's#./##g')
do :
    NAME=${FILE%%.f*q} # file name without extension

    # Submit task
    QSUB_ID=$(qsub << ENDINPUT
#!/bin/sh
#PBS -N fastqc_${NAME}
#PBS -l nodes=1:ppn=1,walltime=2:00:00,vmem=4gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${NAME}_fastqc.log

# Loading modules
module load fastqc

# This is necessary because qsub default working dir is user home
cd ${WORK_DIR}
fastqc ${FILE}
ENDINPUT
)
    echo "FILE: ${FILE}; TASK: ${QSUB_ID}"
    TASKS="$TASKS $QSUB_ID"
done

wait_complete ${TASKS}
check_logs

echo "Processing multiqc"
mkdir -p ${WORK_DIR}/fastqc
mv *_fastqc.* ${WORK_DIR}/fastqc
multiqc ${WORK_DIR}/fastqc
echo "Done. Batch Fastqc: $WORK_DIR"
