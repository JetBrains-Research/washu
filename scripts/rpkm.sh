#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

which bamCoverage &>/dev/null || { echo "Deeptools bamCoverage not found! Install: <http://deeptools.readthedocs.io/en/latest/content/installation.html>"; exit 1; }

# Load technical stuff, not available in qsub emulation
if [ -f "$(dirname $0)/util.sh" ]; then
    source "$(dirname $0)/util.sh"
fi

if [ $# -lt 1 ]; then
    echo "Need 1 parameter! <WORK_DIR>"
    exit 1
fi
WORK_DIR=$1


echo "Batch RPKM: ${WORK_DIR}"
cd ${WORK_DIR}

TASKS=""
for FILE in $(find . -name '*.bam' | sed 's#./##g' | grep -vE ".tr")
do :
    NAME=${FILE%%.bam} # file name without extension

    # Submit task
    QSUB_ID=$(qsub << ENDINPUT
#!/bin/sh
#PBS -N rpkm_${NAME}
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=64gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${NAME}_rpkm.log

# This is necessary because qsub default working dir is user home
cd ${WORK_DIR}

# Check index
module load samtools
if [[ ! -f "${FILE}.bai" ]]; then
    samtools index ${FILE}
fi

bamCoverage --bam ${FILE} --outFileName ${NAME}_rpkm.bw --outFileFormat bigwig --ignoreDuplicates --normalizeUsingRPKM
ENDINPUT
)
    echo "FILE: ${FILE}; TASK: ${QSUB_ID}"
    TASKS="$TASKS $QSUB_ID"
done
wait_complete ${TASKS}
check_logs

echo "Done. Batch RPKM: ${WORK_DIR}"
