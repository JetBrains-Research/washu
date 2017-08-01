#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

which bamCoverage &>/dev/null || {
    echo "deeptools not found! You can install it using:"
    echo "  conda install -c bioconda deeptools"
    echo "For further details see http://deeptools.readthedocs.io/en/latest/content/installation.html"
    exit 1;
}

# Load technical stuff, not available in qsub emulation
if [ -f "$(dirname $0)/util.sh" ]; then
    source "$(dirname $0)/util.sh"
fi

if [ $# -lt 1 ]; then
    echo "Need at least 1 parameter! <WORK_DIR> [<WORK_DIR>]*"
    exit 1
fi
WORK_DIRS=$@


echo "Batch RPKM: ${WORK_DIRS}"

TASKS=""
for WORK_DIR in ${WORK_DIRS}; do
    WORK_DIR_NAME=${WORK_DIR##*/}
    cd ${WORK_DIR}

    for FILE in $(find . -name '*.bam' | sed 's#\./##g' | grep -vE ".tr")
    do :
        NAME=${FILE%%.bam} # file name without extension

        # Submit task
        QSUB_ID=$(qsub << ENDINPUT
#!/bin/sh
#PBS -N rpkm_${WORK_DIR_NAME}_${NAME}
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
        echo "FILE: ${WORK_DIR_NAME}/${FILE}; TASK: ${QSUB_ID}"
        TASKS="$TASKS $QSUB_ID"
    done
done
wait_complete ${TASKS}
check_logs

echo "Done. Batch RPKM: ${WORK_DIRS}"
