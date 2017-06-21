#!/usr/bin/env bash
# author zayats1812@mail.ru

# Load technical stuff, not available in qsub emulation
if [ -f "$(dirname $0)/util.sh" ]; then
    source "$(dirname $0)/util.sh"
fi

if [ $# -lt 1 ]; then
    echo "Need 1 parameter! <WORK_DIR>"
    exit 1
fi
WORK_DIR=$1

echo "Batch fragments: ${WORK_DIR}"
cd ${WORK_DIR}

TASKS=""
for FILE in $(find . -type f -name '*.bam' | sed 's#./##g')
do :
    NAME=${FILE%%.bam} # file name without extension

    # Submit task
    QSUB_ID=$(qsub << ENDINPUT
#!/bin/sh
#PBS -N fragments_${NAME}
#PBS -l nodes=1:ppn=1,walltime=2:00:00,vmem=8gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${NAME}_fragments.log

# Loading modules
module load samtools
module load R

cd ${WORK_DIR}
samtools view -f66 $FILE | cut -f 9 | sed 's/^-//' > ${NAME}_metrics.txt
Rscript $(dirname $0)/../R/fragments.R ${NAME}_metrics.txt ${NAME}_fragments.png
ENDINPUT
)
    echo "FILE: ${FILE}; TASK: ${QSUB_ID}"
    TASKS="$TASKS $QSUB_ID"
done

wait_complete ${TASKS}
check_logs

echo "Done. Batch fragments: $WORK_DIR"