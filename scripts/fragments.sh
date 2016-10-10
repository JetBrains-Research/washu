#!/usr/bin/env bash
# author zayats1812@mail.ru

# Load technical stuff
source ~/work/washu/scripts/util.sh

WORK_DIR=$1
echo "Batch fragments: ${WORK_DIR}"
cd ${WORK_DIR}

TASKS=""
for FILE in $(find . -type f -name '*.bam' -printf '%P\n')
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
samtools view -f66 $NAME.bam | cut -f 9 | sed 's/^-//' > InsertSizeMetrics.txt
Rscript ~/work/washu/R/fragments.R
ENDINPUT
)
    echo "FILE: ${FILE}; JOB: ${QSUB_ID}"
    TASKS="$TASKS $QSUB_ID"
done

wait_complete ${TASKS}
check_logs

echo "Done. Batch fragments: $WORK_DIR"