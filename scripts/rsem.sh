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

REF="$2/human_rsem_100/human_rsem_100"
RSEMPATH="/home/kzaytsev/rna_seq_pipeline/tools/RSEM-1.2.31"


echo "Batch RSEM: ${WORK_DIR} ${REF}"

cd ${WORK_DIR}

TASKS=""
for FILE in $(find . -name '*.tr.bam' | sed 's#\./##g')
do :

    NAME=${FILE%%.tr.bam} # file name without extension

    # Submit task
    QSUB_ID=$(qsub << ENDINPUT
#!/bin/bash
#PBS -N rsem_exp_$NAME
#PBS -j oe
#PBS -l mem=32gb,nodes=1:ppn=8:haswell,walltime=24:00:00

## this is where pre-generated STAR reference genome

cd ${WORK_DIR}
$RSEMPATH/rsem-calculate-expression -p 8 --paired-end --bam --estimate-rspd --no-bam-output $FILE $REF $NAME
ENDINPUT
)
    echo "FILE: ${FILE}; TASK: ${QSUB_ID}"
    TASKS="$TASKS $QSUB_ID"
done
wait_complete ${TASKS}
module load R
Rscript $(dirname $0)/../R/gather.R

echo "Done. Batch RSEM: ${WORK_DIR} ${REF}"
