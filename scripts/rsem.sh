#!/usr/bin/env bash
# author zayats1812@mail.ru

# Load technical stuff
source ~/work/washu/scripts/util.sh

WORK_DIR=$1
REF="$2/human_rsem_100/human_rsem_100"
RSEMPATH="/home/kzaytsev/rna_seq_pipeline/tools/RSEM-1.2.31"


echo "Batch RSEM: ${WORK_DIR} ${REF}"

cd ${WORK_DIR}

TASKS=""
for FILE in $(find . -type f -name '*.tr.bam' -printf '%P\n')
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
    echo "FILE: ${FILE}; JOB: ${QSUB_ID}"
    TASKS="$TASKS $QSUB_ID"
done
wait_complete ${TASKS}
module load R
Rscript ~/work/washu/R/gather.R

echo "Done. Batch RSEM: ${WORK_DIR} ${REF}"
