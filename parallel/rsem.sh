#!/usr/bin/env bash
# author zayats1812@mail.ru
# TODO: fix hardcoded!

# Check configuration
[[ ! -z ${WASHU_ROOT} ]] || { echo "ERROR: WASHU_ROOT not configured"; exit 1; }
source ${WASHU_ROOT}/parallel/util/util.sh

>&2 echo "Batch rsem $@"
if [ $# -lt 1 ]; then
    echo "Need 1 parameter! <WORK_DIR>"
    exit 1
fi
WORK_DIR=$1

REF="$2/human_rsem_100/human_rsem_100"
RSEMPATH="/home/kzaytsev/rna_seq_pipeline/tools/RSEM-1.2.31"


cd ${WORK_DIR}

TASKS=()
for FILE in $(find . -name '*.tr.bam' | sed 's#\./##g')
do :

    NAME=${FILE%%.tr.bam} # file name without extension

    # Submit task
    run_parallel << SCRIPT
#!/bin/bash
#PBS -N rsem_exp_$NAME
#PBS -j oe
#PBS -l mem=32gb,nodes=1:ppn=8:haswell,walltime=24:00:00

## this is where pre-generated STAR reference genome

cd ${WORK_DIR}
$RSEMPATH/rsem-calculate-expression -p 8 --paired-end --bam --estimate-rspd --no-bam-output $FILE $REF $NAME
SCRIPT
    echo "FILE: ${FILE}; TASK: ${QSUB_ID}"
    TASKS+=("$QSUB_ID")
done
wait_complete ${TASKS[@]}
module load R
Rscript ${WASHU_ROOT}/parallel/util/gather.R

>&2 echo "Done. Batch rsem $@"
