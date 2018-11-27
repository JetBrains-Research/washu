#!/usr/bin/env bash
# author zayats1812@mail.ru
# author oleg.shpynov@jetbrains.com

which rsem-calculate-expression &>/dev/null || { echo "RSEM not found! Download RSEM: <https://github.com/deweylab/RSEM>"; exit 1; }

# Check configuration
[[ ! -z ${WASHU_ROOT} ]] || { echo "ERROR: WASHU_ROOT not configured"; exit 1; }
source ${WASHU_ROOT}/parallel/util/util.sh

>&2 echo "Batch rsem $@"
if [ $# -lt 2 ]; then
    echo "Need 2 parameter! <WORK_DIR> <STAR_REF>"
    exit 1
fi
WORK_DIR=$1
REF=$2

cd ${WORK_DIR}

TASKS=()
for FILE in $(find . -name '*.tr.bam' | sed 's#\./##g')
do :

    NAME=${FILE%%.tr.bam} # file name without extension

    # Submit task
    run_parallel << SCRIPT
#!/bin/bash
#PBS -N rsem_${NAME}
#PBS -j oe
#PBS -l mem=32gb,nodes=1:ppn=8:haswell,walltime=24:00:00
#PBS -o ${NAME}_rsem.log

## this is where pre-generated STAR reference genome

cd ${WORK_DIR}
rsem-calculate-expression -p 8 --paired-end --bam --estimate-rspd --no-bam-output ${FILE} ${REF} ${NAME}
SCRIPT
    echo "FILE: ${FILE}; TASK: ${QSUB_ID}"
    TASKS+=("$QSUB_ID")
done
wait_complete ${TASKS[@]}

module load R
Rscript ${WASHU_ROOT}/parallel/util/gather.R

>&2 echo "Done. Batch rsem $@"
