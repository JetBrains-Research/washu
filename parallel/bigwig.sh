#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

# Check configuration
[[ ! -z ${WASHU_ROOT} ]] || { echo "ERROR: WASHU_ROOT not configured"; exit 1; }
source ${WASHU_ROOT}/parallel/util.sh

>&2 echo "Batch bigwig $@"
if [[ $# -lt 2 ]]; then
    echo "Need at least 2 parameters! <CHROM_SIZES> [<genes.gtf>] <WORK_DIR> [<WORK_DIR>]*"
    exit 1
fi
CHROM_SIZES=$1
if [[ -f $2 && $(echo $2 | grep -n '.*\.gtf$') ]]; then
    GENES_GTF=$2
    WORK_DIRS=${@:3}
else
    GENES_GFT=""
    WORK_DIRS=${@:2}
fi

TASKS=()
for WORK_DIR in ${WORK_DIRS}; do
    WORK_DIR_NAME=${WORK_DIR##*/}
    cd ${WORK_DIR}

    for FILE in $(find . -name '*.bam' | sed 's#\./##g' | grep -vE ".tr")
    do :
        NAME=${FILE%%.bam} # file name without extension
        if [[ ! -f ${FILE/bam/bw} ]]; then
            # Submit task
            run_parallel << SCRIPT
#!/bin/sh
#PBS -N bw_${WORK_DIR_NAME}_${NAME}
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${NAME}_bw.log

# This is necessary because qsub default working dir is user home
cd ${WORK_DIR}

module load bedtools2
if [[ -f ${GENES_GFT} ]];
  bash ${WASHU_ROOT}/scripts/exome2bw.sh ${FILE} ${CHROM_SIZES} ${NAME}.bw
else
  bash ${WASHU_ROOT}/scripts/reads2bw.sh ${FILE} ${CHROM_SIZES} ${NAME}.bw
fi
SCRIPT
            echo "FILE: ${WORK_DIR_NAME}/${FILE}; TASK: ${QSUB_ID}"
            TASKS+=("$QSUB_ID")
        fi
    done
done
wait_complete ${TASKS[@]}
check_logs

>&2 echo "Done. Batch bigwig $@"
