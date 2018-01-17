#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

# Load technical stuff
source ${WASHU_ROOT}/parallel/util/util.sh

>&2 echo "Batch tags_bigwig $@"
if [ $# -lt 2 ]; then
    echo "Need at least 2 parameters! <CHROM_SIZES> <INSERT_SIZE> <WORK_DIR> [<WORK_DIR>]*"
    exit 1
fi
CHROM_SIZES=$1
INSERT_SIZE=$2
WORK_DIRS=${@:3}

export TMPDIR=$(type job_tmp_dir &>/dev/null && echo "$(job_tmp_dir)" || echo "/tmp")
mkdir -p $TMPDIR

TASKS=()
for WORK_DIR in ${WORK_DIRS}; do
    WORK_DIR_NAME=${WORK_DIR##*/}
    cd ${WORK_DIR}

    for FILE in $(find . -name '*.bam' | sed 's#\./##g' | grep -vE ".tr")
    do :
        NAME=${FILE%%.bam} # file name without extension
        if [[ ! -f ${WORK_DIR}/${NAME}_tags.bw ]]; then
            # Submit task
            run_parallel << SCRIPT
#!/bin/sh
#PBS -N tags_bw_${WORK_DIR_NAME}_${NAME}
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${NAME}_tags_bw.log

# This is necessary because qsub default working dir is user home
cd ${WORK_DIR}

module load bedtools2
bash ${WASHU_ROOT}/scripts/reads2tagsbw.sh ${FILE} ${INSERT_SIZE} ${CHROM_SIZES}
SCRIPT
            echo "FILE: ${WORK_DIR_NAME}/${FILE}; TASK: ${QSUB_ID}"
            TASKS+=("$QSUB_ID")
        fi
    done
done
wait_complete ${TASKS[@]}
check_logs

# TMP dir cleanup:
type clean_job_tmp_dir &>/dev/null && clean_job_tmp_dir

>&2 echo "Done. Batch tags_bigwig $@"
