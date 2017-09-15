#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

# Load technical stuff
source $(dirname $0)/../parallel/util.sh
SCRIPT_DIR="$(project_root_dir)"

>&2 echo "Batch bigwig $@"
if [ $# -lt 2 ]; then
    echo "Need at least 2 parameters! <CHROM_SIZES> <WORK_DIR> [<WORK_DIR>]*"
    exit 1
fi
CHROM_SIZES=$1
WORK_DIRS=${@:2}

TASKS=""
for WORK_DIR in ${WORK_DIRS}; do
    WORK_DIR_NAME=${WORK_DIR##*/}
    cd ${WORK_DIR}

    for FILE in $(find . -name '*.bam' | sed 's#\./##g' | grep -vE ".tr")
    do :
        NAME=${FILE%%.bam} # file name without extension

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
bash ${SCRIPT_DIR}/scripts/bam2bw.sh ${FILE} ${CHROM_SIZES}
SCRIPT
        echo "FILE: ${WORK_DIR_NAME}/${FILE}; TASK: ${QSUB_ID}"
        TASKS="$TASKS $QSUB_ID"
    done
done
wait_complete ${TASKS}
check_logs

>&2 echo "Done. Batch bigwig $@"
