#!/usr/bin/env bash

###########################################################################
# Batch merge bams:
#    Accepts list on one or many bam containing directories.
#    In each <WORK_DIR> script runs 'samtools merge' for all its *.bam files
#    saves results to <WORK_DIR>/<WORK_DIR>_<GENOME>.bam.
###########################################################################
# author roman.chernyatchik@jetbrains.com

# Check configuration
[[ ! -z ${WASHU_ROOT} ]] || { echo "ERROR: WASHU_ROOT not configured"; exit 1; }
source ${WASHU_ROOT}/parallel/util/util.sh

>&2 echo "Batch samtools-merge $@"
if [ $# -lt 1 ]; then
    echo "Need at least one parameter! <GENOME> <WORK_DIR> [<WORK_DIR>]*"
    exit 1
fi
GENOME=$1
WORK_DIRS=${@:2}


TASKS=()
for WORK_DIR in ${WORK_DIRS}; do
    cd ${WORK_DIR}
    WORK_DIR_NAME=${WORK_DIR##*/}

    BAM_FILES=$(find . -name '*.bam' | sort)
    if [ -z "$BAM_FILES" ]; then
        # No files found
        exit -1
    else
        MERGED_FILE="${WORK_DIR_NAME}_${GENOME}.bam"

        N_BAM_FILES=$(echo "${BAM_FILES}" | wc -l);
        if [ "$N_BAM_FILES" -eq "1" ]; then
            # rename file if only one bam
            echo "rename ${FILE} -> ${MERGED_FILE}"
            for FILE in ${BAM_FILES}; do
                mv "${FILE}" "${MERGED_FILE}"
            done
        else
            # merge bams:
            BAM_FILES_ARG=$(echo "${BAM_FILES}" | tr '\n' ' ')

            # Submit task
            NAME="${MERGED_FILE%%.bam}"
            run_parallel << SCRIPT
#!/bin/sh
#PBS -N samtools_merge_${NAME}
#PBS -l nodes=1:ppn=1,walltime=2:00:00,vmem=4gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${NAME}_samtools_merge.log

# Loading modules
module load samtools

# This is necessary because qsub default working dir is user home
# and our FILE is relative
cd ${WORK_DIR}

samtools merge ${MERGED_FILE} ${BAM_FILES_ARG}

# Cleanup
rm ${BAM_FILES_ARG}

SCRIPT
            echo "FILE: ${FILE}; TASK: ${QSUB_ID}"
            TASKS+=("$QSUB_ID")
        fi
    fi
done

wait_complete ${TASKS[@]}
check_logs

>&2 echo "Done. Batch samtools-merge $@"