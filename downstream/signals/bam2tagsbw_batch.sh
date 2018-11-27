#!/usr/bin/env bash
# Script to batch convert BAM files to BigWig files for tags
# author roman.chernyatchik@jetbrains.com

# Check configuration
[[ ! -z ${WASHU_ROOT} ]] || { echo "ERROR: WASHU_ROOT not configured"; exit 1; }
source ${WASHU_ROOT}/parallel/util.sh

>&2 echo "Batch signals $@"
if [ $# -lt 4 ]; then
    echo "Need at least 4 parameters! <WORK_DIR> <FRAGMENT> <CHROM.SIZES> <BAMS_DIR1> [<BAMS_DIR2>, ...]"
    exit 1
fi

WORK_DIR=$(expand_path $1)
FRAGMENT=$2
CHROM_SIZES=$(expand_path $3)
BAMS_DIRS="${@:4}"

echo "WORK_DIR: $WORK_DIR"
echo "FRAGMENT: $FRAGMENT"
echo "CHROM_SIZES: $CHROM_SIZES"
echo "BAMS_DIRS: $BAMS_DIRS"

########################################################################################################################
# Prepare BW files
########################################################################################################################
echo "Prepare tags BW in ${WORK_DIR}"
TAGS_BW_LOGS="${WORK_DIR}/tags_bw_logs"
if [[ ! -d "${TAGS_BW_LOGS}" ]]; then
    mkdir -p ${TAGS_BW_LOGS}
fi
TASKS=()

cd ${WORK_DIR}
for BAMS_DIR in ${BAMS_DIRS}
do
    for BAM in $(find $(expand_path ${BAMS_DIR}) -name '*.bam')
    do
        # TODO add conversion: YD18_k27ac_hg19_unique_tags.bw -> YD_YD18_H3K27ac_tags.bw
        FILE_NAME=${BAM##*/}
        NAME=${FILE_NAME%%.bam} # file name without extension
        TAGS_BW_NAME=${NAME}_tags_${FRAGMENT}
        RESULT=${WORK_DIR}/${TAGS_BW_NAME}.bw
        if [[ ! -f ${RESULT} ]]; then
            # Submit task
            run_parallel << SCRIPT
#!/bin/sh
#PBS -N ${RESULT}
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${TAGS_BW_LOGS}/${TAGS_BW_NAME}.log

# This is necessary because qsub default working dir is user home
cd ${WORK_DIR}

module load bedtools2
bash ${WASHU_ROOT}/downstream/signals/bam2tagsbw.sh ${BAM} ${FRAGMENT} ${CHROM_SIZES} ${RESULT}
SCRIPT
            echo "FILE: ${FILE_NAME}; TASK: ${QSUB_ID}"
            TASKS+=("$QSUB_ID")
        fi
    done
done

wait_complete ${TASKS[@]}
cd ${WORK_DIR}
check_logs

export TMPDIR=$(type job_tmp_dir &>/dev/null && echo "$(job_tmp_dir)" || echo "/tmp")
mkdir -p $TMPDIR

# Cleanup
type clean_job_tmp_dir &>/dev/null && clean_job_tmp_dir

>&2 echo "Done. Batch signals $@"