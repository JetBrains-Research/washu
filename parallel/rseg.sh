#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

which rseg &>/dev/null || { echo "rseg not found! Download rseg: <http://smithlabresearch.org/software/rseg/>"; exit 1; }

# Load technical stuff
source $(dirname $0)/../parallel/util.sh
SCRIPT_DIR="$(project_root_dir)"

>&2 echo "Batch rseg $@"
if [ $# -lt 3 ]; then
    echo "Need 3 parameters! <work_dir> <genome> <chrom.sizes>"
    exit 1
fi

WORK_DIR=$1
GENOME=$2
CHROM_SIZES=$3

cd ${WORK_DIR}

echo "Prepare chrom_sized.bed"
cat ${CHROM_SIZES} | awk '{print $1, 1, $2}' > ${GENOME}_chrom_sizes.bed

echo "Obtain deadzones file for ${GENOME}"
if [[ ${GENOME} =~ ^hg19$ ]]; then
    DEADZONES="deadzones-k36-hg19.bed"
fi
if [[ ${GENOME} =~ ^mm9$ ]]; then
    DEADZONES="deadzones-k36-mm9.bed"
fi
# Check that DEADZONES are with us
if [[ -z "${DEADZONES}" ]]; then
    echo "Unknown species for rseg: ${GENOME}"
    exit 1
fi
# Download DEADZONES file
if [[ ! -f ${DEADZONES} ]]; then
    wget "http://smithlabresearch.org/data/${DEADZONES}"
fi

TASKS=""
for FILE in $(find . -name '*.bam' | sed 's#\./##g' | grep -v 'input')
do :
    NAME=${FILE%%.bam} # file name without extension
    FILE_BED=${NAME}.bed

    INPUT=$(python ${SCRIPT_DIR}/scripts/util.py find_input ${WORK_DIR}/${FILE})
    echo "${FILE} input: ${INPUT}"
    INPUT_BED=${INPUT%%.bam}.bed

    TMP_FOLDER=${WORK_DIR}/rseg_tmp/${NAME}
    mkdir -p ${TMP_FOLDER}

    # Submit task
    run_parallel << ENDINPUT
#!/bin/sh
#PBS -N rseg_${NAME}
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${NAME}_rseg.log

module load bedtools2

# This is necessary because qsub default working dir is user home
cd ${WORK_DIR}

# RSEG works with BED only
# See sort instructions at http://smithlabresearch.org/wp-content/uploads/rseg_manual_v0.4.4.pdf
export LC_ALL=C
bedtools bamtobed -i ${FILE} | sort -k1,1 -k3,3n -k2,2n -k6,6 > ${TMP_FOLDER}/${FILE_BED}

if [ -f "${INPUT}" ]; then
    # Use tmp files to reduced async problems with same input parallel processing
    echo "${FILE}: control file found: ${INPUT}"
    if [ ! -f ${INPUT_BED} ]; then
        bedtools bamtobed -i ${INPUT} | sort -k1,1 -k3,3n -k2,2n -k6,6 > ${TMP_FOLDER}/${INPUT_BED}
        # Check that we are the first in async calls, not 100% safe
        if [ ! -f ${INPUT_BED} ]; then
            cp ${TMP_FOLDER}/${INPUT_BED} ${WORK_DIR}
        fi
    fi

    rseg-diff -c ${GENOME}_chrom_sizes.bed -o ${NAME}_domains.bed -i 20 -v -d ${DEADZONES} -mode 2 ${TMP_FOLDER}/${FILE_BED} ${INPUT_BED}
else
    echo "${FILE}: no control file"
    rseg -c ${GENOME}_chrom_sizes.bed -o ${NAME}_domains.bed -i 20 -v -d ${DEADZONES} ${TMP_FOLDER}/${FILE_BED}
fi

# Compute Reads in Peaks
bash ${SCRIPT_DIR}/reports/rip.sh ${FILE} ${NAME}_domains.bed
ENDINPUT

    echo "FILE: ${FILE}; TASK: ${QSUB_ID}"
    TASKS="$TASKS $QSUB_ID"
done
wait_complete ${TASKS}
check_logs

# Cleanup
rm -r ${WORK_DIR}/rseg_tmp
>&2 echo "Done. Batch rseg $@"