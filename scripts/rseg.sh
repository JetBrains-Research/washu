#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

which rseg &>/dev/null || { echo "rseg not found! Download rseg: <http://smithlabresearch.org/software/rseg/>"; exit 1; }

# Load technical stuff
source ~/work/washu/scripts/util.sh

if [ $# -lt 3 ]; then
    echo "Need 3 parameters! <work_dir> <genome> <chrom.sizes>"
    exit 1
fi

WORK_DIR=$1
GENOME=$2
CHROM_SIZES=$3

echo "Batch rseg: ${WORK_DIR} ${GENOME} ${CHROM_SIZES}"
cd ${WORK_DIR}

echo "Prepare chrom_sized.bed"
cat ${CHROM_SIZES} | awk '{print $1, 1, $2}' > ${GENOME}_chrom_sizes.bed

echo "Obtain deadzones file for ${GENOME}"
[[ ${GENOME} =~ ^hg18$ ]] && { wget http://smithlabresearch.org/data/deadzones-k36-hg18.bed; DEADZONES="deadzones-k36-hg18.bed" }
[[ ${GENOME} =~ ^hg19$ ]] && { wget http://smithlabresearch.org/data/deadzones-k36-hg19.bed; DEADZONES="deadzones-k36-hg19.bed" }
[[ ${GENOME} =~ ^mm9$ ]] && { wget http://smithlabresearch.org/data/deadzones-k36-mm9.bed; DEADZONES="deadzones-k36-mm9.bed" }
[[ -z "${DEADZONES}" ]] && echo "Unknown species for macs: ${GENOME}" && exit 1


TASKS=""
for FILE in $(find . -name '*.bam' -printf '%P\n')
do :
    INPUT=$(python ~/work/washu/scripts/macs2_find_input.py ${FILE})
    echo "${FILE} input: ${INPUT}"

    NAME=${FILE%%.bam} # file name without extension
    ID=${NAME}

    # Submit task
    QSUB_ID=$(qsub << ENDINPUT
#!/bin/sh
#PBS -N rseg_${ID}
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${ID}_rseg.log

# This is necessary because qsub default working dir is user home
cd ${WORK_DIR}

if [ -f "${INPUT}" ]; then
    echo "${FILE}: control file found: ${INPUT}"
    rseg-diff -c ${GENOME}_chrom_sizes.bed -o ${ID}.bed -i 20 -v -d ${DEADZONES} -mode 2 ${FILE} ${INPUT}
else
    echo "${FILE}: no control file"
    rseg -c ${GENOME}_chrom_sizes.bed -o ${ID}.bed -i 20 -v -d ${DEADZONES} ${FILE}
fi
ENDINPUT
)
    echo "FILE: ${FILE}; JOB: ${QSUB_ID}"
    TASKS="$TASKS $QSUB_ID"
done
wait_complete ${TASKS}
check_logs

echo "Done. Batch rseg: ${WORK_DIR} ${GENOME} ${CHROM_SIZES}"