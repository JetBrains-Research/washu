#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

which rseg &>/dev/null || { echo "rseg not found! Download rseg: <http://smithlabresearch.org/software/rseg/>"; exit 1; }

# Load technical stuff, not available in qsub emulation
if [ -f "$(dirname $0)/util.sh" ]; then
    source "$(dirname $0)/util.sh"
fi

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
if [[ ${GENOME} =~ ^hg19$ ]]; then
    wget http://smithlabresearch.org/data/deadzones-k36-hg19.bed
    DEADZONES="deadzones-k36-hg19.bed"
fi
if [[ ${GENOME} =~ ^mm9$ ]]; then
    wget http://smithlabresearch.org/data/deadzones-k36-mm9.bed
    DEADZONES="deadzones-k36-mm9.bed"
fi
# Check that DEADZONES are with us
if [[ -z "${DEADZONES}" ]]; then
    echo "Unknown species for macs: ${GENOME}"
    exit 1
fi


TASKS=""
for FILE in $(find . -name '*.bam' | sed 's#./##g' | grep -v 'input')
do :
    INPUT=$(python $(dirname $0)/macs_util.py find_input ${FILE})
    echo "${FILE} input: ${INPUT}"

    NAME=${FILE%%.bam} # file name without extension

    # Submit task
    QSUB_ID=$(qsub << ENDINPUT
#!/bin/sh
#PBS -N rseg_${NAME}
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${NAME}_rseg.log

# This is necessary because qsub default working dir is user home
cd ${WORK_DIR}

if [ -f "${INPUT}" ]; then
    echo "${FILE}: control file found: ${INPUT}"
    rseg-diff -c ${GENOME}_chrom_sizes.bed -o ${NAME}.bed -i 20 -v -d ${DEADZONES} -mode 2 ${FILE} ${INPUT}
else
    echo "${FILE}: no control file"
    rseg -c ${GENOME}_chrom_sizes.bed -o ${NAME}.bed -i 20 -v -d ${DEADZONES} ${FILE}
fi
ENDINPUT
)
    echo "FILE: ${FILE}; JOB: ${QSUB_ID}"
    TASKS="$TASKS $QSUB_ID"
done
wait_complete ${TASKS}
check_logs

echo "Done. Batch rseg: ${WORK_DIR} ${GENOME} ${CHROM_SIZES}"