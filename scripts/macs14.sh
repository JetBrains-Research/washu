#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

which macs14 &>/dev/null || { echo "MACS14 not found! Download MACS14: <http://liulab.dfci.harvard.edu/MACS/00README.html>"; exit 1; }

# Load technical stuff
source ~/work/washu/scripts/util.sh

if [ $# -lt 3 ]; then
    echo "Need 3 parameters! <work_dir> <genome> <p>"
    exit 1
fi

WORK_DIR=$1
GENOME=$2
P=$3

echo "Batch macs14: ${WORK_DIR} ${GENOME} ${P}"

SPECIES=$(macs_species $GENOME)

cd ${WORK_DIR}

TASKS=""
for FILE in $(find . -name '*.bam' -printf '%P\n')
do :
    INPUT=$(python ~/work/washu/scripts/find_input.py ${WORK_DIR}/${FILE})
    echo "${FILE} input: ${INPUT}"

    NAME=${FILE%%.bam} # file name without extension
    ID=${NAME}_${P}

    # Submit task
    QSUB_ID=$(qsub << ENDINPUT
#!/bin/sh
#PBS -N macs2_${ID}
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${ID}_macs2.log

# This is necessary because qsub default working dir is user home
cd ${WORK_DIR}
module load bedtools2

if [ -f "${INPUT}" ]; then
    echo "${FILE}: control file found: ${INPUT}"
    macs14 -t ${FILE} -c ${INPUT} -f BAM -g ${SPECIES} -n ${ID} -p ${P}
else
    echo "${FILE}: no control file"
    macs14 -t ${FILE} -f BAM -g ${SPECIES} -n ${ID} -p ${P}
fi

# Compute Reads in Peaks
bash ~/work/washu/logs/rip.sh ${FILE} ${ID}*.narrowPeak
ENDINPUT
)
    echo "FILE: ${FILE}; JOB: ${QSUB_ID}"
    TASKS="$TASKS $QSUB_ID"
done
wait_complete ${TASKS}
check_logs

# Create pdf reports
MODELS=$(ls *.r); for M in ${MODELS[@]}; do Rscript $M; done

echo "Batch macs14: ${WORK_DIR} ${GENOME} ${P}"