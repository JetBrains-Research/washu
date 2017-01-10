#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

# Load technical stuff
source ~/work/washu/scripts/util.sh

if [ $# -lt 4 ]; then
    echo "Need 4 parameters! <work_dir> <genome> <q> <chrom.sizes>"
    exit 1
fi

WORK_DIR=$1
GENOME=$2
Q=$3
CHROM_SIZES=$4

echo "Batch macs2: ${WORK_DIR} ${GENOME} ${Q} ${CHROM_SIZES}"

SPECIES=$(macs2_species $GENOME)

cd ${WORK_DIR}

TASKS=""
for FILE in $(find . -name '*.bam' -printf '%P\n')
do :
    INPUT=$(python ~/work/washu/scripts/find_input.py ${WORK_DIR}/${FILE})
    echo "${FILE} input: ${INPUT}"

    NAME=${FILE%%.bam} # file name without extension
    ID=${NAME}_${Q}

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
    macs2 callpeak -t ${FILE} -c ${INPUT} -f BAM -g ${SPECIES} -n ${ID} -B -q ${Q}
else
    echo "${FILE}: no control file"
    macs2 callpeak -t ${FILE} -f BAM -g ${SPECIES} -n ${ID} -B -q ${Q}
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

echo "Batch macs2: ${WORK_DIR} ${GENOME} ${Q} ${CHROM_SIZES}"