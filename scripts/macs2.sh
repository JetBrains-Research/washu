#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

which macs2 &>/dev/null || { echo "MACS2 not found! Download MACS2: <https://github.com/taoliu/MACS/wiki/Install-macs2>"; exit 1; }

# Load technical stuff, not available in qsub emulation
if [ -f "$(dirname $0)/util.sh" ]; then
    source "$(dirname $0)/util.sh"
fi

if [ $# -lt 5 ]; then
    echo "Need 5 parameters! <work_dir> <genome> <chrom.sizes> <suffix> <params>"
    echo "if <chrom.sizes> file not specified (NONE), no signal will be created"
    exit 1
fi

WORK_DIR=$1
GENOME=$2
CHROM_SIZES=$3
SUFFIX=$4
shift 4
PARAMS=$@

echo "Batch macs2: ${WORK_DIR} ${GENOME} ${CHROM_SIZES} ${SUFFIX} ${PARAMS}"
if [ ! -f ${CHROM_SIZES} ]; then
    echo "chrom.sizes file not specified, no signal"
fi

SPECIES=$(macs_species $GENOME)

cd ${WORK_DIR}

TASKS=""
for FILE in $(find . -name '*.bam' | sed 's#./##g' | grep -v 'input')
do :
    INPUT=$(python $(dirname $0)/util.py find_input ${WORK_DIR}/${FILE})
    echo "${FILE}: control file: ${INPUT}"

    NAME=${FILE%%.bam} # file name without extension
    ID=${NAME}_${SUFFIX}

    # Submit task
    QSUB_ID=$(qsub << ENDINPUT
#!/bin/sh
#PBS -N macs2_${ID}
#PBS -l nodes=1:ppn=4,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${ID}_macs2.log

# This is necessary because qsub default working dir is user home
cd ${WORK_DIR}
# Required for signal track processing
module load bedtools2

if [ -f "${INPUT}" ]; then
    echo "${FILE}: control file found: ${INPUT}"
    macs2 callpeak -t ${FILE} -c ${INPUT} -f BAM -g ${SPECIES} -n ${ID} ${PARAMS}

    if [ -f "${CHROM_SIZES}" ]; then
        echo "Create fold enrichment signal track for ${FILE} and ${INPUT}"
        macs2 bdgcmp -t ${ID}_treat_pileup.bdg -c ${ID}_control_lambda.bdg -o ${NAME}_signal.bdg -m FE
        bash $(dirname $0)/../bdg2bw.sh ${NAME}_signal.bdg ${CHROM_SIZES}
    fi
else
    echo "${FILE}: no control file"
    macs2 callpeak -t ${FILE} -f BAM -g ${SPECIES} -n ${ID} ${PARAMS}
fi

# Compute Reads in Peaks
bash $(dirname $0)/../reports/rip.sh ${FILE} ${ID}*.*Peak
ENDINPUT
)
    echo "FILE: ${FILE}; JOB: ${QSUB_ID}"
    TASKS="$TASKS $QSUB_ID"
done
wait_complete ${TASKS}
check_logs

# Create pdf reports
module load R
MODELS=$(ls *.r); for M in ${MODELS[@]}; do Rscript $M; done

echo "DONE. Batch macs2: ${WORK_DIR} ${GENOME} ${CHROM_SIZES} ${SUFFIX} ${PARAMS}"