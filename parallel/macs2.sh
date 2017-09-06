#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

which macs2 &>/dev/null || { echo "MACS2 not found! Download MACS2: <https://github.com/taoliu/MACS/wiki/Install-macs2>"; exit 1; }

# Load technical stuff, not available in qsub emulation
if [ -f "$(dirname $0)/util.sh" ]; then
    source "$(dirname $0)/util.sh"
fi

>&2 echo "Batch macs2 $@"
if [ $# -lt 5 ]; then
    echo "Need 5 parameters! <genome> <chrom.sizes> <suffix> <params_str> <work_dir> [<work_dir>]*"
    echo "if <chrom.sizes> file not specified (NONE), no signal will be created"
    exit 1
fi

GENOME=$1
CHROM_SIZES=$2
SUFFIX=$3
PARAMS=$4
WORK_DIRS=${@:5}

if [ ! -f ${CHROM_SIZES} ]; then
    echo "chrom.sizes file not specified, no signal"
fi

SPECIES=$(python $(dirname $0)/../scripts/util.py macs_species ${GENOME})

TASKS=""
for WORK_DIR in ${WORK_DIRS}; do :
    echo "${WORK_DIR}"
    WORK_DIR_NAME=${WORK_DIR##*/}
    cd ${WORK_DIR}

    for FILE in $(find . -name '*.bam' | sed 's#\./##g' | grep -v 'input')
    do :
        INPUT=$(python $(dirname $0)/../scripts/util.py find_input ${WORK_DIR}/${FILE})
        echo "${FILE}: control file: ${INPUT}"

        NAME=${FILE%%.bam} # file name without extension
        ID=${NAME}_${SUFFIX}

        # Submit task
        run_parallel << ENDINPUT
#!/bin/sh
#PBS -N macs2_${WORK_DIR_NAME}_${ID}
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=16gb
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
        bash $(dirname $0)/../scripts/bdg2bw.sh ${NAME}_signal.bdg ${CHROM_SIZES}
    fi
else
    echo "${FILE}: no control file"
    macs2 callpeak -t ${FILE} -f BAM -g ${SPECIES} -n ${ID} ${PARAMS}
fi

# Compute Reads in Peaks
bash $(dirname $0)/../reports/rip.sh ${FILE} ${ID}*.*Peak
ENDINPUT

        echo "FILE: ${WORK_DIR_NAME}/${FILE}; TASK: ${QSUB_ID}"
        TASKS="$TASKS $QSUB_ID"
    done
done
wait_complete ${TASKS}
check_logs

for WORK_DIR in ${WORK_DIRS}; do :
    cd ${WORK_DIR}

    # Cleanup BedGraph files
    rm *.bdg

    # Create pdf reports
    module load R
    MODELS=$(ls *.r); for M in ${MODELS[@]}; do Rscript $M; done
done

>&2 echo "Done. Batch macs2 $@"