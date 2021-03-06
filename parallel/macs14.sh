#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

which macs14 &>/dev/null || { echo "ERROR: MACS14 not found! Download MACS14: <http://liulab.dfci.harvard.edu/MACS/00README.html>"; exit 1; }

# Check configuration
[[ ! -z ${WASHU_ROOT} ]] || { echo "ERROR: WASHU_ROOT not configured"; exit 1; }
source ${WASHU_ROOT}/parallel/util.sh

>&2 echo "Batch macs14 $@"
if [[ $# -lt 3 ]]; then
    echo "Need 3 parameters! <work_dir> <genome> <p>"
    exit 1
fi

WORK_DIR=$1
GENOME=$2
P=$3

SPECIES=$(python ${WASHU_ROOT}/scripts/util.py macs_species ${GENOME})

cd ${WORK_DIR}

TASKS=()
for FILE in $(find . -name '*.bam' | sed 's#\./##g' | grep -v 'input')
do :
    INPUT=$(python ${WASHU_ROOT}/scripts/util.py find_input ${WORK_DIR}/${FILE})
    echo "${FILE}: control file: ${INPUT}"

    NAME=${FILE%%.bam} # file name without extension
    ID=${NAME}_${P}

    # Submit task
    run_parallel << SCRIPT
#!/bin/sh
#PBS -N macs2_${ID}
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${ID}_macs2.log

# This is necessary because qsub default working dir is user home
cd ${WORK_DIR}
module load bedtools2

if [[ -f "${INPUT}" ]]; then
    echo "${FILE}: control file found: ${INPUT}"
    macs14 -t ${FILE} -c ${INPUT} -f BAM -g ${SPECIES} -n ${ID} -p ${P}
else
    echo "${FILE}: no control file"
    macs14 -t ${FILE} -f BAM -g ${SPECIES} -n ${ID} -p ${P}
fi
SCRIPT
    echo "FILE: ${FILE}; TASK: ${QSUB_ID}"
    TASKS+=("$QSUB_ID")
done
wait_complete ${TASKS[@]}
check_logs

# Create pdf reports
MODELS=$(ls *.r); for M in ${MODELS[@]}; do Rscript $M; done

>&2 echo "Done. Batch macs14 $@"