#!/usr/bin/env bash
# Script to run diffreps tool differential chip-seq analysis

# Check tools
which bedtools &>/dev/null || { echo "bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>"; exit 1; }

# Load technical stuff
source $(dirname $0)/../parallel/util.sh
SCRIPT_DIR="$(project_root_dir)"
export TMPDIR=$(type job_tmp_dir &>/dev/null && echo "$(job_tmp_dir)" || echo "/tmp")
mkdir -p "${TMPDIR}"

################################################################################
# Configuration start ##########################################################
################################################################################

>&2 echo "diffreps: $@"
if [ $# -lt 3 ]; then
    echo "Need 3 parameters! <NAME> <CHROM_SIZES> <BAM_DIR>"
    exit 1
fi

NAME=$1
CHROM_SIZES=$2
BAM_DIR=$3

FOLDER=$(pwd)
echo "FOLDER"
echo $FOLDER

PREFIX="$FOLDER/$NAME"
echo "PREFIX"
echo $PREFIX

FOLDER=$(pwd)
echo "FOLDER"
echo $FOLDER

for BAM_FILE in $(find ${BAM_DIR} -name '*.bam' | sed 's#.*/##g' | grep -v 'input');
do :
    NAME=${BAM_FILE%%.bam} # file name without extension

    echo "NAME"
    echo $NAME

    # Submit task
    run_parallel << SCRIPT
#!/bin/sh
#PBS -N bamtobed_${BAM_FILE}
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${FOLDER}/${NAME}_bamtobed.log

# This is necessary because qsub default working dir is user home
cd ${FOLDER}
module load bedtools2

bedtools bamtobed -i ${BAM_DIR}/${BAM_FILE} >${NAME}.bed

SCRIPT
    echo "FILE: ${BAM_FILE}; TASK: ${QSUB_ID}"
    TASKS="$TASKS $QSUB_ID"
done
wait_complete ${TASKS}

run_parallel << SCRIPT
#!/bin/sh
#PBS -N diffReps
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${FOLDER}/diffReps.log

# This is necessary because qsub default working dir is user home
cd ${FOLDER}
module load bedtools2

diffReps.pl -co YD*.bed \
            -tr OD*.bed \
            --chrlen ${CHROM_SIZES} \
            -re diff.nb.txt -me nb
SCRIPT

wait_complete ${TASKS}


