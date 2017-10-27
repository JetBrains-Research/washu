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
if [ $# -lt 4 ]; then
    echo "Need 4 parameters! <NAME> <CHROM_SIZES> <BAM_DIR> <DIFFREPS_PARAMS>"
    exit 1
fi

NAME=$1
CHROM_SIZES=$2
BAM_DIR=$3
DIFFREPS_PARAMS=$4

FOLDER=$(pwd)
echo "FOLDER"
echo $FOLDER

PREFIX="$FOLDER/$NAME"
echo "PREFIX"
echo $PREFIX


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
            -re diff.nb.txt \
            ${DIFFREPS_PARAMS} \
            --nproc 8
SCRIPT

wait_complete ${TASKS}

# Cleanup
rm *.bed

cat diff.nb.txt.hotspot | tail -n +4 | cut -f 1-3 >hotspot.bed

cat diff.nb.txt | grep "^chr" | cut -f1-3,11-14 | sort -g -k7 | awk '{ if($7 <= 0.02) { print }}' >enriched.txt

cat enriched.txt | cut -f1-3 >enriched.bed
grep Up enriched.txt | cut -f1-3 >enriched_up.bed
grep Down enriched.txt | cut -f1-3 >enriched_down.bed
