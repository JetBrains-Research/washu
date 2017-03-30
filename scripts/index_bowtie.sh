#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

if [ $# -lt 2 ]; then
    echo "Need 2 parameters! <GENOME> <FOLDER>"
    exit 1
fi
GENOME=$1
FOLDER=$2

# Load technical stuff
source ~/work/washu/scripts/util.sh

echo "Check bowtie indexes ${GENOME} ${FOLDER}"
cd ${FOLDER}
# Check both 32 and 64 large indexes
if ([[ ! -f "$GENOME.1.ebwt" ]] && [[ ! -f "$GENOME.1.ebwtl" ]]); then
    QSUB_ID=$(qsub << ENDINPUT
#!/bin/sh
#PBS -N bowtie_indexes_${GENOME}
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${FOLDER}/${GENOME}_bowtie_indexes.log

# Load module
module load bowtie

# This is necessary because qsub default working dir is user home
cd ${FOLDER}
bowtie-build $(find . -type f -name "*.fa" -printf '%P\n' | paste -sd "," -) ${GENOME}
ENDINPUT
)
    wait_complete ${QSUB_ID}
    check_logs
fi