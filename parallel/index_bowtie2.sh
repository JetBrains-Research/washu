#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

if [ $# -lt 2 ]; then
    echo "Need 2 parameters! <GENOME> <FOLDER>"
    exit 1
fi
GENOME=$1
FOLDER=$2

# Load technical stuff, not available in qsub emulation
if [ -f "$(dirname $0)/util.sh" ]; then
    source "$(dirname $0)/util.sh"
fi

echo "Check bowtie2 indexes ${GENOME}"
cd ${FOLDER}
if ([[ ! -f "$GENOME.1.bt2" ]] && [[ ! -f "$GENOME.1.bt2l" ]]); then
    QSUB_ID=$(qsub << ENDINPUT
#!/bin/sh
#PBS -N bowtie2_indexes_${GENOME}
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=32gb
#PBS -j oe
#PBS -o ${FOLDER}/${GENOME}_bowtie2_indexes.log

# Load module
module load bowtie2

# This is necessary because qsub default working dir is user home
cd ${FOLDER}
bowtie2-build $(find . -type f -name "*.fa" | sed 's#\./##g' | paste -sd "," -) ${GENOME}
ENDINPUT
)
    wait_complete ${QSUB_ID}
    check_logs
fi
