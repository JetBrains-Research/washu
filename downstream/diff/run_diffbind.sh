#!/usr/bin/env bash

# Check tools
which bedtools &>/dev/null || { echo "bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>"; exit 1; }

# Load technical stuff
source $(dirname $0)/../../parallel/util/util.sh
SCRIPT_DIR="$(project_root_dir)"

################################################################################
# Configuration start ##########################################################
################################################################################

>&2 echo "run_diffbind: $@"
if [ $# -lt 2 ]; then
    echo "Need 2 parameters! <NAME> <DIFFBIND_CSV>"
    exit 1
fi

NAME=$1
DIFFBIND_CSV=$2

DIFFBIND=$(pwd)
echo "DIFFBIND"
echo $DIFFBIND

################################################################################
# Configuration end ############################################################
################################################################################

echo
echo "Processing diffbind"

run_parallel << SCRIPT
#!/bin/sh
#PBS -N diffbind_${NAME}
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=64gb
#PBS -j oe
#PBS -o ${NAME}_diffbind.log
# This is necessary because qsub default working dir is user home
cd ${DIFFBIND}
module load R
Rscript ${SCRIPT_DIR}/diffbind.R ${NAME}.csv
SCRIPT

wait_complete "$QSUB_ID"
