#!/usr/bin/env bash

# Check tools
which bedtools &>/dev/null || { echo "bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>"; exit 1; }
which macs2 &>/dev/null || { echo "macs2 not found! Install macs2: <https://github.com/taoliu/MACS/wiki/Install-macs2>"; exit 1; }

# Load technical stuff
source $(dirname $0)/../parallel/util.sh
SCRIPT_DIR="$(project_root_dir)"
export TMPDIR=$(type job_tmp_dir &>/dev/null && echo "$(job_tmp_dir)" || echo "/tmp")
mkdir -p "${TMPDIR}"

>&2 echo "run_diff_macs_pooled: $@"
if [ $# -ne 2 ]; then
    echo "Need 2 parameters! <NAME> <WORKDIR>"
    exit 1
fi

NAME=$1

BROAD_CUTOFF=0.1

MACS_POOLED_Y_VS_O=$2

READS_Y=$(find . -name 'YD*.bam' | sed 's#\./##g' | grep -v 'input' | tr '\n' ' ')

INPUTS_Y=YD_input.bam

READS_O=$(find . -name 'OD*.bam' | sed 's#\./##g' | grep -v 'input' | tr '\n' ' ')

INPUTS_O=OD_input.bam

cd ${MACS_POOLED_Y_VS_O}

echo "MACS_POOLED_Y_VS_O"
echo ${MACS_POOLED_Y_VS_O}

echo
echo "Processing MACS2 pooled Y vs O as control and vice versa"

run_parallel << SCRIPT
#!/bin/sh
#PBS -N ${NAME}_Y_vs_O_macs2_broad
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${MACS_POOLED_Y_VS_O}/${NAME}_Y_vs_O_macs2_broad.log
# This is necessary because qsub default working dir is user home
cd ${MACS_POOLED_Y_VS_O}
macs2 callpeak -t $READS_Y -c $READS_O -f BAM -g hs -n ${NAME}_Y_vs_O_${BROAD_CUTOFF} -B --broad --broad-cutoff ${BROAD_CUTOFF}
SCRIPT
QSUB_ID_Y_vs_O=$QSUB_ID

run_parallel << SCRIPT
#!/bin/sh
#PBS -N ${NAME}_O_vs_Y_macs2_broad
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${MACS_POOLED_Y_VS_O}/${NAME}_O_vs_Y_macs2_broad.log
# This is necessary because qsub default working dir is user home
cd ${MACS_POOLED_Y_VS_O}
macs2 callpeak -t $READS_O -c $READS_Y -f BAM -g hs -n ${NAME}_O_vs_Y_${BROAD_CUTOFF} -B --broad --broad-cutoff ${BROAD_CUTOFF}
SCRIPT

QSUB_ID_O_vs_Y=$QSUB_ID

wait_complete "$QSUB_ID_Y_vs_O $QSUB_ID_O_vs_Y"

