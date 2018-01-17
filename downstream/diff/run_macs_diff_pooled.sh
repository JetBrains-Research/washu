#!/usr/bin/env bash

# Check tools
which bedtools &>/dev/null || { echo "bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>"; exit 1; }
which macs2 &>/dev/null || { echo "macs2 not found! Install macs2: <https://github.com/taoliu/MACS/wiki/Install-macs2>"; exit 1; }

# Load technical stuff
source $(dirname $0)/../../parallel/util/util.sh
PROJECT_ROOT=$(project_root_dir)
export TMPDIR=$(type job_tmp_dir &>/dev/null && echo "$(job_tmp_dir)" || echo "/tmp")
mkdir -p "${TMPDIR}"

>&2 echo "run_diff_macs_pooled: $@"
if [ $# -ne 2 ]; then
    echo "Need 2 parameters! <NAME> <WORKDIR>"
    exit 1
fi

NAME=$1

BROAD_CUTOFF=0.1

DIFF_MACS_POOLED=$2

READS_Y=$(find . -name 'YD*.bam' | sed 's#\./##g' | grep -v 'input' | tr '\n' ' ')

INPUTS_Y=YD_input.bam

READS_O=$(find . -name 'OD*.bam' | sed 's#\./##g' | grep -v 'input' | tr '\n' ' ')

INPUTS_O=OD_input.bam

cd $DIFF_MACS_POOLED

echo "DIFF_MACS_POOLED"
echo $DIFF_MACS_POOLED

echo "Processing MACS2 pooled peaks and compare them";

run_parallel << SCRIPT
#!/bin/sh
#PBS -N ${NAME}_Y_macs2_broad_${BROAD_CUTOFF}
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${DIFF_MACS_POOLED}/${NAME}_Y_macs2_broad_${BROAD_CUTOFF}.log
# This is necessary because qsub default working dir is user home
cd ${DIFF_MACS_POOLED}
macs2 callpeak -t $READS_Y -c $INPUTS_Y -f BAM -g hs -n Y_${BROAD_CUTOFF} -B --broad --broad-cutoff ${BROAD_CUTOFF}
SCRIPT
QSUB_ID1=$QSUB_ID

run_parallel << SCRIPT
#!/bin/sh
#PBS -N ${NAME}_O_macs2_broad_${BROAD_CUTOFF}
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${DIFF_MACS_POOLED}/${NAME}_O_macs2_broad_${BROAD_CUTOFF}.log
# This is necessary because qsub default working dir is user home
cd ${DIFF_MACS_POOLED}
macs2 callpeak -t $READS_O -c $INPUTS_O -f BAM -g hs -n O_${BROAD_CUTOFF} -B --broad --broad-cutoff ${BROAD_CUTOFF}
SCRIPT

QSUB_ID2=$QSUB_ID
wait_complete "$QSUB_ID1 $QSUB_ID2"

check_logs
bash ${PROJECT_ROOT}/bed/compare.sh Y_${BROAD_CUTOFF}_peaks.broadPeak O_${BROAD_CUTOFF}_peaks.broadPeak ${NAME}_${BROAD_CUTOFF}

