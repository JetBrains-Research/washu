#!/usr/bin/env bash

# Check tools
which bedtools &>/dev/null || { echo "ERROR: bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>"; exit 1; }
which macs2 &>/dev/null || { echo "ERROR: macs2 not found! Install macs2: <https://github.com/taoliu/MACS/wiki/Install-macs2>"; exit 1; }

# Check configuration
[[ ! -z ${WASHU_ROOT} ]] || { echo "ERROR: WASHU_ROOT not configured"; exit 1; }
source ${WASHU_ROOT}/parallel/util/util.sh

################################################################################
# Configuration start ##########################################################
################################################################################

>&2 echo "macs_bdgdiff: $@"
if [ $# -lt 3 ]; then
    echo "Need 3 parameters! <NAME> <WORKDIR> <DIFF_MACS_POOLED>"
    exit 1
fi

NAME=$1
MACS_BDGDIFF=$2
DIFF_MACS_POOLED=$3
BROAD_CUTOFF=0.1


macs2_total_tags_control() {
    echo $(cat $1 | grep "total tags in control" | sed 's/.*total tags in control: //g')
}


echo
echo "Processing $MACS_BDGDIFF"
cd ${MACS_BDGDIFF}

echo "Use MACS2 pooled peaks as input for MACS2 bdgdiff"

CONTROL_Y=$(macs2_total_tags_control ${DIFF_MACS_POOLED}/Y_${BROAD_CUTOFF}_peaks.xls)
echo "Control Y: $CONTROL_Y"
CONTROL_O=$(macs2_total_tags_control ${DIFF_MACS_POOLED}/O_${BROAD_CUTOFF}_peaks.xls)
echo "Control O: $CONTROL_O"

run_parallel << SCRIPT
#!/bin/sh
#PBS -N ${NAME}_macs2_broad_bdgdiff
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=48gb
#PBS -j oe
#PBS -o ${MACS_BDGDIFF}/${NAME}_macs2_broad_bdgdiff.log
# This is necessary because qsub default working dir is user home
cd ${MACS_BDGDIFF}
macs2 bdgdiff\
 --t1 ${DIFF_MACS_POOLED}/Y_${BROAD_CUTOFF}_treat_pileup.bdg --c1 ${DIFF_MACS_POOLED}/Y_${BROAD_CUTOFF}_control_lambda.bdg\
 --t2 ${DIFF_MACS_POOLED}/O_${BROAD_CUTOFF}_treat_pileup.bdg --c2 ${DIFF_MACS_POOLED}/O_${BROAD_CUTOFF}_control_lambda.bdg\
  --d1 ${CONTROL_Y} --d2 ${CONTROL_O} --o-prefix ${NAME}_${BROAD_CUTOFF}
SCRIPT

wait_complete $QSUB_ID

