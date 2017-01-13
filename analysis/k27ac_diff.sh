#!/usr/bin/env bash
# Script to process analysis analysis of K27ac OD vs YD
# author Oleg Shpynov (oleg.shpynov@jetbrains.com)

# Check tools
which bedtools &>/dev/null || { echo "bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>"; exit 1; }
which macs2 &>/dev/null || { echo "macs2 not found! Install macs2: <https://github.com/taoliu/MACS/wiki/Install-macs2>"; exit 1; }

# QSUB mock replacement
which qsub &>/dev/null || {
    echo "QSUB system not found, using mock replacement"
    qsub() {
        while read -r line; do CMD+=$line; CMD+=$'\n'; done;
        echo "$CMD" > /tmp/qsub.sh
        LOG=$(echo "$CMD" | grep "#PBS -o" | sed 's/#PBS -o //g')
        bash /tmp/qsub.sh | tee "$LOG"
    }
    qstat() {
        echo ""
    }
    module() {
        echo ""
    }
}

# Load technical stuff
source ~/work/washu/scripts/util.sh

# Configure folder
FOLDER="/scratch/artyomov_lab_aging/Y10OD10/chipseq/processed/3vs3_2"
if [ ! -d $FOLDER ]; then
    FOLDER="/mnt/stripe/bio/raw-data/aging/Y10OD10/chipseq/processed/3vs3_2"
fi
if [ ! -d $FOLDER ]; then
    FOLDER="/Volumes/WD/scratch/artyomov_lab_aging/Y10OD10/chipseq/processed/3vs3_2"
fi
echo "Processing k27ac difference: $FOLDER"

# Hillbilly strategy
DIFF_HB="${FOLDER}/k27ac_diff_hb"
mkdir $DIFF_HB
cd $DIFF_HB
# * Call peaks independently
# * Intersect individual peaks to get common peaks
# * Analyze OD and YD common peaks
QS=( 0.01 )
for Q in "${QS[@]}"; do
    echo "Processing HillBilly $Q";
    # Get aggregated peaks k27ac YD
    bash ~/work/washu/bed/intersect.sh ${FOLDER}/k27ac_bams_macs_broad_${Q}/YD_ac*.broadPeak > YD_peaks_${Q}.bed
    # Get aggregated peaks k27ac OD
    bash ~/work/washu/bed/intersect.sh ${FOLDER}/k27ac_bams_macs_broad_${Q}/OD_ac*.broadPeak > OD_peaks_${Q}.bed

    bash ~/work/washu/bed/compare.sh YD_peaks_${Q}.bed OD_peaks_${Q}.bed diff_YD_OD_${Q}
done


# MACS2 pooled peaks
# * Call pooled peaks for OD and YD
DIFF_MACS_POOLED="${FOLDER}/k27ac_diff_macs_pooled"
mkdir $DIFF_MACS_POOLED
cd $DIFF_MACS_POOLED
WORK_DIR=$DIFF_MACS_POOLED

for Q in "${QS[@]}"; do
    echo "Processing MACS2 pooled $Q";
    # Get aggregated peaks k27ac YD
    QSUB_ID_YD=$(qsub << ENDINPUT
#!/bin/sh
#PBS -N macs2_broad_k27ac_YD
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${WORK_DIR}/macs2_broad_k27ac_YD_${Q}.log
# This is necessary because qsub default working dir is user home
cd ${WORK_DIR}
module load bedtools2
macs2 callpeak -t ${FOLDER}/k27ac_bams/YD_ac*.bam -c ${FOLDER}/k27ac_bams/YD_input.bam\
 -f BAM -g hs -n YD_peaks_${Q} -B --broad --broad-cutoff ${Q}
ENDINPUT
)

    QSUB_ID_OD=$(qsub << ENDINPUT
#!/bin/sh
#PBS -N macs2_broad_k27ac_OD
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${WORK_DIR}/macs2_broad_k27ac_OD_${Q}.log
# This is necessary because qsub default working dir is user home
cd ${WORK_DIR}
module load bedtools2
macs2 callpeak -t ${FOLDER}/k27ac_bams/OD_ac*.bam -c ${FOLDER}/k27ac_bams/OD_input.bam\
 -f BAM -g hs -n OD_peaks_${Q} -B --broad --broad-cutoff ${Q}
ENDINPUT
)
    wait_complete "$QSUB_ID_YD $QSUB_ID_OD"
    check_logs
    bash ~/work/washu/bed/compare.sh YD_peaks_${Q}_peaks.broadPeak OD_peaks_${Q}_peaks.broadPeak diff_YD_OD_${Q}
done



# MACS2 bdgdiff
# See https://github.com/taoliu/MACS/wiki/Call-differential-binding-events
# * Call pooled peaks for OD and YD
# * Use macs2 bdgdiff
DIFF_MACS_BDGDIFF="${FOLDER}/k27ac_diff_macs_bdgdiff"
mkdir $DIFF_MACS_BDGDIFF
cd $DIFF_MACS_BDGDIFF
WORK_DIR=$DIFF_MACS_BDGDIFF

for Q in "${QS[@]}"; do
    echo "Processing MACS2 bdgdiff $Q";
    CONTROL_OD=$(cat ${DIFF_MACS_POOLED}/macs2_broad_k27ac_OD_${Q}.log |\
     grep "total tags in control" | sed 's/.*total tags in control: //g')
    echo "Control OD: $CONTROL_OD"
    CONTROL_YD=$(cat ${DIFF_MACS_POOLED}/macs2_broad_k27ac_YD_${Q}.log |\
     grep "total tags in control" | sed 's/.*total tags in control: //g')
    echo "Control YD: $CONTROL_YD"

    QSUB_ID=$(qsub << ENDINPUT
#!/bin/sh
#PBS -N macs2_broad_k27ac_bdgdiff
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${WORK_DIR}/macs2_bdgdiff_k27ac_${Q}.log
# This is necessary because qsub default working dir is user home
cd ${WORK_DIR}
module load bedtools2
macs2 bdgdiff\
 --t1 ${DIFF_MACS_POOLED}/YD_peaks_${Q}_treat_pileup.bdg --c1 ${DIFF_MACS_POOLED}/YD_peaks_${Q}_control_lambda.bdg\
 --t2 ${DIFF_MACS_POOLED}/OD_peaks_${Q}_treat_pileup.bdg --c2 ${DIFF_MACS_POOLED}/OD_peaks_${Q}_control_lambda.bdg\
  --d1 ${CONTROL_YD} --d2 ${CONTROL_OD} --o-prefix diff_YD_OD_${Q}
ENDINPUT
)
    wait_complete "$QSUB_ID"
    check_logs
done