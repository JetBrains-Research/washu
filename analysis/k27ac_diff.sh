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
module load bedtools

compare_peaks()
# - Two peaks aver overlapping if they share at least one nucleotide
# - Prints only peaks that overlap in all files (merged)
{
    >&2 echo "compare_peaks $@"
    PEAKS_FILE_1=$1
    NAME_1=$2
    PEAKS_FILE_2=$3
    NAME_2=$4
    TMP_FILE=intersection_${NAME_1}_${NAME_2}.txt
    # Compute common and exclusive peaks
    multiIntersectBed -i ${PEAKS_FILE_1} ${PEAKS_FILE_2} |\
    bedtools merge -c 6,7 -o max |\
    # Reproducible with 2 args: max of '0' is 2.225073859e-308 - known floating point issue in bedtools merge
    # See https://groups.google.com/forum/#!topic/bedtools-discuss/RN2U64Y5Z6Q
     awk -v OFS="\t" '{print $1,$2,$3,int($4),int($5)}' > ${TMP_FILE}

    cat ${TMP_FILE} | awk '/\t1\t0/' |\
     awk -v OFS="\t" '{for (i=1; i<=3; i++) printf("%s%s", $i, (i==3) ? "\n" : OFS)}' | sort > ${NAME_1}_exclusive.bed
    cat ${TMP_FILE} | awk '/\t0\t1/' |\
     awk -v OFS="\t" '{for (i=1; i<=3; i++) printf("%s%s", $i, (i==3) ? "\n" : OFS)}' | sort > ${NAME_2}_exclusive.bed
    cat ${TMP_FILE} | awk '/\t1\t1/' |\
     awk -v OFS="\t" '{for (i=1; i<=3; i++) printf("%s%s", $i, (i==3) ? "\n" : OFS)}' | sort > ${NAME_1}_${NAME_2}_common.bed
    # Cleanup
    rm ${TMP_FILE}
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
    bash ~/work/washu/intersect.sh ${FOLDER}/k27ac_bams_macs_broad_${Q}/YD_ac*.broadPeak > YD_peaks_${Q}.bed
    # Get aggregated peaks k27ac OD
    bash ~/work/washu/intersect.sh ${FOLDER}/k27ac_bams_macs_broad_${Q}/OD_ac*.broadPeak > OD_peaks_${Q}.bed

    compare_peaks YD_peaks_${Q}.bed YD_peaks_${Q} OD_peaks_${Q}.bed OD_peaks_${Q}
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
    compare_peaks YD_peaks_${Q}_peaks.broadPeak YD_peaks_${Q} OD_peaks_${Q}_peaks.broadPeak OD_peaks_${Q}
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
    echo "Processing MACS2 pooled $Q";
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
  --d1 ${CONTROL_YD} --d2 ${CONTROL_OD} --o-prefix diff_OD_YD_${Q}
ENDINPUT
)
    wait_complete "$QSUB_ID"
    check_logs
done