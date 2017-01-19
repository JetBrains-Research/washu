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
BASE="/scratch/artyomov_lab_aging/Y10OD10"
if [ ! -d $BASE ]; then
    BASE="/mnt/stripe/bio/raw-data/aging/Y10OD10"
fi
if [ ! -d $BASE ]; then
    BASE="/Volumes/WD/scratch/artyomov_lab_aging/Y10OD10"
fi
FOLDER="$BASE/chipseq/processed/3vs3_2"
echo "Processing k27ac difference: $FOLDER"
CHROM_SIZES="$BASE/../indexes/hg19/hg19.chrom.sizes"

# Hillbilly strategy
DIFF_HB="${FOLDER}/k27ac_diff_hb"
if [ ! -d ${DIFF_HB} ]; then
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
fi

# MACS2 pooled peaks
# * Call pooled peaks for OD and YD
DIFF_MACS_POOLED="${FOLDER}/k27ac_diff_macs_pooled"
if [ ! -d $DIFF_MACS_POOLED ]; then
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
 -f BAM -g hs -n YD_${Q} -B --broad --broad-cutoff ${Q}
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
 -f BAM -g hs -n OD_${Q} -B --broad --broad-cutoff ${Q}
ENDINPUT
)
        wait_complete "$QSUB_ID_YD $QSUB_ID_OD"
        check_logs
        bash ~/work/washu/bed/compare.sh YD_${Q}_peaks.broadPeak OD_${Q}_peaks.broadPeak diff_YD_OD_${Q}
    done
fi

macs2_total_tags_control() {
    echo $(cat $1 | grep "total tags in control" | sed 's/.*total tags in control: //g')
}

# MACS2 bdgdiff
# See https://github.com/taoliu/MACS/wiki/Call-differential-binding-events
# * Call pooled peaks for OD and YD
# * Use macs2 bdgdiff
DIFF_MACS_BDGDIFF="${FOLDER}/k27ac_diff_macs_bdgdiff"
if [ ! -d $DIFF_MACS_BDGDIFF ]; then
    mkdir $DIFF_MACS_BDGDIFF
    cd $DIFF_MACS_BDGDIFF
    WORK_DIR=$DIFF_MACS_BDGDIFF

    for Q in "${QS[@]}"; do
        echo "Processing MACS2 bdgdiff $Q";
        CONTROL_OD=$(macs2_total_tags_control ${DIFF_MACS_POOLED}/macs2_broad_k27ac_OD_${Q}.log)
        echo "Control OD: $CONTROL_OD"
        CONTROL_YD=$(macs2_total_tags_control ${DIFF_MACS_POOLED}/macs2_broad_k27ac_YD_${Q}.log)
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
 --t1 ${DIFF_MACS_POOLED}/YD_${Q}_treat_pileup.bdg --c1 ${DIFF_MACS_POOLED}/YD_${Q}_control_lambda.bdg\
 --t2 ${DIFF_MACS_POOLED}/OD_${Q}_treat_pileup.bdg --c2 ${DIFF_MACS_POOLED}/OD_${Q}_control_lambda.bdg\
  --d1 ${CONTROL_YD} --d2 ${CONTROL_OD} --o-prefix diff_YD_OD_${Q}
ENDINPUT
)
        wait_complete "$QSUB_ID"
        check_logs
    done
fi


bams_to_tags() {
    OUT=$1
    # Shift arguments
    shift 1
    for F in $@; do
        >&2 echo $F
        bedtools bamtobed -i ${F} | grep -E "chr[0-9]+|chrX" | awk '{print $1, $2, $6}' >> $OUT
    done
}

# Pooled ChIPDiff
CHIPDIFF="${FOLDER}/k27ac_diff_chipdiff"
if [ ! -d $CHIPDIFF ]; then
    mkdir $CHIPDIFF
    cd $CHIPDIFF
    WORK_DIR=$CHIPDIFF

    cat >config.txt <<CONFIG
maxIterationNum  500
minP             0.95
maxTrainingSeqNum 10000
minFoldChange    3
minRegionDist    1000
CONFIG

    >&2 echo "Processing OD Tags";
    bams_to_tags OD.tag $(find ${FOLDER}/k27ac_bams/ -name 'OD_ac*.bam')
    >&2 echo "Processing YD Tags";
    bams_to_tags YD.tag $(find ${FOLDER}/k27ac_bams/ -name 'YD_ac*.bam')

    QSUB_ID=$(qsub << ENDINPUT
#!/bin/sh
#PBS -N chipdiff_k27ac
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${WORK_DIR}/chipdiff_k27ac_3.log
# This is necessary because qsub default working dir is user home
cd ${WORK_DIR}

sort -k1,1 -k2,2n -o YD_tags_sorted.tag YD_tags.tag
sort -k1,1 -k2,2n -o OD_tags_sorted.tag OD_tags.tag

ChIPDiff YD_tags_sorted.tag OD_tags_sorted.tag $CHROM_SIZES config.txt diff_YD_OD_3
ENDINPUT
)
    wait_complete "$QSUB_ID"
    check_logs
fi


bams_to_reads() {
    OUT=$1
    # Shift arguments
    shift 1
    for F in $@; do
        >&2 echo $F
        bedtools bamtobed -i ${F} | grep -E "chr[0-9]+|chrX" | awk '{print $1, $2, $3, $6}' >> $OUT
    done
}

macs2_shift() {
    echo $(cat $1 | grep "# d =" | sed 's/.*# d = //g')
}

# MANorm
MANORM="${FOLDER}/k27ac_diff_manorm"
if [ ! -d $MANORM ]; then
    mkdir $MANORM
    cd $MANORM
    WORK_DIR=$MANORM
    Q=0.01
# README.txt
# Create a folder and place in the folder MAnorm.sh, MAnorm.r, and all 4 bed files to be analyzed.
# run command:   ./MAnorm.sh    sample1_peakfile[BED]     sample2_peakfile[BED] \
#                               sample1_readfile[BED]     sample2_readfile[BED]  \
#                               sample1_readshift_lentgh[INT]      sample2_readshift_length[INT]

    cp ~/MAnorm_Linux_R_Package/MAnorm.* ${WORK_DIR}
    cp ${DIFF_MACS_POOLED}/YD_${Q}_peaks.broadPeak ${WORK_DIR}/YD_peaks.bed
    cp ${DIFF_MACS_POOLED}/OD_${Q}_peaks.broadPeak ${WORK_DIR}/OD_peaks.bed

    >&2 echo "Processing OD Pooled Reads";
    bams_to_reads OD_reads.bed $(find ${FOLDER}/k27ac_bams/ -name 'OD_ac*.bam')
    >&2 echo "Processing YD Pooled Reads";
    bams_to_reads YD_reads.bed $(find ${FOLDER}/k27ac_bams/ -name 'YD_ac*.bam')

    # Check MACS2 for shift values
    SHIFT_OD=$(macs2_shift ${DIFF_MACS_POOLED}/OD_${Q}_peaks.xls)
    echo "SHIFT OD: $SHIFT_OD"
    SHIFT_YD=$(macs2_shift ${DIFF_MACS_POOLED}/YD_${Q}_peaks.xls)
    echo "SHIFT YD: $SHIFT_YD"

    QSUB_ID=$(qsub << ENDINPUT
#!/bin/sh
#PBS -N manorm_k27ac
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${WORK_DIR}/manorm_k27ac_3.log
# This is necessary because qsub default working dir is user home
cd ${WORK_DIR}

sort -k1,1 -k2,2n -o YD_reads_sorted.bed YD_reads.bed
sort -k1,1 -k2,2n -o OD_reads_sorted.bed OD_reads.bed

# Load required R module
module load R
bash ${WORK_DIR}/MAnorm.sh YD_peaks.bed OD_peaks.bed YD_reads_sorted.bed OD_reads_sorted.bed $SHIFT_YD $SHIFT_OD
ENDINPUT
)
    wait_complete "$QSUB_ID"
    check_logs
fi