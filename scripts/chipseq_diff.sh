#!/usr/bin/env bash
# Script to process differential chip-seq analysis
# author Oleg Shpynov (oleg.shpynov@jetbrains.com)
#
# Example of ChIP-Seq differential analysis:
# cd /scratch/artyomov_lab_aging/Y20O20/chipseq/processed/k4me3
# mkdir -p k4me3_diff
# bash ~/work/washu/analysis/diff_config.sh k4me3_20vs20_bams k4me3_20vs20_bams_macs2_broad_0.1 > k4me3_diff/k4me3.csv
# cd k4me3_diff
# bash ~/work/washu/analysis/chipseq_diff.sh k4me3 /scratch/artyomov_lab_aging/Y20O20/chipseq/indexes/hg19/hg19.chrom.sizes k4me3.csv

# Check tools
which bedtools &>/dev/null || { echo "bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>"; exit 1; }
which macs2 &>/dev/null || { echo "macs2 not found! Install macs2: <https://github.com/taoliu/MACS/wiki/Install-macs2>"; exit 1; }

# Load technical stuff
source $(dirname $0)/../parallel/util.sh
SCRIPT_DIR="$(project_root_dir)"

################################################################################
# Configuration start ##########################################################
################################################################################

>&2 echo "chipseq_diff: $@"
if [ $# -lt 3 ]; then
    echo "Need 3 parameters! <NAME> <CHROM_SIZES> <DIFFBIND_CSV>"
    exit 1
fi

NAME=$1
CHROM_SIZES=$2
DIFFBIND_CSV=$3

TMP_DIR=~/tmp
mkdir -p "${TMP_DIR}"

FOLDER=$(pwd)
echo "FOLDER"
echo $FOLDER

PREFIX="$FOLDER/$NAME"
echo "PREFIX"
echo $PREFIX

Q=0.01
echo "Q $Q"
BROAD_CUTOFF=0.1
echo "BROAD_CUTOFF $BROAD_CUTOFF"

READS_Y=$(awk -v FS=',' '{ if ($4 == "Y") print $6 }' $DIFFBIND_CSV | sort  -T ${TMP_DIR} --unique | tr '\n' ' ')
echo "READS Y"
echo "$READS_Y"
READS_O=$(awk -v FS=',' '{ if ($4 == "O") print $6 }' $DIFFBIND_CSV | sort -T ${TMP_DIR} --unique | tr '\n' ' ')
echo "READS O"
echo "$READS_O"

INPUTS_Y=$(awk -v FS=',' '{ if ($4 == "Y") print $8 }' $DIFFBIND_CSV | sort -T ${TMP_DIR} --unique | tr '\n' ' ')
echo "INPUT_READS Y"
echo "$INPUTS_Y"
INPUTS_O=$(awk -v FS=',' '{ if ($4 == "O") print $8 }' $DIFFBIND_CSV | sort -T ${TMP_DIR} --unique | tr '\n' ' ')
echo "INPUT_READS O"
echo "$INPUTS_O"

PEAKS_Y=$(awk -v FS=',' '{ if ($4 == "Y") print $9 }' $DIFFBIND_CSV | sed 's#xls#broadPeak#g' | sort -T ${TMP_DIR} --unique | tr '\n' ' ')
echo "PEAKS Y"
echo "$PEAKS_Y"
PEAKS_O=$(awk -v FS=',' '{ if ($4 == "O") print $9 }' $DIFFBIND_CSV | sed 's#xls#broadPeak#g' | sort -T ${TMP_DIR} --unique | tr '\n' ' ')
echo "PEAKS O"
echo "$PEAKS_O"

################################################################################
# Configuration end ############################################################
################################################################################

DIFFBIND="${PREFIX}_diffbind"
echo
echo "Processing $DIFFBIND"
if [ ! -d $DIFFBIND ]; then
    mkdir -p ${DIFFBIND}
    cp ${DIFFBIND_CSV} ${DIFFBIND}/${NAME}.csv
    cd ${DIFFBIND}

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
Rscript ${SCRIPT_DIR}/R/diffbind.R ${NAME}.csv
SCRIPT
    wait_complete "$QSUB_ID"

    # Filter out old and young donors and sort by Q-Value
    cat ${NAME}_result.csv | awk 'NR > 1 {print $0}' | sed 's#"##g' | tr ',' '\t' |\
        awk -v OFS='\t' '$9 > 0 {print $0}' | sort -T ${TMP_DIR} -k10,10g > ${NAME}_cond1.bed
    cat ${NAME}_result.csv | awk 'NR > 1 {print $0}' | sed 's#"##g' | tr ',' '\t' |\
        awk -v OFS='\t' '$9 < 0 {print $0}' | sort -T ${TMP_DIR} -k10,10g > ${NAME}_cond2.bed
    # Save ${NAME} results to simple BED3 format
    awk -v OFS='\t' '{print $1,$2,$3}' ${NAME}_cond1.bed > ${NAME}_cond1.bed3
    awk -v OFS='\t' '{print $1,$2,$3}' ${NAME}_cond2.bed > ${NAME}_cond2.bed3
fi

MACS_POOLED_Y_VS_O="${PREFIX}_macs_pooled_Y_vs_O"
echo
echo "Processing $MACS_POOLED_Y_VS_O"
if [ ! -d $MACS_POOLED_Y_VS_O ]; then
    mkdir -p ${MACS_POOLED_Y_VS_O}
    cd ${MACS_POOLED_Y_VS_O}

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
fi


DIFF_MACS_POOLED="${PREFIX}_macs_pooled"
echo
echo "Processing $DIFF_MACS_POOLED"
if [ ! -d $DIFF_MACS_POOLED ]; then
    mkdir -p ${DIFF_MACS_POOLED}
    cd ${DIFF_MACS_POOLED}
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
    bash ${SCRIPT_DIR}/bed/compare.sh Y_${BROAD_CUTOFF}_peaks.broadPeak O_${BROAD_CUTOFF}_peaks.broadPeak ${NAME}_${BROAD_CUTOFF}
fi

macs2_total_tags_control() {
    echo $(cat $1 | grep "total tags in control" | sed 's/.*total tags in control: //g')
}

MACS_BDGDIFF="${PREFIX}_macs_bdgdiff"
echo
echo "Processing $MACS_BDGDIFF"
if [ ! -d $MACS_BDGDIFF ]; then
    mkdir -p ${MACS_BDGDIFF}
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
fi


# Check MACS2 for shift values
macs2_shift() {
    echo $(cat $1 | grep "# d =" | sed 's/.*# d = //g')
}

CHIPDIFF="${PREFIX}_chipdiff"
echo
echo "Processing $CHIPDIFF"
if [ ! -d $CHIPDIFF ]; then
    mkdir -p ${CHIPDIFF}
    cd ${CHIPDIFF}
    echo "Processing chipdiff as on pooled tags (reads)"

    cat >config.txt <<CONFIG
maxIterationNum  500
minP             0.95
maxTrainingSeqNum 10000
minFoldChange    3
minRegionDist    1000
CONFIG
    >&2 echo "Processing Y Tags"
    SHIFT_Y=$(macs2_shift ${DIFF_MACS_POOLED}/Y_${BROAD_CUTOFF}_peaks.xls)
    echo "SHIFT Y: $SHIFT_Y"
    run_parallel << SCRIPT
#!/bin/sh
#PBS -N Y_bam2tags
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=8gb
#PBS -j oe
#PBS -o ${CHIPDIFF}/Y_bams2tags.log
# This is necessary because qsub default working dir is user home
cd ${CHIPDIFF}
module load bedtools2
for F in ${READS_Y}; do
    >&2 echo \$F
    bash ${SCRIPT_DIR}/scripts/bam2tags.sh \$F $SHIFT_Y >> Y_tags.tag
done
SCRIPT
    QSUB_ID1=$QSUB_ID

    >&2 echo "Processing O Tags"
    SHIFT_O=$(macs2_shift ${DIFF_MACS_POOLED}/O_${BROAD_CUTOFF}_peaks.xls)
    echo "SHIFT O: $SHIFT_O"
    run_parallel << SCRIPT
#!/bin/sh
#PBS -N O_bam2tags
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=8gb
#PBS -j oe
#PBS -o ${CHIPDIFF}/O_bams2tags.log
# This is necessary because qsub default working dir is user home
cd ${CHIPDIFF}
module load bedtools2
for F in ${READS_O}; do
    >&2 echo \$F
    bash ${SCRIPT_DIR}/scripts/bam2tags.sh \$F $SHIFT_O >> O_tags.tag
done
SCRIPT
    QSUB_ID2=$QSUB_ID

    wait_complete "$QSUB_ID1 $QSUB_ID2"

    run_parallel << SCRIPT
#!/bin/sh
#PBS -N ${NAME}_chipdiff_3
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${CHIPDIFF}/${NAME}_chipdiff_3.log
# This is necessary because qsub default working dir is user home
cd ${CHIPDIFF}

# Inplace sort
sort -T ${TMP_DIR} -k1,1 -k2,2n -o Y_tags.tag Y_tags.tag
sort -T ${TMP_DIR} -k1,1 -k2,2n -o O_tags.tag O_tags.tag

ChIPDiff Y_tags.tag O_tags.tag $CHROM_SIZES config.txt ${NAME}_3
cat ${NAME}_3.region | awk -v OFS='\t' '\$4=="-" {print \$1,\$2,\$3}' > ${NAME}_3_cond1.bed
cat ${NAME}_3.region | awk -v OFS='\t' '\$4=="+" {print \$1,\$2,\$3}' < ${NAME}_3_cond2.bed
SCRIPT
    wait_complete $QSUB_ID
fi


# MANorm
MANORM="${PREFIX}_manorm"
echo
echo "Processing $MANORM"
if [ ! -d $MANORM ]; then
    mkdir -p ${MANORM}

    echo "Processing MAnorm using pooled MACS2 peaks as peakfile and pooled reads as readfiles"
# README.txt
# Create a folder and place in the folder MAnorm.sh, MAnorm.r, and all 4 bed files to be analyzed.
# run command:   ./MAnorm.sh    sample1_peakfile[BED]     sample2_peakfile[BED] \
#                               sample1_readfile[BED]     sample2_readfile[BED]  \
#                               sample1_readshift_lentgh[INT]      sample2_readshift_length[INT]
    MANORM_SH=$(which MAnorm.sh)
    echo "Found MAnorm.sh: ${MANORM_SH}"
    cp ${MANORM_SH} ${MANORM_SH%%.sh}.r ${MANORM}

    cp ${DIFF_MACS_POOLED}/Y_${BROAD_CUTOFF}_peaks.broadPeak ${MANORM}/Y_peaks.bed
    cp ${DIFF_MACS_POOLED}/O_${BROAD_CUTOFF}_peaks.broadPeak ${MANORM}/O_peaks.bed

    >&2 echo "Processing Y Pooled Reads"
    run_parallel << SCRIPT
#!/bin/sh
#PBS -N Y_bam2reads
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=8gb
#PBS -j oe
#PBS -o ${MANORM}/Y_bams2reads.log
# This is necessary because qsub default working dir is user home
cd ${MANORM}

module load bedtools2
for F in ${READS_Y}; do
    >&2 echo \$F
    bash ${SCRIPT_DIR}/scripts/bam2tags.sh \$F >> Y_reads.tag
done
SCRIPT
    QSUB_ID1=$QSUB_ID

    >&2 echo "Processing O Pooled Reads"
    run_parallel << SCRIPT
#!/bin/sh
#PBS -N O_bam2reads
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=8gb
#PBS -j oe
#PBS -o ${MANORM}/O_bams2reads.log
# This is necessary because qsub default working dir is user home
cd ${MANORM}

module load bedtools2
for F in ${READS_O}; do
    >&2 echo \$F
    bash ${SCRIPT_DIR}/scripts/bam2tags.sh \$F >> O_reads.tag
done
SCRIPT
    QSUB_ID2=$QSUB_ID

    wait_complete "$QSUB_ID1 $QSUB_ID2"

    run_parallel << SCRIPT
#!/bin/sh
#PBS -N manorm_${Q}
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${MANORM}/manorm_${Q}.log
# This is necessary because qsub default working dir is user home
cd ${MANORM}

# Sort inplace
sort -T ${TMP_DIR} -k1,1 -k2,2n -o Y_reads.bed Y_reads.bed
sort -T ${TMP_DIR} -k1,1 -k2,2n -o Y_peaks.bed Y_peaks.bed

sort -T ${TMP_DIR} -k1,1 -k2,2n -o O_reads.bed O_reads.bed
sort -T ${TMP_DIR} -k1,1 -k2,2n -o O_peaks.bed O_peaks.bed

# Load required modules
module load R
module load bedtools2

bash ${MANORM}/MAnorm.sh Y_peaks.bed O_peaks.bed Y_reads.bed O_reads.bed $SHIFT_Y $SHIFT_O
SCRIPT
    wait_complete $QSUB_ID
fi