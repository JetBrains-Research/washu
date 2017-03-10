#!/usr/bin/env bash
# Script to process differential chip-seq analysis
# author Oleg Shpynov (oleg.shpynov@jetbrains.com)

# Check tools
which bedtools &>/dev/null || { echo "bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>"; exit 1; }
which macs2 &>/dev/null || { echo "macs2 not found! Install macs2: <https://github.com/taoliu/MACS/wiki/Install-macs2>"; exit 1; }

# Load cluster stuff
source ~/work/washu/scripts/util.sh

################################################################################
# Configuration start ##########################################################
################################################################################

>&2 echo "chipseq_diff: $@"
if [ $# -lt 4 ]; then
    echo "Need 4 parameters! <NAME> <CHROM_SIZES> <GENES_GTF> <DIFFBIND_CSV>"
    exit 1
fi

NAME=$1
CHROM_SIZES=$2
GENES_GTF=$3
DIFFBIND_CSV=$4

FOLDER=$(pwd)
echo "FOLDER"
echo $FOLDER

PREFIX="$FOLDER/$NAME"
echo "PREFIX"
echo $PREFIX

# Base Q value threshold for all the experiments
Q=0.01
echo "Q"
echo $Q

READS_Y=$(awk '{ if ($4 == "Y") printf("%s ", $6); }' $DIFFBIND_CSV)
echo "READS Y"
echo $READS_Y
READS_O=$(awk '{ if ($4 == "O") printf("%s ", $6); }' $DIFFBIND_CSV)
echo "READS O"
echo $READS_O

INPUTS_Y=$(awk '{ if ($4 == "Y") printf("%s ", $8); }' $DIFFBIND_CSV)
echo "INPUT_READS Y"
echo $INPUTS_Y
INPUTS_O=$(awk '{ if ($4 == "O") printf("%s ", $8); }' $DIFFBIND_CSV)
echo "INPUT_READS O"
echo $INPUTS_O

PEAKS_Y=$(awk '{ if ($4 == "Y") printf("%s ", $9); }' $DIFFBIND_CSV | sed 's#xls#broadPeak#g')
echo "INDIVIDUAL_PEAKS Y"
echo $PEAKS_Y
PEAKS_O=$(awk '{ if ($4 == "O") printf("%s ", $9); }' $DIFFBIND_CSV | sed 's#xls#broadPeak#g')
echo "INDIVIDUAL_PEAKS O"
echo $PEAKS_O

################################################################################
# Configuration end ############################################################
################################################################################

#DIFF_MACS_POOLED="${PREFIX}_macs_pooled"
#echo
#echo "Processing $DIFF_MACS_POOLED"
#if [ ! -d $DIFF_MACS_POOLED ]; then
#    mkdir ${DIFF_MACS_POOLED}
#    cd ${DIFF_MACS_POOLED}
#    echo "Processing MACS2 pooled peaks and compare them";
#
#    QSUB_ID1=$(qsub << ENDINPUT
##!/bin/sh
##PBS -N ${NAME}_Y_macs2_broad_${Q}
##PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=16gb
##PBS -j oe
##PBS -o ${DIFF_MACS_POOLED}/${NAME}_Y_macs2_broad_${Q}.log
## This is necessary because qsub default working dir is user home
#cd ${DIFF_MACS_POOLED}
#macs2 callpeak -t $READS_Y -c $INPUTS_Y -f BAM -g hs -n Y_${Q} -B --broad --broad-cutoff ${Q}
#ENDINPUT
#)
#
#    QSUB_ID2=$(qsub << ENDINPUT
##!/bin/sh
##PBS -N ${NAME}_O_macs2_broad_${Q}
##PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=16gb
##PBS -j oe
##PBS -o ${DIFF_MACS_POOLED}/${NAME}_O_macs2_broad_${Q}.log
## This is necessary because qsub default working dir is user home
#cd ${DIFF_MACS_POOLED}
#macs2 callpeak -t $READS_O -c $INPUTS_O -f BAM -g hs -n O_${Q} -B --broad --broad-cutoff ${Q}
#ENDINPUT
#)
#    wait_complete "$QSUB_ID1 $QSUB_ID2"
#    check_logs
#    bash ~/work/washu/bed/compare.sh Y_${Q}_peaks.broadPeak O_${Q}_peaks.broadPeak ${NAME}_${Q}
#fi


MACS_POOLED_Y_VS_O="${PREFIX}_macs_pooled_Y_vs_O"
echo
echo "Processing $MACS_POOLED_Y_VS_O"
if [ ! -d $MACS_POOLED_Y_VS_O ]; then
    mkdir ${MACS_POOLED_Y_VS_O}
    cd ${MACS_POOLED_Y_VS_O}

    echo "Processing MACS2 pooled Y vs O as control and vice versa"
    
    QSUB_ID_Y_vs_O=$(qsub << ENDINPUT
#!/bin/sh
#PBS -N ${NAME}_Y_vs_O_macs2_broad
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${MACS_POOLED_Y_VS_O}/${NAME}_Y_vs_O_macs2_broad.log
# This is necessary because qsub default working dir is user home
cd ${MACS_POOLED_Y_VS_O}
macs2 callpeak -t $READS_Y -c $READS_O -f BAM -g hs -n ${NAME}_Y_vs_O_${Q} -B --broad --broad-cutoff ${Q}
ENDINPUT
)
    QSUB_ID_O_vs_Y=$(qsub << ENDINPUT
#!/bin/sh
#PBS -N ${NAME}_O_vs_Y_macs2_broad
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${MACS_POOLED_Y_VS_O}/${NAME}_O_vs_Y_macs2_broad.log
# This is necessary because qsub default working dir is user home
cd ${MACS_POOLED_Y_VS_O}
macs2 callpeak -t $READS_O -c $READS_Y -f BAM -g hs -n ${NAME}_O_vs_Y_${Q} -B --broad --broad-cutoff ${Q}
ENDINPUT
)
    wait_complete "$QSUB_ID_Y_vs_O $QSUB_ID_O_vs_Y"
    check_logs
fi

#macs2_total_tags_control() {
#    echo $(cat $1 | grep "total tags in control" | sed 's/.*total tags in control: //g')
#}
#
#MACS_BDGDIFF="${PREFIX}_macs_bdgdiff"
#echo
#echo "Processing $MACS_BDGDIFF"
#if [ ! -d $MACS_BDGDIFF ]; then
#    mkdir ${MACS_BDGDIFF}
#    cd ${MACS_BDGDIFF}
#
#    echo "Use MACS2 pooled peaks as input for MACS2 bdgdiff"
#
#    CONTROL_Y=$(macs2_total_tags_control ${DIFF_MACS_POOLED}/Y_${Q}_peaks.xls)
#    echo "Control Y: $CONTROL_Y"
#    CONTROL_O=$(macs2_total_tags_control ${DIFF_MACS_POOLED}/O_${Q}_peaks.xls)
#    echo "Control O: $CONTROL_O"
#
#    QSUB_ID=$(qsub << ENDINPUT
##!/bin/sh
##PBS -N ${NAME}_macs2_broad_bdgdiff
##PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=16gb
##PBS -j oe
##PBS -o ${MACS_BDGDIFF}/${NAME}_macs2_broad_bdgdiff.log
## This is necessary because qsub default working dir is user home
#cd ${MACS_BDGDIFF}
#macs2 bdgdiff\
# --t1 ${DIFF_MACS_POOLED}/Y_${Q}_treat_pileup.bdg --c1 ${DIFF_MACS_POOLED}/Y_${Q}_control_lambda.bdg\
# --t2 ${DIFF_MACS_POOLED}/O_${Q}_treat_pileup.bdg --c2 ${DIFF_MACS_POOLED}/O_${Q}_control_lambda.bdg\
#  --d1 ${CONTROL_Y} --d2 ${CONTROL_O} --o-prefix ${NAME}_${Q}
#ENDINPUT
#)
#    wait_complete "$QSUB_ID"
#    check_logs
#fi


bams_to_tags() {
    OUT=$1
    # Shift arguments
    shift 1
    for F in $@; do
        >&2 echo $F
        bedtools bamtobed -i ${F} | grep -E "chr[0-9]+|chrX" | awk '{print $1, $2, $6}' >> $OUT
    done
}

CHIPDIFF="${PREFIX}_chipdiff"
echo
echo "Processing $CHIPDIFF"
if [ ! -d $CHIPDIFF ]; then
    mkdir ${CHIPDIFF}
    cd ${CHIPDIFF}
    echo "Processing chipdiff as on pooled tags (reads)"

    cat >config.txt <<CONFIG
maxIterationNum  500
minP             0.95
maxTrainingSeqNum 10000
minFoldChange    3
minRegionDist    1000
CONFIG
    >&2 echo "Processing Y Tags";
    bams_to_tags Y_tags.tag $READS_Y

    >&2 echo "Processing O Tags";
    bams_to_tags O_tags.tag $READS_O

    QSUB_ID=$(qsub << ENDINPUT
#!/bin/sh
#PBS -N ${NAME}_chipdiff_3
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${CHIPDIFF}/${NAME}_chipdiff_3.log
# This is necessary because qsub default working dir is user home
cd ${CHIPDIFF}

sort -k1,1 -k2,2n -o Y_tags_sorted.tag Y_tags.tag
sort -k1,1 -k2,2n -o O_tags_sorted.tag O_tags.tag

ChIPDiff Y_tags_sorted.tag O_tags_sorted.tag $CHROM_SIZES config.txt ${NAME}_3
cat ${NAME}_3.region | awk -v OFS='\t' '$4=="-" {print $1,$2,$3}' > ${NAME}_3_cond1.bed
cat ${NAME}_3.region | awk -v OFS='\t' '$4=="+" {print $1,$2,$3}' < ${NAME}_3_cond2.bed
ENDINPUT
)
    wait_complete "$QSUB_ID"
    check_logs
fi


#bams_to_reads() {
#    OUT=$1
#    # Shift arguments
#    shift 1
#    for F in $@; do
#        >&2 echo $F
#        bedtools bamtobed -i ${F} | grep -E "chr[0-9]+|chrX" | awk '{print $1, $2, $3, $6}' >> $OUT
#    done
#}
#
#macs2_shift() {
#    echo $(cat $1 | grep "# d =" | sed 's/.*# d = //g')
#}
#
## MANorm
#MANORM="${PREFIX}_manorm"
#echo
#echo "Processing $MANORM"
#if [ ! -d $MANORM ]; then
#    mkdir ${MANORM}
#    mkdir ${MANORM}/${Q}
#    cd ${MANORM}/${Q}
#
#    echo "Processing MAnorm using pooled MACS2 peaks as peakfile and pooled reads as readfiles"
## README.txt
## Create a folder and place in the folder MAnorm.sh, MAnorm.r, and all 4 bed files to be analyzed.
## run command:   ./MAnorm.sh    sample1_peakfile[BED]     sample2_peakfile[BED] \
##                               sample1_readfile[BED]     sample2_readfile[BED]  \
##                               sample1_readshift_lentgh[INT]      sample2_readshift_length[INT]
#    MANORM_SH=$(which MAnorm.sh)
#    echo "Found MAnorm.sh: ${MANORM_SH}"
#    cp ${MANORM_SH} ${MANORM_SH%%.sh}.r ${MANORM}/${Q}
#
#    cp ${DIFF_MACS_POOLED}/Y_${Q}_peaks.broadPeak ${MANORM}/${Q}/Y_peaks.bed
#    cp ${DIFF_MACS_POOLED}/O_${Q}_peaks.broadPeak ${MANORM}/${Q}/O_peaks.bed
#
#    >&2 echo "Processing Y Pooled Reads";
#    bams_to_reads Y_reads.bed $READS_Y
#    >&2 echo "Processing O Pooled Reads";
#    bams_to_reads O_reads.bed $READS_O
#
#    # Check MACS2 for shift values
#    SHIFT_Y=$(macs2_shift ${DIFF_MACS_POOLED}/Y_${Q}_peaks.xls)
#    echo "SHIFT Y: $SHIFT_Y"
#    SHIFT_O=$(macs2_shift ${DIFF_MACS_POOLED}/O_${Q}_peaks.xls)
#    echo "SHIFT O: $SHIFT_O"
#
#    QSUB_ID=$(qsub << ENDINPUT
##!/bin/sh
##PBS -N manorm_k27ac
##PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=16gb
##PBS -j oe
##PBS -o ${MANORM}/${Q}/manorm_k27ac_3.log
## This is necessary because qsub default working dir is user home
#cd ${MANORM}/${Q}
#
#sort -k1,1 -k2,2n -o Y_reads_sorted.bed Y_reads.bed
#sort -k1,1 -k2,2n -o O_reads_sorted.bed O_reads.bed
#
#sort -k1,1 -k2,2n -o Y_peaks_sorted.bed Y_peaks.bed
#sort -k1,1 -k2,2n -o O_peaks_sorted.bed O_peaks.bed
#
## Load required modules
#module load R
#module load bedtools2
#
#bash ${MANORM}/${Q}/MAnorm.sh Y_peaks_sorted.bed O_peaks_sorted.bed \
#Y_reads_sorted.bed O_reads_sorted.bed $SHIFT_Y $SHIFT_O
#ENDINPUT
#)
#    wait_complete "$QSUB_ID"
#    check_logs
#fi


DIFFBIND="${PREFIX}_diffbind"
echo
echo "Processing $DIFFBIND"
if [ ! -d $DIFFBIND ]; then
    mkdir ${DIFFBIND}
    cd ${DIFFBIND}
    echo "Processing diffbind"
    CSV_NAME=${DIFFBIND_CSV%%.csv}
    
    QSUB_ID=$(qsub << ENDINPUT
#!/bin/sh
#PBS -N ${NAME}_diffbind
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${NAME}_diffbind.log
# This is necessary because qsub default working dir is user home
cd ${DIFFBIND}
module load R
Rscript ~/work/washu/analysis/diffbind.R ${DIFFBIND_CSV}
ENDINPUT
)
    wait_complete "$QSUB_ID"
    check_logs

    # Filter out old and young donors and sort by Q-Value
    cat ${CSV_NAME}_result.csv | awk 'NR > 1 {print $0}' | sed 's#"##g' | tr ',' '\t' |\
        awk -v OFS='\t' '$9 > 0 {print $0}' | sort -k1,1 -k2,2n > ${CSV_NAME}_cond1.bed
    cat ${CSV_NAME}_result.csv | awk 'NR > 1 {print $0}' | sed 's#"##g' | tr ',' '\t' |\
        awk -v OFS='\t' '$9 < 0 {print $0}' | sort -k1,1 -k2,2n > ${CSV_NAME}_cond2.bed

    # Save diffbind results to simple BED3 format
    awk -v OFS='\t' '{ print $1,$2,$3}' ${CSV_NAME}_cond1.bed > ${CSV_NAME}_cond1.bed3
    awk -v OFS='\t' '{ print $1,$2,$3}' ${CSV_NAME}_cond2.bed > ${CSV_NAME}_cond2.bed3

    ~/work/washu/bed/closest_gene.sh ${GENES_GTF} \
        ${CSV_NAME}_cond1.bed3 > ${CSV_NAME}_cond1_closest_genes.tsv
    ~/work/washu/bed/closest_gene.sh ${GENES_GTF} \
        ${CSV_NAME}_cond2.bed3 > ${CSV_NAME}_cond2_closest_genes.tsv
fi