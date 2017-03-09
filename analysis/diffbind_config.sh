#!/usr/bin/env bash
# Script to prepare configuration csv for diffbind
# author Oleg Shpynov (oleg.shpynov@jetbrains.com)

>&2 echo "diffbind: $@"
if [ $# -lt 2 ]; then
    echo "Need 2 parameters! <WORKING_FOLDER> <NAME>"
    echo "Example: scratch/artyomov_lab_aging/Y10OD10/chipseq/processed k27ac_10vs10"
    echo "Necessary folders: \${NAME}_bams, \${NAME}_bams_macs_broad_0.01"
    exit 1
fi
WORK_DIR=$1
NAME=$2

cd $WORK_DIR
READS_DIR=${NAME}_bams
PEAKS_DIR=${NAME}_bams_macs_broad_0.01
RESULT=${NAME}_diffbind.csv
T=$'\t'

# Don't use -printf "%P\n", because it doesn't work on MacOS
PEAKS_FILES=$(find $PEAKS_DIR -name "*.xls" | grep -v input | sort)
echo "SampleID${T}Tissue${T}Factor${T}Condition${T}Replicate${T}bamReads${T}ControlID${T}bamControl${T}Peaks${T}PeakCaller" > $RESULT
for P in $PEAKS_FILES; do
    F=${P##*1/};
    SAMPLE=${F%%_R*};
    CONDITION=${SAMPLE%%D*};
    REPLICATE=${SAMPLE##*D};
    READ=$(ls $READS_DIR/${SAMPLE}*.bam);
    CONTROL=$(ls $READS_DIR/${CONDITION}*input*.bam);
    PEAK=$(ls $PEAKS_DIR/${SAMPLE}*.xls);
    echo "$SAMPLE${T}CD14${T}Age${T}$CONDITION${T}$REPLICATE${T}$(pwd)/$READ${T}${CONDITION}_pooled${T}$(pwd)/$CONTROL${T}$(pwd)/${PEAK}${T}macs" >> $RESULT;
done

echo "DONE: ${RESULT}"