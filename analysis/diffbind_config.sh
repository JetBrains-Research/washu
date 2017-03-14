#!/usr/bin/env bash
# Script to prepare configuration csv for diffbind
# author Oleg Shpynov (oleg.shpynov@jetbrains.com)

>&2 echo "diffbind_config: $@"
if [ $# -lt 2 ]; then
    echo "Need 2 parameters! <WORKING_FOLDER> <NAME> <Q>"
    echo "Example: scratch/artyomov_lab_aging/Y10OD10/chipseq/processed k27ac_10vs10"
    echo "Necessary folders: \${NAME}_bams, \${NAME}_bams_macs_broad_\$Q"
    exit 1
fi
WORK_DIR=$1
NAME=$2
Q=$3

cd $WORK_DIR
READS_DIR=${NAME}_bams
PEAKS_DIR=${NAME}_bams_macs_broad_${Q}

# Start with reads, so that we can move outliers to separate folders and process only valid data
READS_FILES=$(find $READS_DIR -name "*.bam" | grep -v input | sort) # Don't use -printf "%P\n" - doesn't work on MacOS
echo "SampleID,Tissue,Factor,Condition,Replicate,bamReads,ControlID,bamControl,Peaks,PeakCaller"
for R in $READS_FILES; do
    F=${R##*1/}
    SAMPLE=${F%%_R*}
    CONDITION=${SAMPLE%%D*}
    REPLICATE=${SAMPLE##*D}
    READ=$(ls $READS_DIR/${SAMPLE}*.bam)
    CONTROL=$(ls $READS_DIR/${CONDITION}*input*.bam)
    PEAK=$(ls $PEAKS_DIR/${SAMPLE}*.xls)
    echo "$SAMPLE,CD14,Age,$CONDITION,$REPLICATE,$WORK_DIR/$READ,${CONDITION}_pooled,$WORK_DIR/$CONTROL,$WORK_DIR/${PEAK},macs"
done