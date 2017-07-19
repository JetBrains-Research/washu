#!/usr/bin/env bash
# Script to prepare configuration csv for differential peak calling in diffbind format
# author Oleg Shpynov (oleg.shpynov@jetbrains.com)

>&2 echo "diff_config: $@"
if [ $# -lt 2 ]; then
    echo "Need 2 parameters! <READS_DIR> <PEAKS_DIR>"
    echo "Example: "
    echo "  bash diffbind_config.sh k4me3_20vs20_bams k4me3_20vs20_bams_macs2_broad_0.1"
    exit 1
fi

READS_DIR=$1
>&2 echo "READS_DIR: $READS_DIR"
if [[ ! -d ${READS_DIR} ]]; then
    echo "Missing folder ${READS_DIR}"
    exit 1
fi

PEAKS_DIR=$2
>&2 echo "PEAKS_DIR: $PEAKS_DIR"
if [[ ! -d ${PEAKS_DIR} ]]; then
    echo "Missing folder ${PEAKS_DIR}"
    exit 1
fi

# http://stackoverflow.com/questions/1055671/how-can-i-get-the-behavior-of-gnus-readlink-f-on-a-mac
readlink_(){
    TARGET_FILE=$1

    cd `dirname $TARGET_FILE`
    TARGET_FILE=`basename $TARGET_FILE`

    # Iterate down a (possible) chain of symlinks
    while [ -L "$TARGET_FILE" ]
    do
        TARGET_FILE=`readlink $TARGET_FILE`
        cd `dirname $TARGET_FILE`
        TARGET_FILE=`basename $TARGET_FILE`
    done

    # Compute the canonicalized name by finding the physical path
    # for the directory we're in and appending the target file.
    PHYS_DIR=`pwd -P`
    RESULT=$PHYS_DIR/$TARGET_FILE
    echo $RESULT
}

# To absolute paths
READS_DIR=$(readlink_ $READS_DIR)
PEAKS_DIR=$(readlink_ $PEAKS_DIR)

# Start with reads, so that we can move outliers to separate folders and process only valid data
READS_FILES=$(find $READS_DIR -name "*.bam" | grep -v input | sort)
echo "SampleID,Tissue,Factor,Condition,Replicate,bamReads,ControlID,bamControl,Peaks,PeakCaller"
for R in $READS_FILES; do
    >&2 echo "READ: $R"
    FNAME=${R##*/}
    # Assume name pattern .D[0-9]+_.*
    SAMPLE=${FNAME%%_*.bam}
    >&2 echo "SAMPLE: $SAMPLE"
    CONDITION=${SAMPLE%%D*}
    >&2 echo "CONDITION: $CONDITION"
    REPLICATE=${SAMPLE##*D}
    >&2 echo "REPLICATE: $REPLICATE"
    READ=$(ls $READS_DIR/${SAMPLE}*.bam)
    >&2 echo "READ: $READ"
    INPUT=$(python $(dirname $0)/../scripts/util.py find_input ${R})
    >&2 echo "INPUT: $INPUT"
    PEAK=$(ls $PEAKS_DIR/${SAMPLE}*.xls)
    >&2 echo "PEAK: $PEAK"
    echo "$SAMPLE,CD14,Age,$CONDITION,$REPLICATE,$READ,${CONDITION}_input,$INPUT,${PEAK},macs"
done