#!/usr/bin/env bash
# Script to prepare configuration csv for differential peak calling in diffbind format
# author Oleg Shpynov (oleg.shpynov@jetbrains.com)

# Load technical stuff
source $(dirname $0)/../parallel/util.sh
SCRIPT_DIR="$(project_root_dir)"

>&2 echo "diff_config: $@"
if [ $# -lt 2 ]; then
    echo "Need 2 parameters! <READS_DIR> <PEAKS_DIR>"
    echo "Example params:"
    echo "  k4me3_20vs20_bams k4me3_20vs20_bams_macs2_broad_0.1"
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


# To absolute paths
READS_DIR="$(expand_path "${READS_DIR}")"
PEAKS_DIR="$(expand_path "${PEAKS_DIR}")"

# Start with reads, so that we can move outliers to separate folders and process only valid data
READS_FILES=$(find "${READS_DIR}" -name "*.bam" | grep -v input | sort)
echo "SampleID,Tissue,Factor,Condition,Replicate,bamReads,ControlID,bamControl,Peaks,PeakCaller"
for R in ${READS_FILES}; do
    FNAME=${R##*/}
    # Assume name pattern .D[0-9]+_.*
    SAMPLE=${FNAME%%_*.bam}
    CONDITION=${SAMPLE%%D*}
    REPLICATE=${SAMPLE##*D}
    READ=$(ls ${READS_DIR}/${SAMPLE}*.bam)
    INPUT=$(python ${SCRIPT_DIR}/scripts/util.py find_input ${R})
    PEAK=$(ls ${PEAKS_DIR}/${SAMPLE}*.xls)
    echo "$SAMPLE,CD14,Age,$CONDITION,$REPLICATE,$READ,${CONDITION}_input,$INPUT,${PEAK},macs"
done