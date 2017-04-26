#!/usr/bin/env bash
# Script to prepare configuration csv for diffbind
# author Oleg Shpynov (oleg.shpynov@jetbrains.com)

>&2 echo "diffbind_config: $@"
if [ $# -lt 2 ]; then
    echo "Need 2 parameters! <READS_DIR> <PEAKS_DIR>"
    echo "Example: "
    echo "  bash diffbind_config.sh scratch/artyomov_lab_aging/Y10OD10/chipseq/processed/k4me1_10vs10_reseq_bams\\"
    echo "    scratch/artyomov_lab_aging/Y10OD10/chipseq/processed/k4me1_10vs10_reseq_bams_macs_broad_0.1"
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
READS_DIR=$(readlink -f $READS_DIR)
PEAKS_DIR=$(readlink -f $PEAKS_DIR)

# Start with reads, so that we can move outliers to separate folders and process only valid data
READS_FILES=$(find $READS_DIR -name "*.bam" | grep -v input | sort)
echo "SampleID,Tissue,Factor,Condition,Replicate,bamReads,ControlID,bamControl,Peaks,PeakCaller"
for R in $READS_FILES; do
    >&2 echo "READ: $R"
    FNAME=${R##*/}
    # Should be changed for particular naming scheme
    SAMPLE=${FNAME%%_K9me3.bam}
    >&2 echo "SAMPLE: $SAMPLE"
    CONDITION=${SAMPLE%%D*}
    >&2 echo "CONDITION: $CONDITION"
    REPLICATE=${SAMPLE##*D}
    >&2 echo "REPLICATE: $REPLICATE"
    READ=$(ls $READS_DIR/${SAMPLE}*.bam)
    >&2 echo "READ: $READ"
    CONTROL=$(ls $READS_DIR/${CONDITION}*input*.bam)
    >&2 echo "CONTROL: $CONTROL"
    PEAK=$(ls $PEAKS_DIR/${SAMPLE}*.xls)
    >&2 echo "PEAK: $PEAK"
    echo "$SAMPLE,CD14,Age,$CONDITION,$REPLICATE,$READ,${CONDITION}_pooled,$CONTROL,${PEAK},macs"
done