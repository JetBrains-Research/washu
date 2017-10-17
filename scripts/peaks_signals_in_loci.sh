#!/usr/bin/env bash
# Counts coverage of signal peaks which belongs to loci of interest.
#
# Arguments:
#   <WORK_DIR_WITH_BAMS> Bam files with signal
#   <INSERT_LENGTH>      Signal fragment size
#   <SIGNAL_PEAKS>       Desired peaks called from signal
#   <LOCI>               Loci of interest
#   <ID>                 Result "id" used to save results
#
#   Result will be saved to <WORK_DIR_WITH_BAMS>/coverages/<ID> folder
#
#
# Example
#   ./peaks_signals_in_loci.sh . 150 YO_macs_broad_0.1.bed cpgIslands.bed CGI
#

# Load technical stuff
source $(dirname $0)/../parallel/util.sh
SCRIPT_DIR="$(project_root_dir)"

if [ $# -lt 5 ]; then
    echo "Need 5 parameters! <WORK_DIR_WITH_BAMS> <INSERT_LENGTH> <SIGNAL_PEAKS> <LOCI> <ID>"
    exit 1
fi

echo "Batch loci_signals: $@"

WORK_DIR=$1
INSERT_LENGTH=$2
SIGNAL_PEAKS=$3
LOCI=$4
ID=$5

echo "WORK_DIR: $WORK_DIR"
echo "INSERT_LENGTH: $INSERT_LENGTH"
echo "SIGNAL_PEAKS: $SIGNAL_PEAKS"
echo "LOCI: $LOCI"
echo "ID: $ID"

COVERAGES_FOLDER=${WORK_DIR}/coverages/${ID}
echo "RESULTS FOLDER: $COVERAGES_FOLDER"
mkdir -p ${COVERAGES_FOLDER}

# Intesect loci with peaks
SIGNAL_PEAKS_FNAME=${SIGNAL_PEAKS##*/}
PEAKS_IN_LOCI="$COVERAGES_FOLDER/${SIGNAL_PEAKS_FNAME%.*}.${ID}.bed"
echo "PEAKS_IN_LOCI: $PEAKS_IN_LOCI"

if [[ ! -f ${PEAKS_IN_LOCI} ]]; then
    "${SCRIPT_DIR}/bed/intersect.sh" ${LOCI} ${SIGNAL_PEAKS} > ${PEAKS_IN_LOCI}
fi

echo "Signal total peaks: $(cat ${SIGNAL_PEAKS} | wc -l)"
echo "Peaks in loci: $(cat ${PEAKS_IN_LOCI} | wc -l)"
echo "Loci: $(cat ${LOCI} | wc -l)"

# Get coverage:
"${SCRIPT_DIR}/scripts/peaks_signals.sh" ${WORK_DIR} ${INSERT_LENGTH} ${PEAKS_IN_LOCI} ${ID}


echo "Done. Batch peaks_signals_in_loci: $@"
