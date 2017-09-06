#!/usr/bin/env bash
# This script is used to filter peaks by given FDR from given MACS2 peaks folder
# author Oleg Shpynov (oleg.shpynov@jetbrains.com)

which bedtools &>/dev/null || { echo "bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>"; exit 1; }

# Load technical stuff
source $(dirname $0)/../parallel/util.sh
SCRIPT_DIR="$(project_root_dir)"

>&2 echo "peaks_frip.sh: $@"

if [ $# -lt 2 ]; then
    echo "Need 2 parameters! <PEAKS_FOLDER> <READS_FOLDER>"
    exit 1
fi

PEAKS_FOLDER=$1
READS_FOLDER=$2

echo "Compute FRIPs for READS_FOLDER: $READS_FOLDER; PEAKS_FOLDER: $PEAKS_FOLDER"
cd ${PEAKS_FOLDER}
for F in $(ls *.*Peak | grep -v gapped); do
    NAME=${F%%_broad*} # filter if broad peaks
    NAME=${NAME%%_q0.*}   # filter suffix if narrow peaks
	BAM=${READS_FOLDER}/${NAME}*.bam
	bash ${SCRIPT_DIR}/reports/rip.sh ${BAM} ${F}
done
python ${SCRIPT_DIR}/reports/peaks_logs.py ${PEAKS_FOLDER}
