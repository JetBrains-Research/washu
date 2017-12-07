#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

which bedtools &>/dev/null || { echo "bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>"; exit 1; }
>&2 echo "bam_qc: $@"

# Load technical stuff
source $(dirname $0)/../parallel/util.sh
SCRIPT_DIR="$(project_root_dir)"

>&2 echo "Batch bam_qc $@"
if [ $# -lt 2 ]; then
    echo "Need 2 parameters! <phantompeakqualtools> <work_dir> [<work_dir>]*"
    exit 1
fi

PHANTOMPEAKQUALTOOLS=$1
WORK_DIRS=${@:2}

TASKS=()
for WORK_DIR in ${WORK_DIRS}; do :
    echo "${WORK_DIR}"
    WORK_DIR_NAME=${WORK_DIR##*/}
    cd ${WORK_DIR}

    for FILE in $(find . -name '*.bam' | sed 's#\./##g')
    do :
        NAME=${FILE%%.bam} # file name without extension

        # Submit task
        run_parallel << SCRIPT
#!/bin/sh
#PBS -N bam_qc_${WORK_DIR_NAME}_${NAME}
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=8gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${NAME}_bam_qc.log

source "${SCRIPT_DIR}/parallel/util.sh"
export TMPDIR=\$(type job_tmp_dir &>/dev/null && echo "\$(job_tmp_dir)" || echo "/tmp")

# This is necessary because qsub default working dir is user home
cd ${WORK_DIR}

module load bedtools2
module load R
module load samtools

#col.	abbreviation	description
#1	Filename	tagAlign/BAM filename
#2	numReads	effective sequencing depth i.e. total number of mapped reads in input file
#3	estFragLen	comma separated strand cross-correlation peak(s) in decreasing order of correlation.
#4	corr_estFragLen	comma separated strand cross-correlation value(s) in decreasing order (COL2 follows the same order)
#5	phantomPeak	Read length/phantom peak strand shift
#6	corr_phantomPeak	Correlation value at phantom peak
#7	argmin_corr	strand shift at which cross-correlation is lowest
#8	min_corr	minimum value of cross-correlation
#9	NSC	Normalized strand cross-correlation coefficient (NSC) = COL4 / COL8
#10	RSC	Relative strand cross-correlation coefficient (RSC) = (COL4 - COL8) / (COL6 - COL8)
#11	QualityTag	Quality tag based on thresholded RSC (codes= -2:veryLow, -1:Low, 0:Medium, 1:High, 2:veryHigh)
Rscript ${PHANTOMPEAKQUALTOOLS}/run_spp.R -c=${FILE} -savp -out=${WORK_DIR}/${NAME}.phantom.txt

#col.	abbreviation
#1   TotalReadPairs
#2   DistinctReadPairs
#3   OneReadPair
#4   TwoReadPairs
#5   NRF=Distinct/Total
#6   PBC1=OnePair/Distinct
#7   PBC2=OnePair/TwoPair
bash ${SCRIPT_DIR}/reports/pbc_nrf.sh ${FILE} > ${WORK_DIR}/${NAME}.pbc_nfr.txt
SCRIPT
        echo "FILE: ${WORK_DIR_NAME}/${FILE}; TASK: ${QSUB_ID}"
        TASKS+=("$QSUB_ID")
    done
done
wait_complete ${TASKS[@]}
check_logs


>&2 echo "Done. Batch bams2_qc $@"