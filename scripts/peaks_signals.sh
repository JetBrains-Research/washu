#!/usr/bin/env bash
# Script to create _pileup.bed files for BAM alignment to compute reads coverage.
# author oleg.shpynov@jetbrains.com

# Load technical stuff, not available in qsub emulation
if [ -f "$(dirname $0)/util.sh" ]; then
    source "$(dirname $0)/util.sh"
fi

if [ $# -lt 4 ]; then
    echo "Need 4 parameters! <WORK_DIR_WITH_BAMS> <FRAGMENT> <REGIONS> <ID>"
    exit 1
fi
echo "Batch peaks_signals: $@"

WORK_DIR=$1
FRAGMENT=$2
REGIONS=$3
ID=$4

echo "WORK_DIR: $WORK_DIR"
echo "FRAGMENT: $FRAGMENT"
echo "REGIONS: $REGIONS"
echo "ID: $ID"

TAGS_FOLDER=${WORK_DIR}/tags_${FRAGMENT}
COVERAGES_FOLDER=${WORK_DIR}/coverages/${ID}
mkdir -p ${TAGS_FOLDER}
mkdir -p ${COVERAGES_FOLDER}

echo "TAGS SHIFTED TO FRAGMENT/2: $TAGS_FOLDER"
echo "RESULTS FOLDER: $COVERAGES_FOLDER"

PROCESSED=""
TASKS=""

REGIONS3=${COVERAGES_FOLDER}/regions.bed3
if [[ ! -f ${REGIONS3} ]]; then
    echo "Create BED3 regions file"
    cat ${REGIONS} | awk -v OFS='\t' '{print($1,$2,$3)}' | sort -k1,1 -k3,3n -k2,2n > ${REGIONS3}
fi

cd ${WORK_DIR}
for FILE in $(find . -name '*.bam' | sed 's#./##g' | sort)
do :
    NAME=${FILE%%.bam}
    TAGS=${TAGS_FOLDER}/${NAME}.tag
    COVERAGE_CSV=${COVERAGES_FOLDER}/${NAME}.csv
    if [[ ! -f ${COVERAGE_CSV} ]]; then
        # Submit task
        QSUB_ID=$(qsub << ENDINPUT
#!/bin/sh
#PBS -N peaks_coverage_${NAME}
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=8gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${NAME}_peaks_coverage.log

cd ${WORK_DIR}
module load bedtools2

if [[ ! -f ${TAGS} ]]; then
    bedtools bamtobed -i ${FILE} |
        awk -v OFS='\t' -v F=${FRAGMENT} '{if (\$6=="-") {print(\$1, \$2+F/2, \$2+F/2+1)} else {if (\$3-F/2>=1) {print(\$1, \$3-F/2, \$3-F/2+1)}}}' |
        sort -k1,1 -k3,3n -k2,2n > ${TAGS}
fi

bedtools intersect -wa -c -a ${REGIONS3} -b ${TAGS} -sorted |\
    awk -v OFS=',' -v NAME=${NAME} '{ print(\$1,\$2,\$3,\$4,NAME) }' > ${COVERAGE_CSV}

ENDINPUT
)
        echo "FILE: ${FILE}; JOB: ${QSUB_ID}"
        TASKS="$TASKS $QSUB_ID"
    fi
done

wait_complete ${TASKS}
check_logs

# Merge all the coverages files into a single file for further python processing
cd $COVERAGES_FOLDER
cat $(ls *.csv | grep -v ${ID})  > ${ID}_coverage.csv
# Cleanup
rm $(ls *.csv | grep -v ${ID})

# Process libraries sizes
cd ${TAGS_FOLDER}
if [[ ! -f sizes.csv ]]; then
    for FILE in $(find . -name '*.tag' | sed 's#./##g' | sort)
    do :
        NAME=${FILE%%.tag}
        LINES=$(cat $FILE | wc -l)
        echo "$NAME,$LINES" >> sizes.csv
    done
fi

# Create resulting tables
QSUB_ID=$(qsub << ENDINPUT
#!/bin/sh
#PBS -N peaks_signal_${ID}
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${ID}_peaks_signal.log

cd $COVERAGES_FOLDER
source activate py3.5
python $(dirname $0)/peaks_signals.py ${COVERAGES_FOLDER}/${ID}_coverage.csv ${TAGS_FOLDER}/sizes.csv $ID

ENDINPUT
)
echo "CREATE TABLES JOB: ${QSUB_ID}"
wait_complete ${QSUB_ID}
check_logs

echo "Done. Batch peaks_signals: $@"