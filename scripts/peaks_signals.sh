#!/usr/bin/env bash
# Script to create _pileup.bed files for BAM alignment to compute reads coverage.
# author oleg.shpynov@jetbrains.com

# Load technical stuff, not available in qsub emulation
if [ -f "$(dirname $0)/util.sh" ]; then
    source "$(dirname $0)/util.sh"
fi

if [ $# -lt 4 ]; then
    echo "Need 4 parameters! <WORK_DIR_WITH_BAMS> <INSERT_LENGTH> <REGIONS> <ID>"
    exit 1
fi
echo "Batch peaks_signals: $@"

WORK_DIR=$1
INSERT_LENGTH=$2
REGIONS=$3
ID=$4

echo "WORK_DIR: $WORK_DIR"
echo "INSERT_LENGTH: $INSERT_LENGTH"
echo "REGIONS: $REGIONS"
echo "ID: $ID"

TAGS_FOLDER=${WORK_DIR}/tags_${INSERT_LENGTH}
COVERAGES_FOLDER=${WORK_DIR}/coverages_${INSERT_LENGTH}/${ID}
mkdir -p ${TAGS_FOLDER}
mkdir -p ${COVERAGES_FOLDER}

echo "TAGS INSERT_LENGTH PROCESSED: $TAGS_FOLDER"
SHIFT=$(($INSERT_LENGTH / 2))
echo "SHIFT FOR TAGS USED: $SHIFT"
echo "RESULTS FOLDER: $COVERAGES_FOLDER"

PROCESSED=""
TASKS=""

REGIONS3=${COVERAGES_FOLDER}/${ID}.bed3
echo "Create BED3 regions file ${REGIONS3}"
cat ${REGIONS} | awk -v OFS='\t' '{print($1,$2,$3)}' | sort -k1,1 -k3,3n -k2,2n > ${REGIONS3}

cd ${WORK_DIR}
for FILE in $(find . -name '*.bam' | sed 's#\./##g' | sort)
do :
    NAME=${FILE%%.bam}
    TAGS=${TAGS_FOLDER}/${NAME}.tag
    COVERAGE_TSV=${COVERAGES_FOLDER}/${NAME}.tsv
    if [[ ! -f ${COVERAGE_TSV} ]]; then
        # Submit task
        QSUB_ID=$(qsub << ENDINPUT
#!/bin/sh
#PBS -N peaks_coverage_${NAME}
#PBS -l nodes=1:ppn=1,walltime=4:00:00,vmem=8gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${NAME}_peaks_coverage.log

cd ${WORK_DIR}
module load bedtools2

if [[ ! -f ${TAGS} ]]; then
    bedtools bamtobed -i ${FILE} |
        awk -v OFS='\t' -v S=${SHIFT} '{if (\$6 != "-") {print(\$1, \$2+S, \$2+S+1)} else {if (\$3-S>=1) {print(\$1, \$3-S, \$3-S+1)}}}' |
        sort -k1,1 -k3,3n -k2,2n > ${TAGS}
fi

if [[ ! -f ${COVERAGE_TSV} ]]; then
    bedtools intersect -wa -c -a ${REGIONS3} -b ${TAGS} -sorted > ${COVERAGE_TSV}
fi

ENDINPUT
)
        echo "FILE: ${FILE}; TASK: ${QSUB_ID}"
        TASKS="$TASKS $QSUB_ID"
    fi
done

wait_complete ${TASKS}
check_logs

cd $COVERAGES_FOLDER
if [[ ! -f ${ID}_coverage.tsv ]]; then
    for FILE in $(ls *.tsv | grep -v ${ID}); do
        NAME=${FILE%%.tsv}
        cat ${FILE} | awk -v OFS='\t' -v NAME=${NAME} '{print $1,$2,$3,$4,NAME}' >> ${ID}_coverage.tsv
    done
fi

# Process libraries sizes
cd ${TAGS_FOLDER}
if [[ ! -f sizes.tsv ]]; then
    for FILE in $(find . -name '*.tag' | sed 's#\./##g' | sort)
    do :
        NAME=${FILE%%.tag}
        LINES=$(cat $FILE | wc -l)
        echo "${NAME}"$'\t'"${LINES}" >> sizes.tsv
    done
fi

QSUB_ID=$(qsub << ENDINPUT
#!/bin/sh
#PBS -N peaks_signal_${ID}
#PBS -l nodes=1:ppn=1,walltime=4:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${ID}_peaks_signal.log

cd $COVERAGES_FOLDER
PY_MAJOR_VERS=\$(python -c 'import sys; print(sys.version_info[0])')
if [[ \$PY_MAJOR_VERS != "3" ]]
then
    source activate py3.5
fi

python $(dirname $0)/peaks_signals.py ${COVERAGES_FOLDER}/${ID}_coverage.tsv ${TAGS_FOLDER}/sizes.tsv $ID

ENDINPUT
)
echo "CREATE TABLES TASK: ${QSUB_ID}"
wait_complete ${QSUB_ID}
check_logs

echo "Done. Batch peaks_signals: $@"