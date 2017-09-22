#!/usr/bin/env bash
# Script to create _pileup.bed files for BAM alignment to compute reads coverage.
# author oleg.shpynov@jetbrains.com

# Load technical stuff
source $(dirname $0)/../parallel/util.sh
SCRIPT_DIR="$(project_root_dir)"

>&2 echo "Batch peaks_signals $@"
if [ $# -lt 4 ]; then
    echo "Need 4 parameters! <WORK_DIR_WITH_BAMS> <INSERT_LENGTH> <REGIONS> <ID>"
    exit 1
fi

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

echo "TAGS PROCESSED: $TAGS_FOLDER"
echo "INSERT SIZE FOR TAGS USED: $INSERT_LENGTH"
echo "RESULTS FOLDER: $COVERAGES_FOLDER"

PROCESSED=""
TASKS=""

TMPDIR=$(type job_tmp_dir &>/dev/null && echo "$(job_tmp_dir)" || echo "/tmp")
mkdir -p "${TMPDIR}"

REGIONS3=${COVERAGES_FOLDER}/${ID}.bed3
echo "Create BED3 regions file ${REGIONS3}"
cat ${REGIONS} | awk -v OFS='\t' '{print($1,$2,$3)}' | sort -k1,1 -k3,3n -k2,2n -T ${TMPDIR} > ${REGIONS3}

cd ${WORK_DIR}
for FILE in $(find . -name '*.bam' | sed 's#\./##g' | sort)
do :
    NAME=${FILE%%.bam}
    TAGS=${TAGS_FOLDER}/${NAME}.tag
    COVERAGE_TSV=${COVERAGES_FOLDER}/${NAME}.tsv
    if [[ ! -f ${COVERAGE_TSV} ]]; then
        # Submit task
        run_parallel << SCRIPT
#!/bin/sh
#PBS -N peaks_coverage_${NAME}
#PBS -l nodes=1:ppn=1,walltime=4:00:00,vmem=8gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${NAME}_peaks_coverage.log

cd ${WORK_DIR}
module load bedtools2

if [[ ! -f ${TAGS} ]]; then
    bash ${SCRIPT_DIR}/scripts/bam2tags.sh ${FILE} $INSERT_LENGTH > ${TAGS}
fi

if [[ ! -f ${COVERAGE_TSV} ]]; then
    bedtools intersect -wa -c -a ${REGIONS3} -b ${TAGS} -sorted > ${COVERAGE_TSV}
fi

SCRIPT
        echo "FILE: ${FILE}; TASK: ${QSUB_ID}"
        TASKS="$TASKS $QSUB_ID"
    fi
done

wait_complete ${TASKS}
check_logs

cd $COVERAGES_FOLDER
if [[ ! -f ${ID}_coverage.tsv ]]; then
    for FILE in $(ls *.tsv | grep -v ${ID}_coverage.tsv); do
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

run_parallel << SCRIPT
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

python ${SCRIPT_DIR}/scripts/peaks_signals.py ${COVERAGES_FOLDER}/${ID}_coverage.tsv ${TAGS_FOLDER}/sizes.tsv $ID

SCRIPT
wait_complete $QSUB_ID
check_logs

# TMP dir cleanup:
type clean_job_tmp_dir &>/dev/null && clean_job_tmp_dir

>&2 echo "Done. Batch peaks_signals $@"