#!/usr/bin/env bash
# Script to create _pileup.bed files for BAM alignment to compute reads coverage.
# author oleg.shpynov@jetbrains.com

# Load technical stuff, not available in qsub emulation
if [ -f "$(dirname $0)/util.sh" ]; then
    source "$(dirname $0)/util.sh"
fi

if [ $# -lt 3 ]; then
    echo "Need 3 parameters! <WORK_DIR_WITH_BAMS> <REGIONS> <ID>"
    exit 1
fi

WORK_DIR=$1
REGIONS=$2
ID=$3

echo "Batch peaks_signals: $@"
cd ${WORK_DIR}

BEDS_FOLDER=${WORK_DIR}/beds
COVERAGES_FOLDER=${WORK_DIR}/coverages/${ID}
mkdir -p ${BEDS_FOLDER}
mkdir -p ${COVERAGES_FOLDER}

PROCESSED=""
TASKS=""

REGIONS3=${REGIONS}.bed3
if [[ ! -f ${REGIONS3} ]]; then
    echo "Create BED regions file"
    cat ${REGIONS} | awk -v OFS='\t' '{print($1,$2,$3)}' > ${REGIONS3}
fi

cd ${WORK_DIR}
for FILE in $(find . -name '*.bam' | sed 's#./##g' | sort)
do :
    NAME=${FILE%%.bam}
    FILE_BED=${BEDS_FOLDER}/${NAME}.bed
    COVERAGE_BED=${COVERAGES_FOLDER}/${NAME}.bed
    if [[ ! -f ${COVERAGE_BED} ]]; then
        # Submit task
        QSUB_ID=$(qsub << ENDINPUT
#!/bin/sh
#PBS -N peaks_coverage_${NAME}
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=8gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${NAME}_peaks_coverage.log

cd ${WORK_DIR}
module load bedtools2

if [[ ! -f ${FILE_BED} ]]; then
    bedtools bamtobed -i ${FILE} |
            grep -v 'chrU' | grep -v 'random' |
            awk '{if (\$6=="-") {print(\$1, \$3-1, \$3)} else {print(\$1, \$2, \$2+1)}}' |
            sort -k1,1 -k3,3n -k2,2n > ${FILE_BED}
fi

bedtools intersect -wa -wb -a ${REGIONS3} -b ${FILE_BED} -sorted |
awk -v OFS='\t' 'BEGIN{c="";s=0;e=0;x=0}\
{   if (\$1!=c||\$2!=s||\$3!=e) {\
        if (x!=0) print(\$1,\$2,\$3,x);c=\$1;s=\$2;e=\$3;x=1\
    } else {x+=1}\
}\
END{print(\$1,\$2,\$3,x)}' > ${COVERAGE_BED}

ENDINPUT
)
        echo "FILE: ${FILE}; JOB: ${QSUB_ID}"
        TASKS="$TASKS $QSUB_ID"
    fi
done

wait_complete ${TASKS}
check_logs

# Process BED sizes
cd ${BEDS_FOLDER}
if [[ ! -f sizes.csv ]]; then
    for FILE in $(find . -name '*.bed' | sed 's#./##g' | sort)
    do :
        NAME=${FILE%%.bed}
        LINES=$(cat $FILE | wc -l)
        echo "$NAME,$LINES" >> sizes.csv
    done
fi

# Merge all the coverages files into a single file for further python processing
cd $COVERAGES_FOLDER
for FILE in $(find . -name '*.bed' | sed 's#./##g' | sort)
do :
    NAME=${FILE%%.bed}
    cat $FILE | awk -v OFS=',' -v NAME=$NAME '{print($1,$2,$3,NAME,$4)}' >> coverage.csv
done

# Finally create tables
python $(dirname $0)/peaks_signals.sh ${COVERAGES_FOLDER}/coverage.csv ${BEDS_FOLDER}/sizes.csv $ID
echo "Done. Batch peaks_signals: $@"