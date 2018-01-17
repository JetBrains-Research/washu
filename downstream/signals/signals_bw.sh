#!/usr/bin/env bash
# Script to compute signal for BW files at given regions
# author oleg.shpynov@jetbrains.com

which bigWigAverageOverBed &>/dev/null || {
    echo "bigWigAverageOverBed not found!"
    echo "Download: <http://hgdownload.cse.ucsc.edu/admin/exe/> or"
    echo "  conda install -c bioconda ucsc-bigwigaverageoverbed"
    exit 1
}
# Load technical stuff
source $(dirname $0)/../../parallel/util/util.sh
PROJECT_ROOT=$(project_root_dir)

>&2 echo "Batch bw_signal $@"
if [ $# -lt 3 ]; then
    echo "Need at least 3 parameters! <WORK_DIR_WITH_BWS> <REGIONS.BED> <CHROM.SIZES> [<PEAKS_FILE.BED>]"
    exit 1
fi

WORK_DIR=$1
REGIONS_BED=$2
CHROM_SIZES=$3
PEAKS_FILE_BED=$4

echo "WORK_DIR: $WORK_DIR"
echo "REGIONS: $REGIONS_BED"
echo "CHROM_SIZES: $CHROM_SIZES"
echo "PEAKS_FILE: $PEAKS_FILE_BED"

ID=${REGIONS_BED%%.bed};
ID=${ID##*/};
RESULTS_FOLDER=${WORK_DIR}/${ID}
echo "RESULTS FOLDER: $RESULTS_FOLDER"

# Skip if folder already processed
if [[ -d "${RESULTS_FOLDER}" ]]; then
    echo "FOLDER exists, skip";
    exit 0
else
    mkdir -p ${RESULTS_FOLDER}
fi

export TMPDIR=$(type job_tmp_dir &>/dev/null && echo "$(job_tmp_dir)" || echo "/tmp")
mkdir -p $TMPDIR

cd ${WORK_DIR}

# Function to process all the summary coverages for given file
process_coverage()
{
    _BED4=$1
    _ID=$2
    _RESULT=$3
    _LOGS_DIR=$4
    if [[ -f $_RESULT ]]; then
        rm -r $_RESULT
    fi

    TASKS=()
    TSVS=()
    for FILE in $(find . -name '*.bw' | sed 's#\./##g' | sort)
    do :
        NAME=${FILE%%.bw}
        TSV=$TMPDIR/$NAME.tsv
        TSVS+=("$TSV")
        # Submit task
        run_parallel << SCRIPT
#!/bin/sh
#PBS -N bw_signals_${ID}_${NAME}
#PBS -l nodes=1:ppn=1,walltime=4:00:00,vmem=8gb
#PBS -j oe
#PBS -o ${_LOGS_DIR}/bw_signals_${ID}_${NAME}.log

# Process regions coverage
#   sum - sum of values over all bases covered
#   mean0 - average over bases with non-covered bases counting as zeroes
#   mean - average over just covered bases
# Fields \$6 \$7 \$8 - sum, mean0, mean values, after chr#start#end split by #
bigWigAverageOverBed ${FILE} ${_BED4} ${TSV}.tmp
cat ${TSV}.tmp | tr '#' '\t' | awk -v NAME=${NAME} -v OFS='\t' '{print \$1,\$2,\$3,\$6,\$7,\$8,NAME}' > ${TSV}
rm ${TSV}.tmp
SCRIPT
        echo "FILE: ${FILE}; TASK: ${QSUB_ID}"
        TASKS+=("$QSUB_ID")
    done
    wait_complete ${TASKS[@]}

    for FILE in ${TSVS[@]}; do
        cat ${FILE} >> $_RESULT
        rm ${FILE}
    done
}

LIBRARIES_SIZES=${WORK_DIR}/${CHROM_SIZES##*/}.tsv
echo "Compute libraries size ${LIBRARIES_SIZES}"
if [[ ! -f ${LIBRARIES_SIZES} ]]; then
    # Prepare BED4 region
    cat ${CHROM_SIZES} | awk '{printf("%s\t1\t%s\t%s#1#%s\n",$1,$2,$1,$2)}' > ${TMPDIR}/chrom.sizes.bed4

    process_coverage ${TMPDIR}/chrom.sizes.bed4 "chrom.sizes" ${TMPDIR}/chrom.sizes.tsv ${TMPDIR}

    for FILE in $(find . -name '*.bw' | sed 's#\./##g' | sort)
    do :
        NAME=${FILE%%.bw}
        SIZE=$(cat ${TMPDIR}/chrom.sizes.tsv | grep ${NAME} | awk 'BEGIN{S=0} {S+=$4} END{print(S)}')
        echo "${NAME}"$'\t'"${SIZE}" >> ${LIBRARIES_SIZES}
    done
fi

if [[ -f ${PEAKS_FILE_BED} ]]; then
    LIBRARIES_PEAKS_SIZES=${WORK_DIR}/${PEAKS_FILE_BED##*/}.tsv
    echo "Compute libraries peaks size ${LIBRARIES_PEAKS_SIZES}"
    if [[ ! -f ${LIBRARIES_PEAKS_SIZES} ]]; then
        cat ${PEAKS_FILE_BED} | awk '{printf("%s\t%s\t%s\t%s#%s#%s\n",$1,$2,$3,$1,$2,$3)}' |\
            sort -k1,1 -k3,3n -k2,2n --unique -T $TMPDIR > ${TMPDIR}/peaks.sizes.bed4

        process_coverage ${TMPDIR}/peaks.sizes.bed4 "peaks.sizes" ${TMPDIR}/peaks.sizes.tsv ${WORK_DIR}

        for FILE in $(find . -name '*.bw' | sed 's#\./##g' | sort)
        do :
            NAME=${FILE%%.bw}
            SIZE=$(cat ${TMPDIR}/peaks.sizes.tsv | grep ${NAME} | awk 'BEGIN{S=0} {S+=$4} END{print(S)}')
            echo "${NAME}"$'\t'"${SIZE}" >> ${LIBRARIES_PEAKS_SIZES}
        done
    fi
else
    echo "NO peaks file given"
fi

echo "Compute regions coverage ${REGIONS_BED}"
cat $REGIONS_BED | awk '{printf("%s\t%s\t%s\t%s#%s#%s\n",$1,$2,$3,$1,$2,$3)}' |\
    sort -k1,1 -k3,3n -k2,2n --unique -T $TMPDIR > ${TMPDIR}/regions.bed4

process_coverage ${TMPDIR}/regions.bed4 ${ID} ${RESULTS_FOLDER}/${ID}.tsv ${RESULTS_FOLDER}

echo "Processing data, rpm, rpkm, and rpm_peaks for ${ID}.tsv"
run_parallel << SCRIPT
#!/bin/sh
#PBS -N peaks_signal_${ID}
#PBS -l nodes=1:ppn=1,walltime=4:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${RESULTS_FOLDER}/${ID}_signal.log
PY_MAJOR_VERS=\$(python -c 'import sys; print(sys.version_info[0])')
if [[ \$PY_MAJOR_VERS != "3" ]]
then
    source activate py35 || source activate py3.5
fi
cd $RESULTS_FOLDER
python ${PROJECT_ROOT}/downstream/signals/signals.py ${RESULTS_FOLDER}/${ID}.tsv ${LIBRARIES_SIZES} ${LIBRARIES_PEAKS_SIZES}
SCRIPT
wait_complete $QSUB_ID

# TMP dir cleanup:
type clean_job_tmp_dir &>/dev/null && clean_job_tmp_dir

>&2 echo "Done. Batch bw_signal $@"