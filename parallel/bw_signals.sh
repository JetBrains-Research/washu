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
source $(dirname $0)/../parallel/util.sh
SCRIPT_DIR="$(project_root_dir)"

>&2 echo "Batch bw_signal $@"
if [ $# -lt 4 ]; then
    echo "Need 4 parameters! <WORK_DIR_WITH_BWS> <REGIONS> <ID> <CHROM.SIZES>"
    exit 1
fi

WORK_DIR=$1
REGIONS=$2
ID=$3
CHROM_SIZES=$4

echo "WORK_DIR: $WORK_DIR"
echo "REGIONS: $REGIONS"
echo "ID: $REGIONS"
echo "CHROM_SIZES: $CHROM_SIZES"

RESULTS_FOLDER=${WORK_DIR}/${ID}
mkdir -p ${RESULTS_FOLDER}
echo "RESULTS FOLDER: $RESULTS_FOLDER"

export TMPDIR=$(type job_tmp_dir &>/dev/null && echo "$(job_tmp_dir)" || echo "/tmp")
mkdir -p "${TMPDIR}"

echo "Process libraries sizes ${RESULTS_FOLDER}/sizes.tsv"
cd ${WORK_DIR}
if [[ ! -f ${RESULTS_FOLDER}/sizes.tsv ]]; then
    cat ${CHROM_SIZES} | awk -v OFS='\t' '{print $1,1,$2,$1$2}' > ${RESULTS_FOLDER}/chrom.sizes.bed
    for FILE in $(find . -name '*.bw' | sed 's#\./##g' | sort)
    do :
        NAME=${FILE%%.bw}
        bigWigAverageOverBed ${FILE} ${RESULTS_FOLDER}/chrom.sizes.bed ${RESULTS_FOLDER}/${NAME}.size.tab
        SIZE=$(cat ${RESULTS_FOLDER}/${NAME}.size.tab | awk 'BEGIN{S=0} {S+=$4} END{print(S)}')
        echo "${NAME}"$'\t'"${SIZE}" >> ${RESULTS_FOLDER}/sizes.tsv
    done
    rm ${RESULTS_FOLDER}/chrom.sizes.bed
    rm ${RESULTS_FOLDER}/*size.tab
fi

REGIONS4=${RESULTS_FOLDER}/${ID}.bed4
echo "Create BED4 regions file ${REGIONS4}"
cat ${REGIONS} | awk '{printf("%s\t%s\t%s\t%s_%s_%s\n",$1,$2,$3,$1,$2,$3)}' |\
    sort -k1,1 -k3,3n -k2,2n -T ${TMPDIR} > ${REGIONS4}

echo "Batch bw processing"
TASKS=""
for FILE in $(find . -name '*.bw' | sed 's#\./##g' | sort)
do :
    NAME=${FILE%%.bw}
    TSV=${RESULTS_FOLDER}/${NAME}.tsv
    # Submit task
    run_parallel << SCRIPT
#!/bin/sh
#PBS -N bw_signals_${NAME}
#PBS -l nodes=1:ppn=1,walltime=4:00:00,vmem=8gb
#PBS -j oe
#PBS -o ${RESULTS_FOLDER}/${NAME}_bw_signals.log

cd ${WORK_DIR}
bigWigAverageOverBed ${FILE} ${REGIONS4} ${TSV}.tmp
cat ${TSV}.tmp | tr '_' '\t' | awk -v NAME=${NAME} -v OFS='\t' '{print \$1,\$2,\$3,\$6,NAME}' > ${TSV}
rm ${TSV}.tmp
SCRIPT
        echo "FILE: ${FILE}; TASK: ${QSUB_ID}"
        TASKS="$TASKS $QSUB_ID"
done

wait_complete ${TASKS}
check_logs


echo "Merge all the tsv files into  ${ID}.tsv"
cd $RESULTS_FOLDER
if [[ ! -f ${ID}.tsv ]]; then
    for FILE in $(ls *.tsv | grep -v ${ID}.tsv | grep -v sizes.tsv); do
        cat ${FILE} >> ${ID}.tsv
        rm ${FILE}
    done
fi

echo "Processing coverage, rpm, rpkm, and rpm_peaks for ${ID}.tsv"
run_parallel << SCRIPT
#!/bin/sh
#PBS -N peaks_signal_${ID}
#PBS -l nodes=1:ppn=1,walltime=4:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${ID}_peaks_signal.log

cd $RESULTS_FOLDER
PY_MAJOR_VERS=\$(python -c 'import sys; print(sys.version_info[0])')
if [[ \$PY_MAJOR_VERS != "3" ]]
then
    source activate py3.5
fi

python ${SCRIPT_DIR}/scripts/peaks_signals.py ${RESULTS_FOLDER}/${ID}.tsv ${RESULTS_FOLDER}/sizes.tsv $ID

SCRIPT
wait_complete $QSUB_ID
check_logs

# TMP dir cleanup:
type clean_job_tmp_dir &>/dev/null && clean_job_tmp_dir

# Cleanup
rm ${REGIONS4}

>&2 echo "Done. Batch bw_signal $@"