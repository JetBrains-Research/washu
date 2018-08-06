#!/usr/bin/env bash
# Script to compute signal for BAM files at given regions
# author oleg.shpynov@jetbrains.com

which bigWigAverageOverBed &>/dev/null || {
    echo "bigWigAverageOverBed not found!"
    echo "Download: <http://hgdownload.cse.ucsc.edu/admin/exe/> or"
    echo "  conda install -c bioconda ucsc-bigwigaverageoverbed"
    exit 1
}
# Check configuration
[[ ! -z ${WASHU_ROOT} ]] || { echo "ERROR: WASHU_ROOT not configured"; exit 1; }
source ${WASHU_ROOT}/parallel/util/util.sh

>&2 echo "Batch signals $@"
if [ $# -lt 5 ]; then
    echo "Need at least 5 parameters! <WORK_DIR> <TAGS_BW> <FRAGMENT> <FOLDER or REGIONS.BED> <CHROM.SIZES> [<PEAKS_FILE.BED>]"
    exit 1
fi

WORK_DIR=$(expand_path $1)
TAGS_BW=$(expand_path $2)
FRAGMENT=$3
REGIONS_ARG=$(expand_path $4)
CHROM_SIZES=$(expand_path $5)
if [[ -f $6 ]]; then
    PEAKS_FILE_BED=$(expand_path $6)
fi
echo "WORK_DIR: $WORK_DIR"
echo "TAGS_BW: $TAGS_BW"
echo "FRAGMENT: $FRAGMENT"
echo "REGIONS_ARG: $REGIONS_ARG"
echo "CHROM_SIZES: $CHROM_SIZES"
echo "PEAKS_FILE: $PEAKS_FILE_BED"

export TMPDIR=$(type job_tmp_dir &>/dev/null && echo "$(job_tmp_dir)" || echo "/tmp")
mkdir -p $TMPDIR

# Work in dedicated folder now
WORK_DIR_FRAGMENT=${WORK_DIR}/${FRAGMENT}
mkdir -p ${WORK_DIR_FRAGMENT}
cd ${WORK_DIR_FRAGMENT}

tag_file_2_donor_name() {
    FILE=$1
    NAME=${FILE##*/}
    echo ${NAME%%_tags.bw}
}
TAGS_BW_FILES=$(find ${TAGS_BW} -name '*.bw' | sed 's#\./##g' | sort)
########################################################################################################################
# Function to process all the summary coverages for given file
########################################################################################################################
process_coverage()
{
    _BED4=$1
    _ID=$2
    _RESULT=$3
    _LOGS_DIR=$4
    if [[ -f $_RESULT ]]; then
        rm -r $_RESULT
    fi

    NAMES=()
    TASKS=()
    TSVS=()
    LOGS=()
    cd ${WORK_DIR_FRAGMENT}
    for FILE in ${TAGS_BW_FILES}
    do :
        NAME=$(tag_file_2_donor_name ${FILE})
        NAMES+=("$NAME")
        TSV=$TMPDIR/${_ID}_${NAME}.tsv
        TSVS+=("$TSV")
        LOG=${_LOGS_DIR}/${_ID}_${NAME}_tsv.log
        LOGS+=("$LOG")
        # Submit task
        run_parallel << SCRIPT
#!/bin/sh
#PBS -N ${_ID}_${NAME}_tsv
#PBS -l nodes=1:ppn=1,walltime=4:00:00,vmem=8gb
#PBS -j oe
#PBS -o ${LOG}
cd ${WORK_DIR_FRAGMENT}
bigWigAverageOverBed ${FILE} ${_BED4} ${TSV}
SCRIPT
        echo "FILE: ${FILE}; TASK: ${QSUB_ID}"
        TASKS+=("$QSUB_ID")
    done
    wait_complete ${TASKS[@]}

    # Master log
    MASTER_LOG=${_LOGS_DIR}/${_ID}_tsv.log
    for I in $(seq 0 $((${#NAMES[@]} - 1))); do
        # Merge tsv result
        NAME=${NAMES[$I]}
        TSV=${TSVS[$I]}
        # Process regions coverage
        #   sum - sum of values over all bases covered
        #   mean0 - average over bases with non-covered bases counting as zeroes
        #   mean - average over just covered bases
        # Fields in ${TSV}:  $6 $7 $8 - sum, mean0, mean values, after chr#start#end split by #
        cat ${TSV} | tr '#' '\t' | awk -v N=${NAME} -v OFS='\t' '{print $1,$2,$3,$6,$7,$8,N}' >> ${_RESULT}
        rm ${TSV}

        # Merge log
        LOG=${LOGS[$I]}
        echo "$LOG" >> ${MASTER_LOG}
        cat ${LOG} >> ${MASTER_LOG}
        rm ${LOG}
    done
}


########################################################################################################################
# Process libraries
########################################################################################################################
LIBRARIES_SIZES=${WORK_DIR_FRAGMENT}/${CHROM_SIZES##*/}.tsv
echo "Compute FULL ${CHROM_SIZES} libraries size ${LIBRARIES_SIZES}"
if [[ ! -f ${LIBRARIES_SIZES} ]]; then
    # Prepare BED4 region
    cat ${CHROM_SIZES} | awk '{printf("%s\t1\t%s\t%s#1#%s\n",$1,$2,$1,$2)}' > ${TMPDIR}/chrom.sizes.bed4

    process_coverage ${TMPDIR}/chrom.sizes.bed4 "chrom.sizes" ${TMPDIR}/chrom.sizes.tsv ${TMPDIR}

    for FILE in ${TAGS_BW_FILES}
    do :
        NAME=$(tag_file_2_donor_name ${FILE})
        SIZE=$(cat ${TMPDIR}/chrom.sizes.tsv | grep ${NAME} | awk 'BEGIN{S=0} {S+=$4} END{print(S)}')
        echo "${NAME}"$'\t'"${SIZE}" >> ${LIBRARIES_SIZES}
    done
fi

########################################################################################################################
# Processing Peaks files
########################################################################################################################
if [[ -f ${PEAKS_FILE_BED} ]]; then
    LIBRARIES_PEAKS_SIZES=${WORK_DIR_FRAGMENT}/${PEAKS_FILE_BED##*/}.tsv
    echo "Compute peaks ${PEAKS_FILE_BED} size ${LIBRARIES_PEAKS_SIZES}"
    if [[ ! -f ${LIBRARIES_PEAKS_SIZES} ]]; then
        cat ${PEAKS_FILE_BED} | awk '{printf("%s\t%s\t%s\t%s#%s#%s\n",$1,$2,$3,$1,$2,$3)}' |\
            sort -k1,1 -k3,3n -k2,2n --unique -T $TMPDIR > ${TMPDIR}/peaks.sizes.bed4

        process_coverage ${TMPDIR}/peaks.sizes.bed4 "peaks.sizes" ${TMPDIR}/peaks.sizes.tsv ${WORK_DIR_FRAGMENT}

        for FILE in ${TAGS_BW_FILES}
        do :
            NAME=$(tag_file_2_donor_name ${FILE})
            SIZE=$(cat ${TMPDIR}/peaks.sizes.tsv | grep ${NAME} | awk 'BEGIN{S=0} {S+=$4} END{print(S)}')
            echo "${NAME}"$'\t'"${SIZE}" >> ${LIBRARIES_PEAKS_SIZES}
        done
    fi
else
    echo "NO peaks file given"
fi


########################################################################################################################
# Processing BED files
########################################################################################################################
if [[ -d ${REGIONS_ARG} ]]; then
    echo "Folder with BED files provided: ${REGIONS_ARG}"
    REGIONS=$(find ${REGIONS_ARG} -name '*.bed')
    REGIONS_ROOT=${REGIONS_ARG}
    echo "REGIONS_ROOT: ${REGIONS_ROOT}"
else
    echo "Single BED file provided"
    REGIONS=( "${REGIONS_ARG}" )
    REGIONS_ROOT=$(dirname ${REGIONS_ARG})
    echo "REGIONS_ROOT: ${REGIONS_ROOT}"
fi


for REGIONS_BED in ${REGIONS}; do
    echo "REGIONS: ${REGIONS_BED}"
    SUB_PATH=$(dirname ${REGIONS_BED} | sed "s#$REGIONS_ROOT##g")
    echo "SUB_PATH: ${SUB_PATH}"
    ID=${REGIONS_BED%%.bed}
    ID=${ID##*/}
    RESULTS_FOLDER=${WORK_DIR_FRAGMENT}/${SUB_PATH}/${ID}
    echo "RESULTS FOLDER: $RESULTS_FOLDER"
    if [[ ! -d "${RESULTS_FOLDER}" ]]; then
        mkdir -p ${RESULTS_FOLDER}
    fi

    REGIONS_SIZES=${RESULTS_FOLDER}/${ID}.tsv
    if [[ ! -f ${REGIONS_SIZES} ]]; then
        echo "Compute regions ${REGIONS_BED} coverage ${REGIONS_SIZES}"
        SHIFT=$(($FRAGMENT / 2))
        # Here we extend all the regions left and right to FRAGMENT / 2
        # to ensure that we counted all the intersections between fragments and regions.
        cat ${REGIONS_BED} |\
            awk -v S=${SHIFT} '{if ($2-S<1){L=1}else{L=$2-S};R=$3+S;printf("%s\t%s\t%s\t%s#%s#%s\n",$1,L,R,$1,$2,$3)}' |\
         sort -k1,1 -k3,3n -k2,2n --unique -T $TMPDIR > ${TMPDIR}/${ID}.bed4

        process_coverage ${TMPDIR}/${ID}.bed4 ${ID} ${REGIONS_SIZES} ${RESULTS_FOLDER}
    fi

    echo "Processing signals normalization and visualization for ${ID}.tsv"
    run_parallel << SCRIPT
#!/bin/sh
#PBS -N signal_${ID}
#PBS -l nodes=1:ppn=1,walltime=4:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${RESULTS_FOLDER}/${ID}_signal.log
PY_MAJOR_VERS=\$(python -c 'import sys; print(sys.version_info[0])')
if [[ \$PY_MAJOR_VERS != "3" ]]
then
    source activate py35 || source activate py3.5
fi
cd ${RESULTS_FOLDER}
python ${WASHU_ROOT}/downstream/signals/signals_normalize.py ${REGIONS_SIZES} ${LIBRARIES_SIZES} ${LIBRARIES_PEAKS_SIZES}
SCRIPT
    wait_complete $QSUB_ID

    cd ${RESULTS_FOLDER}
    check_logs
done


# Cleanup
type clean_job_tmp_dir &>/dev/null && clean_job_tmp_dir

>&2 echo "Done. Batch signals $@"