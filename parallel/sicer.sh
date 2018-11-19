#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

which SICER.sh &>/dev/null || {
    echo "SICER not found! Download SICER: <http://home.gwu.edu/~wpeng/Software.htm>"
    echo "Please refer to README for installation instructions, modify scripts, i.e."
    echo "sed -i 's#/home/data/SICER1.1#<YOUR_INSTALLATION_FOLDER>#g' SICER.sh"
    echo "sed -i 's#/home/data/SICER1.1#<YOUR_INSTALLATION_FOLDER>#g' SICER-rb.sh"
    echo "SICER is python2 library, force it!"
    echo "sed -i 's#python#python2#g' SICER.sh"
    echo "sed -i 's#python#python2#g' SICER-rb.sh"
    exit 1
}

# Check configuration
[[ ! -z ${WASHU_ROOT} ]] || { echo "ERROR: WASHU_ROOT not configured"; exit 1; }
source ${WASHU_ROOT}/parallel/util/util.sh

>&2 echo "Batch sicer $@"
if [[ $# -lt 4 ]]; then
    echo "Need at least 4 parameters! <work_dir> <genome> <chrom.sizes> <FDR> [window size (bp)] [fragment size] [gap size (bp)]"
    exit 1
fi

WORK_DIR=$1
GENOME=$2
CHROM_SIZES=$3
FDR=$4

WINDOW_SIZE=200
if [[ $# -ge 5 ]]; then WINDOW_SIZE=$5 ; fi

FRAGMENT_SIZE=150
if [[ $# -ge 6 ]]; then FRAGMENT_SIZE=$6 ; fi

GAP_SIZE=600
if [[ $# -ge 7 ]]; then GAP_SIZE=$7 ; fi

EFFECTIVE_GENOME_FRACTION=$(python ${WASHU_ROOT}/scripts/util.py effective_genome_fraction ${GENOME} ${CHROM_SIZES})
echo "EFFECTIVE_GENOME_FRACTION: ${EFFECTIVE_GENOME_FRACTION}"

if [[ -z "${EFFECTIVE_GENOME_FRACTION}" ]]; then
    echo "EFFECTIVE_GENOME_FRACTION is not determined"
    exit 1
fi

export TMPDIR=$(type job_tmp_dir &>/dev/null && echo "$(job_tmp_dir)" || echo "/tmp")
mkdir -p ${TMPDIR}

cd ${WORK_DIR}

TASKS=()
for FILE in $(find . -name '*.bam' | sed 's#\./##g' | grep -v 'input')
do :
    NAME=${FILE%%.bam} # file name without extension
    ID=${NAME}_F${FRAGMENT_SIZE}_W${WINDOW_SIZE}_G${GAP_SIZE}_FDR${FDR}

    # Check if file already processed
    # Naming example: OD_OD10_H3K27me3-W200-G0-FDR0.01-island.bed
    ISLAND_BED="${NAME}-W${WINDOW_SIZE}-G${GAP_SIZE}-FDR${FDR}-island.bed"
    if [[ ! -f "${ISLAND_BED}" ]]; then
        FILE_BED=${NAME}.bed # It is used for results naming
        INPUT=$(python ${WASHU_ROOT}/scripts/util.py find_input ${WORK_DIR}/${FILE})
        echo "${FILE} input: ${INPUT}"
        INPUT_BED=${INPUT/.bam/.bed}

        # Submit task
        run_parallel << SCRIPT
#!/bin/sh
#PBS -N sicer_${ID}
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${ID}_sicer.log

source ${WASHU_ROOT}/parallel/util/util.sh

export TMPDIR=\$(type job_tmp_dir &>/dev/null && echo "\$(job_tmp_dir)" || echo "/tmp")
SICER_FOLDER=\${TMPDIR}/${ID}
SICER_OUT_FOLDER=\${SICER_FOLDER}/out
# Create folders
mkdir -p \${SICER_FOLDER}
mkdir -p \${SICER_OUT_FOLDER}

# This is necessary because qsub default working dir is user home
cd ${WORK_DIR}
module load bedtools2

# SICER works with BED only, reuse _pileup.bed if possible
PILEUP_BED=\$(pileup ${WORK_DIR}/${FILE})
ln -s \${PILEUP_BED} \${SICER_FOLDER}/${FILE_BED}

if [ -f "${INPUT}" ]; then
    echo "${FILE}: control file found: ${INPUT}"
    INPUT_PILEUP_BED=\$(pileup ${WORK_DIR}/${INPUT})
    ln -sf \${INPUT_PILEUP_BED} \${SICER_FOLDER}/${INPUT_BED}
fi

cd \${SICER_FOLDER}
# Usage:
# SICER.sh    ["InputDir"] ["bed file"] ["control file"] ["OutputDir"] ["Species"] ["redundancy threshold"] \
#   ["window size (bp)"] ["fragment size"] ["effective genome fraction"] ["gap size (bp)"] ["FDR"]
# SICER-rb.sh ["InputDir"] ["bed file"]                  ["OutputDir"] ["species"] ["redundancy threshold"] \
#   ["window size (bp)"] ["fragment size"] ["effective genome fraction"] ["gap size (bp)"] ["E-value"]
#
# Defaults:
#   redundancy threshold    = 1
#   window size (bp)        = 200
#   fragment size           = 150
#   gap size (bp)           = 600

if [ -f \${SICER_FOLDER}/${INPUT_BED} ]; then
    SICER.sh    \${SICER_FOLDER} ${FILE_BED} ${INPUT_BED} \${SICER_OUT_FOLDER} ${GENOME} 1 ${WINDOW_SIZE} \
        ${FRAGMENT_SIZE} ${EFFECTIVE_GENOME_FRACTION} ${GAP_SIZE} ${FDR}
else
    SICER-rb.sh \${SICER_FOLDER} ${FILE_BED}              \${SICER_OUT_FOLDER} ${GENOME} 1 ${WINDOW_SIZE} \
        ${FRAGMENT_SIZE} ${EFFECTIVE_GENOME_FRACTION} ${GAP_SIZE} ${FDR}
fi

# SICER generates lots of output, ignore it: resulting BED only.
# See https://github.com/JetBrains-Research/washu/issues/27

# Cleanup everything else
rm -r \${SICER_FOLDER}
SCRIPT

        echo "FILE: ${FILE}; TASK: ${QSUB_ID}"
        TASKS+=("$QSUB_ID")
    fi
done
wait_complete ${TASKS[@]}
check_logs

type clean_job_tmp_dir &>/dev/null && clean_job_tmp_dir
>&2 echo "Done. Batch sicer $@"