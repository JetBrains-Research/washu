#!/usr/bin/env bash
# author zayats1812@mail.ru
# author oleg.shpynov@jetbrains.com

which STAR &>/dev/null || { echo "STAR not found! Download STAR: <https://github.com/alexdobin/STAR>"; exit 1; }

# Check configuration
[[ ! -z ${WASHU_ROOT} ]] || { echo "ERROR: WASHU_ROOT not configured"; exit 1; }
source ${WASHU_ROOT}/parallel/util.sh

>&2 echo "Batch star $@"
if [ $# -lt 3 ]; then
    echo "Need 3 parameter! <WORK_DIR> <GENOME> <STAR_REF>"
    exit 1
fi
WORK_DIR=$1
GENOME=$2
REF=$3


PROCESSED=()
TASKS=()

cd ${WORK_DIR}
for FILE in $(find . -name '*.f*q' | sed 's#\./##g' | sort)
do :
    if $(echo "${PROCESSED[@]}"  | fgrep -q "${FILE}");
    then
        echo "$FILE was already processed"
        continue
    fi

    # Assumption: the only difference between paired-end read files is _1 and _2 / _R1 and _R2
    FILE_PAIRED=""
    if $(echo "${FILE}"  | fgrep -q "_1"); then
        FILE_PAIRED=$(echo "${FILE}" | sed 's#_1#_2#g')
        NAME=$(echo ${FILE} | sed 's#_1##g' | sed -r 's#^.*/##g' | sed -r 's#\.f.*q$##g')
    elif $(echo "${FILE}"  | fgrep -q "_R1"); then
        FILE_PAIRED=$(echo "${FILE}" | sed 's#_R1#_R2#g')
        NAME=$(echo ${FILE} | sed 's#_R1##g' | sed -r 's#^.*/##g'| sed -r 's#\.f.*q$##g')
    fi
    ID=${NAME}_${GENOME}

    # Check correct paired name
    if [[ -f "${FILE_PAIRED}" ]]; then
        echo "PAIRED END reads detected: $FILE and $FILE_PAIRED"
        # Mark it as already processed
        PROCESSED+=("${FILE_PAIRED}")
    fi

    if [[ -f "${ID}.bam" ]]; then
        echo "   [Skipped] ${WORK_DIR}/${ID}.bam already exists."
        continue
    fi

    # Submit task
    run_parallel << SCRIPT
#!/bin/bash
#PBS -N star_${ID}
#PBS -j oe
#PBS -l mem=64gb,nodes=1:ppn=8,walltime=24:00:00
#PBS -o ${WORK_DIR}/${ID}_star.log

# This is necessary because qsub default working dir is user home
cd ${WORK_DIR}

# Create tmp folder to avoid any conflicts
export TMPDIR=\$(type job_tmp_dir &>/dev/null && echo "\$(job_tmp_dir)" || echo "/tmp")
STAR_FOLDER=\${TMPDIR}/${NAME}
mkdir -p \${STAR_FOLDER}

ln -sf ${WORK_DIR}/${FILE} \${STAR_FOLDER}/${FILE}
if [[ -f "${FILE_PAIRED}" ]]; then
    ln -sf ${WORK_DIR}/${FILE_PAIRED} \${STAR_FOLDER}/${FILE_PAIRED}
fi

cd \${STAR_FOLDER}

# LoadAndKeep option allows to save lots of time on genome loading
# this command also will generate a transcripts BAM file which can be used for RSEM or other similar tools
if [[ -f "${FILE_PAIRED}" ]]; then
    STAR --genomeDir ${REF} --genomeLoad LoadAndKeep \
        --readFilesIn ${FILE} ${FILE_PAIRED} --runThreadN 8 \
        --outFilterMultimapNmax 15 --outFilterMismatchNmax 6 \
        --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 30000000000 \
        --outFileNamePrefix "${ID}_" \
        --quantMode TranscriptomeSAM
else
    STAR --genomeDir ${REF} --genomeLoad LoadAndKeep \
        --readFilesIn ${FILE} --runThreadN 8 \
        --outFilterMultimapNmax 15 --outFilterMismatchNmax 6 \
        --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 30000000000 \
        --outFileNamePrefix "${ID}_" \
        --quantMode TranscriptomeSAM
fi
mv \${STAR_FOLDER}/${ID}_Aligned.sortedByCoord.out.bam ${WORK_DIR}/${ID}.bam
mv \${STAR_FOLDER}/${ID}_Aligned.toTranscriptome.out.bam ${WORK_DIR}/${ID}.tr.bam
mv \${STAR_FOLDER}/${ID}_Log.out ${WORK_DIR}/${ID}.star_run.log
mv \${STAR_FOLDER}/${ID}_Log.final.out ${WORK_DIR}/${ID}.star_final.log

# Cleanup everything else
rm -r \${STAR_FOLDER}
SCRIPT
    echo "FILE: $NAME; TASK: ${QSUB_ID}"
    TASKS+=("$QSUB_ID")
done
wait_complete ${TASKS[@]}

# Purge the genome from RAM
STAR --genomeDir ${REF} --genomeLoad Remove
check_logs

>&2 echo "Done. Batch star $@"