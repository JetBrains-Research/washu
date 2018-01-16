#!/usr/bin/env bash
# author zayats1812@mail.ru

# Load technical stuff
source $(dirname $0)/../parallel/util/util.sh
SCRIPT_DIR="$(project_root_dir)"

>&2 echo "Batch fragments $@"
if [ $# -lt 1 ]; then
    echo "Need at least 1 parameter! <WORK_DIR> [<WORK_DIR>]*"
    exit 1
fi

WORK_DIRS="$@"

TASKS=()
for WORK_DIR in ${WORK_DIRS}; do :
    cd ${WORK_DIR}

    for FILE in $(find . -type f -name '*.bam' | sed 's#\./##g')
    do :
        WORK_DIR_NAME=${WORK_DIR##*/}
        NAME=${FILE%%.bam} # file name without extension

        # Submit task
        run_parallel << SCRIPT
#!/bin/sh
#PBS -N fragments_${WORK_DIR_NAME}_${NAME}
#PBS -l nodes=1:ppn=1,walltime=2:00:00,vmem=8gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${NAME}_fragments.log

# Loading modules
module load samtools

cd ${WORK_DIR}
samtools view -f66 $FILE | cut -f 9 | sed 's/^-//' > ${NAME}_metrics.txt

module unload samtools # unload samtools, because it conflicts with R at the moment
module load R
Rscript ${SCRIPT_DIR}/parallel/util/fragments.R ${NAME}_metrics.txt ${NAME}_fragments.png
SCRIPT
        echo "FILE: ${WORK_DIR_NAME}:${FILE}; TASK: ${QSUB_ID}"
        TASKS+=("$QSUB_ID")
    done
done

wait_complete ${TASKS[@]}
check_logs

>&2 echo "Done. Batch fragments $@"
