#!/usr/bin/env bash

echo "Pipeline script1"
echo "Working directory: `pwd`"

echo "Submitting sra to fastq.gz tasks"
find . -type f -name "*.sra" | xargs -n1 ~/work/washu/scripts/sra2fastq.sh

# See sra2fastq script for tasks naming convention
echo "Collecting tasks: sra2fastq"
SRA_FILES=$(find . -type f -name "*.sra")
SRA_TASKS=""
for FILE in ${SRA_FILES}
do :
    if [[ -z "$SRA_TASKS" ]]; then
        $SRA_TASKS="sra2fastq_$FILE"
  else
      $SRA_TASKS="$SRA_TASKS,sra2fastq_$FILE"
  fi
done

echo "Waiting for tasks $SRA_TASKS to finish"
qsub -hold_jid $SRA_TASKS -cwd ~/work/washu/scripts/makemehappy2.sh

