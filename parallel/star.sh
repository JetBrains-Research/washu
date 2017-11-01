#!/usr/bin/env bash
# author zayats1812@mail.ru

# Load technical stuff
source $(dirname $0)/../parallel/util.sh

>&2 echo "Batch star $@"
if [ $# -lt 3 ]; then
    echo "Need 3 parameter! <WORK_DIR> <STAR_REF> <READ_FILES>"
    exit 1
fi
WORK_DIR=$1
REF="$2/human_star_100"
READ_FILES=($3)
cd ${WORK_DIR}

TASKS=()
for (( i=0; i<${#READ_FILES[@]} ; i+=2 )) ;
do :
    NAME=${READ_FILES[i]%%_1.f*q.gz} # file name without extension

    # Submit task
    run_parallel << SCRIPT
#!/bin/bash
#PBS -N star_align_$NAME
#PBS -j oe
#PBS -l mem=64gb,nodes=1:ppn=8:haswell,walltime=24:00:00

cd $WORK_DIR
module load star
## paired-end allignment
STAR --genomeDir $REF --genomeLoad LoadAndKeep \
 --readFilesIn ${READ_FILES[i]} ${READ_FILES[i+1]} --runThreadN 8 --readFilesCommand zcat \
 --outFilterMultimapNmax 15 --outFilterMismatchNmax 6  --outSAMstrandField All \
 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 30000000000 \
 --outFileNamePrefix "./${NAME}_" \
 --quantMode TranscriptomeSAM &> ${NAME}_star_stdout.log

## purge the genome from RAM and remove temporary files
STAR --genomeDir $REF --genomeLoad Remove

mv ${NAME}_Aligned.sortedByCoord.out.bam $NAME.bam
mv ${NAME}_Aligned.toTranscriptome.out.bam $NAME.tr.bam
mv ${NAME}_Log.out $NAME.star_run.log
mv ${NAME}_Log.final.out $NAME.star_final.log
SCRIPT
    echo "FILE: $NAME; TASK: ${QSUB_ID}"
    TASKS+=("$QSUB_ID")
done

wait_complete ${TASKS[@]}

rm Aligned.out.sam
rm *_Log.progress.out
rm Log.out
rm *_SJ.out.tab
rm -rf _STARtmp

>&2 echo "Done. Batch star $@"