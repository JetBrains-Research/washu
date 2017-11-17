#!/usr/bin/env bash

# Create module command alias
# We need this for "which module" command
ln -s /bin/echo /usr/bin/module
module() { source /opt/module.sh $@; }
export -f module

module load bedtools2
module load R
export IS_TEST=TRUE
export PYTHONPATH="/washu:$PYTHONPATH"

#############
# Fast tests#
#############
python -m pytest test/*.py
python -m pytest --pep8 -m pep8

##################
# Pipeline tests #
##################
echo "Check tools available by module load"
which R &>/dev/null && { echo "R: HIDE until module load R"; exit 1; }
which Rscript &>/dev/null && { echo "Rscript: HIDE until module load R"; exit 1; }
which bedtools &>/dev/null && { echo "bedtools: HIDE until module load bedtools2"; exit 1; }
which bowtie &>/dev/null && { echo "bowtie: HIDE until module load bowtie"; exit 1; }
which bowtie2 &>/dev/null && { echo "bowtie2: HIDE until module load bowtie2"; exit 1; }
which samtools &>/dev/null && { echo "samtools: HIDE until module load samtools"; exit 1; }
which fastq-dump &>/dev/null && { echo "fastq-dump: HIDE until module load sratoolkit"; exit 1; }
which fastqc &>/dev/null && { echo "fastqc: HIDE until module load fastqc"; exit 1; }
which java &>/dev/null && { echo "java: HIDE until module load java"; exit 1; }

echo "Check tools required by other scripts"
which macs2 &>/dev/null || { echo "macs2 not found"; exit 1; }
which SICER.sh &>/dev/null || { echo "SICER.sh not found"; exit 1; }
which rseg &>/dev/null || { echo "rseg not found"; exit 1; }
which bedGraphToBigWig &>/dev/null || { echo "bedGraphToBigWig not found"; exit 1; }
which bedClip &>/dev/null || { echo "bedClip not found"; exit 1; }
which bigWigAverageOverBed &>/dev/null || { echo "bigWigAverageOverBed not found"; exit 1; }
which multiqc &>/dev/null || { echo "multiqc not found"; exit 1; }
which bamCoverage &>/dev/null || { echo "bamCoverage not found"; exit 1; }

echo "Check downloaded files"
if [ ! -f ~/picard.jar ]; then
    echo "Picard tools not found! Download Picard: <http://broadinstitute.github.io/picard/>"
fi
if [ ! -f ~/zinbra.jar ]; then
    echo "Zinbra not found! Download ZINBRA: <https://github.com/JetBrains-Research/zinbra>"
fi
python -m pytest test/pipeline/*.py