#!/usr/bin/env bash

# Create module command alias
# We need this for "which module" command
ln -s /bin/echo /usr/bin/module
module() { source /opt/module.sh $@; }
export -f module

export IS_TEST=TRUE
export PYTHONPATH="/washu:$PYTHONPATH"

##################
# Pipeline tests #
##################
echo "Check tools available by module load"
module load bedtools2
#which bedtools &>/dev/null && { echo "bedtools: HIDE until module load bedtools2"; exit 1; }
which R &>/dev/null && { echo "FAILED R: HIDE until module load R"; exit 1; }
which Rscript &>/dev/null && { echo "FAILED Rscript: HIDE until module load R"; exit 1; }
which bowtie &>/dev/null && { echo "FAILED bowtie: HIDE until module load bowtie"; exit 1; }
which bowtie2 &>/dev/null && { echo "FAILED bowtie2: HIDE until module load bowtie2"; exit 1; }
which samtools &>/dev/null && { echo "FAILED samtools: HIDE until module load samtools"; exit 1; }
which fastq-dump &>/dev/null && { echo "FAILED fastq-dump: HIDE until module load sratoolkit"; exit 1; }
which fastqc &>/dev/null && { echo "FAILED fastqc: HIDE until module load fastqc"; exit 1; }
which java &>/dev/null && { echo "FAILED java: HIDE until module load java"; exit 1; }

echo "Check tools required by other scripts"
which macs2 &>/dev/null || { echo "FAILED: macs2 not found"; exit 1; }
which SICER.sh &>/dev/null || { echo "FAILED: SICER.sh not found"; exit 1; }
which rseg &>/dev/null || { echo "FAILED: rseg not found"; exit 1; }
which bedGraphToBigWig &>/dev/null || { echo "FAILED: bedGraphToBigWig not found"; exit 1; }
which bedClip &>/dev/null || { echo "FAILED: bedClip not found"; exit 1; }
which bigWigAverageOverBed &>/dev/null || { echo "FAILED: bigWigAverageOverBed not found"; exit 1; }
which multiqc &>/dev/null || { echo "FAILED: multiqc not found"; exit 1; }
which bamCoverage &>/dev/null || { echo "FAILED: bamCoverage not found"; exit 1; }

echo "Check downloaded files"
if [ ! -f ~/picard.jar ]; then
    echo "FAILED: Picard tools not found! Download Picard: <http://broadinstitute.github.io/picard/>"
fi
if [ ! -f ~/zinbra.jar ]; then
    echo "FAILED: Zinbra not found! Download ZINBRA: <https://github.com/JetBrains-Research/zinbra>"
fi

# Limit parallelism level for tests
export WASHU_PARALLELISM=2

python -m pytest test/pipeline/*.py

# For the ease of troubleshooting
mkdir -p /washu/out
cp -r /root/fastq* /washu/out