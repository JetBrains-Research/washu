#!/usr/bin/env bash

# Create module command alias
module() { source /opt/module.sh $@; }
export -f module

export IS_TEST=TRUE
export WASHU_ROOT="/washu"
export PYTHONPATH="$WASHU_ROOT:$PYTHONPATH"

##################
# Pipeline tests #
##################
module load bedtools2
module load samtools
#which bedtools &>/dev/null && { echo "bedtools: HIDE until module load bedtools2"; exit 1; }
#which samtools &>/dev/null && { echo "ERROR samtools: HIDE until module load samtools"; exit 1; }

echo "Check tools available by module load"
which R &>/dev/null && { echo "ERROR R: HIDE until module load R"; exit 1; }
which Rscript &>/dev/null && { echo "ERROR Rscript: HIDE until module load R"; exit 1; }
which bowtie &>/dev/null && { echo "ERROR bowtie: HIDE until module load bowtie"; exit 1; }
which bowtie2 &>/dev/null && { echo "ERROR bowtie2: HIDE until module load bowtie2"; exit 1; }
which fastq-dump &>/dev/null && { echo "ERROR fastq-dump: HIDE until module load sratoolkit"; exit 1; }
which fastqc &>/dev/null && { echo "ERROR fastqc: HIDE until module load fastqc"; exit 1; }
which java &>/dev/null && { echo "ERROR java: HIDE until module load java"; exit 1; }

echo "Check tools required by other scripts"
which macs2 &>/dev/null || { echo "ERROR: macs2 not found"; exit 1; }
which SICER.sh &>/dev/null || { echo "ERROR: SICER.sh not found"; exit 1; }
which rseg &>/dev/null || { echo "ERROR: rseg not found"; exit 1; }
which bedGraphToBigWig &>/dev/null || { echo "ERROR: bedGraphToBigWig not found"; exit 1; }
which bedClip &>/dev/null || { echo "ERROR: bedClip not found"; exit 1; }
which bigWigAverageOverBed &>/dev/null || { echo "ERROR: bigWigAverageOverBed not found"; exit 1; }
which multiqc &>/dev/null || { echo "ERROR: multiqc not found"; exit 1; }
which bamCoverage &>/dev/null || { echo "ERROR: bamCoverage not found"; exit 1; }

echo "Check downloaded files"
if [ ! -f ~/picard.jar ]; then
    echo "ERROR: Picard tools not found! Download Picard: <http://broadinstitute.github.io/picard/>"
fi
if [ ! -f ~/zinbra.jar ]; then
    echo "ERROR: Zinbra not found! Download ZINBRA: <https://github.com/JetBrains-Research/zinbra>"
fi

echo "Prepare data"

if [[ ! -d ~/washu_test_data ]]; then
    tar -xf ~/washu_test_data.tar.gz --directory ~
fi

cp -r ~/washu_test_data/fastq ~/
cp -r ~/washu_test_data/index ~/
cp -r ~/washu_test_data/data ~/

# Prepare regions for RSEG
if [[ ! -d ~/fastq_bams ]]; then
    mkdir -p ~/fastq_bams
    ln -sf ~/index/hg19/deadzones-k36-hg19.bed ~/fastq_bams/deadzones-k36-hg19.bed
fi

# Limit parallelism level for tests
export WASHU_PARALLELISM=2

python -m pytest test/pipeline/*.py

# For the ease of troubleshooting
mkdir -p /washu/out
cp -r ~/fastq* /washu/out
cp -r ~/index* /washu/out
cp -r ~/data* /washu/out
cp -r ~/pileup* /washu/out
cp -r ~/signals* /washu/out