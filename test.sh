#!/usr/bin/env bash

# All the modules loaded
# module load R
# module load bedtools2
# module load bowtie
# module load bowtie2
# module load samtools
# module load sratoolkit
# module load fastqc
# module load java

echo

R=$(which R)
echo "R: $R"
echo "HIDE until module load R"

RSCRIPT=$(which Rscript)
echo "Rscript: $RSCRIPT"
echo "TODO: HIDE until module load R"

BEDTOOLS=$(which bedtools)
echo "bedtools: $BEDTOOLS"
echo "TODO: HIDE until module load bedtools2"

BOWTIE=$(which bowtie)
echo "bowtie: $BOWTIE"
echo "TODO: HIDE until module load bowtie"

BOWTIE2=$(which bowtie2)
echo "bowtie2: $BOWTIE2"
echo "TODO: HIDE until module load bowtie2"

SAMTOOLS=$(which samtools)
echo "samtools: $SAMTOOLS"
echo "TODO: HIDE until module load samtools"

FASTQDUMP=$(which fastq-dump)
echo "fastq-dump: $FASTQDUMP"
echo "TODO: HIDE until module load sratoolkit"

FASTQC=$(which fastqc)
echo "fastqc: $FASTQC"
echo "TODO: HIDE until module load fastqc"

JAVA=$(which java)
echo "java: $JAVA"
echo "TODO: HIDE until module load java"

# Tools required by other scripts
MACS2=$(which macs2)
echo "MACS2: $MACS2"
SICER=$(which SICER.sh)
echo "SICER: $SICER"
RSEG=$(which rseg)
echo "rseg: $RSEG"
BDGTOBW=$(which bedGraphToBigWig)
echo "bedGraphToBigWig: $BDGTOBW"
BEDCLIP=$(which bedClip)
echo "bedClip: $BEDCLIP"
BWAVGOVERBED=$(which bigWigAverageOverBed)
echo "bigWigAverageOverBed: $BWAVGOVERBED"
MULTIQC=$(which multiqc)
echo "multiqc: $MULTIQC"
BAMCOVERAGE=$(which bamCoverage)
echo "bamcoverage: $BAMCOVERAGE"
if [ ! -f ~/picard.jar ]; then
    echo "Picard tools not found! Download Picard: <http://broadinstitute.github.io/picard/>"
fi
if [ ! -f ~/zinbra.jar ]; then
    echo "Zinbra not found! Download ZINBRA: <https://github.com/JetBrains-Research/zinbra>"
fi
echo

ln -s /bin/echo /usr/bin/module

module() { source /opt/module.sh $@; }
export -f module

export IS_TEST=TRUE

# Launch all the tests
python -m pytest test/
python -m pytest --pep8 -m pep8
