#!/usr/bin/env bash

# Mock of module command for tests
>&2 echo "MOCK module $@"
if [ $# -ne 2 ]; then
    echo "Need 2 parameters! module load [module_name]"
elif [ "load" != $1 ]; then
    echo "Only load command supported"
elif [ "samtools" == $2 ]; then
    which samtools &>/dev/null  || export PATH=$PATH:/opt/conda/envs/samtools/bin/
elif [ "sratoolkit" == $2 ]; then
    which fastq-dump &>/dev/null || export PATH=$PATH:/opt/conda/envs/sratoolkit/bin/
elif [ "bedtools2" == $2 ]; then
    which bedtools &>/dev/null || export PATH=$PATH:/opt/conda/envs/bedtools/bin/
elif [ "R" == $2 ]; then
    which R &>/dev/null || export PATH=$PATH:/opt/conda/envs/r/bin/
elif [ "bowtie" == $2 ]; then
    which bowtie &>/dev/null || export PATH=$PATH:/opt/conda/envs/bowtie/bin/
elif [ "bowtie2" == $2 ]; then
    which bowtie2 &>/dev/null || export PATH=$PATH:/opt/conda/envs/bowtie2/bin/
elif [ "fastqc" == $2 ]; then
    which java &>/dev/null || {
        export PATH=$PATH:/opt/fastqc/:/opt/conda/envs/java/bin/ ;
        export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/conda/envs/java/lib/:/opt/conda/envs/java/jre/lib/amd64/
    }
elif [ "java" == $2 ]; then
    which java &>/dev/null || {
        export PATH=$PATH:/opt/conda/envs/java/bin/
        export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/conda/envs/java/lib/:/opt/conda/envs/java/jre/lib/amd64/
    }
else
    echo "ERROR: Unsupported module $2"
    exit 1
fi

