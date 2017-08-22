#!/usr/bin/env bash

# Mock of module command for tests

if [ $# -ne 2 ]; then
    echo "Need 2 parameters! module load [module_name]"
elif [ "load" != $1 ]; then
    echo "Only load command supported"
elif [ "samtools" == $2 ]; then
    which samtools &>/dev/null || export PATH=$PATH:/opt/conda/envs/samtools/bin/
elif [ "bedtools2" == $2 ]; then
    which bedtools &>/dev/null || export PATH=$PATH:/opt/conda/envs/bedtools/bin/
elif [ "R" == $2 ]; then
    which R &>/dev/null || export PATH=$PATH:/opt/conda/envs/r/bin/
elif [ "bowtie" == $2 ]; then
    which bowtie &>/dev/null || export PATH=$PATH:/opt/bowtie-1.2.1.1/
else
    echo "ERROR: Unsupported module $2"
fi

