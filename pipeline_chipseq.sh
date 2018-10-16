#!/usr/bin/env bash

# Create module command alias
module() { source /opt/module.sh $@; }
export -f module

source activate py3.5

python $WASHU_ROOT/pipeline_chipseq.py /data/fastq /data/indexes $@