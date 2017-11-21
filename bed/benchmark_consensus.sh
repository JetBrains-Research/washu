#!/usr/bin/env bash

# Load technical stuff
source $(dirname $0)/../parallel/util.sh

>&2 echo "benchmark_consensus $@"
if [ $# -lt 2 ]; then
    echo "Need 2 parameter! <BENCHMARK_ROOT_DIR> <OUTPUT_DIR>"
    exit 1
fi

BENCHMARK_ROOT=$1
LOCI_ROOT=$2

cd ${BENCHMARK_ROOT}
for HIST_DIR in $(find . -type d -name "H*" -maxdepth 1); do
    cd ${HIST_DIR}
    HIST_NAME=${HIST_DIR##*/}
    # We decided to exclude MACS2 narrow peaks
    # Sicer not ready yet, let's exclude it too
    for DIR in $(find . -type d -maxdepth 1 ! -path . ! -path "./macs_narrow" ! -path "./sicer"); do
        /home/user/work/tsurinov/washu/bed/consensus.sh -p 50 ${DIR} ${HIST_NAME}
        # "$(project_root_dir)/bed/consensus.sh" -p 50  ${DIR} ${HIST_NAME}
        cp "${DIR}/*_consensus.bed" ${LOCI_ROOT} ${HIST_NAME}
    done
done
