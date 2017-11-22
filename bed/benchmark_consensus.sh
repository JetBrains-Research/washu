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

mkdir -p ${LOCI_ROOT}
for HIST_DIR in $(find ${BENCHMARK_ROOT} -maxdepth 1 -type d -name "H*"); do
    HIST_NAME=$(basename ${HIST_DIR})

    # We decided to exclude MACS2 narrow peaks
    # Sicer not ready yet, let's exclude it too
    cd ${HIST_DIR}
    # exclude from find results: same dir ('.'), and some tools
    for DIR in $(find . -maxdepth 1 -type d ! -path . ! -path "./macs_narrow" ! -path "./sicer"); do
        echo "Processing: $(expand_path ${DIR})"
        bash /home/user/work/tsurinov/washu/bed/consensus.sh -p 50 ${DIR} ${HIST_NAME}
        # bash "$(project_root_dir)/bed/consensus.sh" -p 50  ${DIR} ${HIST_NAME}
        find ${DIR} -maxdepth 1 -name "*_consensus.bed" | xargs -I fname cp fname ${LOCI_ROOT}
    done
done

>&2 echo "Done. benchmark_consensus $@"