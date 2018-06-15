import os
from pipeline_utils import run_bash, PROJECT_ROOT_PATH

import pytest
from test.fixtures import test_data, tmp_dir


# 0 10 20 30 40 50 60 70 80 90    110 120 130 140 150 160 170 180 190 200 210 220 230
# OD1_peaks.bed
#      |---------------|
# OD2_peaks.bed
#            |--------------------------|
# OD3_peaks.bed
#         |---------------|
# OD4_peaks.bed
#      |--------------------------------|
# YD1_peaks.bed
# |--------------------|                          |---------------------|
# YD2_peaks.bed
#   |---------------|                                 |-------------|
# YD3_peaks.bed
#      |------------------|                               |---------------------|
# YD4_peaks.bed
#         |------------------|                |-----------------------------|
# YD5_peaks.bed
#            |------------|                       |-------------|
# YD6_peaks.bed
#               |-------------------|     |-----------------------------------------|
# 0 10 20 30 40 50 60 70 80 90    110 120 130 140 150 160 170 180 190 200 210 220 230
# 290 300   330   350   370   390   420   440   460 470   490   520   540
# OD1_peaks.bed
#     |-------------|
# OD2_peaks.bed
#
# OD3_peaks.bed
# |-----------|
# OD4_peaks.bed
#
# YD1_peaks.bed
#                                                         |-------------|
# YD2_peaks.bed
#                       |-------------------------|
# YD3_peaks.bed
#                                                   |-------------|
# YD4_peaks.bed
#                       |-------------------|
# YD5_peaks.bed
#
# YD6_peaks.bed
#                             |-------|
# 290 300   330   350   370   390   420   440   460 470   490   520   540

@pytest.mark.parametrize("folder,consensus_file,percent,count", [
    ("tracks_for_consensus", "consensus/consensus_weak_span_consensus.bed", False, 2),
    ("tracks_for_consensus", "consensus/consensus_median_span_consensus.bed", 50, False),
    ("tracks_for_consensus", "consensus/consensus_strong_span_consensus.bed", 100, False)
])
def test_consensus(capfd, test_data, folder, consensus_file, percent, count):
    files = []
    for file in os.listdir(test_data("bed/" + folder)):
        files.append(test_data("bed/" + folder) + "/" + file)

    args = [' '.join(files)]
    if percent:
        args.insert(0, "-p")
        args.insert(1, percent)
    if count:
        args.insert(0, "-c")
        args.insert(1, count)

    run_bash("bed/consensus.sh", *args)

    out, _err = capfd.readouterr()
    res = out.replace(test_data("bed/"), "").replace(PROJECT_ROOT_PATH, ".")
    with open(test_data("bed/" + consensus_file), 'r') as consensus:
        assert "bash ./bed/consensus.sh {}{}{}\n{}".format(
            ("-p %d " % percent) if percent else "",
            ("-c %d " % count) if count else "",
            ' '.join(files).replace(test_data("bed/"), "").replace(PROJECT_ROOT_PATH, "."),
            consensus.read()
        ) == res
