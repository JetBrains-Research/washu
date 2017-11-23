import os
import shutil
import subprocess
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

@pytest.mark.parametrize("folder,consensus_folder,percent,count", [
    ("tracks_for_consensus", "consensus_weak", False, 2),
    ("tracks_for_consensus", "consensus_median", 50, False),
    ("tracks_for_consensus", "consensus_strong", 100, False)
])
def test_consensus(capfd, tmp_dir, test_data, folder, consensus_folder, percent, count):
    os.makedirs(tmp_dir + "/bed")
    for file in os.listdir(test_data("bed/" + folder)):
        shutil.copy(test_data("bed/" + folder) + "/" + file, tmp_dir + "/bed")

    args = [tmp_dir + "/bed", "H3K4me3"]
    if percent:
        args.insert(0, "-p")
        args.insert(1, percent)
    if count:
        args.insert(0, "-c")
        args.insert(1, count)

    run_bash("bed/consensus.sh", *args)

    files_suff = ["_zinbra_consensus.bed", "_zinbra_ODS_consensus.bed", "_zinbra_YDS_consensus.bed",
                  "_zinbra_ODS_without_YDS_consensus.bed", "_zinbra_YDS_without_ODS_consensus.bed"]

    out, _err = capfd.readouterr()
    res = out.replace(test_data("bed/"), "").replace(PROJECT_ROOT_PATH, ".")
    assert "bash ./bed/consensus.sh {}{}{} {}\n".format(
        ("-p %d " % percent) if percent else "",
        ("-c %d " % count) if count else "",
        tmp_dir + "/bed",
        "H3K4me3"
    ) == res

    for file_suff in files_suff:
        print("Comparing %s" % file_suff)
        with open(tmp_dir + "/bed/H3K4me3" + file_suff) as act_cons_file:
            act_cons = act_cons_file.readlines()
            with open(test_data("bed/" + consensus_folder + "/" + consensus_folder + file_suff)) as\
                    exp_cons_file:
                exp_cons = exp_cons_file.readlines()
                assert act_cons == exp_cons


@pytest.mark.parametrize("folder,percent,count", [
    ("tracks_for_consensus", 50, 2)
])
def test_incorrect_consensus_call(capfd, test_data, folder, percent, count):
    args = [test_data(folder)]
    if percent:
        args.insert(0, "-p")
        args.insert(1, percent)
    if count:
        args.insert(0, "-c")
        args.insert(1, count)
    with pytest.raises(subprocess.CalledProcessError) as e:
        run_bash("bed/consensus.sh", *args)

    out, _err = capfd.readouterr()
    assert "returned non-zero exit status 1" in str(e.value)
    assert """Calculate consensus for peaks in selected folder

Usage: consensus.sh [OPTIONS] folder_path

Options:
  -p number       Consensus percent should be taken as number (cannot be used with -c)
  -c number       Count of tracks for consensus should be taken as number
                          (cannot be used with -p)
  -h|--help       Show help

With no arguments - show help""" in out


@pytest.mark.parametrize("help_arg", [True, False])
def test_consensus_help(capfd, help_arg):
    if help_arg:
        args = ["-s", "foo", "-h"]
    else:
        args = []

    with pytest.raises(subprocess.CalledProcessError) as e:
        run_bash("bed/consensus.sh", *args)

    out, _err = capfd.readouterr()
    assert "returned non-zero exit status 1" in str(e.value)
    assert """Calculate consensus for peaks in selected folder

Usage: consensus.sh [OPTIONS] folder_path

Options:
  -p number       Consensus percent should be taken as number (cannot be used with -c)
  -c number       Count of tracks for consensus should be taken as number
                          (cannot be used with -p)
  -h|--help       Show help

With no arguments - show help""" in out
