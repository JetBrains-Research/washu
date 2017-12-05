from pipeline_utils import run_bash, PROJECT_ROOT_PATH

import pytest
from test.fixtures import test_data


# OD1
#   |--------|
# OD2
#      |-----------|
# YD1
#                |--------------|
# YD2
#                 |-------------|
# 0 10 20 30 40 50 60 70 80 90  100 ..

def test_consensus(capfd, test_data):
    path = test_data("consensus_diffscore/tracks")
    run_bash("bed/consensus_diffscore.sh", path)

    out, err = capfd.readouterr()
    assert "" == err
    assert """bash {}/bed/consensus_diffscore.sh {}
chr	start	end	cons	od_cons	yd_cons	absdiff
chr1 10 20 1 1 0 1
chr1 20 40 2 2 0 2
chr1 40 50 1 1 0 1
chr1 50 55 2 1 1 0
chr1 55 60 3 1 2 1
chr1 60 100 2 0 2 2
chr1 1010 1020 1 1 0 1
chr1 1020 1040 2 2 0 2
chr1 1040 1050 1 1 0 1
chr1 1050 1060 2 1 1 0
chr1 1060 1070 1 0 1 1
chr1 1070 1100 2 0 2 2
""".format(PROJECT_ROOT_PATH, path) == out
