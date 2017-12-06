from pipeline_utils import run_bash, PROJECT_ROOT_PATH

import pytest
from test.fixtures import test_data


# OD1
#   |--------|
# OD2
#      |--------------|
# YD1
#               |--------------|
# YD2
#                  |-------------|
# YD3
#      |-----------------|
# 0 10 20 30 40 50 60 70 80 90  100 ..
@pytest.mark.parametrize("cons,file", [
    (2, "result_cons2.bed"),
    (3, "result_cons3.bed"),
    (4, "result_cons4.bed"),
])
def test_consensus(capfd, test_data, cons, file):
    path = test_data("consensus_diffscore/tracks")
    run_bash("bed/consensus_diffscore.sh", path, cons)

    out, err = capfd.readouterr()
    assert "" == err
    with open(test_data("consensus_diffscore/{}".format(file)), "r") as f:
        exp = "".join(f.readlines())
        assert "bash {}/bed/consensus_diffscore.sh {} {}\n{}".format(
            PROJECT_ROOT_PATH, path, cons, exp
        ) == out
