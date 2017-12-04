from pipeline_utils import run_bash, PROJECT_ROOT_PATH

import pytest
from test.fixtures import test_data


# 0      100  200    300  400    500  600    700
# A.bed
# |---------| |---------| |---------| |---------|
# B.bed
#            |------------------|              |--|
# C.bed
#  |-| |-|        |-|          |-|                  |-|
def test_metapeaks(capfd, test_data):
    run_bash("bed/metapeaks.sh", test_data("bed/A.bed"),
             test_data("bed/B.bed"), test_data("bed/C.bed"))
    out, _err = capfd.readouterr()
    res = out.replace(test_data("bed/"), "").replace(PROJECT_ROOT_PATH, ".")

    assert res == """bash ./bed/metapeaks.sh A.bed B.bed C.bed
PEAKS:\t4\t2\t5
0 0 1	1
1 0 1	1
1 1 0	1
1 1 1	1
1, 1, 1, 1
"""


def test_metapeaks_empty():
    try:
        run_bash("bed/metapeaks.sh")
        # Should fail with error code 1
        assert False
    except:  # nopep8
        pass
