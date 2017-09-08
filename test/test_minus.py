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
@pytest.mark.parametrize("file1,file2,line", [
    ("A.bed", "B.bed", "chr1	0	100"),
    ("A.bed", "C.bed", "chr1	600	700"),
    ("B.bed", "C.bed", "chr1	650	750"),
])
def test_minus(capfd, test_data, file1, file2, line):
    run_bash("bed/minus.sh", test_data("bed/" + file1),
             test_data("bed/" + file2))
    out, _err = capfd.readouterr()
    res = out.replace(test_data("bed/"), "").replace(PROJECT_ROOT_PATH, ".")
    assert res == "bash ./bed/minus.sh {} {}\n{}\n".format(file1, file2, line)
