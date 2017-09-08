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
@pytest.mark.parametrize("files,line", [
    (["A.bed", "B.bed", "C.bed"], "chr1	150	500"),
    (["A.bed", "B.bed"], "chr1	150	500\nchr1	600	750"),
    (["A.bed", "C.bed"], "chr1	0	100\nchr1	200	300\nchr1	400	500"),
    (["B.bed", "C.bed"], "chr1	150	460"),
])
def test_intersect(capfd, test_data, files, line):
    run_bash("bed/intersect.sh", *[test_data("bed/" + f) for f in files])

    out, _err = capfd.readouterr()
    res = out.replace(test_data("bed/"), "").replace(PROJECT_ROOT_PATH, ".")
    assert "bash ./bed/intersect.sh {}\n{}\n".format(" ".join(files),
                                                     line) == res
