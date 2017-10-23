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
# D.bed
#                                                   |-|
@pytest.mark.parametrize("file1,file2,jaccard", [
    ("A.bed", "A.bed", "1.00000000000000000000"),
    ("A.bed", "B.bed", "0.33333333333333333333"),
    ("B.bed", "A.bed", "0.33333333333333333333"),
    ("B.bed", "C.bed", "0.02083333333333333333"),
    ("A.bed", "C.bed", "0.08888888888888888888"),
    ("A.bed", "D.bed", "0"),
    ("A.unsorted.bed", "B.bed", "0.33333333333333333333"),
    ("A.bed", "B.unsorted.bed", "0.33333333333333333333"),
    ("A.unsorted.bed", "B.unsorted.bed", "0.33333333333333333333"),
])
def test_intersect(capfd, test_data, file1, file2, jaccard):
    files = [file1, file2]
    run_bash("bed/jaccard.sh", *[test_data("bed/" + f) for f in files])

    out, _err = capfd.readouterr()
    res = out.replace(test_data("bed/"), "").replace(PROJECT_ROOT_PATH, ".")
    assert "bash ./bed/jaccard.sh {}\n{}\n".format(" ".join(files),
                                                   jaccard) == res
