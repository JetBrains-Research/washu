import subprocess
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
@pytest.mark.parametrize("file1,file2,jaccard,resort", [
    ("A.bed", "A.bed", "1.00000000000000000000", False),
    ("A.bed", "B.bed", "0.33333333333333333333", False),
    ("A.bed", "B.bed", "0.33333333333333333333", True),
    ("B.bed", "A.bed", "0.33333333333333333333", False),
    ("B.bed", "C.bed", "0.02083333333333333333", False),
    ("A.bed", "C.bed", "0.08888888888888888888", False),
    ("A.bed", "D.bed", "0", False),
    ("A.unsorted.bed", "B.bed", "0", False),
    ("A.bed", "B.unsorted.bed", "0", False),
    ("A.unsorted.bed", "B.unsorted.bed", "0", False),
    ("A.unsorted.bed", "B.bed", "0.33333333333333333333", True),
    ("A.bed", "B.unsorted.bed", "0.33333333333333333333", True),
    ("A.unsorted.bed", "B.unsorted.bed", "0.33333333333333333333", True),
    ("E.bed", "F.bed", "0.48611111111111111111", False),
])
def test_intersect(capfd, test_data, file1, file2, jaccard, resort):
    files = [file1, file2]
    args = [test_data("bed/" + f) for f in files]
    if not resort:
        args.insert(0, "-s")
    run_bash("bed/jaccard.sh", *args)

    out, _err = capfd.readouterr()
    res = out.replace(test_data("bed/"), "").replace(PROJECT_ROOT_PATH, ".")
    assert "bash ./bed/jaccard.sh {}{}\n{}\n".format(
        "" if resort else "-s ",
        " ".join(files),
        jaccard
    ) == res


@pytest.mark.parametrize("help_arg", [ True, False ])
def test_intersect_help(capfd, help_arg):
    if help_arg:
        args = ["-s", "foo", "-h"]
    else:
        args = []

    with pytest.raises(subprocess.CalledProcessError) as e:
        run_bash("bed/jaccard.sh", *args)

    out, _err = capfd.readouterr()
    assert "returned non-zero exit status 1" in str(e.value)
    assert """Calculate Jaccard Index for A.bed B.bed

Usage: jaccard.sh [OPTIONS] A.bed B.bed

Options:
  -s              Files already sorted, do not resort them before calculations
  -h|--help       Show help

With no arguments - show help""" in out