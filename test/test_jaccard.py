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
@pytest.mark.parametrize("file1,file2,jaccard,sorted,merged", [
    ("A.bed", "A.bed", "1.00000000000000000000", True, True),
    ("A.bed", "B.bed", "0.33333333333333333333", True, True),
    ("A.bed", "B.bed", "0.33333333333333333333", False, True),
    ("B.bed", "A.bed", "0.33333333333333333333", True, True),
    ("B.bed", "C.bed", "0.02083333333333333333", True, True),
    ("A.bed", "C.bed", "0.08888888888888888888", True, True),
    ("A.bed", "D.bed", "0", True, True),
    ("A.unsorted.bed", "B.bed", "0.33333333333333333333", False, True),
    ("A.bed", "B.unsorted.bed", "0.33333333333333333333", False, True),
    ("A.unsorted.bed", "B.unsorted.bed", "0.33333333333333333333", False,
     True),
    ("E.bed", "F.bed", "0.48611111111111111111", True, False),
    ("E.bed", "F.bed", "0.48611111111111111111", False, False),

    # Examples of not-valid behaviour due to lack of sort/merge:
    # 0.48611111111111111111 expected:
    ("E.bed", "F.bed", "1.22222222222222222222", False, True),
    ("E.bed", "F.bed", "1.22222222222222222222", True, True),
    # 0.33333333333333333333 expected:
    ("A.unsorted.bed", "B.bed", "0", True, True),
    ("A.bed", "B.unsorted.bed", "0", True, True),
    ("A.unsorted.bed", "B.unsorted.bed", "0", True, True),
])
def test_jaccard(capfd, test_data, file1, file2, jaccard, sorted, merged):
    files = [file1, file2]
    args = [test_data("bed/" + f) for f in files]
    if sorted:
        args.insert(0, "-s")
    if merged:
        args.insert(0, "-m")
    run_bash("bed/jaccard.sh", *args)

    out, _err = capfd.readouterr()
    res = out.replace(test_data("bed/"), "").replace(PROJECT_ROOT_PATH, ".")
    assert "bash ./bed/jaccard.sh {}{}{}\n{}\n".format(
        "-m " if merged else "",
        "-s " if sorted else "",
        " ".join(files),
        jaccard
    ) == res


@pytest.mark.parametrize("help_arg", [True, False])
def test_jaccard_help(capfd, help_arg):
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
  -s              Files already sorted, skip resort step
  -m              Files already merged, skip merge step
  -h|--help       Show help

With no arguments - show help""" in out
