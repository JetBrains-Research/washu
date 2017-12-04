from pipeline_utils import run_bash

import pytest
from test.fixtures import test_data


# 0      100  200    300  400    500  600    700
# A.bed
# |---------| |---------| |---------| |---------|
# B.bed
#            |------------------|              |--|
# C.bed
#  |-| |-|        |-|          |-|                  |-|
def test_union(capfd, test_data):
    run_bash("bed/union.sh", test_data('bed/A.bed'), test_data('bed/B.bed'),
             test_data('bed/C.bed'))
    out, _err = capfd.readouterr()
    res = out.replace(test_data("bed/"), "")
    expected = """chr1	0	100	1_A.bed|3_C.bed
chr1	150	500	1_A.bed|2_B.bed|3_C.bed
chr1	600	750	1_A.bed|2_B.bed
chr1	800	850	3_C.bed
"""
    assert expected in res


def test_union_empty():
    try:
        run_bash("bed/union.sh")
        # Should fail with error code 1
        assert False
    except:  # nopep8
        pass
