import os

from pathlib import Path
from pipeline_utils import run_bash

import pytest
from test.fixtures import test_data, tmp_dir


# 0      100  200    300  400    500  600    700
# A.bed
# |---------| |---------| |---------| |---------|
# B.bed
#            |------------------|              |--|
# C.bed
#  |-| |-|        |-|          |-|                  |-|
@pytest.mark.parametrize("cond,expected", [
    ("cond1", "chr1	0	100\n"),
    ("cond2", ""),
    ("common", "chr1	150	500\nchr1	600	750\n"),
])
def test_compare(tmp_dir, test_data, cond, expected):
    os.chdir(tmp_dir)
    run_bash("bed/compare.sh", test_data("bed/A.bed"), test_data("bed/B.bed"),
             "metabeds")

    result_path = Path(tmp_dir) / "metabeds_{}.bed".format(cond)
    assert result_path.read_text() == expected
