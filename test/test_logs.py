import os
import importlib

import pandas as pd

import pytest
from test.fixtures import test_data


@pytest.mark.parametrize("util_name", ["macs2", "peaks", "bowtie", "bowtie2"])
def test_process_logs(tmpdir, test_data, util_name):
    fun = getattr(importlib.import_module("reports.{}_logs".format(util_name)),
                  "process_{}_logs".format(util_name))
    fun(test_data(util_name), tmpdir)

    output_file = os.path.join(tmpdir, util_name + "_report.csv")
    assert os.path.exists(output_file)

    expected = str(pd.read_csv(test_data("{0}/{0}_report_correct.csv".format(util_name)),
                               sep=',')).strip()

    assert str(pd.read_csv(output_file, sep=',')).strip() == expected
