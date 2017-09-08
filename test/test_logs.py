import os
import importlib

import pandas as pd
import pytest
from test.fixtures import test_data, tmp_dir


@pytest.mark.parametrize(
    "util_name", ["macs2", "peaks", "bowtie", "bowtie2"])
def test_process_logs(tmp_dir, test_data, util_name):
    logger_fun = getattr(importlib.import_module("reports.{}_logs"
                                                 .format(util_name)),
                         "process_{}_logs".format(util_name))
    logger_fun(test_data(util_name), tmp_dir)

    output_file = os.path.join(tmp_dir, util_name + "_report.csv")
    assert os.path.exists(output_file)

    expected = str(pd.read_csv(test_data("{0}/{0}_report_correct.csv"
                                         .format(util_name)),
                               sep=',')).strip()

    print(expected)
    actual = str(pd.read_csv(output_file, sep=',')).strip()
    print(actual)
    assert actual == expected
