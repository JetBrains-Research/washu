import os
import importlib

import pandas as pd
import pytest
from test.fixtures import test_data, tmp_dir


@pytest.mark.parametrize(
    "util_name", ["macs2", "peaks", "bowtie", "bowtie2"])
def test_process_logs(tmp_dir, test_data, util_name):
    logger_fun = getattr(importlib.import_module("parallel.util.{}_logs"
                                                 .format(util_name)),
                         "process_{}_logs".format(util_name))
    logger_fun(test_data(util_name), tmp_dir)

    output_file = os.path.join(tmp_dir, util_name + "_report.csv")
    assert os.path.exists(output_file)

    expected = str(pd.read_csv(test_data("{0}/{0}_report_correct.csv"
                                         .format(util_name)),
                               sep=',')).strip()

    actual = str(pd.read_csv(output_file, sep=',')).strip()
    assert actual == expected


@pytest.mark.parametrize("util_name,log_folder,err_type,msg", [
    ("macs2", "macs2_empty_rip", "AssertionError",
     "Expected 6 comma separated values, but was 1: line = '',"
     " file = empty_rip.csv"),
])
def test_problem_report(tmp_dir, test_data, util_name, log_folder,
                        err_type, msg):
    logger_fun = getattr(importlib.import_module("parallel.util.{}_logs"
                                                 .format(util_name)),
                         "process_{}_logs".format(util_name))
    try:
        logger_fun(test_data(log_folder), tmp_dir)
    except Exception as e:
        assert str(e) == msg
        assert type(e).__name__ == err_type
    else:
        assert False, "Should throw exception: {}".format(msg)
