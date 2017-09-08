import os

from bed.bedtrace import run
from pipeline_utils import run_bash, PROJECT_ROOT_PATH

import pytest
from test.fixtures import test_data


@pytest.mark.parametrize("q_target", ["0.01", "0.05"])
def test_macs2_filter_fdr(tmpdir, test_data, q_target):
    output_folder = tmpdir
    run_bash("bed/macs2_filter_fdr.sh", test_data("macs2_filter_fdr"),
             output_folder, "0.1", q_target)

    for _dirpath, _dirs, files in os.walk(output_folder):
        peaks = [f for f in files if q_target in f]
        assert len(peaks) > 0
        assert peaks[0].endswith('.broadPeak')
        with open(os.path.join(output_folder, peaks[0])) as peak:
            fst_line = [l for l in peak][1]
            assert fst_line == "chr1\t241186\t241332\tfoo_{}_peak_2\t" \
                               "20\t.\t2.82218\t3.36676\t2.05718\t" \
                               "\n".format(q_target)
