import os
from pathlib import Path
import pytest
from pipeline_utils import run_bash, PROJECT_ROOT_PATH
from test.fixtures import test_data, tmp_dir


@pytest.mark.parametrize("bam,peaks", [
    ("a.bam", "peaks.bed"),
    ("b.bam", "peaks.bed")
])
def test_rip_sh(test_data, bam, peaks):
    bam = test_data("rip/" + bam)
    peaks = test_data("rip/" + peaks)
    rip_csv = peaks + '_rip.csv'
    if os.path.exists(rip_csv):
        os.remove(rip_csv)
    try:
        run_bash("{}/reports/rip.sh".format(PROJECT_ROOT_PATH),
                 bam, peaks)
        assert Path(rip_csv).read_text() == """file,peaks_file,reads,peaks,rip
{},{},11,1,12
""".format(bam, peaks)
    finally:
        if os.path.exists(rip_csv):
            os.remove(rip_csv)
