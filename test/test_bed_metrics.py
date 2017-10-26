from pathlib import Path
import pandas as pd
import pytest
from test.fixtures import test_data, bedtrace_cleanup

from reports.bed_metrics import intersections_table


@pytest.mark.parametrize("jaccard,fname,swap_ab", [
    (False, "metric1.csv", False),
    (False, "metric1_ba.csv", True),
    (True, "metric_j.csv", False),
    (True, "metric_j_ba.csv", True)
])
def test_intersections_table(test_data, jaccard, fname, swap_ab):
    a_paths \
        = [Path(test_data("metrics/c{}.bed".format(i))) for i in range(1, 4)]
    b_paths = [Path(test_data("metrics/a.bed")), *a_paths]

    df = intersections_table(a_paths if not swap_ab else b_paths,
                             b_paths if not swap_ab else a_paths,
                             jaccard=jaccard)
    expected = pd.DataFrame.from_csv(test_data("metrics/{}".format(fname)))

    assert str(expected) == str(df)
