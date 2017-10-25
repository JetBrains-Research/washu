import os
from multiprocessing import Pool, TimeoutError
from typing import List

from pathlib import Path
import pandas as pd
import numpy as np

from scripts.util import run


def _run_metric_intersection(a, b, *args):
    """
    Metric intersection (#1):
        (wc -l $a) / (bedtools intersect -a $a -b $b -wa | uniq | wc -l)

    :param a: A.bed
    :param b: B.bed
    :param *args: other args, may be useful for parallel exec
    :return: (metric value, *args)
    """
    output = run((["bedtools", "intersect", "-a", str(a),
                   "-b", str(b), "-wa"],
                  ["uniq"], ["wc", "-l"]))
    n_intersecting = int(output[0].decode().strip())

    output = run([["wc", "-l", str(a)]])
    n_total = int(output[0].decode().strip().split()[0])
    return (n_intersecting / n_total, *args)


def _run_metric_jaccard(a, b, *args):
    """
    :param id: Id for parallel execution
    :param a: A.bed
    :param b: B.bed
    :return: Jaccard index by `~/work/washu/bed/jaccard.sh $a $b`
    """
    script = os.path.expanduser("~/work/washu/bed/jaccard.sh")
    output = run([["bash", script, str(a), str(b)]])
    stdout = output[0].decode().strip()
    return (float(stdout), *args)


def intersections_table(a_paths: List[Path], b_paths: List[Path],
                        jaccard=False,
                        threads=4, timeout_hours=10):

    path_pairs = []
    for i, a in enumerate(a_paths, 0):
        for j, b in enumerate(b_paths, 0):
            path_pairs.append((a, b, (i, j)))

    inner_metric = _run_metric_jaccard if jaccard else _run_metric_intersection

    pool = Pool(processes=threads)
    multiple_results = [pool.apply_async(inner_metric,
                                         (a, b, ij)) for a, b, ij in
                        path_pairs]
    values = [res.get(timeout=3600 * timeout_hours) for res in
              multiple_results]

    x = np.zeros((len(a_paths), len(b_paths)), np.float32)
    for value, (i, j) in values:
        x[i, j] = value

    df = pd.DataFrame(x,
                      index=[f.name for f in a_paths],
                      columns=[f.name for f in b_paths])
    # for i, a in enumerate(a_paths, 0):
    #     output = run((["cat", str(a)], ["wc", "-l"],))
    #     x[i, 0] = int(output[0].decode().strip())
    #
    # df = pd.DataFrame(x,
    #                   index=[f.name for f in a_paths],
    #                   columns=["total"] + [f.name for f in b_paths])
    return df
