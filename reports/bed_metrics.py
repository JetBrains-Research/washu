import os
import sys
from multiprocessing import Pool, TimeoutError
from typing import List, Tuple
from collections import defaultdict

from pathlib import Path
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from scripts.util import run
from pipeline_utils import PROJECT_ROOT_PATH


def _try_parse_stdout(func, stdout, stderr):
    try:
        return func(stdout, stderr)
    except ValueError:
        print("Cannot parse stdout value:\n"
              "--- STDOUT -------\n"
              "{}\n"
              "--- STDERR -------\n"
              "{}\n"
              "------------------\n".format(stdout, stderr),
              file=sys.stderr)
        raise


def _run_metric_intersection(a, b, *args, **kw):
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
    n_intersecting = _try_parse_stdout(
        lambda stdout, stderr: int(stdout.decode().strip()),
        *output
    )

    output = run([["wc", "-l", str(a)]])
    n_total = _try_parse_stdout(
        lambda stdout, stderr: int(stdout.decode().strip().split()[0]),
        *output
    )

    if n_total == 0:
        print("Warning: Bed file is empty:", str(a), file=sys.stderr)
        metric = 0
    else:
        metric = n_intersecting / n_total

    return (metric, *args)


def _run_metric_jaccard(a, b, *args, **kw):
    """
    :param id: Id for parallel execution
    :param a: A.bed
    :param b: B.bed
    :return: Jaccard index by `~/work/washu/bed/jaccard.sh $a $b`
    """
    script = os.path.join(PROJECT_ROOT_PATH, "bed/jaccard.sh")
    cmdline = ["bash", script, str(a), str(b)]
    if kw.get("sorted", False):
        cmdline.append("-s")
    if kw.get("merged", False):
        cmdline.append("-m")
    output = run([cmdline])

    return (_try_parse_stdout(
        lambda stdout, stderr: float(stdout.decode().strip()),
        *output
    ), *args)


def bed_metric_table(a_paths: List[Path], b_paths: List[Path],
                     jaccard=False,
                     threads=4, timeout_hours=10, **kw):
    """
    :param a_paths: First paths set
    :param b_paths:  Second paths set
    :param jaccard: If True use Jaccard metric, else Intersection metric
    :param threads: Threads number for parallel execution
    :param timeout_hours: Execution timeout in hours
    :param kw: Additional args, e.g.: "sorted" or "merged" for jaccard
    :return: Dataframe selected metric for A x B
    """
    path_pairs = []
    for i, a in enumerate(a_paths, 0):
        for j, b in enumerate(b_paths, 0):
            path_pairs.append((a, b, (i, j)))

    inner_metric = _run_metric_jaccard if jaccard else _run_metric_intersection

    pool = Pool(processes=threads)
    multiple_results = [pool.apply_async(inner_metric, (a, b, ij), kw)
                        for a, b, ij in path_pairs]
    values = [res.get(timeout=3600 * timeout_hours) for res in
              multiple_results]

    x = np.zeros((len(a_paths), len(b_paths)), np.float32)
    for value, (i, j) in values:
        x[i, j] = value

    df = pd.DataFrame(x,
                      index=[f.name for f in a_paths],
                      columns=[f.name for f in b_paths])
    return df


# def load_intersection_table(beds, loci, path_to_bt, result_path, *kw):
#     if result_path.exists():
#         df = pd.DataFrame.from_csv(result_path)
#         print("Loaded: ", result_path)
#     else:
#         print("Calculating: ", result_path)
#         df = intersections_table(beds, loci, *kw)
#         result_path.parent.mkdir(parents=True, exist_ok=True)
#         df.to_csv(str(result_path))
#         print("  Saved: ", result_path)
#
#     return df


def heatmap_donor_color_fun(label) -> Tuple[Tuple[str, str]]:
    chunks = [ch.lower() for ch in label.split("_")]

    for chunk in chunks:
        ch = chunk.lower()
        if ch.startswith("od"):
            return (("age", "b"),)
        elif ch.startswith("yd"):
            return (("age", "r"),)

    return (("age", "gray"),)


def plot_metric_heatmap(title, df, figsize=(10, 10),
                        save_to=None,
                        vmin=0, vmax=1,
                        col_color_fun=None, row_color_fun=None,
                        col_label_fun=None, row_label_fun=None,
                        col_cluster=False, row_cluster=False):

    """
    :param title: Plot title
    :param df: Dataframe with metric value
    :param figsize: Plot size
    :param save_to: Destination: 1) string path 2) Pdf obj 3) None to
           plot on the screen
    :param vmin: Heatmap min value (use None to infer from data)
    :param vmax: Heatmap max value (use None to infer from data)
    :param col_color_fun: Function which highlights cols
    :param row_color_fun: Function which highlights rows
    :param col_label_fun: Function which modifies cols names
    :param row_label_fun: Function which modifies rows names
    :param col_cluster: see seaborn.clustermap(..) details
    :param row_cluster: see seaborn.clustermap(..) details
    """
    ncol, nrow = df.shape

    if not ncol or not nrow:
        plt.figure(figsize=figsize)
    else:
        if col_label_fun or row_label_fun:
            df = df.copy()

            if col_label_fun:
                df.columns = [col_label_fun(s) for s in df.columns]

            if row_label_fun:
                df.index = [row_label_fun(s) for s in df.index]

        def as_colors_df(color_fun, items):
            if len(items) == 0:
                return None

            data = defaultdict(list)
            for item in items:
                for k, v in color_fun(item):
                    data[k].append(v)
            assert len({len(colors) for col, colors in data.items()}) == 1,\
                "All color list should be equal size:\n{}\n{}".format(
                    items,
                    {col: len(colors) for col, colors in data.items()}
                )

            df = pd.DataFrame.from_dict(data)
            df.index = items
            return df

        c_colors = None if not col_color_fun else as_colors_df(col_color_fun,
                                                               df.columns)
        r_colors = None if not row_color_fun else as_colors_df(row_color_fun,
                                                               df.index)

        # TODO: for jaccard use dist function? matrix could be not square here
        g = sns.clustermap(
            df, figsize=figsize, cmap="rainbow",
            col_cluster=col_cluster, row_cluster=row_cluster, metric="chebyshev",
            col_colors=c_colors, row_colors=r_colors,
            vmin=vmin, vmax=vmax,
            robust=True,  # robust=True: ignore color outliers
        )

        plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)

    plt.title(title)
    if save_to is None:
        plt.show()
    elif isinstance(save_to, str):
        plt.savefig(str(save_to))
        plt.close()
    elif isinstance(save_to, PdfPages):
        save_to.savefig()
        plt.close()
    else:
        raise ValueError("Unsupported value type: {}".format(type(save_to)))
