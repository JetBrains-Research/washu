import os
import sys
import re
from multiprocessing import Pool
from typing import List, Tuple
from collections import defaultdict
from itertools import chain
import argparse

from pathlib import Path
import pandas as pd
import numpy as np
import seaborn as sns

if __name__ == "__main__":
    parent_dir = os.path.dirname(os.path.realpath(__file__))
    project_root = os.path.abspath(os.path.join(parent_dir) + "/..")
    sys.path.insert(0, project_root)

    # Force matplotlib to not use any Xwindows backend.
    import matplotlib
    matplotlib.use('Agg')
    sns.set()

import matplotlib.pyplot as plt  # nopep8
from matplotlib.backends.backend_pdf import PdfPages  # nopep8
from scripts.util import run  # nopep8
from pipeline_utils import PROJECT_ROOT_PATH  # nopep8


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

    with Pool(processes=threads) as pool:
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


def color_annotator_age(label) -> Tuple[Tuple[str, str]]:
    chunks = [ch.lower() for ch in label.split("_")]

    for ch in chunks:
        if ch.startswith("od"):
            return (("age", "b"),)
        elif ch.startswith("yd"):
            return (("age", "r"),)

    return (("age", "gray"),)


def color_annotator_chain(*annotators):
    def inner(label):
        return tuple(chain(*(f(label) for f in annotators)))

    return inner


def color_annotator_outlier(outliers_df, data_type):
    outlier_mapping = outliers_df.loc[:, data_type]

    def inner(label):
        chunks = [ch.upper() for ch in label.split("_") if len(ch) > 2]

        for ch in chunks:
            if ch.startswith("OD") or ch.startswith("YD"):
                if ch in outlier_mapping:
                    value = outlier_mapping[ch]
                    if value == 0:
                        # ok
                        return (("outlier", "g"),)
                    elif value == 1:
                        # outlier
                        return (("outlier", "lightgray"),)

        # unknown
        return (("outlier", "white"),)

    return inner


def label_converter_donor_and_tool(name):
    chunks = []
    suffix_tool_map = {"Peak": "macs2", "island.bed": "sicer", "peaks.bed": "zinbra"}
    match = re.match("(?:^|.*_)([yo]d(?:s|\d+)).*", name, flags=re.IGNORECASE)

    if match:
        chunks.append(match.group(1))
    if "consensus" in name:
        chunks.append("consensus")

    for suffix, tool in suffix_tool_map.items():
        if tool in name or name.endswith(suffix):
            chunks.append(tool)

    return "_".join(chunks) if chunks else name


def plot_metric_heatmap(title, df, *, figsize=(14, 14),
                        save_to=None,
                        vmin=0, vmax=1,
                        col_color_annotator=None, row_color_annotator=None,
                        col_label_converter=None, row_label_converter=None,
                        col_cluster=False, row_cluster=False,
                        adjustments=None,
                        cbar=True,
                        **kw):

    """
    :param title: Plot title
    :param df: Dataframe with metric value
    :param figsize: Plot size
    :param save_to: Destination: 1) string path 2) Pdf obj 3) None to
           plot on the screen
    :param vmin: Heatmap min value (use None to infer from data)
    :param vmax: Heatmap max value (use None to infer from data)
    :param col_color_annotator: Function which highlights cols
    :param row_color_annotator: Function which highlights rows
    :param col_label_converter: Function which modifies cols names
    :param row_label_converter: Function which modifies rows names
    :param col_cluster: see seaborn.clustermap(..) details
    :param row_cluster: see seaborn.clustermap(..) details
    :param adjustments: Right / left / top /  bottom insets dict
    :param kw: extra arguments for easier API usages
    """
    ncol, nrow = df.shape

    if not ncol or not nrow:
        plt.figure(figsize=figsize)
    else:
        if col_label_converter or row_label_converter:
            df = df.copy()

            if col_label_converter:
                df.columns = [col_label_converter(s) for s in df.columns]

            if row_label_converter:
                df.index = [row_label_converter(s) for s in df.index]

        def as_colors_df(color_fun, items):
            if len(items) == 0:
                return None

            data = defaultdict(list)
            for item in items:
                for k, v in color_fun(item):
                    data[k].append(v)
            assert len({len(colors) for col, colors in data.items()}) == 1, \
                "All color list should be equal size:\n{}\n{}".format(
                    items,
                    {col: len(colors) for col, colors in data.items()}
                )

            df = pd.DataFrame.from_dict(data)
            df.index = items
            return df

        c_colors = None if not col_color_annotator \
            else as_colors_df(col_color_annotator, df.columns)
        r_colors = None if not row_color_annotator \
            else as_colors_df(row_color_annotator, df.index)

        # TODO: for jaccard use dist function? matrix could be not square here
        g = sns.clustermap(
            df, figsize=figsize, cmap="rainbow",
            # cbar_kws={"orientation": "horizontal", "label": "Metric"},
            col_cluster=col_cluster, row_cluster=row_cluster,
            metric="chebyshev",
            col_colors=c_colors, row_colors=r_colors,
            vmin=vmin, vmax=vmax,
            # linewidths=0.75,  # separator line width
            robust=True,  # robust=True: ignore color outliers
        )
        # Turn off color bar
        if not cbar:
            g.cax.set_visible(False)

        plt.setp(g.ax_heatmap.get_xticklabels(), rotation=-90)
        if g.row_colors is not None:
            plt.setp(g.ax_row_colors.get_xticklabels(), rotation=-90)
        plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)

    plt.title(title)
    adjustments = adjustments or {}
    plt.subplots_adjust(left=adjustments.get('left', 0.2),
                        right=adjustments.get('right', 0.8),
                        top=adjustments.get('top', 0.8),
                        bottom=adjustments.get('bottom', 0.2))
    save_plot(save_to)


def save_plot(save_to):
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


def _cli():
    parser = argparse.ArgumentParser(
        description="For given two loci paths sets build intersection heatmap pdf",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-a", required=True, action='append', metavar="PATHS",
                        help="First set: comma separated list of files. Could be used multiple "
                             "times")
    parser.add_argument("-b", action='append', required=True, metavar="PATHS",
                        help="Second set: comma separated list of files. Could be used multiple "
                             "times")
    parser.add_argument('-o', '--out', default='result', metavar="PREFIX",
                        help="Output path prefix")

    parser.add_argument('--rowc', help="Rows clustering", action="store_true")
    parser.add_argument('--colc', help="Cols clustering", action="store_true")
    parser.add_argument('-p', '--threads', help="Threads number for parallel processing",
                        type=int, default=4)
    parser.add_argument('--outliers', metavar="PATH",
                        default="/mnt/stripe/bio/experiments/aging/Y20O20.outliers.csv",
                        help="Outliers *.csv path")
    parser.add_argument('--hist', help="Histone modification name (is used for outliers "
                                       "highlighting), e.g. H3K4me3")
    parser.add_argument('--age', help="Highlight donors age", action="store_true")
    parser.add_argument('--jaccard', help="Use Jaccard metric instead intersection",
                        action="store_true")
    parser.add_argument('--size', help="Figure size: width and height separated by space",
                        type=int, nargs=2, metavar="INT",
                        default=[14, 14])

    args = parser.parse_args()

    a_paths = [Path(s.strip()) for s in chain(*(s.split(',') for s in args.a))]
    b_paths = [Path(s.strip()) for s in chain(*(s.split(',') for s in args.b))]
    prefix = args.out.strip()
    plot_path = prefix + ".pdf"
    df_path = Path(prefix + ".df")

    # Df
    df = load_or_build_metrics_table(a_paths, b_paths, df_path, jaccard=args.jaccard)

    anns = []
    # age
    if args.age:
        anns.append(color_annotator_age)
    # outliers
    hist_mod = args.hist
    outliers_df = pd.read_csv(args.outliers, delimiter="\t", skiprows=1, index_col="donor")
    if hist_mod and outliers_df is not None:
        if hist_mod in outliers_df.columns:
            anns.append(color_annotator_outlier(outliers_df, hist_mod))
    annotator = None if not anns else color_annotator_chain(*anns)

    # Do plot:
    plot_metric_heatmap("IM: {}".format(df_path.name), df,
                        save_to=plot_path,
                        row_cluster=args.rowc, col_cluster=args.colc,
                        row_color_annotator=annotator,
                        cow_color_annotator=annotator,
                        figsize=args.size,
                        threads=args.threads)


def load_or_build_metrics_table(a_paths, b_paths, df_path, **kw):
    if df_path.exists():
        df = pd.DataFrame.from_csv(str(df_path))
        print("[Skipped]: Already exists", str(df_path))
    else:
        print("Calculating metrics: ", str(df_path))
        df = bed_metric_table(a_paths, b_paths, **kw)
        df.to_csv(str(df_path))
        print("  [Saved]", str(df_path))
    return df


if __name__ == "__main__":
    _cli()
